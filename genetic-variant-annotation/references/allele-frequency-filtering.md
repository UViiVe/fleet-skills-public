# Allele frequency filtering — gnomAD, TOPMed, and population strategy

Population allele frequencies are the single most powerful filter in any
variant prioritisation pipeline. Almost every clinical, rare-disease, or
target-discovery workflow filters first on "rare in unaffected
populations" — usually with a gnomAD threshold — and only then layers on
consequence and pathogenicity. This guide collects the thresholds, the
population-stratification subtleties, and the handling of missing
frequencies.

---

## Reference databases

| Database | Samples | Ancestry coverage | Notes |
|---|---|---|---|
| gnomAD v4 (exomes + genomes) | ~807 K total | African, Latino, East Asian, Finnish, non-Finnish European, South Asian, Ashkenazi Jewish, Middle Eastern, "remaining" | Default for human variant filtering |
| gnomAD v2 (exomes + genomes) | ~141 K total | Same population set | Older; use only if you need GRCh37 coordinates |
| TOPMed Bravo | ~133 K WGS | Diverse U.S. population | Better African ancestry coverage; rich for genomes |
| 1000 Genomes Phase 3 | ~2.5 K | AFR / AMR / EAS / EUR / SAS | Legacy; superseded by gnomAD for most uses |
| UK10K | ~3.8 K | British | Niche; specific UK population reference |

For most workflows gnomAD v4 GRCh38 is the default. TOPMed is a strong
complement when your cohort has significant African or Latin American
ancestry — gnomAD's coverage of those populations is improved over v2
but still asymmetric.

---

## Standard thresholds by use case

| Use case | Maximum allele frequency | Rationale |
|---|---|---|
| Mendelian dominant disease | < 1e-4 (≤ 0.01 %) | The disease itself is rare; carriers should be rare |
| Mendelian recessive disease | < 0.01 (≤ 1 %) | Carriers can reach percent-level frequencies even for rare disease |
| X-linked dominant | < 1e-4 in the relevant sex-specific population | Sex-specific frequencies matter |
| Cancer somatic | < 0.001 in gnomAD | Hard cutoff to exclude common germline; somatic variants can be at any tumour VAF |
| Variant prioritisation for novel-target discovery | < 0.01 | Standard rare-variant burden test threshold |
| Population genetics (no rarity filter) | none | Do not filter on AF; the AF itself is the signal |

These are starting points; loosen or tighten in either direction with an
explicit decision. Document the threshold in any clinical report.

The `--max_af gnomAD` flag in VEP captures the maximum allele frequency
across all gnomAD populations and is the single most useful aggregate
field. Use it as the primary filter; drop to per-population frequencies
only when there is an ancestry-specific question.

---

## Population stratification

Allele frequencies are not uniform across populations. A variant that is
rare globally can be common in a specific population, and vice versa.
gnomAD reports per-population AFs that you should consult when the cohort
ancestry is known.

### When to filter on per-population AF

- The cohort is ancestry-homogeneous (e.g. a Finnish biobank — use the
  Finnish-specific AF).
- A specific population is known to harbour founder variants (e.g.
  Ashkenazi Jewish BRCA1 founder mutations).
- The disease itself has strong ancestry stratification (sickle cell,
  cystic fibrosis).

### When to use the maximum AF across populations

- Cohort ancestry is heterogeneous, unknown, or globally diverse.
- You want a conservative rarity filter that does not let a variant slip
  through because it happens to be rare in the cohort's home population.
- This is the default in most clinical pipelines (`--max_af gnomAD`).

### Sex-specific frequencies for X chromosome

For X-linked variants in males (hemizygous) use the male-specific AF;
females carrying a single copy use the female-specific AF. gnomAD reports
`gnomAD_AF_male` and `gnomAD_AF_female` for this purpose. Mixing them up
is a common error in X-linked recessive analyses.

---

## Missing frequencies — how to handle

A variant absent from gnomAD has no allele frequency. There are two
common interpretations:

1. **Truly novel** — the variant has never been observed in 807 K humans,
   which is a strong signal that it is either ultra-rare or the result
   of a sequencing artefact in your cohort.
2. **Artefact** — your variant caller produced a false-positive call at
   a position that is hard to sequence; gnomAD's blacklist of low-quality
   regions excludes the same position for QC reasons.

Practical handling:

- `missing_as_rare=True` — treat absent AF as `0` (default in most
  pipelines, including the skill's `filter_by_frequency()`). Keeps novel
  variants in the candidate set.
- `missing_as_rare=False` — drop variants without a recorded AF. Strict;
  useful only when you have strong reason to believe most absent calls
  are artefacts.

The right answer depends on the quality of the upstream variant calling.
If your pipeline produces a lot of low-quality calls (e.g. amplicon
sequencing with low coverage), favour the strict treatment plus a
manual inspection of the novel calls.

---

## Compound heterozygosity and zygosity

For autosomal recessive disease the relevant filter is not just "rare"
but "rare and biallelic in the proband". Two failure modes:

- **Heterozygous AF > 1 %** but biallelic carriers in the patient — keep
  if the gene is known recessive and there is a second hit in trans.
- **Two heterozygous variants in cis** — false compound het. Use phasing
  (read-based, family-based, or population-based) before classifying.

`filter_by_zygosity(mode='recessive', require_compound_het=True)` in
the skill's filter set handles the basic logic; trans-phasing is the
hard part and usually needs parental sequencing.

---

## Population-frequency edge cases

- **Pseudogenes and segmental duplications.** Variant calling in
  paralogous regions produces inflated AFs because reads from multiple
  genomic locations collapse to one. gnomAD's QC removes most but not
  all of these; cross-check with `--custom`-loaded segmental duplication
  tracks if your candidates land in known paralogous regions.
- **Centromeric and telomeric repeats.** Same problem; gnomAD masks the
  worst regions but the boundaries are conservative.
- **Mitochondrial AF.** gnomAD reports mtDNA frequencies separately
  (`gnomAD_AF_mt`) with a different denominator. Do not mix nuclear and
  mitochondrial thresholds.
- **Cohort overlap.** If your discovery cohort is part of gnomAD
  (e.g. UK Biobank participants), the AF is inflated by your own cases.
  This is rarely a problem for rare-disease work but can bias case-
  enriched somatic studies.

---

## Per-population thresholds

When ancestry is known and homogeneous, per-population AF cutoffs are
often stricter than the maximum-across-populations cutoff:

| Population | gnomAD field | Notes |
|---|---|---|
| African / African-American | `AF_afr` | Higher overall genetic diversity; same threshold may admit more variants |
| Latino / Admixed American | `AF_amr` | Moderate diversity |
| East Asian | `AF_eas` | Distinct allele frequency landscape from European |
| Finnish | `AF_fin` | Founder population; many alleles enriched relative to other Europeans |
| Non-Finnish European | `AF_nfe` | The "majority" in most clinical databases |
| South Asian | `AF_sas` | Distinct landscape; under-represented historically |
| Ashkenazi Jewish | `AF_asj` | Strong founder effects; many disease-associated variants enriched |
| Middle Eastern | `AF_mid` | Newer population in gnomAD v4 |
| Remaining | `AF_remaining` | Catch-all for unclassified ancestry |

For founder-population queries, **always** use the founder-specific AF.
A canonical Ashkenazi BRCA1 founder mutation has AF ≈ 1 % among AJ but
is essentially absent in NFE; filtering on the cross-population maximum
will incorrectly include it as common.

---

## Putting it together — recommended default

For a generic rare-disease or rare-variant analysis on a mixed-ancestry
cohort:

```python
filter_by_frequency(
    df,
    max_af=0.01,                  # 1 % across all gnomAD populations
    population="gnomAD_AF",
    missing_as_rare=True,         # keep novel variants
)
```

For a dominant Mendelian disease with strict requirements:

```python
filter_by_frequency(
    df,
    max_af=1e-4,
    population="gnomAD_AF",
    missing_as_rare=True,
)
```

For an ancestry-specific cohort (e.g. all-Finnish):

```python
filter_by_frequency(
    df,
    max_af=1e-4,
    population="gnomAD_AF_fin",
    missing_as_rare=True,
)
```

Document the chosen threshold and population in any report. The
threshold is a tunable parameter, not a constant of nature.
