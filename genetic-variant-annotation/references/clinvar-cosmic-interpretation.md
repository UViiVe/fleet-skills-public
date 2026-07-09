# ClinVar, COSMIC, and clinical-significance interpretation

This guide explains how to translate clinical and somatic databases into a
defensible classification per variant. ClinVar carries the germline
clinical-significance landscape and is the entry point for ACMG / AMP
interpretation; COSMIC carries the somatic cancer mutation landscape with
recurrence counts and oncogenicity hotspots. The two databases overlap
but their interpretation frameworks are distinct.

---

## ClinVar — clinical-significance categories

ClinVar (NCBI) reports a clinical-significance assertion per variant per
submitter, plus an aggregated assertion when assertions disagree. The
controlled vocabulary is:

| Term | Abbreviation | Meaning |
|---|---|---|
| Pathogenic | P | Disease-causing |
| Likely pathogenic | LP | Disease-causing with > 90 % posterior |
| Uncertain significance | VUS | Insufficient evidence for either class |
| Likely benign | LB | Benign with > 90 % posterior |
| Benign | B | Not disease-causing |
| Drug response | — | Affects drug response (pharmacogenetic) |
| Risk factor | — | Modifies risk; not Mendelian causation |
| Conflicting interpretations | — | Submitters disagree |
| Not provided | — | No assertion |

Read the **review status** alongside the assertion. The "stars" indicate
the strength of curation:

| Stars | Meaning |
|---|---|
| ★★★★ | Practice guideline (e.g. ACMG SF list) |
| ★★★ | Expert panel reviewed |
| ★★ | Multiple submitters with no conflicts |
| ★ | Single submitter, criteria provided |
| 0 | No assertion criteria provided |

For clinical reporting, evidence at ★★ or above is generally usable
without further investigation; ★ or 0 requires you to inspect the
underlying submission and may need orthogonal evidence.

---

## ACMG/AMP framework (germline)

The American College of Medical Genetics and Association for Molecular
Pathology framework (Richards et al. 2015) is the canonical scaffolding
for translating raw evidence into a clinical classification. It defines
21 evidence codes that combine into a five-tier classification.

### Pathogenic evidence

| Code | Strength | Description |
|---|---|---|
| PVS1 | Very strong | Null variant in a gene where LoF is a known disease mechanism |
| PS1–PS4 | Strong | Same amino acid as a known pathogenic variant (PS1), de novo (PS2), functional studies (PS3), prevalence in cases (PS4) |
| PM1–PM6 | Moderate | Mutational hotspot or critical domain (PM1), absent from controls (PM2), in trans with pathogenic for recessive (PM3), protein length change (PM4), novel missense at known-pathogenic position (PM5), assumed de novo (PM6) |
| PP1–PP5 | Supporting | Co-segregation (PP1), missense in a gene with low benign-missense rate (PP2), computational evidence (PP3), specific phenotype (PP4), reputable source asserts pathogenicity (PP5) |

### Benign evidence

| Code | Strength | Description |
|---|---|---|
| BA1 | Stand-alone | Allele frequency > 5 % in any population |
| BS1–BS4 | Strong | Frequency higher than expected (BS1), seen in healthy individuals (BS2), functional studies (BS3), lack of segregation (BS4) |
| BP1–BP7 | Supporting | Missense in a gene where missense is rarely pathogenic (BP1), in trans with dominant pathogenic (BP2), in-frame indel in repeat (BP3), computational predicts benign (BP4), other-gene explanation exists (BP5), reputable source asserts benign (BP6), silent variant with no splice impact (BP7) |

### Classification combinations

- **Pathogenic** — `1× PVS1 + ≥1× PS` OR `≥2× PS` OR `1× PS + ≥3× PM` OR `1× PS + 2× PM + ≥2× PP` OR `1× PS + 1× PM + ≥4× PP` OR `≥3× PM` OR `2× PM + ≥3× PP` OR `1× PM + ≥4× PP`
- **Likely pathogenic** — `1× PVS1 + 1× PM` OR `1× PS + 1–2× PM` OR `1× PS + ≥2× PP` OR `≥2× PM` OR `1× PM + ≥2× PP` OR `≥4× PP`
- **Benign** — `1× BA1` OR `≥2× BS`
- **Likely benign** — `1× BS + 1× BP` OR `≥2× BP`
- **VUS** — anything that does not satisfy the criteria above

The skill's `prioritize_variants.py` runs a simplified version of this
logic using the populated evidence columns. Full ACMG curation requires
expert review and access to family / clinical data that VCF annotation
cannot supply.

---

## COSMIC — somatic mutation interpretation

COSMIC (Catalogue Of Somatic Mutations In Cancer) is structured around
the *recurrence* of somatic mutations across tumour samples. Two
identifiers matter:

- **COSMIC mutation ID** (`COSV`...) — a stable identifier for a single
  protein-level change at a specific genomic position.
- **COSMIC sample count** — how many tumour samples carry the mutation.
  High counts (≥ 50) flag well-established cancer mutations; low counts
  (1–5) are weak evidence.

### Driver vs passenger

A handful of summary fields condense the recurrence signal:

| Field | What it tells you |
|---|---|
| `Tier 1` | Genes with extensive evidence for driver role (the Cancer Gene Census top tier) |
| `Tier 2` | Genes with suggestive driver evidence |
| `mutation_count` | Number of independent tumour samples |
| `hotspot` | Boolean — does the variant sit at a recurrent hotspot |
| `oncogenicity` | Curated assertion (Oncogenic / Likely oncogenic / Uncertain / Benign) — only populated for a subset of variants |

For somatic interpretation use **mutation_count ≥ 5 within Tier 1 or
Tier 2 genes** as a practical first-pass driver filter. Combine with:

- VAF (variant allele fraction) > 5 % in the tumour
- Population frequency (gnomAD) < 0.001 to exclude germline contamination
- ACMG/AMP for cancer or AMP/ASCO/CAP-2017 tier when you need a clinical
  oncology classification

### COSMIC vs ClinVar — they are not interchangeable

Use ClinVar for germline variants. Use COSMIC for tumour somatic
mutations. Cross-checking is fine — a known pathogenic ClinVar variant
that also appears as a COSMIC hotspot in cancer is unusual but not
impossible (e.g. `BRCA1`/`BRCA2` germline + somatic). But do not apply
ClinVar's "Pathogenic" label to a somatic call, and do not apply COSMIC
recurrence to a germline call.

---

## Common interpretation pitfalls

- **VUS does not mean benign.** It means "we don't know". The variant
  should remain on the candidate list pending additional evidence.
- **Single-submitter Pathogenic with 0 stars** is weaker than VUS with
  expert-panel review. Read the stars, not just the term.
- **Conflicting interpretations** require human review. Inspect each
  submitter's assertion and review-status.
- **Multiple ClinVar IDs per genomic position.** A position can carry
  several variants; do not collapse without checking REF / ALT.
- **Treating recurrence counts as evidence of pathogenicity.** A
  somatic recurrence count is evidence the mutation is *selected* in
  cancer, not that it causes disease in germline carriers.
- **Pharmacogenetic ClinVar entries** ("Drug response") are about
  drug metabolism, not disease causation. They should not appear in a
  Mendelian disease report unless the gene is also a Mendelian gene.

---

## Workflow snippet

A pragmatic pipeline for a clinical report:

1. **Annotate** with VEP `--everything` + `--custom ClinVar.vcf.gz` +
   relevant pathogenicity plugins.
2. **Filter** to coding consequences (HIGH + MODERATE) and rare
   population frequency (gnomAD AF < 0.01).
3. **Stratify by ClinVar significance** — separate piles for
   Pathogenic / Likely pathogenic, VUS, Likely benign / Benign.
4. **Annotate ACMG codes** computationally for the P/LP and VUS piles
   (PVS1 from consequence + LoF-mechanism gene list, PM2 from gnomAD,
   PP3 from pathogenicity predictors).
5. **Manual review** of the P/LP/VUS pile against clinical context and
   family / segregation data.
6. **For cancer reports**, add a parallel COSMIC stratification on
   the tumour-only somatic call set and apply AMP/ASCO/CAP tiers
   rather than ACMG.
