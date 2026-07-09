# MIQE 2.0 Guidelines for qPCR

MIQE stands for *Minimum Information for Publication of Quantitative Real-Time
PCR Experiments*. The 2025 revision (MIQE 2.0) is the authoritative checklist
for any qPCR result that goes into a manuscript, regulatory submission, or
clinical assay.

This file translates MIQE 2.0 into the design and validation decisions you
actually need to make. The Standard Workflow in [SKILL.md](../SKILL.md) flags
which steps are MIQE-essential.

Citation: Bustin SA et al. *Clinical Chemistry* 71(6):634–660 (2025).

← Back to [SKILL.md](../SKILL.md)
Related: [Parameter Ranges](parameter_ranges.md) · [Best Practices](primer_design_best_practices.md) · [Code Examples](code_examples.md)

---

## Essential vs Desirable

MIQE marks each item as Essential (E) or Desirable (D).

- **Essential** items belong in the Methods of every publication. If your
  reviewer is paying attention, missing essentials get the paper rejected.
- **Desirable** items strengthen the work but are not strictly required.

This document tags each item with **[E]** or **[D]** where relevant.

---

## 1. Experimental Design

### Purpose **[E]**

State the biological question and the qPCR application explicitly. Gene
expression quantification, copy-number variation, microbial load, and
genotyping each have different validation expectations.

### Replication **[E]**

- **Biological replicates:** minimum 3 independent samples per condition
- **Technical replicates:** 2 minimum, 3 is the practical default
- **Controls:** positive control (known template), negative control (no
  primer), no-template control (water instead of template); also a no-RT
  control for RT-qPCR
- **Sample size justification:** power calculation or pilot data

### Sample description **[E]**

Document the source: cell line and passage, tissue and dissection method,
organism, treatment conditions, collection method, storage temperature and
duration, processing time (sample → flash-freeze).

---

## 2. Sample Quality

### RNA quality **[E for RT-qPCR]**

- **Integrity** — measure RIN (Bioanalyzer) or RQN (Fragment Analyzer):
  - RIN ≥ 8 — excellent
  - RIN 7–8 — good
  - RIN 5–7 — marginal; document and justify
  - RIN < 5 — generally not usable; consider amplification before qPCR
- **Purity** — 260/280 between 1.8 and 2.0; 260/230 between 2.0 and 2.2
- **Quantification** — Qubit (preferred) or NanoDrop (spectrophotometric)

### DNA quality **[E for qPCR on gDNA]**

- Gel: a single high-molecular-weight band
- 260/280 between 1.8 and 2.0
- Qubit dsDNA quantification (not absorbance)

### Reverse transcription **[E for RT-qPCR]**

Disclose enzyme name and supplier; priming strategy (oligo-dT vs random
hexamer vs gene-specific); reaction temperature, time, and volume; RNA input
amount (commonly 100 ng – 1 µg).

Run RT controls every batch:
- **+RT** — the standard RT reaction
- **−RT** — same reagents, no enzyme; catches gDNA contamination
- **NTC** — water instead of RNA

---

## 3. Assay Design

### Primer design **[E]**

Publication-quality methods sections must include:
- Forward and reverse primer sequences (5'→3')
- Gene symbol and RefSeq accession with exon positions
- Amplicon length (target 70–140 bp)
- Primer length (18–22 nt)
- Tm range (58–62°C, matched within 2°C)
- GC content (40–60%)

Pre-filling these from the design output is exactly what
`scripts/export_results.py:_export_miqe_checklist` produces.

### Exon-spanning **[E for RT-qPCR]**

The amplicon should cross an exon–exon junction, **or** the intervening intron
should be > 1 kb so that any contaminating gDNA simply won't amplify within
the qPCR cycle time. Document which strategy you used and the exon positions.

### TaqMan probe **[E for TaqMan]**

- Sequence disclosed
- Length 18–30 nt
- Tm 65–70°C (5–10°C above primers)
- Fluorophore and quencher named (FAM/TAMRA, FAM/BHQ1, etc.)
- Probe positioned between primers, no 5' G

### In-silico specificity **[E]**

A BLAST search against the target species' RefSeq must show the target as the
top hit with no other amplifiable products. Save the search parameters
(database, mismatch limit) and the top off-target if any.

---

## 4. Assay Validation

### Specificity **[E]**

SYBR Green: single sharp melt-curve peak. Multiple peaks or broad shoulders
signal non-specific products or dimers.

Agarose gel: single band at the expected size.

Sequencing (Desirable): Sanger the qPCR product to confirm identity. Required
in some journals for novel assays.

### Amplification efficiency **[E]**

Build a standard curve:
- 5–7 dilution points
- 5- to 10-fold dilution steps
- Cover the range of expected sample concentrations

Compute efficiency:
```
E = 10^(−1/slope) − 1
```

Accept 0.90–1.10 (90–110%). Outside that range:
- E < 0.90 — inhibition, secondary structure, poor design
- E > 1.10 — primer dimers, contamination, or pipetting variance

Linearity (R²) must be ≥ 0.98. Anything lower means the dilutions or
pipetting weren't consistent.

### Linear dynamic range **[E]**

Report the lowest and highest concentrations where E is within the acceptable
window. Minimum range: 5 orders of magnitude.

### Limit of detection **[D]**

The lowest template input that still amplifies reliably. Report alongside Cq
variance at the LOD.

### Cq variation **[E]**

- Intra-assay (technical replicates within a plate): CV < 5%
- Inter-assay (between plates/runs): CV < 10%

---

## 5. Protocol

### Reaction setup **[E]**

- Reaction volume (10, 20, or 25 µL — state explicitly)
- Master mix: vendor, catalog number, lot
- Primer concentration (200–500 nM)
- Probe concentration (100–250 nM, TaqMan only)
- Template amount: equivalent of 10–100 ng RNA for cDNA reactions

### Cycling **[E]**

A canonical SYBR program:
```
95°C, 10 min initial denaturation
[ 95°C 15 s | 60°C 60 s ] × 40 cycles
Melt curve: 60–95°C, 0.5°C increments
```

Disclose all temperatures, times, ramp rates, and cycle counts.

### Instrument **[E]**

Manufacturer, model, software version, last calibration date (D).

---

## 6. Normalization

### Reference genes **[E]**

Use multiple validated reference genes — never one alone.

- **Minimum:** 2–3 reference genes
- **Recommended:** 3–4 for high-stakes experiments

Common candidates: GAPDH, ACTB, RPL13A, HPRT1, TBP, B2M, PPIA, YWHAZ. Validate
stability in your specific samples using geNorm M-value, NormFinder stability
score, or BestKeeper SD. Report the metric used.

### Normalization method **[E]**

State whether you used:
- ΔΔCq (relative quantification, requires equal efficiencies)
- Pfaffl method (efficiency-corrected relative)
- Absolute quantification with a standard curve in copies

---

## 7. Data Analysis

### Cq determination **[E]**

State the method: automatic threshold (default in most software), manual
threshold, second-derivative maximum (LightCycler default), or sigmoid fit.
Document baseline correction.

### Quality controls **[E]**

- NTC Cq > 35, or undetected. If anywhere lower, treat as contamination.
- −RT control Cq > target Cq + 5, or undetected. Lower means gDNA leaks
  through.
- Positive control should hit its expected Cq within ± 0.5.

### Outlier handling **[E]**

State the method (Grubbs test, ROUT, ± 0.5 Cq from median, etc.) and report
how many points were removed.

### Statistics **[E]**

- Test used: paired/unpaired t-test, ANOVA, Mann-Whitney, etc.
- Multiple-testing correction (Bonferroni, Benjamini-Hochberg) if comparing
  many genes
- Sample size per group
- Error bars: clearly label as SD, SEM, or 95% CI

---

## 8. Reporting Templates

### Methods

```
Total RNA was extracted with [kit], quantified by Qubit, and quality-checked
on a Bioanalyzer (RIN ≥ X). cDNA was synthesized from Y ng RNA using
[RT enzyme] with [priming method] for Z min at T°C. qPCR was performed on a
[instrument] using [master mix] in 20 µL reactions containing 200 nM each
primer (sequences in Supplementary Table SX). Cycling: 95°C 10 min,
followed by 40 cycles of 95°C 15 s / 60°C 60 s, with a final melt curve
60–95°C at 0.5°C increments. Each sample was run in technical triplicate
across at least three biological replicates. Reference genes
([gene list]) were validated using geNorm; the geometric mean Cq was used
for ΔΔCq normalization. NTC and −RT controls were included on every plate.
```

### Supplementary

Include a primer sequence table, the standard curve and computed efficiency
for each assay, raw Cq values, and a completed MIQE checklist. The MIQE
checklist Excel produced by `export_primers(..., format="miqe_checklist")`
already includes the fields the journal expects.

---

## 9. Pre-Order Checklist

Before sending a primer order:

- [ ] Forward Tm 58–62°C, reverse within 2°C of forward
- [ ] Primer lengths 18–22 nt, GC content 40–60%
- [ ] Amplicon 70–140 bp
- [ ] Spans exon–exon junction or intron > 1 kb
- [ ] Heterodimer ΔG > −5 kcal/mol
- [ ] Hairpin ΔG > −2 kcal/mol
- [ ] No 3' self-complementarity ≥ 3 bases
- [ ] In-silico specificity passes (BLAST or local in-silico PCR)

## 10. Post-Ordering Validation

After primers arrive:

- [ ] Gradient PCR finds a clear optimum
- [ ] Standard curve: 5–7 points, R² ≥ 0.98, E in 0.90–1.10
- [ ] Linear dynamic range ≥ 5 logs
- [ ] Single peak on melt curve
- [ ] NTC Cq > 35 or undetected
- [ ] −RT control Cq > target + 5 or undetected
- [ ] Technical replicate CV < 5%

## 11. Publication Submission Checklist

- [ ] ≥ 3 biological replicates per condition
- [ ] ≥ 2 technical replicates per sample
- [ ] ≥ 2 validated reference genes
- [ ] Complete Methods section using the template above
- [ ] Standard curve data in supplementary
- [ ] Raw Cq table in supplementary
- [ ] MIQE checklist completed and attached

---

## Frequently Cited Violations

| Violation | Why reviewers reject |
|-----------|----------------------|
| Single reference gene | Statistically unsound; normalize using a geometric mean of 2+ genes |
| No efficiency reported | Without E, fold-changes are not interpretable |
| Amplicon > 150 bp for qPCR | Longer amplicons amplify with lower efficiency, especially on degraded RNA |
| Missing NTC results | Reviewer can't rule out reagent contamination |
| No melt curve (SYBR) | Reviewer can't rule out non-specific amplification |
| Insufficient replicates | Statistical claims unsupported |
| Missing primer sequences | The work is not reproducible |
| Outlier rule unstated | Looks like cherry-picking |
| No reference-gene validation | Reviewer can't assess normalization soundness |

---

## Resources

- MIQE 2.0 (2025): Bustin SA et al. Clin Chem 71(6):634–660
- Original MIQE (2009): Bustin SA et al. Clin Chem 55(4):611–622
- RDML data format: https://rdml.org/
- geNorm: https://genorm.cmgg.be/
- RefFinder: https://www.heartcure.com.au/reffinder/

---

Last updated: 2026-05-20.
