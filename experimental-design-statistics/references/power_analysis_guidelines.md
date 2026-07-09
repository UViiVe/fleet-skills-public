# Power Analysis — Practical Guide

How to compute statistical power for genomics experiments without lying to yourself about effect sizes.

---

## What power actually is

Power is the probability of rejecting the null when the alternative is true — the probability of finding a real effect when it exists. It depends on four quantities:

1. **Effect size (δ):** the difference you want to detect, expressed as fold-change or standardized effect.
2. **Sample size (n):** biological replicates per group.
3. **α:** the significance threshold, typically 0.05.
4. **Variability (σ):** standard deviation of the measurements; for count data, dispersion.

Holding the other three fixed, power increases with larger effect sizes, larger n, smaller α (no — smaller α *decreases* power), and lower variability.

A power calculation is solved one of two ways: fix n and compute power, or fix the target power and solve for n. The skill's scripts support both directions.

---

## Prospective vs. retrospective

**Prospective (a priori):** done before data collection, to size the experiment. This is the only kind worth doing.

**Retrospective (post-hoc):** done after the experiment using observed effect sizes. Mathematically equivalent to the p-value and provides no new information; commonly misused to "explain" a null result.

If a reviewer asks you for a post-hoc power calculation, push back and offer a sensitivity analysis instead (what effect sizes would have been detectable at your achieved n?).

---

## Workflow

### Step 1 — Pick the effect size

Three sources, in order of preference:

1. **Pilot data.** Extract observed fold-changes for genes of interest, use them as targets, recognize that pilot estimates with small n have wide intervals and lean toward the smaller realistic value.
2. **Published studies.** Search comparable experiments. Be conservative; published effect sizes are subject to winner's-curse inflation.
3. **Minimum biologically meaningful effect.** What's the smallest difference that would mean anything? For RNA-seq, 1.5–2× is typical; smaller and you start asking whether the effect, even if real, matters.

Typical realistic effects:

| Assay | Strong treatment | Subtle effect |
|---|---|---|
| Bulk RNA-seq | 2–4× | 1.2–1.5× |
| scRNA-seq cell-type markers | 3–10× | 1.5–2× |
| ATAC-seq | 1.5–3× | 1.2–1.5× |
| ChIP-seq | 2–10× enrichment | varies |

### Step 2 — Estimate variability

**From pilot data (preferred):** use DESeq2 dispersion. The dispersion parameter accounts for the mean-variance relationship of count data, which a single CV does not.

```r
library(DESeq2)
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
median_disp <- median(dispersions(dds), na.rm = TRUE)
```

**From literature CV:** see [cv_tissue_database.csv](cv_tissue_database.csv) for tissue-specific values. Reasonable ranges:

| Source | CV |
|---|---|
| Cell lines | 0.1–0.2 |
| Inbred mice | 0.2–0.3 |
| Primary cells, outbred animals | 0.3–0.4 |
| Human cohorts | 0.4–0.5 |
| Heterogeneous disease populations | 0.6–1.0+ |

**Conservative fallback:** if uncertain, use the upper end of the range. Over-estimating variability gives a larger n and a safer study.

### Step 3 — Set α and target power

Standard:
- α = 0.05.
- Power = 0.80 for most studies.
- Power = 0.90 for grant proposals, expensive validations, clinical translation.

Lower power (0.60–0.70) can be acceptable for explicitly exploratory pilots — but say so.

### Step 4 — Run the calculation

For bulk RNA-seq, the simplest path:

```r
source("scripts/power_rnaseq.R")
calc_power_rnaseq(depth = 20, n_per_group = 6, cv = 0.4, fold_change = 2, alpha = 0.05)
```

To solve for n at a target power:

```r
calc_sample_size_rnaseq(depth = 20, cv = 0.4, fold_change = 1.5, target_power = 0.8)
```

For pilot-driven power:

```r
source("scripts/power_pilot_based.R")
calc_power_from_pilot(pilot_dds, n_per_group = 6, fold_change = 1.5)
```

For scRNA-seq, use `scripts/sample_size_scrna.R` and remember: cells per sample doesn't substitute for samples.

### Step 5 — Interpret the result

| Power | Verdict |
|---|---|
| ≥ 0.90 | Excellent — proceed |
| 0.80 – 0.89 | Adequate — proceed, document |
| 0.60 – 0.79 | Borderline — increase n or relax effect size |
| < 0.60 | Underpowered — redesign or label exploratory |

---

## Assay-specific notes

### Bulk RNA-seq

- Tools: `RNASeqPower` (fast, Hart et al. 2013), `ssizeRNA` (accounts for multiple testing), `PROPER` (simulation), `DESeq2` + pilot (most accurate).
- Power varies by expression level: lowly-expressed genes have lower power.
- ~15–20 K genes tested → multiple-testing burden material.

### Single-cell RNA-seq

Two levels of variation interact: sample-to-sample (biological) and cell-to-cell (technical + biological). Most subtle bugs in scRNA-seq power calculations come from conflating the two.

Rules of thumb:

| Goal | Samples/group | Cells/sample |
|---|---|---|
| Cell-type discovery | 3–4 | 3000–10000 |
| Differential expression | 5–6 | 1000–2000 |
| Trajectory inference | 4–6 | 3000–5000 |
| Rare cell type (< 5%) | 8–10 | 10000–50000 |

**Cells per sample ≠ replicates.** Five samples with 1000 cells each is statistically more powerful than one sample with 5000 cells for population-level inference.

Tools: `powsimR` for simulation-based power; `scPower` for DE-specific calculations.

### ATAC-seq

Higher technical variability than RNA-seq (CV ~ 0.3–0.5). Plan ~50% more replicates than the equivalent RNA-seq for a given effect.

Standard parameters:
- 50–150 K peaks tested.
- Dispersion 0.2–0.5.
- 4–6 replicates per group for moderate effects.

### ChIP-seq

Power is dominated by antibody quality, not n. ENCODE accepts n = 2 for peak calling but requires concordance metrics; for differential binding, plan n ≥ 3–4.

Narrow peaks (TFs): 10–50 K features. Broad peaks (histone marks): 50–200 K features.

---

## Multiple-testing burden in the power calculation

Testing 20 K genes shifts the effective threshold from 0.05 to roughly 10⁻³ for BH-FDR at 0.05. Three ways to handle this:

1. **Use a tool that builds it in:** `ssizeRNA` accepts FDR directly.
2. **Use a corrected α in the power calculation:** approximate effective α as 0.005–0.01 for genome-wide BH-FDR at 0.05.
3. **Frame power as a discovery rate:** "We have 80% power to detect 50% of genuinely 2-fold-changed genes."

The third framing is honest and informative, but requires assumptions about the true positive rate.

---

## Common scenarios

### Two-group case-control

```r
calc_power_rnaseq(depth = 20, n_per_group = 5, cv = 0.4, fold_change = 2)
```

Most common scenario; the default workflow.

### Paired design (before/after on same subject)

Within-subject correlation reduces variance, typically by 30–50%. Use `paired_samplesize()` from `sample_size_de.R`, which adjusts dispersion downward. Plan ~30% fewer pairs than independent samples.

### Multi-group (k conditions)

Pairwise comparisons share the α budget. A rough rule: multiply the two-group n by √k. For three conditions, plan n × 1.7 per group if all pairwise comparisons matter.

If only one specific contrast matters, design for that contrast directly.

### Multi-factor (interaction)

Interaction effects are typically smaller than main effects. Plan n = 6–8 per cell to detect interactions reliably. Many factorial designs end up underpowered for interactions even when well-powered for main effects.

### Pilot → main study staged design

```r
source("scripts/power_pilot_based.R")
pilot_dds <- readRDS("pilot_deseq2.rds")
calc_samplesize_from_pilot(pilot_dds, fold_change = 1.5, target_power = 0.8)
```

Pilot-driven estimates routinely differ from literature CV by 30–50%. The pilot pays for itself if it changes n by even 1 sample at $1000 per sample.

---

## Reporting

### Methods section template

> Sample size was determined by power analysis using the RNASeqPower R package (v1.36.0) implementing the method of Hart et al. (2013). Assuming a coefficient of variation of 0.40 estimated from pilot data (n = 3 per group), 20 million reads per sample, α = 0.05, and a target fold-change of 1.5, n = 8 per group was required to achieve 80% power. Multiple testing was controlled with the Benjamini-Hochberg procedure at FDR < 0.05.

### Grant section template

> The proposed sample size of n = 8 per group provides ≥ 80% statistical power to detect 1.5-fold expression changes at FDR < 0.05, based on tissue-matched dispersion estimates from preliminary data (Aim 1 pilot, supplementary Figure S2). Power was computed using RNASeqPower (Hart et al. 2013) and validated with simulation in ssizeRNA. Sensitivity analysis is shown in supplementary Table S3.

### What not to write

- "We performed a post-hoc power analysis and found power was 0.45." (Don't.)
- "A trend toward significance was observed (p = 0.08)." (Either it's significant or it isn't.)
- "We chose n = 3 as is conventional in the field." (Convention is not a power calculation.)

---

## Decision rules summary

| Question | Default | Override when |
|---|---|---|
| Pilot or literature CV? | Pilot if available | Pilot has n < 3 per group — blend |
| Add samples or add depth? | Samples | Below depth floor for assay |
| α | 0.05 | GWAS (5×10⁻⁸); confirmatory (0.01) |
| Power target | 0.80 | Grants (0.90); pilots (0.60–0.70) |
| Multiple-testing in power calc | Build in via `ssizeRNA` | Quick first pass: use corrected α |

---

## Citations

- Hart SN, Therneau TM, Zhang Y, Poland GA, Kocher J-P. *Calculating sample size estimates for RNA sequencing data.* J Comput Biol 2013;20:970–978.
- Schurch NJ, Schofield P, Gierliński M, et al. *How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?* RNA 2016;22:839–851.
- Ching T, Huang S, Garmire LX. *Power analysis and sample size estimation for RNA-Seq differential expression.* RNA 2014;20:1684–1696.
- Vieth B, Ziegenhain C, Parekh S, Enard W, Hellmann I. *powsimR: a flexible tool for power analysis and sample size estimation of bulk and single-cell RNA-seq experiments.* Bioinformatics 2017;33:3486–3488.

---

## Cross-references

- Design principles — [experimental_design_best_practices.md](experimental_design_best_practices.md).
- Multiple-testing — [multiple_testing_guide.md](multiple_testing_guide.md).
- Tissue CV — [cv_tissue_database.csv](cv_tissue_database.csv).
- QC checklists — [qc_guidelines.md](qc_guidelines.md).
