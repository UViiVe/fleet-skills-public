# Experimental Design — Best Practices

A working reference for the design questions that actually determine whether a genomics experiment yields interpretable results. Where the SKILL.md tells you *what* to run, this document explains *why*.

---

## 1. Replication

### Biological vs. technical replicates

Biological replicates are independent biological units (different individuals, different mice, different organoid wells started from different passages). Technical replicates are the same biological material processed more than once. They answer different questions:

- Biological replication estimates the variance you'll need to overcome to claim a population-level effect. Without it, you cannot generalize.
- Technical replication estimates pipeline noise. With modern protocols this is small and usually not worth the budget.

**Default replication tiers (per condition):**

| Tier | n | When to use |
|---|---|---|
| Absolute minimum | 3 | Pilots, irreplaceable samples, exploratory |
| Standard | 5–6 | Most published bulk RNA-seq, ATAC-seq |
| High confidence | 10+ | Detecting < 1.5-fold effects, clinical/translational |

n = 2 cannot estimate variance. Don't run statistical tests on it; treat it as a case report.

### Pseudoreplication — the most common silent error

Pseudoreplication is treating non-independent measurements as independent. Examples that look fine in the methods section but are wrong:

- Sequencing the same library on two lanes and calling it n = 2.
- Splitting one mouse's liver into three tubes and calling it n = 3.
- Treating 5000 cells from one donor as n = 5000 in a scRNA-seq DE test (correct n is 1 — the donor).

For single-cell experiments, the unit of replication is the sample (donor / mouse / well), not the cell. Use pseudobulk aggregation or a mixed-effects model that accounts for sample as a random effect.

---

## 2. Sample size vs. sequencing depth

When the budget is fixed, the trade-off is: more samples at shallower depth, or fewer samples at deeper depth?

**Default answer: more samples.** Power gains from sample size dominate gains from depth, except below assay-specific depth floors.

| Assay | Depth floor (below this, add depth first) | Standard | Deeper helps when |
|---|---|---|---|
| Bulk RNA-seq | 10 M reads | 20–30 M | Isoform-level, lowly-expressed transcripts |
| scRNA-seq | 20 K reads/cell | 50 K | Rare transcripts, allele-specific expression |
| ATAC-seq | 15 M reads | 30–50 M | Footprinting, TF motif analysis |
| ChIP-seq (narrow) | 10 M | 20 M | Diffuse / weak binding |
| ChIP-seq (broad) | 20 M | 40 M | Broad histone marks (H3K27me3, H3K9me3) |

The empirical rule from Schurch et al. (2016) for RNA-seq: doubling depth gains ~10–20% power, doubling n gains ~40–50%. Spend the marginal dollar on samples.

---

## 3. Randomization

What to randomize:

- Assignment of biological units to conditions, when applicable.
- Sample collection order, processing order within a batch, plate well positions.
- Sequencing lane and run assignment.
- Operator assignment if multiple technicians are involved.

How to randomize: use a computer-generated random assignment (set a seed for reproducibility), stratified by major covariates. Never assign "by eye" or in arrival order.

A useful diagnostic: if a reviewer asked "which sample is in well A01?", could you give them a defensible reason that isn't "it was first on the list"?

---

## 4. Covariate balance

Balance the variables that have large known effects, in priority order:

1. **Must balance:** condition (the primary variable), sex, disease subtype, treatment site.
2. **Should balance:** age band, genetic background, collection date or season, time of day.
3. **Nice to balance:** technician, storage time, sample volume.

When perfect balance is impossible (e.g., 22 samples and batch size of 8), prioritize the must-balance layer and include unbalanced lower-priority covariates as model terms in the analysis. Don't sacrifice condition balance to balance sex.

---

## 5. Pilot studies pay for themselves

A pilot of 2–4 samples per condition can:

- Provide a real dispersion estimate (often 1.5–4× different from a literature CV).
- Validate the entire collection-through-sequencing pipeline.
- Surface protocol issues while there's still time to fix them.
- Detect contamination, batch artifacts, or sample misidentification before scaling.

Run pilots at the **same depth** as the planned main study. A shallow pilot gives misleading dispersion estimates.

Using pilot data:
- Build a `DESeqDataSet` from the pilot counts, run `DESeq()` to estimate dispersion.
- Feed the dispersion into `samplesize_from_pilot()` from `scripts/sample_size_de.R`.
- Compare the recommended n to your budget; if it exceeds, escalate to the budget-optimization workflow below.

---

## 6. Common pitfalls

### Pitfall: confounding batch with condition

Processing all controls in week 1 and all treatments in week 2 makes the experiment uninterpretable. Even strong statistical adjustment cannot rescue this. Run `check_confounding()` before any sample is touched, and refuse to proceed if it returns FAIL.

### Pitfall: pooling

Pooling several individuals into one tube to save sequencing cost destroys the variance estimate. n = 1, no inference. If the budget forces pooling, reduce depth instead.

### Pitfall: under-powering then explaining it away

If power < 0.6 at the maximum feasible n, the honest paths are: redesign for larger effects, escalate to a paired design, run a different experiment, or label this explicitly as exploratory / hypothesis-generating. Running anyway and reporting "trends" is not a path; reviewers see through it.

### Pitfall: post-hoc power analysis

Computing power after the experiment, using the observed effect size, gives no additional information. It tracks the p-value monotonically. Don't do it; don't request it.

### Pitfall: ignoring multiple-testing burden in the power calculation

A power calculation that uses α = 0.05 for an RNA-seq experiment is wrong by orders of magnitude. After BH-FDR at 0.05, the effective per-test threshold is ~10⁻³ for typical genome-wide settings. Either use a tool that builds the correction in (`ssizeRNA`) or use a corrected α.

---

## 7. Sample size quick reference

Per-condition n for detecting moderate effects with standard power (0.80) and CV from literature:

| Assay | Minimum | Standard | High confidence |
|---|---|---|---|
| Bulk RNA-seq | 3 | 5–6 | 10+ |
| scRNA-seq (samples) | 3 | 5–6 | 8–10 |
| ATAC-seq | 4 | 6–8 | 10+ |
| ChIP-seq (narrow) | 2 (ENCODE) | 4–5 | 6+ |
| ChIP-seq (broad) | 3 | 4–5 | 6+ |
| Methylation array | 5 | 8–10 | 15+ |
| Proteomics (TMT) | 3 | 5–6 | 8+ |

Multiply by these factors for smaller effects:

| Target fold-change | Multiplier |
|---|---|
| ≥ 2.0 | × 1 |
| 1.5 – 2.0 | × 1.5–2 |
| 1.2 – 1.5 | × 3–4 |
| < 1.2 | × 20–50 (often impractical) |

These are starting points, not pre-registration values. Always confirm with a proper power calculation using your CV or pilot dispersion.

---

## 8. Budget optimization

**Strategy 1 — depth ↔ replicates sweep.** Compute power for a grid of (depth, n) combinations within the budget. Take the cheapest combination meeting power ≥ 0.80. Almost always favors more samples at moderate depth.

**Strategy 2 — staged design.** Run a small pilot, refit the dispersion estimate, then run the main study with the corrected n. This costs ~15–25% of total budget for the pilot but routinely saves 30–50% on the main study by removing the over-conservative literature CV padding.

**Strategy 3 — adaptive sample size.** Pre-specify stopping rules (e.g., interim analysis after n = 4 per group; stop if a clear signal emerges, else continue to n = 8). Requires proper group-sequential boundaries — see O'Brien-Fleming or Lan-DeMets — and a pre-registered protocol.

**Strategy 4 — drop secondary comparisons.** If a multi-condition design has primary and exploratory contrasts, allocate the power budget to the primary contrast. Document which comparisons are confirmatory vs. hypothesis-generating.

---

## 9. Documentation requirements

What goes in the design artifact you hand to collaborators / reviewers / future-you:

**Mandatory:**
- n per condition, with the source of the number (pilot dispersion, literature CV, both).
- Sequencing depth target per sample.
- Batch structure and the rule used to assign samples to batches.
- Randomization seed.
- Covariate balance plan.
- Power calculation: effect size, variability, α, target power, method, software version.
- Multiple-testing method and threshold.
- Statistical analysis plan covering the primary contrast.

**Strongly recommended:**
- Pilot data location and key statistics.
- QC pass/fail criteria for sample inclusion.
- Failure replacement policy.
- Processing timeline and batch dates.
- Software environment (R + Bioconductor versions).

**For grants and publications:**
- Sample size justification with the calculation shown.
- A statement that batch will not be confounded with condition.
- The pre-specified multiple-testing method (no swapping after seeing data).

---

## 10. When to consult a statistician

Consider escalating to a biostatistician for: factorial designs with interaction terms of interest, repeated measures, hierarchical/multi-site studies, adaptive designs, very rare-disease cohorts, regulatory submissions, clinical trials. Earlier is cheaper than later — bring them in at design, not at analysis.

---

## Cross-references

- Power calculation mechanics — [power_analysis_guidelines.md](power_analysis_guidelines.md).
- Batch-effect prevention and validation — [batch_effect_mitigation.md](batch_effect_mitigation.md).
- Multiple-testing method selection — [multiple_testing_guide.md](multiple_testing_guide.md).
- Pre-run and post-run checklists — [qc_guidelines.md](qc_guidelines.md).
- CV by tissue — [cv_tissue_database.csv](cv_tissue_database.csv).
