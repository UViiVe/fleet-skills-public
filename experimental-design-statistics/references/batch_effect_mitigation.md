# Batch Effect Mitigation

How to design experiments so batch effects don't masquerade as biology. The cheapest batch correction is the one you didn't need.

---

## What batch effects are

Batch effects are systematic non-biological variation introduced when samples are processed in distinct groups — separate days, separate plates, separate kits, separate operators, separate sequencing runs. They are technical artifacts. They can be larger than the biological effect you're trying to find.

Sources, in roughly the order they bite:

- **Sample processing:** day-to-day, technician-to-technician, reagent lot.
- **Library prep:** kit lot, barcoding plate, prep day.
- **Sequencing:** flow cell, lane, instrument.
- **Computational:** alignment software version, reference build, pipeline parameters.

Magnitude: batch effects routinely affect 20–50% of features in RNA-seq data. Variance attributable to batch can rival or exceed variance from biological condition.

---

## The cardinal rule

**Every batch must contain samples from every experimental condition.** If condition correlates with batch, you cannot separate the two — no statistical method, no surrogate variable, nothing.

Confounded (bad):

```
Batch 1: control_1, control_2, control_3
Batch 2: treated_1, treated_2, treated_3
```

Balanced (good):

```
Batch 1: control_1, control_2, treated_1, treated_2
Batch 2: control_3, control_4, treated_3, treated_4
```

The script `check_confounding()` in `batch_validation.R` enforces this with a chi-square test. If it returns FAIL, regenerate the design. Do not proceed.

---

## Design strategies

### Strategy 1 — Complete randomized assignment

For a single experimental factor with no other covariates to balance, randomly assign samples to batches so each batch has an equal number of each condition.

```r
batch <- assign_samples_to_batches(metadata, batch_size = 8, balance_vars = "condition")
```

### Strategy 2 — Stratified randomization

When important covariates need to be balanced (sex, age, site), stratify the random assignment within combinations of those variables.

```r
batch <- assign_samples_to_batches(metadata, batch_size = 8,
                                   balance_vars = c("condition", "sex", "age_group"))
```

OSAT optimization tries thousands of random assignments and picks the one with the best balance score.

### Strategy 3 — Blocked / paired design

When samples come in natural pairs (before/after on the same subject; matched twins; case-control pairs), keep pairs in the same batch when possible. This removes between-batch noise from the within-pair comparison.

```
Batch 1: subj1_pre, subj1_post, subj2_pre, subj2_post
Batch 2: subj3_pre, subj3_post, subj4_pre, subj4_post
```

Statistical model: include the pairing factor as a blocking variable in the design formula.

### Strategy 4 — Reference samples across batches

For multi-batch studies that span weeks or months, include one or two reference aliquots (a pool, a cell line) in every batch. Reference samples enable empirical estimation of batch shifts. Standard in GTEx, TCGA, and most large cohort studies.

---

## Batch size

Trade-off: smaller batches mean more batches (more between-batch variance) but easier to fill. Larger batches mean fewer batches and tighter balance but harder logistics.

| Assay | Typical batch size | Driven by |
|---|---|---|
| Bulk RNA-seq | 8–24 | Library prep plate |
| scRNA-seq | 4–8 | 10x Chromium lane |
| ATAC-seq | 8–12 | Tagmentation reaction |
| ChIP-seq | 4–8 | Antibody / ChIP success rate |
| Methylation array | 96 | 96-well array plate |
| Proteomics TMT | 6, 10, 16, 18 | TMT plex chemistry |

Rule: choose batch size as a multiple of `n_conditions × 2` so each condition has at least 2 samples per batch. If you can't, use an incomplete block design (`incomplete_block_design()`).

---

## Covariate balancing

### Priority tiers

1. **Must balance:** condition. Non-negotiable. Failure here = unrecoverable confounding.
2. **High priority:** sex (huge effects in many tissues), disease subtype, collection site.
3. **Medium priority:** age band, genetic background, time of collection.
4. **Nice to balance:** technician, storage time, instrument.

When you can't balance everything, hold the higher tier and include lower-tier variables as model terms in the analysis (`~ batch + age + condition`).

### Diagnostic tests

For each covariate, run a chi-square test of independence with batch:

| chi-square p | Interpretation |
|---|---|
| > 0.05 | Well balanced |
| 0.01–0.05 | Moderate imbalance — include as model term |
| < 0.01 | Poorly balanced — regenerate if covariate is high-priority |

---

## Temporal batching

When samples accrue over time (rolling clinical enrollment, sequential cohorts), three patterns:

- **Fixed-interval batches.** Wait until enough samples accumulate, process a full balanced batch, repeat. Default for clinical studies.
- **Continuous processing with reference.** Process samples as they arrive in small batches, with a shared reference in every batch. Suits diagnostics or surveillance.
- **Drop-and-recollect.** If a sample arrives too late and would unbalance the design, exclude it and document.

If processing takes longer than expected (months instead of weeks), record exact processing dates as a covariate. "Date as batch" is a real model term.

---

## Validation, end-to-end

### Before processing

```r
source("scripts/batch_validation.R")
check_confounding(batch, "condition")          # MUST PASS
check_balance(batch, c("sex", "age_group"))    # SHOULD show p > 0.05
visualize_batch_design(batch, "condition", output_file = "out/batch_design")
```

If `check_confounding()` returns FAIL, do not proceed. Regenerate the assignment.

### During processing

Record for each batch (the lab protocol checklist enforces this):

- Date and start/end time.
- Operator(s).
- Reagent lots used (every lot, especially library prep kit).
- Instrument/equipment IDs.
- Plate barcodes and well positions.
- Any protocol deviations.

### After processing, before DE analysis

- PCA colored by batch: ideally samples cluster by biology, not by batch.
- Per-sample QC metrics by batch: are library sizes, alignment rates, duplication rates similar across batches?
- Variance partition: how much variance is explained by batch vs. condition?

If batch explains > 10% of variance, statistical correction is warranted. If > 30%, examine whether the experiment is salvageable.

---

## Statistical correction (when prevention isn't enough)

If batch effects remain after good design, four common approaches:

| Method | Strength | Weakness |
|---|---|---|
| Linear model with batch covariate (`~ batch + condition` in DESeq2) | Simple, transparent, well understood | Limited if batch effects are non-linear |
| ComBat (`sva` package) | Strong batch removal | Can over-correct, remove real biology |
| RUVSeq | Data-driven, uses control genes | Needs negative controls or replicates |
| SVA (surrogate variable analysis) | Discovers hidden factors | Can be aggressive; hard to interpret |

When statistical correction fails:

- Complete confounding (batch = condition) — unrecoverable.
- Very strong batch effects (> 50% variance) — even good methods struggle.
- Non-linear or batch-specific biology — model misspecification.

In these cases the honest options are dropping problematic batches or repeating the experiment with a better design.

---

## Common mistakes

- **Sequential processing by condition.** Processing all controls week 1 and all treatments week 2 creates complete confounding. Always interleave.
- **Ignoring sex.** Sex effects can hit thousands of genes. Balance it; if you can't, include it in the model.
- **Filling batches as samples arrive.** Early arrivals end up in batch 1, late arrivals in batch 2. If arrival timing correlates with condition (e.g., disease severity), you've created a hidden confounder.
- **Variable batch sizes.** A batch of 12 and a batch of 4 are harder to combine analytically. Keep them within 20% of each other.
- **No metadata.** Process samples without recording batch IDs and you cannot correct for batch later. Capture metadata as you process.

---

## Special cases

### One sample arrives late

Options, in rough order of preference:

1. Drop it. If you're already well-powered (n ≥ 6), one sample is unlikely to swing results.
2. Process in a follow-up batch with paired controls from the freezer; include batch in the model.
3. Save for a future replication study.

Don't replace it with a sample from another condition — that creates confounding.

### Sample count doesn't divide cleanly

22 samples with batch size 8: use 8 + 8 + 6, balance conditions within each. Variable batch sizes are fine as long as no batch is dramatically smaller.

### Library-prep batch and sequencing-run batch differ

Record both in metadata. Library-prep batch is usually the bigger effect. Use `~ seq_run + prep_batch + condition`, or drop the smaller-effect term if it's tiny.

---

## Checklist

Before processing:
- [ ] All conditions present in every batch.
- [ ] Covariates (priority 1–2) balanced or included in model plan.
- [ ] `check_confounding()` returns PASS.
- [ ] `check_balance()` reviewed for each covariate.
- [ ] Batch layout visualized, manually inspected.
- [ ] Batch CSV exported, lab team has the file.

During processing:
- [ ] Date, operator, reagent lots recorded for each batch.
- [ ] Randomized order within batch.
- [ ] Reference samples included if multi-week run.

After processing:
- [ ] All samples have batch column in metadata.
- [ ] PCA colored by batch reviewed.
- [ ] Variance-by-batch quantified.
- [ ] Batch handling documented in methods section.

---

## Citations

- Leek JT, Scharpf RB, Bravo HC, et al. *Tackling the widespread and critical impact of batch effects in high-throughput data.* Nat Rev Genet 2010;11:733–739.
- Goh WWB, Wang W, Wong L. *Why batch effects matter in omics data, and how to avoid them.* Trends Biotechnol 2017;35:498–507.
- Johnson WE, Li C, Rabinovic A. *Adjusting batch effects in microarray expression data using empirical Bayes methods.* Biostatistics 2007;8:118–127.
- Yasrebi H. *OSAT: a tool for sample-to-batch allocations in genomics experiments.* BMC Genomics 2012;13:689.

---

## Cross-references

- Design principles — [experimental_design_best_practices.md](experimental_design_best_practices.md).
- QC checklists — [qc_guidelines.md](qc_guidelines.md).
- Scripts — `batch_assignment.R`, `batch_validation.R`.
