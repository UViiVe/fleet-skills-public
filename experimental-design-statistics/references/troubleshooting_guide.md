# Troubleshooting — Common Failures and Fixes

A symptom-driven guide for when the workflow doesn't behave.

---

## Power analysis

### Computed power is < 0.60

Likely cause: effect size too small for the available CV and n.

Fixes, in order:

1. Increase n. Use `calc_sample_size_rnaseq()` to find the n that hits target power.
2. Increase sequencing depth, but only if currently below the assay's depth floor.
3. Revise the target effect size upward. Smaller effects, even if real, may not be biologically meaningful.
4. Reduce biological variability — age-matched samples, controlled diet, standardized collection.
5. Switch to a paired design if biologically appropriate.

### Required n exceeds budget

Fixes:

1. Sweep depth-vs-n combinations within budget; usually more samples wins.
2. Accept detecting only larger effects; document the limitation in methods.
3. Run a pilot first — often reveals lower CV than literature assumed.
4. Switch to a paired design (within-subject correlation lowers required n by 20–30%).

### No pilot data available

Use literature CV from [cv_tissue_database.csv](cv_tissue_database.csv). If unsure, use the upper end of the range (conservative). When confidence is low, plan a small pilot before the main study.

### Pilot CV higher than expected

Don't dismiss it. Investigate:

- Outlier samples? Check QC metrics.
- Hidden batch effects in the pilot? Recompute CV after correction.
- Genuine biology with subgroups? Stratify analysis.

Then size the main study using the pilot estimate — your samples will behave like your pilot, not like a textbook.

---

## Batch design

### `check_confounding()` returns FAIL 🚨

Do not proceed. Regenerate the assignment:

```r
batch <- assign_samples_to_batches(metadata, batch_size = 8,
                                   balance_vars = c("condition", "sex"))
check_confounding(batch, "condition")
```

If confounding persists:

- Verify `condition` is in `balance_vars`.
- Ensure batch_size allows balance (multiple of n_conditions × 2 is ideal).
- Make sure n per condition is reasonably divisible by the number of batches.
- If only one timepoint per batch is feasible (e.g., longitudinal collection), use shared reference samples across batches and acknowledge the limitation.

### Cannot balance all conditions in batches

For three conditions and batch size 8, use:

```r
batch <- incomplete_block_design(metadata, batch_size = 8,
                                 primary_var = "condition",
                                 secondary_vars = c("sex"))
```

Or adjust the batch size: 9, 12, or 15 divides cleanly. A small residual imbalance is acceptable if no condition is missing from any batch.

### Sample failure mid-processing

Document immediately: which sample, which batch, what failed.

- Don't replace with a sample from a different condition.
- Add the replacement to the next batch with paired controls.
- Or accept the imbalance and include batch as a model term.

### Processing dragged on for months

- Record exact processing dates as a covariate.
- Add reference samples in remaining batches.
- Be prepared to model "time" as a batch variable.
- Consider stopping and re-collecting if the spread becomes pathological.

---

## Sample size

### n = 2 per group is all we can do

Strongly recommended: don't run a confirmatory experiment at n = 2. Variance can't be estimated. Treat it as a pilot for hypothesis generation only. Plan a properly-powered follow-up.

If absolutely forced: report results descriptively, don't run statistical tests, do not claim significance.

### Paired design — is it worth the constraint?

A paired design helps if within-subject correlation is > 0.5. Estimate from pilot:

```r
# In DESeq2 with paired design
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ pair + condition)
```

If dispersion is markedly lower than an unpaired version of the same data, pairing helps. Rough rule: CV_paired ≈ CV_unpaired × 0.7.

---

## Multiple testing

### Many raw p < 0.05 but nothing after BH

Check:

- Power: was the study adequately powered for the observed effect sizes?
- Try IHW with mean expression as covariate (often 10–20% more discoveries).
- Inspect QC: hidden batch effects, library size issues, sample swaps.
- For genuinely exploratory work, q < 0.10 can be acceptable — but pre-specify, don't shop.

### Too many significant (> 50%)

Probably an unresolved technical effect (batch not corrected, library size normalization off, sample mix-up). Diagnose with PCA, sample-sample correlation, before/after batch correction.

Alternatively, the comparison is genuinely large (gross tissue type, strong drug). Filter by effect size to focus on biologically meaningful hits (`|log2FC| > 1` AND q < 0.05).

### Don't know which correction to use

When in doubt: BH-FDR. It's the field default.

```r
recommend_method(n_tests = 20000, effect_heterogeneity = "moderate", sample_size = 6)
```

returns a justified recommendation given the test characteristics.

---

## Visualization

### `svglite` not installed warning

Ignored. The script falls back to base R's `svg()` device. Do not install `svglite` — many failures originate from trying to fix this warning. PNG is always saved.

### Plot doesn't save

Check the working directory and write permissions. Use absolute paths if uncertain:

```r
ggplot2::ggsave("/abs/path/out.png", plot = p, width = 8, height = 6, dpi = 300)
```

Verify the plot object exists before saving with `print(p)`.

---

## Installation

### BiocManager fails to install

```r
install.packages("BiocManager", repos = "https://cran.r-project.org")
```

### Specific Bioconductor package fails

```r
BiocManager::install(version = "3.18")   # update to current release first
BiocManager::install("RNASeqPower")
```

### `powsimR` install is very slow

Expected — large package with many dependencies. Plan for 10–30 minutes. Run once and cache the environment.

### Compilation errors

macOS:
```bash
xcode-select --install
```

Linux (Debian/Ubuntu):
```bash
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

---

## When statistical correction won't save the experiment

Some failure modes cannot be corrected post-hoc:

- Complete confounding between batch and condition.
- Batch effects > 50% of variance.
- Non-linear batch effects (e.g., dose-dependent kit performance).
- Batch-specific biology (different protocols across batches).

Honest options:

- Drop the problematic batch(es) and accept reduced n.
- Repeat the experiment with a corrected design.
- Acknowledge the limitation in the manuscript and downgrade claims accordingly.

Don't try to rescue a confounded experiment with increasingly aggressive computational correction. The result won't be trustworthy and reviewers will see through it.

---

## Cross-references

- Design principles — [experimental_design_best_practices.md](experimental_design_best_practices.md).
- Software install — [software_requirements.md](software_requirements.md).
- QC — [qc_guidelines.md](qc_guidelines.md).
