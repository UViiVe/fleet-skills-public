# DESeq2 Comprehensive Reference

Full code patterns for DESeq2 differential expression. Adapt these to the experimental design you actually have — copying without thinking is how confounded analyses happen.

Related: [decision-guide.md](decision-guide.md) for choosing methods, [troubleshooting.md](troubleshooting.md) for failure modes.

---

## End-to-end skeleton

```r
library(DESeq2)
library(apeglm)

# 1. Construct
dds <- DESeqDataSetFromMatrix(
  countData = counts,           # integer matrix, genes x samples
  colData   = coldata,          # row.names match colnames(counts)
  design    = ~ condition
)

# 2. Pre-filter
dds <- dds[rowSums(counts(dds)) >= 10, ]

# 3. Reference level (be explicit, never trust alphabetical)
dds$condition <- relevel(dds$condition, ref = "control")

# 4. Fit
dds <- DESeq(dds)                 # size factors + dispersions + Wald

# 5. Extract
res       <- results(dds, name = "condition_treated_vs_control")
res_shrnk <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")

# 6. Filter
sig <- subset(res, padj < 0.05 & abs(log2FoldChange) >= 1)
```

This is the entire backbone. Almost every variation is a different design formula in step 1 or a different extraction in step 5.

---

## Design formulas

DESeq2 reads its design as an R formula; the *last* term is the one the Wald test reports by default. Covariates go before condition.

### Two-group

```r
design = ~ condition
```
Single factor, no batch concerns. The most common — and the place to start every analysis before you add complexity.

### Adjusting for a known technical covariate

```r
design = ~ batch + condition
```
Use when PCA shows batch separation. Hard requirement: every condition must appear in every batch. Check it: `table(coldata$batch, coldata$condition)`. If any cell is zero you are confounded — the design matrix will not be full rank and `DESeq()` will error.

### Paired samples

```r
design = ~ subject + condition
```
Tumor/normal from the same patient, pre/post treatment, longitudinal. Treats subject as a blocking factor — controls between-subject variance and gains power.

### Interaction

```r
design = ~ genotype + treatment + genotype:treatment
# equivalent to:
design = ~ genotype * treatment
```
Tests whether the treatment effect differs by genotype. Pull each effect:
```r
resultsNames(dds)
# "Intercept"
# "genotype_KO_vs_WT"
# "treatment_drug_vs_vehicle"
# "genotypeKO.treatmentdrug"     # the interaction
res_int <- results(dds, name = "genotypeKO.treatmentdrug")
```
Interaction tests need more samples — aim for n ≥ 4 per genotype × treatment cell.

### Multi-factor with several covariates

```r
design = ~ sex + age_bin + batch + condition
```
Order doesn't change the model fit but does change the default `results()` extraction (last term). Be sure every coefficient has at least three samples contributing, or `DESeq()` becomes very unstable.

### No-intercept ("group means") parameterization

```r
coldata$group <- factor(paste(coldata$genotype, coldata$treatment, sep = "_"))
design        = ~ 0 + group
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "KO_drug", "WT_vehicle"))
```
Worth using when (a) you want arbitrary pairwise contrasts among many cells, or (b) you need to dissolve a confounded design into one composite factor.

---

## Extracting results

### By coefficient name

```r
resultsNames(dds)
# [1] "Intercept" "condition_treated_vs_control"
res <- results(dds, name = "condition_treated_vs_control")
```
Format: `factor_level_vs_reference`. This is what `apeglm` shrinkage needs.

### By contrast (works for any pair of levels)

```r
res <- results(dds, contrast = c("condition", "treated", "control"))
# c(factor_name, numerator, denominator) — positive LFC = higher in numerator
```
Use this when the comparison you want isn't a stored coefficient (e.g., treatmentA vs treatmentB when both are non-reference levels).

### By numeric contrast vector (advanced)

```r
resultsNames(dds)  # know your column order
res <- results(dds, contrast = c(0, 0.5, -0.5))   # average vs average etc.
```
Reserved for sums or weighted contrasts (e.g., "average of two treatments vs control").

---

## LFC shrinkage

Why shrink: with low counts, the negative-binomial MLE for log2FC has high variance. Genes with `count = 0 vs 1` show up as |log2FC| > 8 even though their effect is unstable. Shrinkage pulls these toward zero proportionally to their standard error, while leaving high-count, high-confidence genes essentially untouched.

| Use case | Use shrunken | Use unshrunken |
|----------|--------------|----------------|
| Ranked GSEA / fgsea | ✅ | ❌ |
| Volcano / MA plots | ✅ | ❌ |
| Top-N "interesting" genes | ✅ | ❌ |
| Hypothesis testing / padj | ❌ | ✅ |
| Reporting effect size + CI | ❌ | ✅ |

### apeglm (default — recommended)

```r
res_shrnk <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")
```
Heavy-tailed prior, preserves large effects best. Requires `coef =`, not `contrast =`.

### ashr (when you need contrast-form shrinkage)

```r
res_shrnk <- lfcShrink(
  dds,
  contrast = c("condition", "treated", "control"),
  type     = "ashr"
)
```
Slightly more conservative on large effects but works with arbitrary contrasts and is fast.

### normal (legacy — avoid)

`type = "normal"` uses the original DESeq2 prior. Inferior to both apeglm and ashr; only use for backward compatibility.

### Building a ranked list for GSEA

```r
ranks <- res_shrnk$log2FoldChange
names(ranks) <- rownames(res_shrnk)
ranks <- na.omit(ranks)
ranks <- sort(ranks, decreasing = TRUE)
writeLines(
  paste(names(ranks), ranks, sep = "\t"),
  con = "ranked_for_gsea.rnk"
)
```

---

## Transformations

For *visualization and clustering* only — never feed transformed values back into DESeq's testing path.

```r
vsd <- vst(dds, blind = FALSE)    # >30 samples, fast
rld <- rlog(dds, blind = FALSE)   # <30 samples, slower, more accurate at low counts
mat <- assay(vsd)                 # plain numeric matrix
```

`blind = FALSE` lets the transformation know about the design (recommended once you have a design you trust). `blind = TRUE` for fully exploratory PCA where you don't want condition labels leaking into the variance estimate.

### Normalized counts (no log transform)

```r
counts(dds, normalized = TRUE)    # size-factor-normalized
sizeFactors(dds)                  # one per sample; ratio between actual and expected library size
```

---

## Pre-filtering

```r
keep <- rowSums(counts(dds)) >= 10        # standard
keep <- rowSums(counts(dds) >= 10) >= 3   # require 10 counts in at least 3 samples
keep <- rowMeans(counts(dds)) >= 5        # mean expression cutoff
dds  <- dds[keep, ]
```

Heavier filtering speeds the fit and shrinks memory; it does *not* change p-values for kept genes much because DESeq2's independent filtering downstream removes low-mean genes from BH adjustment anyway. The main benefit is operational: 60k-row dispersion fits take real time.

---

## Custom alpha and LFC thresholds

```r
res <- results(dds, alpha = 0.05)                          # for summary() / independent filtering
res <- results(dds, lfcThreshold = 1, altHypothesis = "greaterAbs")   # test |LFC| > 1 directly
```

Testing against a non-zero null (`lfcThreshold`) is statistically the right thing to do when you care about biological magnitude — though most people filter post-hoc on `abs(log2FoldChange) > 1` for simplicity.

---

## Likelihood ratio test (LRT)

Use LRT when the question is "does *any* term in the full model matter" rather than a single pairwise comparison — most often for testing a factor with many levels or a time-course.

```r
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)         # all-vs-null
res_lrt <- results(dds_lrt)

# Or test a specific term: is "treatment" needed once we have batch + genotype?
dds_lrt <- DESeq(dds, test = "LRT",
                 full    = ~ batch + genotype + treatment,
                 reduced = ~ batch + genotype)
```

LRT yields per-gene p-values for "is the more complex model better"; log2FoldChange in the results table refers to the last term of the full design (often not what you want — pull individual coefficients with `results(dds, name = ...)`).

---

## Alternative inputs

### From SummarizedExperiment

```r
dds <- DESeqDataSet(se, design = ~ condition)
```

### From tximport (Salmon / Kallisto / RSEM)

```r
library(tximport)
files       <- file.path("salmon_out", samples$sample_id, "quant.sf")
names(files) <- samples$sample_id
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
```
This is the only correct way to bring transcript-level quantifications into DESeq2 — do not round to integers and use `FromMatrix`.

### From featureCounts (Rsubread)

```r
library(Rsubread)
fc  <- featureCounts(files = bam_files, annot.ext = gtf, isGTFAnnotationFile = TRUE)
dds <- DESeqDataSetFromMatrix(fc$counts, coldata, design = ~ condition)
```

---

## Working with the dds object

```r
design(dds) <- ~ batch + condition                       # change design
dds <- DESeq(dds)                                        # refit

dds_sub <- dds[, dds$cell_type == "B"]                   # subset samples
dds_sub <- DESeq(dds_sub)

dds_genes <- dds[rownames(dds) %in% panel_genes, ]       # subset genes (no refit needed)
```

---

## Accessing and summarizing results

```r
summary(res, alpha = 0.05)

resOrdered <- res[order(res$padj), ]                     # by significance
resOrdered <- res[order(abs(res$log2FoldChange), TRUE), ]# by absolute effect

res_df         <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

sum(res$padj < 0.05, na.rm = TRUE)
sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)   # up
sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)   # down
```

---

## Saving and reloading

```r
saveRDS(dds, "dds.rds")
saveRDS(res, "res.rds")

dds <- readRDS("dds.rds")
res <- results(dds, name = "condition_treated_vs_control")  # re-extract is cheap
```

Save the `dds`. A week later when someone asks for a different contrast you can re-extract in seconds instead of rerunning the whole pipeline.

---

## Reporting checklist

Always include in the methods section:
1. DESeq2 version (`packageVersion("DESeq2")`) and R version.
2. Pre-filter rule (e.g., "genes with total counts ≥ 10 across all samples").
3. Full design formula and reference level.
4. Significance threshold (padj cutoff) and any LFC threshold.
5. Shrinkage method (`apeglm` / `ashr`) and what plots/rankings used shrunken vs unshrunken LFC.
6. Independent filtering settings (default unless you changed them).

This pair (design + thresholds) is what reviewers will ask about first.
