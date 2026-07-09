# DESeq2 Troubleshooting Guide

DESeq2 errors fall into four categories: data shape problems, design problems, expected-but-unfamiliar behavior, and performance. Each section below names the symptom, explains what's happening, and shows the fix.

---

## Data shape problems

### `error in checkFullRank(modelMatrix): the model matrix is not full rank`

Your design has perfectly correlated terms. Most common case:

```r
coldata <- data.frame(
  condition = c("control","control","treated","treated"),
  batch     = c("A","A","B","B")
)
design = ~ batch + condition
```
Every "A" sample is a control and every "B" is treated — batch and condition are the same variable. The fix depends on what you can do:

1. **Drop the confounded covariate.** `design = ~ condition`. You lose batch correction but at least you get a valid fit.
2. **Combine into one composite factor and use contrast extraction:**
   ```r
   coldata$grp <- factor(paste(coldata$batch, coldata$condition, sep="_"))
   design <- ~ 0 + grp
   ```
3. **Get more samples.** This is the only fix that actually solves the problem if batch genuinely confounds condition.

Verify any new design before refitting:
```r
mm <- model.matrix(design, coldata)
Matrix::rankMatrix(mm) == ncol(mm)   # TRUE = full rank
```

---

### `error: some values in assay are not integers`

Either you passed normalized data (TPM/FPKM/RPKM) or you imported transcript-level quantifications without the right wrapper.

```r
# Wrong
counts <- read.csv("tpm.csv")
dds    <- DESeqDataSetFromMatrix(counts, ...)

# Right
library(tximport)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
```

`DESeqDataSetFromTximport` uses the transcript abundance estimates as offsets to recover gene-level counts that DESeq2 can model — much better than rounding the abundances yourself.

If your counts have minor fractional artifacts from preprocessing (rare), `round()` is a last-ditch option. Don't round TPMs and pretend.

---

### `error: every gene contains at least one zero`

DESeq2 needs to compute log-geometric-means and your matrix has zero rows everywhere. Almost always: the matrix is transposed.

```r
dim(counts)         # should be (n_genes, n_samples), usually genes >> samples
head(counts)        # row names should look like genes; columns like samples
counts <- t(counts) # if dim was (samples, genes)
```

Less commonly: every gene actually does have a zero in some sample, which happens for ultra-shallow sequencing. Filter more aggressively (`rowSums >= 50`) or revisit whether you have enough depth to test.

---

### `factor levels not in colData` / `subscript out of bounds`

Sample IDs don't line up between `counts` and `coldata`. The pasilla example dataset trips this — its metadata uses a `fb` suffix not present in the count columns.

```r
colnames(counts)        # e.g. "untreated1" "untreated2" ...
rownames(coldata)       # e.g. "untreated1fb" "untreated2fb" ...

all(colnames(counts) == rownames(coldata))     # FALSE

# Fix: align names then reorder
rownames(coldata) <- gsub("fb$", "", rownames(coldata))
coldata           <- coldata[colnames(counts), , drop = FALSE]
```

The `validate_input_data()` helper in `scripts/load_example_data.R` does these checks for you.

---

## Design problems

### Zero significant genes

The full menu of causes, roughly in order of frequency:

1. **Wrong reference level.** `levels(dds$condition)` puts the alphabetically first level as reference. If that wasn't your intent, `relevel()` and refit.
2. **Batch dominates.** Run `plotPCA(vst(dds, blind=FALSE), intgroup = c("condition","batch"))`. If samples cluster by batch, add it: `design(dds) <- ~ batch + condition; dds <- DESeq(dds)`.
3. **n too small.** Two-vs-two has very little power unless effects are dramatic. Dispersion estimates are unreliable below n = 3 per group.
4. **There is no signal.** Sometimes the treatment didn't do anything detectable at the bulk level. Check that the positive controls (housekeepers stable, known responders changing) behave as expected.
5. **Too-aggressive pre-filter** discarded the genes of interest.

Diagnose:
```r
plotDispEsts(dds)        # genes should cluster around the red fitted line
summary(res, alpha=0.1)  # relax the alpha to see if there is any signal at all
```

If you have a candidate gene that "should" change:
```r
plotCounts(dds, gene = "ENSG00000123456", intgroup = "condition")
```
A flat plot tells you something different from a noisy one.

---

### Unexpected fold-change direction

`log2FoldChange > 0` means the condition is up *relative to the reference level*. If the result looks backwards, the reference is set wrong.

```r
levels(dds$condition)                 # first element is the reference
dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)                     # refit
```

Or flip the contrast direction on extraction:
```r
res <- results(dds, contrast = c("condition", "control", "treated"))
```
No need to refit — `contrast =` just inverts the comparison.

---

### Suspiciously large LFCs (|log2FC| > 8)

Almost always low-count genes. A gene with 0 vs 2 reads gives an enormous unshrunken LFC that means nothing. Apply shrinkage:

```r
res_shrnk <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")
```

If you have *many* such genes, pre-filter harder before fitting:
```r
keep <- rowSums(counts(dds) >= 5) >= 3
dds  <- dds[keep, ]
```

Inspect any specific gene that still looks suspicious:
```r
plotCounts(dds, gene = "weird_gene", intgroup = "condition")
```

---

## Expected behavior that looks wrong

### Many `NA` values in `padj`

DESeq2 sets padj to NA in three cases:
1. **Independent filtering** — low-mean genes are excluded from BH adjustment to gain power. This is intentional and increases the number of significant calls overall.
2. **Cook's distance outliers** — genes with extreme outliers in at least one sample are flagged.
3. **All-zero rows** that survived pre-filtering somehow.

```r
mcols(res)$filterThreshold      # the mean-cutoff DESeq2 chose
table(is.na(res$padj), is.na(res$pvalue))   # pvalue NA but padj NA → independent filtering
```

If you want padj on every gene that has a pvalue:
```r
res <- results(dds, independentFiltering = FALSE)
```
This will give you padj for more genes but with somewhat less power overall.

---

### Cook's distance flagging a lot of genes

```r
boxplot(log10(assays(dds)[["cooks"]]), las = 2)
```
If one sample is an outlier across many genes, it likely has a library prep or labeling issue. Decide whether to drop it:
```r
dds <- dds[, colnames(dds) != "bad_sample"]
dds <- DESeq(dds)
```
Or, if you're confident the data is fine but want the test anyway:
```r
res <- results(dds, cooksCutoff = FALSE)
```

---

### High dispersion estimates

`plotDispEsts(dds)` shows points scattered far from the fitted red trend.

Possible reasons:
- High biological variability (e.g., primary tumors). Real, accept lower power.
- Unmodeled batch effect — add to design.
- Outlier samples — inspect with PCA.

You can't "fix" high dispersion in the data; you can only fit a better model or accept the power cost.

---

## Performance problems

### `DESeq()` is slow

Symptoms: minutes for n=20 samples is normal; tens of minutes is a smell.

- Pre-filter more aggressively: `dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]`.
- Use parallel processing:
  ```r
  library(BiocParallel)
  register(MulticoreParam(4))
  dds <- DESeq(dds, parallel = TRUE)
  ```
- Drop unnecessary covariates from the design — every extra term costs degrees of freedom and fit iterations.

### `rlog()` is unbearably slow

It will be on > 100 samples. Switch to `vst()` — they produce similar results for downstream PCA and clustering.

---

## Quick error → fix lookup

| Error message | Likely cause | One-line fix |
|---------------|--------------|--------------|
| `model matrix not full rank` | Confounded design | Drop covariate or merge factors |
| `should be integers` | TPM/FPKM data | Use raw counts or `FromTximport` |
| `every gene contains at least one zero` | Transposed counts | `counts <- t(counts)` |
| `factor levels not in colData` | Name/ordering mismatch | `coldata <- coldata[colnames(counts), ]` |
| `subscript out of bounds` | Same as above | Same |
| Many `NA` in padj | Independent filtering | Expected; set `independentFiltering = FALSE` to disable |
| Zero significant genes | Reference / batch / no signal | Check PCA + reference level first |
| Huge \|LFC\| > 10 | Low counts | `lfcShrink()` + heavier pre-filter |
| Slow `DESeq()` | Big matrix | Pre-filter + `parallel = TRUE` |
| Slow `rlog()` | Large n | Switch to `vst()` |

---

## When all else fails

1. `sessionInfo()` and `packageVersion("DESeq2")` — version drift is real.
2. Reduce to a minimal reproducible example: 100 random genes, 4-6 samples. If it works, the issue is data; if it still fails, the issue is the design or the install.
3. The Bioconductor support site at `support.bioconductor.org` has nearly every common failure already answered, usually by Mike Love.
