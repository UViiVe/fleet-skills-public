# DESeq2 Usage Guide

Quick-orientation document for first-time users. For the full reference see [comprehensive-reference.md](comprehensive-reference.md); for decision logic see [decision-guide.md](decision-guide.md).

## What DESeq2 does

DESeq2 models bulk RNA-seq counts with a negative binomial GLM. The output you care about per gene:

- `baseMean` — average normalized expression across all samples
- `log2FoldChange` — effect direction and magnitude (numerator vs denominator)
- `lfcSE` — standard error of the LFC
- `stat` — Wald statistic
- `pvalue` — raw p-value
- `padj` — Benjamini-Hochberg-adjusted p-value (use this for significance)

## Minimum pipeline

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds$condition <- relevel(dds$condition, ref = "control")
dds <- DESeq(dds)
res <- results(dds, name = "condition_treated_vs_control")
res_shrnk <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")
```

That's the whole thing. Every additional script in this skill is a wrapper around this skeleton with validation, plots, or export logic added.

## Inputs you'll be asked for

| Input | Format | Example |
|-------|--------|---------|
| Count matrix | Integer matrix, genes × samples | CSV with row.names = gene IDs |
| Metadata | data.frame, row.names = sample IDs | One row per sample, columns for condition + covariates |
| Design formula | R formula referring to metadata columns | `~ batch + condition` |
| Reference level | Character string in the condition column | `"control"` |

## Example prompts that work

- "Run DESeq2 on these counts using treated vs control as the comparison, controlling for batch."
- "I have Salmon outputs in `salmon_out/`. Build a DESeqDataSet from them and run the standard pipeline."
- "Re-extract results for treatmentA vs treatmentB from this dds object."
- "Apply LFC shrinkage and export a ranked list ready for GSEA."

## What the agent will do

1. Validate the count matrix is integer and that sample IDs line up.
2. Build the `DESeqDataSet` and pre-filter low-count genes.
3. Set the reference level explicitly.
4. Run `DESeq()`, extract Wald results, apply `lfcShrink()`.
5. Produce QC plots (dispersion, PCA, MA, volcano).
6. Write a results table, a significant-only table, normalized counts, and the saved `dds` object.

## Pitfalls to avoid

- Do not pass TPM/FPKM/RPKM. DESeq2 needs raw integer counts.
- Do not "remove" the batch effect from counts before fitting — include batch in the design instead.
- Do not use raw p-values for filtering — use `padj`.
- Do not use unshrunken LFC for plotting or ranking — shrinkage matters most for low-count genes that would otherwise dominate.
- Do not assume alphabetical ordering picks the right reference. Use `relevel()`.

## When to use which script

| Goal | Script |
|------|--------|
| Get example data and validate user data | `scripts/load_example_data.R` |
| Run the full pipeline end-to-end | `scripts/basic_workflow.R` |
| Pull a specific contrast or apply shrinkage manually | `scripts/extract_results.R` |
| Generate QC plots | `scripts/qc_plots.R` |
| Compute / export VST or rlog matrices | `scripts/transformations.R` |
| Save everything to disk | `scripts/export_results.R` |
| Try a multi-level factor pipeline | `scripts/multi_condition.R` |
| Try a batch-corrected pipeline on simulated data | `scripts/batch_correction.R` |
