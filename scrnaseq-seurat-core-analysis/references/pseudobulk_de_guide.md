# Pseudobulk Differential Expression (Seurat → DESeq2)

The core rule for condition contrasts in multi-sample scRNA-seq:

> **Cell-level statistics are exploratory only. For inference, aggregate to pseudobulk and use DESeq2 (or edgeR / limma-voom).**

This is true whether you're working in Python (scanpy) or R (Seurat). The R toolchain is `muscat` + `DESeq2` (or just `DESeq2` directly on aggregated counts).

## The pseudoreplication problem

In multi-sample experiments, cells from the same sample are not independent observations. Two CD4 T cells from the same donor share donor genotype, cell isolation, library prep, sequencing run. Treating them as 60 000 independent observations when you have 6 samples is statistically illegitimate.

Empirically (Squair et al. 2021): cell-level Wilcoxon for condition contrasts can produce false discovery rates of 50–90% instead of the nominal 5%. Findings don't replicate.

## The pseudobulk fix

Sum the raw counts of all cells in each `sample × cell_type` combination. This produces a `n_samples × n_genes` count matrix per cell type — exactly what bulk RNA-seq tools expect. You test condition on biological replicates (samples), not on cells.

`aggregate_to_pseudobulk()` does this with sane defaults: minimum 10 cells per combination, raw counts from `RNA$counts`, sum aggregation.

## When pseudobulk is the right answer

| Use case                                               | Method                       |
| ------------------------------------------------------ | ---------------------------- |
| Treated vs control across replicates                   | Pseudobulk + DESeq2          |
| Disease vs healthy across patients                     | Pseudobulk + DESeq2          |
| Timepoint A vs B (paired)                              | Pseudobulk + DESeq2 (paired) |
| "What are the markers of cluster 5?"                   | `FindAllMarkers` (Wilcoxon)  |
| "Which genes define this cluster vs all others?"       | `FindAllMarkers` (Wilcoxon)  |

`FindMarkers(ident.1 = "treated", ident.2 = "control")` is statistically dishonest for condition contrasts. Use pseudobulk.

## Workflow

### 1. Aggregate

`aggregate_to_pseudobulk()` returns a named list keyed on cell type, each entry holding a counts data frame and a sample metadata data frame:

```r
pb <- aggregate_to_pseudobulk(
  seurat_obj,
  sample_col = "sample_id",
  celltype_col = "cell_type",
  min_cells = 10
)
```

Or use `muscat`:

```r
library(muscat)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(seurat_obj)
pb <- aggregateData(
  sce,
  assay = "counts",
  fun = "sum",                       # critical: sum, not mean
  by = c("sample_id", "cluster_id")
)
```

Critical: use **sum**, not mean or median; use **raw counts**, not normalized.

### 2. Run DESeq2

`run_pseudobulk_deseq2()` calls DESeq2 per cell type and returns a list of result data frames:

```r
de_results <- run_pseudobulk_deseq2(
  pb,
  formula = "~ batch + condition",
  contrast = c("condition", "treated", "control")
)
```

Raw `DESeq2` usage:

```r
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = counts_matrix,
  colData   = sample_meta,
  design    = ~ batch + condition
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "treated", "control"))
res_shrunk <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "ashr")
```

`muscat` style:

```r
res_DS <- pbDS(
  pb,
  design   = ~ condition,
  contrast = c("condition", "treated", "control"),
  method   = "DESeq2"
)
res_table <- resDS(sce, res_DS)
```

### 3. Filter and interpret

```r
library(dplyr)
sig <- res_shrunk |>
  as.data.frame() |>
  filter(padj < 0.05, abs(log2FoldChange) > 0.5, baseMean > 10) |>
  arrange(padj)
```

Filtering:

- `padj < 0.05`: FDR-controlled significance.
- `abs(log2FoldChange) > 0.5`: ~ 1.4× fold change, biologically interpretable.
- `baseMean > 10`: drops lowly expressed genes that dominate the noise.

## Design formulas

| Goal                          | Formula                       |
| ----------------------------- | ----------------------------- |
| Two conditions, no covariates | `~ condition`                 |
| Account for batch             | `~ batch + condition`         |
| Paired (donor-matched)        | `~ donor + condition`         |
| Multiple covariates           | `~ batch + sex + condition`   |

The variable of interest goes last in the formula — DESeq2's default contrast is on the last term.

Always check that batch is not perfectly confounded with condition. If batch 1 = control and batch 2 = treated, the model can't separate them.

## Common pitfalls

- **Cell-level Wilcoxon for condition contrasts.** Pseudoreplication. Never use for inference.
- **Mean aggregation instead of sum.** Breaks DESeq2's dispersion model.
- **Normalized counts into DESeq2.** DESeq2 expects raw counts and does its own normalization.
- **Ignoring batch in design.** Inflates dispersion estimates and depresses power.
- **Reporting cell-level p-values.** Even when the effect is real, the p-value is meaningless under pseudoreplication.

## Minimum sample size

| Comparison           | Hard minimum         | Recommended            | Note                                                |
| -------------------- | -------------------- | ---------------------- | --------------------------------------------------- |
| Two conditions       | 2 per group          | 5–6 per group          | N=1 in any group ⇒ DESeq2 can't estimate dispersion |
| Paired               | 3 pairs              | 5–6 pairs              |                                                     |
| Multi-group (k > 2)  | 3 per group          | 4–5 per group          |                                                     |

If N=1 is unavoidable (rare disease, pilot data):

1. **Do not run pseudobulk DE.** Report this as a study limitation.
2. **Use `FindMarkers` as candidate-gene exploration only.** Label it explicitly.
3. **Report descriptive statistics on composition.** "Condition X has 12% more T cells" is fine; "p = 0.003" is not.
4. **Recommend follow-up samples.** This is the right answer.

Power expectations:

- N = 2: low power, only huge effects.
- N = 3: large effects (> 2×) with moderate power.
- N = 5–6: good power for 1.5–2×.
- N ≥ 8: robust for smaller effects.

## Reporting

For each cell type:

1. Number of samples per condition.
2. Number of cells aggregated per sample (median, range).
3. Design formula.
4. Filtering criteria.
5. Multiple-testing correction (default: BH).

Example output:

| Gene  | baseMean | log2FC | lfcSE | stat | pvalue   | padj    | n_samples | n_cells_mean |
| ----- | -------- | ------ | ----- | ---- | -------- | ------- | --------- | ------------ |
| GENE1 | 245.3    | 2.14   | 0.31  | 6.9  | 1.2e-11  | 3.4e-08 | 6         | 421          |
| GENE2 | 189.7    | -1.87  | 0.28  | -6.7 | 2.1e-11  | 3.4e-08 | 6         | 421          |

## Cell-level vs pseudobulk: quick reference

| Goal                              | Method                  | Valid for inference?  |
| --------------------------------- | ----------------------- | --------------------- |
| Find cluster markers              | Wilcoxon (cell-level)   | No (exploratory)      |
| Compare conditions                | Pseudobulk + DESeq2     | Yes                   |
| Identify cell types (within data) | Wilcoxon (cell-level)   | Yes (within-dataset)  |
| Test treatment effect             | Pseudobulk + DESeq2     | Yes                   |

## Citations

1. Squair JW, Gautier M, Kathe C, et al. *Confronting false discoveries in single-cell differential expression.* Nat Commun 12:5692 (2021).
2. Crowell HL, Soneson C, Germain P-L, et al. *muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data.* Nat Commun 11:6077 (2020).
3. Zimmerman KD, Espeland MA, Langefeld CD. *A practical solution to pseudoreplication bias in single-cell studies.* Nat Commun 12:738 (2021).
4. Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology 15:550 (2014).
5. Murphy AE, Skene NG. *A balanced measure shows superior performance of pseudobulk methods.* Nat Commun 13:7851 (2022).
