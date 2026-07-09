# Pseudobulk Differential Expression

The core rule for any condition contrast in multi-sample scRNA-seq:

> **Cell-level statistics are exploratory only. For inference, aggregate to pseudobulk and use a bulk RNA-seq tool (DESeq2, edgeR, limma-voom).**

This document explains why and how.

## The pseudoreplication problem

In a multi-sample experiment, cells from the same sample are not independent observations. Two CD4 T cells captured from the same donor share an enormous amount of state: donor genotype, cell isolation, library prep, sequencing run. Treating them as 60 000 independent observations when you have 6 samples is statistically illegitimate.

Empirically (Squair et al. 2021): when cell-level Wilcoxon is used for condition contrasts in multi-sample data, the false discovery rate can rise to 50–90% instead of the nominal 5%. Findings don't replicate. Reviewers know this now and increasingly reject papers that do it.

## The pseudobulk fix

Sum the raw counts of all cells in each `sample × cell_type` combination. This produces a `n_samples × n_genes` count matrix per cell type — exactly what a bulk RNA-seq tool expects. You then test the condition contrast on biological replicates (samples), not on cells.

```
Individual cells (invalid)         Pseudobulk (valid)
S1_cell1                          S1_CellTypeA = sum(all CellTypeA in S1)
S1_cell2                  →       S2_CellTypeA = sum(all CellTypeA in S2)
...                               S3_CellTypeA = sum(all CellTypeA in S3)
```

This:

- Respects the biological replicate structure.
- Produces p-values you can actually report.
- Works with the dispersion-shrinkage machinery in DESeq2/edgeR.
- Lets you include batch / donor / sex as covariates in the model.

## When pseudobulk is the right answer

| Use case                                               | Method                       |
| ------------------------------------------------------ | ---------------------------- |
| Treated vs control across replicates                   | Pseudobulk + DESeq2          |
| Disease vs healthy across patients                     | Pseudobulk + DESeq2          |
| Timepoint A vs timepoint B (paired)                    | Pseudobulk + DESeq2 (paired) |
| "What are the markers of cluster 5?"                   | Cell-level Wilcoxon          |
| "Which genes define this cluster against all others?"  | Cell-level Wilcoxon          |

Cell-level Wilcoxon is fine for cluster characterization within a single dataset; it's the wrong tool for any contrast across samples.

## Workflow

### 1. Aggregate

Sum raw counts per `sample × cell_type`, filtering small groups:

```python
def aggregate_to_pseudobulk(adata, sample_key, celltype_key, min_cells=10):
    pseudobulk_obs = []
    counts_rows = []
    for sample in adata.obs[sample_key].unique():
        for ct in adata.obs[celltype_key].unique():
            mask = (adata.obs[sample_key] == sample) & (adata.obs[celltype_key] == ct)
            n_cells = mask.sum()
            if n_cells < min_cells:
                continue
            X = adata[mask, :].layers["counts"]
            counts_rows.append(np.asarray(X.sum(axis=0)).ravel())
            pseudobulk_obs.append({"sample": sample, "cell_type": ct, "n_cells": n_cells})
    return pd.DataFrame(counts_rows), pd.DataFrame(pseudobulk_obs)
```

R equivalent (muscat):

```r
library(muscat)
pb <- aggregateData(
  sce,
  assay = "counts",
  fun = "sum",                     # critical: sum, not mean
  by = c("sample_id", "cluster_id")
)
```

Critical: use **sum**, not mean or median, and use **raw counts**, not normalized.

### 2. Build sample metadata

```python
sample_meta = pd.DataFrame({
    "sample_id": ["S1", "S2", "S3", "S4", "S5", "S6"],
    "condition": ["control", "control", "control", "treated", "treated", "treated"],
    "batch":     ["b1",      "b1",      "b2",      "b1",      "b2",      "b2"],
    "donor":     ["D1",      "D2",      "D3",      "D4",      "D5",      "D6"],
})
```

### 3. Run DESeq2

Python (pydeseq2):

```python
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

dds = DeseqDataSet(counts=counts_df, metadata=sample_meta, design_factors=["batch", "condition"])
dds.deseq2()
stats = DeseqStats(dds, contrast=["condition", "treated", "control"])
stats.summary()
results = stats.results_df
```

R native:

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

### 4. Filter and interpret

```r
sig <- subset(res_shrunk,
              padj < 0.05 &
              abs(log2FoldChange) > 0.5 &
              baseMean > 10)
sig <- sig[order(sig$padj), ]
```

Filtering criteria:

- `padj < 0.05`: FDR-controlled significance.
- `|log2FC| > 0.5`: ~1.4× fold change, biologically interpretable.
- `baseMean > 10`: drops lowly expressed genes that dominate the noise.

## Design formulas

| Goal                          | Formula                       |
| ----------------------------- | ----------------------------- |
| Two conditions, no covariates | `~ condition`                 |
| Account for batch             | `~ batch + condition`         |
| Paired design (donor-matched) | `~ donor + condition`         |
| Multiple covariates           | `~ batch + sex + condition`   |

The variable of interest goes last in the formula — that's the convention DESeq2 uses for the default contrast.

Always check that your batch variable is not perfectly confounded with condition. If batch 1 = control and batch 2 = treated, the model can't separate them and you'll see misleading effects.

## muscat (R, integrated workflow)

`muscat` provides aggregation, design, and DE in one package:

```r
library(muscat)
library(SingleCellExperiment)

sce <- as.SingleCellExperiment(seurat_obj)
pb  <- aggregateData(sce, assay = "counts", fun = "sum",
                     by = c("sample_id", "cluster_id"))
res_DS <- pbDS(pb,
               design   = ~ condition,
               contrast = c("condition", "treated", "control"),
               method   = "DESeq2")
res_table <- resDS(sce, res_DS)
```

## Common pitfalls

- **Cell-level Wilcoxon for condition contrasts.** Pseudoreplication. Never do this for inference.
- **Mean aggregation instead of sum.** Loses count scaling — DESeq2's dispersion model breaks.
- **Normalized counts into DESeq2.** DESeq2 assumes raw counts and does its own normalization. Feeding it logCPMs gives garbage.
- **Ignoring batch in design.** If batch is a real source of variation, leaving it out inflates dispersion and depresses power.
- **Reporting cell-level p-values.** Even if the effect is real, the p-value is meaningless. Use pseudobulk for the test, then visualize at cell level.

## Minimum sample size

| Comparison           | Hard minimum         | Recommended            | Note                                            |
| -------------------- | -------------------- | ---------------------- | ----------------------------------------------- |
| Two conditions       | 2 per group          | 5–6 per group          | N=1 in any group ⇒ DESeq2 can't fit dispersion  |
| Paired               | 3 pairs              | 5–6 pairs              |                                                 |
| Multi-group (k > 2)  | 3 per group          | 4–5 per group          |                                                 |

If N=1 in a condition is unavoidable (rare disease, pilot data):

1. **Do not run pseudobulk DE.** Report this as a study limitation.
2. **Use cell-level Wilcoxon as candidate-gene exploration only.** Make this label explicit in the writeup.
3. **Report descriptive statistics on composition.** "Condition X has 12% more T cells" is fine; "p=0.003" is not when N=1.
4. **Recommend follow-up samples.** This is the right answer; report writers should be honest about it.

Power notes:

- N=2: can fit DESeq2 but very low power; only detects huge effects.
- N=3: detects large effects (> 2×) with moderate power.
- N=5–6: good power for 1.5–2× effects.
- N ≥ 8: robust power for smaller effect sizes.

## Reporting

For each cell type:

1. Number of samples per condition.
2. Number of cells aggregated per sample (median, range).
3. Design formula.
4. Filtering criteria (count threshold, presence threshold).
5. Multiple-testing correction method (default: BH).

Example results table:

| Gene  | baseMean | log2FC | lfcSE | stat | pvalue   | padj    | n_samples | n_cells_mean |
| ----- | -------- | ------ | ----- | ---- | -------- | ------- | --------- | ------------ |
| GENE1 | 245.3    | 2.14   | 0.31  | 6.9  | 1.2e-11  | 3.4e-08 | 6         | 421          |
| GENE2 | 189.7    | -1.87  | 0.28  | -6.7 | 2.1e-11  | 3.4e-08 | 6         | 421          |

## Cell-level vs pseudobulk: quick reference

| Goal                              | Method                      | Valid for inference?      |
| --------------------------------- | --------------------------- | ------------------------- |
| Find cluster markers              | Wilcoxon (cell-level)       | No (exploratory only)     |
| Compare conditions                | Pseudobulk + DESeq2         | Yes                       |
| Identify cell types (within data) | Wilcoxon (cell-level)       | Yes (within-dataset)      |
| Test treatment effect             | Pseudobulk + DESeq2         | Yes                       |

## Citations

1. Squair JW, Gautier M, Kathe C, et al. *Confronting false discoveries in single-cell differential expression.* Nat Commun 12:5692 (2021).
2. Crowell HL, Soneson C, Germain P-L, et al. *muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data.* Nat Commun 11:6077 (2020).
3. Zimmerman KD, Espeland MA, Langefeld CD. *A practical solution to pseudoreplication bias in single-cell studies.* Nat Commun 12:738 (2021).
4. Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology 15:550 (2014).
5. Murphy AE, Skene NG. *A balanced measure shows superior performance of pseudobulk methods in single-cell RNA-sequencing analysis.* Nat Commun 13:7851 (2022).
