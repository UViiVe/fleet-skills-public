---
name: bulk-rnaseq-counts-to-de-deseq2
id: "skill_134e37784ea64cd2be9da0822c2eb767"
display-name: "Bulk RNA-seq DE with DESeq2"
short-description: "End-to-end DESeq2 pipeline that turns a raw count matrix into a ranked, shrunken, QC-checked DE result table. Use when you have integer counts and biological replicates and need to identify differe..."
category: transcriptomics
keywords: [rnaseq, deseq2, differential expression, counts, bulk]
last-updated: "2026-07-07"
author: "Fleet Skills contributors"
created: 2026-05-20
description: "End-to-end DESeq2 pipeline that turns a raw count matrix into a ranked, shrunken, QC-checked DE result table. Use when you have integer counts and biological replicates and need to identify differentially expressed genes between two or more conditions."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
version: 1.0.1
provenance:
  wraps: "DESeq2"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Love MI et al., Genome Biology 2014, doi:10.1186/s13059-014-0550-8"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Bulk RNA-seq Differential Expression with DESeq2

DESeq2 is the workhorse for bulk RNA-seq differential expression: a negative-binomial GLM with empirical Bayes dispersion shrinkage. This skill walks a raw count matrix through every step you actually need — validation, normalization, dispersion estimation, hypothesis testing, LFC shrinkage, QC, and export — without reinventing the pieces inline.

## When to Use This Skill

✅ **Use this skill when you have:**
- A **raw integer count matrix** (genes × samples) from featureCounts, HTSeq, Salmon/tximport, or a quantification SummarizedExperiment
- **Biological replicates** (n ≥ 2 per condition; n ≥ 3-4 is where the dispersion estimator starts to behave)
- A **defined comparison**: condition A vs B, treated vs control, paired pre/post, factorial design, etc.
- Need for **shrunken LFCs** so that low-count genes don't dominate volcano plots or gene rankings

❌ **Don't use this skill for:**
- TPM/FPKM/RPKM or otherwise normalized data — DESeq2 needs integer counts. Use limma-voom on a count matrix instead.
- Single-cell RNA-seq — pseudobulk first, or use scran/glmGamPoi/MAST. DESeq2 on per-cell counts will overdisperse.
- Tiny experiments (n = 1 per group) — no replicate, no dispersion estimate, no test.
- ChIP-seq differential binding (use DiffBind), ATAC-seq (csaw), proteomics (limma/DEqMS), or microarrays.

If your data is normalized but you really must use DESeq2, the right move is to go back to upstream alignment/quantification and pull integer counts — not to round TPMs.

## Quick Start

Smoke-test the pipeline in roughly two minutes with the pasilla dataset (Brooks et al. 2011, Drosophila pasilla RNAi knockdown):

```r
source("scripts/load_example_data.R")
ds       <- load_pasilla_data()         # installs pasilla on first run (~50 MB)
counts   <- ds$counts                   # ~14,600 genes × 7 samples
coldata  <- ds$coldata                  # condition: treated / untreated

source("scripts/basic_workflow.R")      # creates dds, res, res_shrunk
source("scripts/qc_plots.R");   run_all_qc(dds, res, output_dir = "qc")
source("scripts/export_results.R"); export_all(dds, res, res_shrunk, output_dir = "results")
```

You should land around 800-1,200 significant genes at padj < 0.1 (the DESeq2 default), with strong PCA separation between the three treated and four untreated replicates.

## Installation

| Package | Version | License | Commercial use | Install |
|---------|---------|---------|----------------|---------|
| DESeq2 | ≥ 1.40 | LGPL-3 | ✅ | `BiocManager::install("DESeq2")` |
| apeglm | ≥ 1.20 | GPL-2 | ✅ | `BiocManager::install("apeglm")` |
| ashr | ≥ 2.2 | GPL-3 | ✅ | `install.packages("ashr")` |
| ggplot2 | ≥ 3.4 | MIT | ✅ | `install.packages("ggplot2")` |
| ggrepel | ≥ 0.9 | GPL-3 | ✅ | `install.packages("ggrepel")` |
| ggprism | ≥ 1.0 | GPL-3 | ✅ | `install.packages("ggprism")` |
| svglite | ≥ 2.1 | GPL-2 | ✅ (optional) | `install.packages("svglite")` |
| pasilla / airway | latest | Artistic-2.0 | ✅ (example data only) | `BiocManager::install(c("pasilla","airway"))` |

Bootstrap install:

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm", "ashr"))
install.packages(c("ggplot2", "ggrepel", "ggprism", "svglite"))
```

## Inputs

**Required:**

1. **Count matrix** — integer counts, genes in rows, samples in columns. Acceptable as `matrix`, `data.frame`, CSV/TSV (with row names = gene IDs), RDS, or `SummarizedExperiment`. Gene IDs can be Ensembl, gene symbols, RefSeq — the pipeline does not care, but stay consistent for downstream annotation.
2. **Sample metadata (`colData`)** — a `data.frame` whose row names match `colnames(counts)` exactly (the loader will fail loudly if they don't). Must contain at least one factor column for the comparison of interest.

**Optional:**

- **Batch / covariate columns** in `colData` (e.g., `batch`, `sex`, `cell_line`, `sequencing_run`). Include in the design formula if PCA shows separation by them.
- **tx2gene** lookup table if importing from Salmon/Kallisto via `tximport` — the pipeline supports `DESeqDataSetFromTximport`.
- **`lfcThreshold`** if you want to test against a fold-change null other than zero (rarely needed; see `scripts/extract_results.R`).

## Outputs

Everything lands under `output_dir/` (default: `results/`). Files:

| File | Contents | Downstream consumer |
|------|----------|---------------------|
| `deseq2_results.csv` | Full DE table: baseMean, log2FoldChange, lfcSE, stat, pvalue, padj | Filtering, reporting |
| `deseq2_results_shrunk.csv` | Same columns but with `apeglm`-shrunken LFC | Volcano/MA plots, ranked GSEA |
| `deseq2_significant.csv` | Rows passing padj + |LFC| cutoffs | Pathway enrichment input |
| `normalized_counts.csv` | Size-factor-normalized counts | Heatmaps, custom plots |
| `vst_transformed.csv` or `rlog_transformed.csv` | Variance-stabilized matrix | PCA, clustering, WGCNA |
| `top100_genes.csv` | Top 100 by padj | Quick scan, supplementary |
| `dds_object.rds` | Full `DESeqDataSet` | Re-extract other contrasts later |
| `qc/dispersion_plot.{png,svg}` | Dispersion-vs-mean fit | Sanity-check normalization |
| `qc/pca_plot.{png,svg}` | PC1/PC2 with sample labels | Detect batch effects, outliers |
| `qc/ma_plot.{png,svg}` | log2FC vs mean, top genes labeled | Visualize bias |
| `qc/volcano_plot.{png,svg}` | -log10(padj) vs log2FC, top genes labeled | Headline figure |

PNG is always produced; SVG is best-effort (`svglite` if available, else base `svg()` device, else skipped with a warning).

## Clarification Questions

**⚠️ CRITICAL: ASK THIS FIRST.**

1. **Is the count data raw integer counts, and where is it?** Either point to a file (CSV/TSV/RDS/SummarizedExperiment) or confirm we should use the bundled pasilla/airway example. If the user provides TPM/FPKM/RPKM, **stop** and explain why DESeq2 won't accept it — you need integer counts from upstream quantification.

2. **What is the primary comparison?** Concretely: which column of `colData` defines the contrast, what is the *reference* level, and what is the *treatment* level? `log2FoldChange > 0` means higher in the non-reference level.

3. **Are there known technical covariates?** Batch, library prep date, sequencing lane, sex, cell line — anything that might co-vary with condition. If the user doesn't know, run PCA first and decide from the plot.

4. **What sample size are we working with per group?** n ≥ 4 is comfortable. n = 2-3 will run but is fragile; flag this and recommend cautious interpretation. n = 1 is a non-starter.

5. **Significance thresholds?** Defaults: `padj < 0.05`, `|log2FC| ≥ 1` (2-fold). If the user has a hypothesis-generating screen, relax to `padj < 0.1`; for a validation list relax LFC but tighten padj.

6. **Multiple contrasts wanted?** A three-level factor exposes A-vs-ref, B-vs-ref directly; A-vs-B needs a `contrast=` extraction. Make this explicit before running.

## Standard Workflow

The scripts under `scripts/` cover the full pipeline. Source them rather than rewriting their logic inline — they encapsulate the validation, reference-level handling, and SVG fallback that you would otherwise have to redo.

```r
# 1. Load and validate
source("scripts/load_example_data.R")
ds      <- load_pasilla_data()              # or load_airway_data()
counts  <- ds$counts
coldata <- ds$coldata

# 2. Run the pipeline (creates dds, res, res_shrunk in the global env)
source("scripts/basic_workflow.R")

# 3. QC plots (dispersion, PCA, MA, volcano)
source("scripts/qc_plots.R")
run_all_qc(dds, res, output_dir = "qc")

# 4. Export everything
source("scripts/export_results.R")
export_all(dds, res, res_shrunk, output_dir = "results")
```

Step-by-step responsibilities:

1. **Validate** sample IDs match between count columns and metadata rows. The loader reorders `coldata` to align — do not skip this, mismatches silently produce nonsensical DE tables.
2. **Construct `DESeqDataSet`** with `design = ~ condition` (or `~ batch + condition`, `~ individual + condition`, `~ genotype * treatment` — see `references/decision-guide.md`).
3. **Pre-filter** genes with `rowSums(counts) < 10`. Cuts memory and speeds dispersion estimation; does not change results materially.
4. **Set reference level** with `relevel()`. Without this, R defaults to alphabetical order — which is rarely what you mean.
5. **`DESeq(dds)`** runs size-factor normalization (median-of-ratios), gene-wise dispersion estimation, dispersion shrinkage toward a fitted trend, and the Wald test.
6. **Extract results** with `results(dds, name = ...)` or `results(dds, contrast = c(factor, num, denom))`. Use `resultsNames(dds)` to discover the coefficient names.
7. **Shrink LFC** with `lfcShrink(..., type = "apeglm")` for plotting and ranking. Unshrunk LFCs are for hypothesis testing and reporting effect sizes with CIs.
8. **QC plots** — if PCA doesn't separate by condition, the DE list is suspect regardless of how many "significant" genes you see.

## When Scripts Fail

Try in this order:

1. **Fix and retry** — most failures are missing packages or a CRAN mirror not set. Install and re-source. The example loaders set the mirror automatically.
2. **Modify the script** — for instance, swap `apeglm` for `ashr` if your contrast needs `contrast=` rather than `coef=`. Edit the file and re-source.
3. **Use as reference** — read the function, write a slimmer adaptation in your context. Cite the source in a comment.
4. **Write from scratch** — only if your environment fundamentally can't run the script (e.g., no internet for installs and packages aren't preinstalled).

Concrete examples:
- `package not available`: BiocManager + retry.
- `model matrix not full rank`: design is confounded — your batch and condition perfectly correlate. Drop the confounded variable or collapse to a single composite factor.
- `every gene contains at least one zero`: count matrix is transposed — `counts <- t(counts)`.
- `lfcShrink ... type='apeglm' requires coef`: pass `coef = resultsNames(dds)[i]` not `contrast=`.

## Decision Points

**Design formula.** Always run a blind VST + PCA (`vst(dds, blind = TRUE); plotPCA(vsd, intgroup = c("condition","batch"))`) before committing. Add a covariate to the design only if it explains visible separation in PCA *and* is not confounded with condition. Check confounding with `table(coldata$batch, coldata$condition)` — every cell should be > 0.

**Transformation.** `vst()` for n > 30 (fast), `rlog()` for n < 30 (better at low counts, slow at scale). The QC script auto-selects.

**LFC shrinkage method.** `apeglm` is the default — best ranking, preserves large effects. `ashr` if you need a non-pairwise `contrast=`. `normal` is legacy; do not use for new work.

**Significance cutoffs.** `padj < 0.05` and `|log2FC| ≥ 1` for a balanced default. Tighten LFC for in-vivo data where effect sizes tend to be modest but real; loosen for an exploratory screen. Always report the cutoff used.

## Common Issues

| Symptom | Root cause | Fix |
|---------|------------|-----|
| `model matrix is not full rank` | Batch and condition perfectly correlated | Drop the confounded variable or merge into a single factor `~ 0 + group` |
| `counts matrix should contain integer values` | You passed TPM/FPKM, or unrounded tximport output | Use raw counts, or `DESeqDataSetFromTximport()` instead of `FromMatrix` |
| `every gene contains at least one zero` | Matrix transposed (samples in rows) | `counts <- t(counts)`; verify `dim(counts)` = (genes, samples) |
| `subscript out of bounds` / `factor levels not in colData` | `colnames(counts)` and `rownames(coldata)` don't match | Reorder: `coldata <- coldata[colnames(counts), ]`; check for stray suffixes (the pasilla loader strips `fb`) |
| All `padj == NA` | Independent filtering removed low-count genes (this is correct) | Expected; if too aggressive, `results(dds, independentFiltering = FALSE)` |
| Zero significant genes | True null, wrong reference level, confounded batch, or n too small | Check PCA, verify `levels(dds$condition)`, look at dispersion plot |
| log2FC values > 10 for some genes | Low counts dividing near zero | Use `lfcShrink()`; pre-filter more aggressively |
| `lfcShrink` says "apeglm requires 'coef'" | You passed `contrast=` with `type="apeglm"` | Switch to `type="ashr"` or convert to a coefficient via `resultsNames(dds)` |
| PCA shows samples cluster by batch, not condition | Batch effect dominates | Add `batch` to design; do not "remove" batch from counts — let the model handle it |
| Top hits dominated by ribosomal / mitochondrial genes | Library composition issue, possibly a normalization artifact | Inspect raw counts; consider filtering rRNA reads upstream |
| SVG files missing, only PNG written | `svglite` not installed | `install.packages("svglite")`; PNG at 300 DPI is fine for almost everything |
| `DESeq()` extremely slow | Tens of thousands of samples or unfiltered low-count genes | Pre-filter harder; `register(MulticoreParam(n))` + `DESeq(dds, parallel = TRUE)` |
| Cook's distance outliers flagged on many genes | One sample with extreme counts | Inspect with `boxplot(log10(assays(dds)[["cooks"]]))`; consider removing the outlier sample and re-running |

## Best Practices

- Use `padj`, never raw `pvalue`, for significance calls. See [bulk-rnaseq-differential-expression](../bulk-rnaseq-differential-expression/SKILL.md) for the methodology rationale.
- Report shrunken LFCs in ranked lists and plots, unshrunk LFCs in headline numerical claims (with the SE).
- Always set the reference level explicitly — alphabetical order will surprise you on the third Friday of a project.
- Pre-filter low-count genes *before* `DESeq()` for speed and stability, not after, and document the cutoff used.
- Save the `DESeqDataSet` (`saveRDS(dds, ...)`) — re-extracting another contrast a week later is much faster than re-running the whole pipeline.
- Record `sessionInfo()` and `packageVersion("DESeq2")` alongside the result tables for reproducibility.

## Suggested Next Steps

1. **[functional-enrichment-from-degs](../functional-enrichment-from-degs/SKILL.md)** — feed the shrunken-LFC ranked list into GSEA (Hallmark + Reactome by default) to identify enriched pathways.
2. **[bulk-rnaseq-differential-expression](../bulk-rnaseq-differential-expression/SKILL.md)** — read this *before* trusting your filter; it covers interpretation pitfalls (FDR vs raw p, LFC thresholds by study type, common errors).
3. **[bulk-omics-clustering](../bulk-omics-clustering/SKILL.md)** — cluster the VST matrix on significant genes to look for subgroup structure (responders/non-responders, latent tissue states).
4. Annotate gene IDs with symbols and descriptions via `AnnotationDbi`/`org.Hs.eg.db` before publication.

## Related Skills

- [bulk-rnaseq-differential-expression](../bulk-rnaseq-differential-expression/SKILL.md) — methodology and interpretation companion to this pipeline skill.
- [functional-enrichment-from-degs](../functional-enrichment-from-degs/SKILL.md) — downstream pathway analysis.
- [bulk-omics-clustering](../bulk-omics-clustering/SKILL.md) — sample/feature clustering on the normalized matrix.

## References

- `references/comprehensive-reference.md` — full code patterns: designs, contrasts, shrinkage, transformations, LRT.
- `references/decision-guide.md` — decision trees for design formula, transformation, shrinkage, thresholds.
- `references/troubleshooting.md` — extended error catalogue with worked fixes.
- `references/usage-guide.md` — quick prompts and example invocations.

**Primary citations:**

- Love MI, Huber W, Anders S (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology 15: 550. DOI: 10.1186/s13059-014-0550-8
- Zhu A, Ibrahim JG, Love MI (2018). *Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.* Bioinformatics 35: 2084-2092. (apeglm) DOI: 10.1093/bioinformatics/bty895
- Stephens M (2017). *False discovery rates: a new deal.* Biostatistics 18: 275-294. (ashr) DOI: 10.1093/biostatistics/kxw041
- Bourgon R, Gentleman R, Huber W (2010). *Independent filtering increases detection power for high-throughput experiments.* PNAS 107: 9546-9551. DOI: 10.1073/pnas.0914005107
