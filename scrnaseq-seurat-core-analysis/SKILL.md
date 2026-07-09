---
name: scrnaseq-seurat-core-analysis
id: "skill_6f14308e7afb492081b12e0e3eedd73a"
display-name: "scRNA-seq Seurat Core Analysis"
short-description: "End-to-end scRNA-seq analysis in R using Seurat v5 — QC, SCTransform, Harmony / CCA / RPCA integration, clustering, annotation, pseudobulk DE. Use when working in R/Bioconductor and the upstream st..."
category: transcriptomics
keywords: [seurat, scRNA-seq, clustering, R, single cell]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "End-to-end scRNA-seq analysis in R using Seurat v5 — QC, SCTransform, Harmony / CCA / RPCA integration, clustering, annotation, pseudobulk DE. Use when working in R/Bioconductor and the upstream stack is already Seurat-based."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "Seurat v5"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Hao Y et al., Nature Biotechnology 2024, doi:10.1038/s41587-023-01767-y"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Single-Cell RNA-seq Core Analysis (Seurat)

A complete R workflow for single-cell RNA-seq using Seurat v5. Takes a raw count matrix through QC, normalization, integration, clustering, annotation, and pseudobulk differential expression. Designed to be drop-in for labs already using Seurat / Bioconductor — the same `seurat_obj` flows through every step.

## When to Use This Skill

Use this skill when:

- You have 10X Chromium output, Drop-seq, Smart-seq2, or other scRNA-seq counts and you want to stay in R.
- You're working alongside an existing Seurat / Bioconductor pipeline (DESeq2, muscat, SingleR, celldex).
- You have multiple batches and want Harmony / CCA / RPCA integration with diagnostic LISI/ASW scores.
- You have biological replicates and a condition contrast and need pseudobulk + DESeq2 to do it correctly.
- You want curated defaults (SCTransform, MAD-based QC, multi-resolution clustering) rather than re-deriving them.

Do NOT use this skill for:

- **Python-based scRNA-seq.** Use `scrnaseq-scanpy-core-analysis` for `scanpy` + `scvi-tools`.
- **Bulk RNA-seq.** Use `bulk-rnaseq-counts-to-de-deseq2`.
- **Spatial transcriptomics.** Use `spatial-transcriptomics` — Visium / Xenium / Slide-seq have different QC and neighbor assumptions.
- **Trajectory inference as the primary question.** Annotate here, then pass the object to `scrna-trajectory-inference` (Monocle3 / Slingshot / RNA velocity).
- **Cell-cell communication.** Use `cell-cell-communication` (CellChat / LIANA / NicheNet) on the annotated object.

## Quick Start

PBMC 3k in ~15 minutes on a laptop:

```r
source("scripts/setup_and_import.R")
source("scripts/load_example_data.R")
source("scripts/qc_metrics.R")
source("scripts/filter_cells.R")
source("scripts/normalize_data.R")
source("scripts/scale_and_pca.R")
source("scripts/cluster_cells.R")
source("scripts/run_umap.R")
source("scripts/find_markers.R")
source("scripts/plot_dimreduction.R")
source("scripts/annotate_celltypes.R")
source("scripts/export_results.R")

setup_seurat_libraries()
seurat_obj <- load_seurat_data("pbmc3k")

seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
seurat_obj <- filter_cells_by_qc(seurat_obj, min_features = 200, max_features = 2500, max_mt_percent = 5)
seurat_obj <- run_sctransform(seurat_obj, vars_to_regress = "percent.mt")
seurat_obj <- run_pca_analysis(seurat_obj, n_pcs = 50)
seurat_obj <- cluster_multiple_resolutions(seurat_obj, dims = 1:30, resolutions = c(0.4, 0.6, 0.8, 1.0))
seurat_obj <- run_umap_reduction(seurat_obj, dims = 1:30)
markers   <- find_all_cluster_markers(seurat_obj, resolution = 0.8)
plot_umap_clusters(seurat_obj, output_dir = "results/umap")
export_all(seurat_obj, output_dir = "results", all_markers = markers)
```

Runtime: ~3 minutes for PBMC 3k on a 4-core laptop. SCTransform on 30k+ cells benefits from a server and 32+ GB of RAM.

## Installation

| Package        | Version  | License | Commercial use | Install command                              |
| -------------- | -------- | ------- | -------------- | -------------------------------------------- |
| R              | ≥ 4.2    | GPL-2   | Permitted      | OS package manager                           |
| Seurat         | ≥ 5.0    | MIT     | Permitted      | `install.packages("Seurat")`                 |
| SeuratObject   | ≥ 5.0    | MIT     | Permitted      | `install.packages("SeuratObject")`           |
| ggplot2        | ≥ 3.4    | MIT     | Permitted      | `install.packages("ggplot2")`                |
| patchwork      | ≥ 1.1    | MIT     | Permitted      | `install.packages("patchwork")`              |
| dplyr          | ≥ 1.1    | MIT     | Permitted      | `install.packages("dplyr")`                  |
| Matrix         | ≥ 1.6    | GPL-2   | Permitted      | `install.packages("Matrix")`                 |
| harmony        | ≥ 1.0    | GPL-3   | Permitted      | `install.packages("harmony")`                |
| DoubletFinder  | ≥ 2.0.3  | CC0     | Permitted      | `remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")` |
| SoupX          | ≥ 1.6    | GPL-2   | Permitted      | `install.packages("SoupX")`                  |
| DESeq2         | ≥ 1.40   | LGPL-3  | Permitted      | `BiocManager::install("DESeq2")`             |
| muscat         | ≥ 1.14   | GPL-3   | Permitted      | `BiocManager::install("muscat")`             |
| SingleR        | ≥ 2.4    | GPL-3   | Permitted      | `BiocManager::install("SingleR")`            |
| celldex        | ≥ 1.12   | GPL-3   | Permitted      | `BiocManager::install("celldex")`            |
| ggrepel        | ≥ 0.9    | GPL-3   | Permitted      | `install.packages("ggrepel")`                |

`scripts/install_packages.R` will install everything in one go. Bioconductor packages require `BiocManager` first: `install.packages("BiocManager")`.

## Inputs

**Required (one of):**

- CellRanger directory: `filtered_feature_bc_matrix/` (or `raw_feature_bc_matrix/` if you'll run SoupX first).
- 10X `.h5` file.
- Saved Seurat object (`.rds`).
- Count matrix (CSV/TSV) — genes × cells with row/column names.
- A SeuratData dataset name (`pbmc3k`, `ifnb`).

**Optional:**

- Sample metadata CSV with `sample_id`, `condition`, `batch`, `donor`, `tissue` etc.
- Custom marker gene table for annotation.
- Reference Seurat object for label transfer (Azimuth) or SingleR reference set.

**Data requirements:** Human (GRCh38) or mouse (GRCm39) work out of the box. Minimum ~500 cells per sample (1000+ recommended). Default QC thresholds assume 10X Chromium v3; adjust for Smart-seq2 (higher gene counts) and Drop-seq (lower).

## Outputs

**Seurat artifact (primary):**

- `seurat_processed.rds` — annotated object with:
  - `Assays$RNA` raw counts and `Assays$RNA$data` log-normalized values.
  - `Assays$SCT` (if SCTransform used) with `counts`, `data`, `scale.data`.
  - `meta.data` with QC metrics (`nFeature_RNA`, `nCount_RNA`, `percent.mt`, `percent.ribo`, `doublet_score`), clustering at multiple resolutions, and cell-type labels.
  - `reductions$pca`, `reductions$umap`, and (when applicable) `reductions$harmony` or `reductions$integrated.cca`.
  - Marker tables under `misc$markers`.

Load with: `seurat_obj <- readRDS("seurat_processed.rds")`.

**CSV exports** (for downstream Python consumers / publication tables):

- `cell_metadata.csv` — `meta.data` table.
- `normalized_counts.csv` / `normalized_counts.rds` — log-normalized expression.
- `pca_coordinates.csv`, `umap_coordinates.csv` — embeddings.
- `cluster_markers_all.csv` — Wilcoxon ranked genes.
- `{celltype}_deseq2_results.csv` — pseudobulk DE per cell type.

**Figures (PNG + SVG at 300 DPI):**

- QC: violin and scatter plots before/after filtering, per-batch.
- Dimensionality reduction: PCA elbow, UMAP by cluster, by batch, by feature.
- Markers: heatmap of top markers per cluster, dot plot, violin plots of canonical markers.
- DE: per-cell-type volcano and MA plots.

**PDF report style rules** (when emitting a report):

- US Letter (8.5 × 11 in); set page size explicitly.
- No Unicode superscripts in numeric output — write `3.4e-06` or `3.4 × 10^(-6)`.
- Multi-panel figures ≥ 80% page width.
- Keep section headings glued to the first paragraph of content.

## Clarification Questions

**⚠️ CRITICAL: ASK THIS FIRST**

**Q1. Where is the data, and is it filtered or raw?**

- Path to CellRanger directory, `.h5`, `.rds`, or SeuratData dataset name?
- Filtered (CellRanger has called cells) or raw (SoupX recommended first for high-soup tissues)?
- Or run the bundled PBMC 3k demo first to verify the pipeline?

**Q2. Species and tissue.**

- Human (default `MT-`, `^RP[SL]`) or mouse (`mt-`, `^Rp[sl]`).
- Tissue dictates QC tolerance — PBMC at 5% MT vs cardiomyocytes naturally ~ 20%.

**Q3. Single batch or multi-batch?**

- Single sample: skip Harmony, fixed thresholds work.
- Multi-batch: batch-aware MAD QC, then Harmony (default), CCA, or RPCA.

**Q4. Comparing conditions?**

- Yes: how many biological replicates per condition? Pseudobulk + DESeq2 needs N ≥ 2 per group (≥ 3 is honest).
- No: stop after annotation.

**Q5. Clustering granularity.**

- Major lineages only: resolution 0.3–0.5.
- Standard (subtypes): 0.6–0.8.
- Fine (substates): 1.0–1.5.
- Default: sweep 0.4 / 0.6 / 0.8 / 1.0 and pick by marker coherence and stability.

## Standard Workflow

Scripts under `scripts/` are the canonical entry points. They source clean — use them rather than open-coding equivalents in your notebook. The failure hierarchy (next section) covers what to do if a script throws.

### Phase 1 — Load and QC

Scripts: `setup_and_import.R`, `load_example_data.R`, `qc_metrics.R`, `plot_qc.R`, `filter_cells.R`.

```r
source("scripts/setup_and_import.R")
source("scripts/qc_metrics.R")
source("scripts/filter_cells.R")

setup_seurat_libraries()
seurat_obj <- import_10x_data("filtered_feature_bc_matrix/")
seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
seurat_obj <- batch_mad_outlier_detection(seurat_obj, batch_col = "batch", nmads = 5)
seurat_obj <- run_doubletfinder(seurat_obj, batch_col = "batch")
seurat_obj <- filter_by_mad_outliers(seurat_obj, remove_doublets = TRUE)
```

For raw matrices in high-soup tissues (brain, lung, tumor), prepend SoupX via `scripts/remove_ambient_rna.R`. Target ≥ 70% cell retention; widen MAD or drop back to tissue-fixed thresholds if you fall under.

### Phase 2 — Normalize and reduce dimensions

Scripts: `normalize_data.R`, `find_variable_features.R`, `scale_and_pca.R`.

```r
source("scripts/normalize_data.R")
source("scripts/scale_and_pca.R")

# SCTransform is the default for UMI data.
seurat_obj <- run_sctransform(seurat_obj, vars_to_regress = "percent.mt")

# Alternative: classic LogNormalize + ScaleData
# seurat_obj <- run_lognormalize(seurat_obj)

seurat_obj <- run_pca_analysis(seurat_obj, n_pcs = 50)
n_pcs <- suggest_n_pcs(seurat_obj)
```

**Never use < 15 PCs.** The standard range is 20–30. `suggest_n_pcs()` is variance-knee–based; trust it over your eye on the elbow plot.

### Phase 3 — Integrate (multi-batch only)

Scripts: `integrate_batches.R`, `integration_diagnostics.R`.

```r
source("scripts/integrate_batches.R")
source("scripts/integration_diagnostics.R")

# Harmony — fast, simple, robust to over-correction (default)
seurat_obj <- run_harmony_integration(seurat_obj, batch_var = "batch", dims_use = 1:30)

# Or Seurat CCA / RPCA when batch effects are complex
# seurat_obj <- run_seurat_cca_integration(seurat_obj, split_by = "batch")
# seurat_obj <- run_seurat_rpca_integration(seurat_obj, split_by = "batch")

lisi <- compute_lisi_scores(seurat_obj, batch_var = "batch", reduction = "harmony")
```

Targets: batch iLISI ≈ `n_batches`, cell-type cLISI ≈ 1. If your iLISI is way below max and cells still cluster by batch, raise `theta` to 3 or try CCA.

### Phase 4 — Cluster, UMAP, find markers

Scripts: `cluster_cells.R`, `run_umap.R`, `find_markers.R`, `plot_dimreduction.R`.

```r
source("scripts/cluster_cells.R")
source("scripts/run_umap.R")
source("scripts/find_markers.R")
source("scripts/plot_dimreduction.R")

reduction <- if ("harmony" %in% names(seurat_obj@reductions)) "harmony" else "pca"
seurat_obj <- cluster_multiple_resolutions(
  seurat_obj, dims = 1:n_pcs, reduction = reduction,
  resolutions = c(0.4, 0.6, 0.8, 1.0)
)
seurat_obj <- run_umap_reduction(seurat_obj, dims = 1:n_pcs, reduction = reduction)

markers <- find_all_cluster_markers(seurat_obj, resolution = 0.8)
plot_umap_clusters(seurat_obj, output_dir = "results/umap")
plot_top_markers_heatmap(seurat_obj, markers, n_top = 10, output_dir = "results/markers")
```

`FindAllMarkers` runs Wilcoxon by default. These statistics are exploratory — they characterize clusters, they do not test condition effects.

### Phase 5 — Annotate cell types

Script: `annotate_celltypes.R`. Reference: `references/marker_gene_database.md`.

```r
source("scripts/annotate_celltypes.R")

annotations <- c(
  "0" = "CD4 T cells",
  "1" = "CD14+ Monocytes",
  "2" = "CD8 T cells",
  "3" = "B cells",
  "4" = "NK cells",
  "5" = "CD16+ Monocytes"
)
seurat_obj <- annotate_clusters_manual(seurat_obj, annotations, resolution = 0.8)

# Optional automated validation
seurat_obj <- annotate_with_singler(seurat_obj, reference = "HPCA")
# or
seurat_obj <- annotate_with_azimuth(seurat_obj, reference = "pbmc")
```

The hybrid pattern — manual labels grounded in `markers`, validated with SingleR — beats either approach alone.

### Phase 6 — Pseudobulk DE (multi-sample, multi-condition)

Script: `pseudobulk_de.R`. Reference: `references/pseudobulk_de_guide.md`.

```r
source("scripts/pseudobulk_de.R")

pb <- aggregate_to_pseudobulk(
  seurat_obj, sample_col = "sample_id", celltype_col = "cell_type", min_cells = 10
)
de_results <- run_pseudobulk_deseq2(
  pb, formula = "~ batch + condition",
  contrast = c("condition", "treated", "control")
)

for (ct in names(de_results)) {
  plot_volcano(de_results[[ct]], celltype = ct, output_dir = "results/de")
}
```

The function refuses to run if any condition has N = 1. That's correct. Cell-level `FindMarkers` with `ident.1 = "treated"` would technically produce p-values, but pseudoreplication makes them dishonest — use Wilcoxon at the cell level only as exploratory candidate-gene mining and label it as such.

### Phase 7 — Export

Script: `export_results.R`.

```r
source("scripts/export_results.R")
export_all(seurat_obj, output_dir = "results", all_markers = markers)
```

Writes `seurat_processed.rds`, normalized count matrices, metadata, embeddings, markers, and `analysis_summary.txt`. The export is idempotent.

## When Scripts Fail

Use this hierarchy, in order. Document the level you used.

1. **Fix and retry (90% of failures):** install the missing package (BiocManager / remotes), re-run. R is fragile about dependency order — DoubletFinder needs Seurat installed first; SingleR pulls in SummarizedExperiment.
2. **Modify the script in place (5%):** edit the function, comment the change, re-run.
3. **Use the script as reference (4%):** read it, write a one-off `sourced` cell that follows the same algorithm, cite the source.
4. **Write from scratch (1%):** only when the algorithm itself doesn't fit your data. Justify in the report.

Never skip past (1) and write inline `Seurat::FindClusters` / `DimPlot` / `FindAllMarkers` without trying the script first.

## Decision Points

| Decision                             | Default                              | Override when                                                                            |
| ------------------------------------ | ------------------------------------ | ---------------------------------------------------------------------------------------- |
| Ambient RNA correction               | Skip (filtered PBMC)                 | Raw matrix or high-soup tissue (brain/lung/tumor) → run SoupX                            |
| QC thresholding                      | MAD with `nmads=5`                   | Single batch with strong tissue priors → `filter_cells_by_tissue(tissue=...)`            |
| DoubletFinder expected rate          | Auto: ~ 0.075 per 1000 cells         | Plate-based (Smart-seq2) → skip; very high cell loading → 10–12%                          |
| Normalization                        | SCTransform (`vst.flavor="v2"`)      | Memory pressure (> 50k cells) → LogNormalize; non-UMI (Smart-seq2) → LogNormalize        |
| `n_pcs`                              | 30                                   | Use `suggest_n_pcs()`; never < 15                                                        |
| Integration                          | Harmony (`theta=2`)                  | Complex batch structure → CCA; > 100k cells with similar comps → RPCA                    |
| Clustering resolution                | Sweep 0.4 / 0.6 / 0.8 / 1.0          | Pick by marker coherence and stability                                                   |
| Annotation                           | Manual from markers                  | Quick-look on standard tissues → SingleR or Azimuth (treat as validation, not truth)     |
| Cell-level DE for condition contrast | **Never**                            | Always pseudobulk + DESeq2 for inference                                                 |

## Common Issues

| Issue                                                              | Likely cause                                                            | Fix                                                                                                |
| ------------------------------------------------------------------ | ----------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------- |
| `there is no package called 'Seurat'`                              | Package not installed                                                   | `install.packages("Seurat")`; ensure R ≥ 4.2                                                       |
| `Error: package 'SeuratObject' could not be loaded`                | Version mismatch with Seurat                                            | `install.packages(c("Seurat", "SeuratObject"))` together                                            |
| Cell retention < 50% after QC                                      | MAD too strict or genuinely degraded sample                             | Raise `nmads` to 7, or visualize `percent.mt` distribution before deciding                         |
| `Cannot find mitochondrial genes`                                  | Wrong species pattern or Ensembl IDs in `rownames`                      | Pass correct `species`; convert Ensembl → symbols with `biomaRt` or `annotables`                    |
| DoubletFinder predicts ~ 100% doublets                             | Expected rate misconfigured                                             | Use the script's auto-rate; if overriding, check that you didn't pass a fraction > 1               |
| `Error: 'harmony' not found in reductions`                         | Tried to cluster on harmony before running integration                  | Run `run_harmony_integration()` first, or fall back to `reduction = "pca"`                          |
| SCTransform crashes with `memory not mapped`                       | RAM exhausted on > 30–50k cells                                         | Switch to `run_lognormalize()`; or run on a server with > 64 GB                                     |
| UMAP looks like one blob                                           | Too few PCs, or wrong reduction                                         | Verify `reduction` arg matches what was integrated; raise `dims` to 1:50                            |
| Top markers are all `RPS*` / `RPL*` / `MT-*`                       | Forgot to exclude them before ranking                                   | `seurat_obj <- subset(seurat_obj, features = setdiff(rownames(seurat_obj), grep("^(MT-|RPS|RPL)", rownames(seurat_obj), value = TRUE)))` |
| `Error in DESeqDataSetFromMatrix(...)`                             | Mean instead of sum aggregation; or normalized counts as input          | `aggregate_to_pseudobulk()` already does it right; use the script                                  |
| Pseudobulk DE refuses to run (N=1)                                 | Script guard fires when any condition has one sample                    | Don't fight it — report N=1 as a design limitation; use cell-level Wilcoxon only as exploratory     |
| `Error: 'svglite' is not installed`                                | Plot script falls back to PNG only                                      | Optional — install `svglite` for SVG, or accept PNG-only output                                    |
| `Error: object 'seurat_obj' not found` after kernel restart        | Session lost the object                                                 | Re-run from `import_*` or `readRDS("results/seurat_processed.rds")`                                 |
| `cluster_multiple_resolutions` very slow on large data             | Default n_neighbors too high                                            | Drop `k.param` to 10; integrate on Harmony before clustering                                       |
| SingleR returns mostly `NA` labels                                 | Reference doesn't match tissue                                          | Pick the right reference: `HPCA` for human PBMC, `BlueprintEncode` for human tissue, `ImmGen` for mouse immune |

Deeper triage: `references/troubleshooting_guide.md`.

## Best Practices

- **Always store raw counts.** Seurat v5 keeps them in `Assays$RNA$counts`. If you create a new assay (SCT), do NOT delete RNA — pseudobulk DE needs the raw counts.
- **Set the seed.** `set.seed(42)` at the top of your script. UMAP, Leiden, DoubletFinder, SCTransform — all are stochastic.
- **Visualize QC before filtering.** A 12% MT peak in cardiomyocytes is biology, not filtering signal.
- **Document integration parameters.** Store them in `seurat_obj@misc$integration` so the next person can reproduce.
- **Don't use SCT-normalized counts for DE inference.** Pseudobulk on `RNA$counts`, with batch as a covariate.
- **Re-cluster after every QC change.** Cluster labels are downstream of filtering; you can't reuse them after re-filtering.
- **Save markers separately.** Don't rely only on `seurat_obj@misc$markers` — write the CSV too.

## Suggested Next Steps

1. **`scrna-trajectory-inference`** — pass the annotated object to Monocle3 / Slingshot / RNA velocity (`velocyto.R`) when the tissue has developmental structure.
2. **`cell-cell-communication`** — once cell types are annotated, run CellChat or LIANA on the same object to score ligand–receptor activity.
3. **`functional-enrichment-from-degs`** — take the pseudobulk DE gene lists into Enrichr / clusterProfiler / fgsea for pathway analysis.
4. **`scrnaseq-scanpy-core-analysis`** — convert to AnnData with `SeuratDisk::SaveH5Seurat` if you need the scverse downstream stack (scVI, scANVI, CellTypist).
5. **`spatial-transcriptomics`** — if you also have Visium / Xenium / Slide-seq from the same tissue, use this object as the dissociated reference for cell-type deconvolution.

## Related Skills

- **Python alternative:** `scrnaseq-scanpy-core-analysis`.
- **Downstream interpretation:** `functional-enrichment-from-degs`, `coexpression-network`, `upstream-regulator-analysis`.
- **Adjacent omics:** `bulk-rnaseq-differential-expression`, `bulk-omics-clustering`.
- **Experiment planning:** `experimental-design-statistics` for sample size and replicate calculations.

## References

Files under `references/`:

- `qc_guidelines.md` — tissue-specific QC thresholds and adaptive filtering recipes.
- `seurat_best_practices.md` — SCTransform, LogNormalize, PCA, clustering, UMAP defaults.
- `integration_methods.md` — Harmony vs CCA vs RPCA, with LISI/ASW interpretation.
- `marker_gene_database.md` — curated cell-type markers across PBMC, brain, lung, gut, stromal, endothelial.
- `pseudobulk_de_guide.md` — pseudoreplication, DESeq2 design formulas, minimum sample size.
- `ambient_rna_correction.md` — SoupX (and CellBender alternative).
- `troubleshooting_guide.md` — error messages and fixes.
- `workflow-details.md` — expanded code per phase.
- `common-patterns.md` — copy-paste recipes.
- `decision-guide.md` — choose between normalization, integration, annotation strategies.

- **[Seurat v5]**: Hao Y, Stuart T, Kowalski MH, et al. "Dictionary learning for integrative, multimodal and scalable single-cell analysis." *Nature Biotechnology* 2024;42:293–304. doi:10.1038/s41587-023-01767-y — [link](https://doi.org/10.1038/s41587-023-01767-y)
- **[Harmony]**: Korsunsky I, Millard N, Fan J, et al. "Fast, sensitive and accurate integration of single-cell data with Harmony." *Nature Methods* 2019;16:1289–1296. doi:10.1038/s41592-019-0619-0 — [link](https://doi.org/10.1038/s41592-019-0619-0)
- **[SCTransform]**: Hafemeister C, Satija R. "Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression." *Genome Biology* 2019;20:296. doi:10.1186/s13059-019-1874-1 — [link](https://doi.org/10.1186/s13059-019-1874-1)
- **[DoubletFinder]**: McGinnis CS, Murrow LM, Gartner ZJ. "DoubletFinder: doublet detection in single-cell RNA sequencing data using artificial nearest neighbors." *Cell Systems* 2019;8:329–337. doi:10.1016/j.cels.2019.03.003 — [link](https://doi.org/10.1016/j.cels.2019.03.003)
- **[Pseudobulk DE]**: Squair JW, Gautier M, Kathe C, et al. "Confronting false discoveries in single-cell differential expression." *Nature Communications* 2021;12:5692. doi:10.1038/s41467-021-22786-0 — [link](https://doi.org/10.1038/s41467-021-22786-0)
- **[SingleR]**: Aran D, Looney AP, Liu L, et al. "Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage." *Nature Immunology* 2019;20:163–172. doi:10.1038/s41590-018-0276-y — [link](https://doi.org/10.1038/s41590-018-0276-y)
- **[muscat]**: Crowell HL, Soneson C, Germain P-L, et al. "muscat detects subpopulation-specific state transitions from multi-sample multi-condition single-cell transcriptomics data." *Nature Communications* 2020;11:6077. doi:10.1038/s41467-020-19894-4 — [link](https://doi.org/10.1038/s41467-020-19894-4)
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)