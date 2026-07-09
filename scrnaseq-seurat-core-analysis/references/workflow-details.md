# Workflow Details (Seurat)

Expanded R code for every step. The `scripts/` directory is the canonical source — this reference shows how the pieces fit together.

## Phase 1: QC and filtering

### Step 1.1 — Optional: SoupX (raw input only)

Skip if starting from `filtered_feature_bc_matrix/`. Otherwise:

```r
source("scripts/remove_ambient_rna.R")
seurat_obj <- run_soupx_correction(
  raw_matrix_dir = "raw_feature_bc_matrix/",
  filtered_matrix_dir = "filtered_feature_bc_matrix/",
  output_dir = "results/ambient"
)
```

### Step 1.2 — Load and QC

```r
source("scripts/setup_and_import.R")
source("scripts/qc_metrics.R")
source("scripts/plot_qc.R")

setup_seurat_libraries()
seurat_obj <- import_10x_data("filtered_feature_bc_matrix/")
seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
seurat_obj <- batch_mad_outlier_detection(
  seurat_obj, batch_col = "batch", nmads = 5
)

plot_qc_violin(seurat_obj, batch_col = "batch", output_dir = "results/qc")
plot_qc_scatter(seurat_obj, output_dir = "results/qc")
```

### Step 1.3 — Doublet detection and filtering

```r
source("scripts/filter_cells.R")

seurat_obj <- run_doubletfinder(seurat_obj, batch_col = "batch")
seurat_obj <- filter_by_mad_outliers(
  seurat_obj, remove_doublets = TRUE
)
message(sprintf(
  "Retained: %d / %d (%.1f%%)",
  ncol(seurat_obj), n_before,
  100 * ncol(seurat_obj) / n_before
))
```

## Phase 2: Normalization and feature selection

### Step 2.1 — Normalize

SCTransform (default):

```r
source("scripts/normalize_data.R")
seurat_obj <- run_sctransform(
  seurat_obj, vars_to_regress = "percent.mt"
)
```

LogNormalize (alternative for very large data or non-UMI):

```r
seurat_obj <- run_lognormalize(seurat_obj)
seurat_obj <- find_hvgs(seurat_obj, n_features = 2000)
```

### Step 2.2 — PCA

```r
source("scripts/scale_and_pca.R")
seurat_obj <- run_pca_analysis(seurat_obj, n_pcs = 50)
plot_elbow(seurat_obj, output_dir = "results/pca")
n_pcs <- suggest_n_pcs(seurat_obj)
```

## Phase 3: Integration (multi-batch only)

### Step 3.1 — Harmony (default)

```r
source("scripts/integrate_batches.R")
seurat_obj <- run_harmony_integration(
  seurat_obj, batch_var = "batch", dims_use = 1:30, theta = 2
)
```

### Step 3.2 — Seurat CCA

```r
seurat_obj <- run_seurat_cca_integration(
  seurat_obj, split_by = "batch", dims = 1:30
)
```

### Step 3.3 — Seurat RPCA

```r
seurat_obj <- run_seurat_rpca_integration(
  seurat_obj, split_by = "batch", dims = 1:30
)
```

### Step 3.4 — Integration diagnostics

```r
source("scripts/integration_diagnostics.R")

lisi <- compute_lisi_scores(
  seurat_obj, batch_var = "batch",
  label_var = "cell_type", reduction = "harmony"
)
asw <- compute_asw_scores(
  seurat_obj, batch_var = "batch",
  label_var = "cell_type", reduction = "harmony"
)
plot_integration_metrics(seurat_obj, lisi, asw, output_dir = "results/integration")
compare_integration_methods(
  seurat_obj, batch_var = "batch",
  reductions = c("pca", "harmony"),
  output_dir = "results/integration"
)
```

Targets: batch iLISI ≈ `n_batches`, cell-type cLISI ≈ 1, batch ASW ≈ 0, cell-type ASW > 0.5.

## Phase 4: Clustering and visualization

### Step 4.1 — Cluster

```r
source("scripts/cluster_cells.R")
reduction <- if ("harmony" %in% names(seurat_obj@reductions)) "harmony" else "pca"
seurat_obj <- cluster_multiple_resolutions(
  seurat_obj, dims = 1:n_pcs, reduction = reduction,
  resolutions = c(0.4, 0.6, 0.8, 1.0)
)
```

### Step 4.2 — UMAP

```r
source("scripts/run_umap.R")
seurat_obj <- run_umap_reduction(
  seurat_obj, dims = 1:n_pcs, reduction = reduction, seed = 42
)
```

### Step 4.3 — UMAP plots

```r
source("scripts/plot_dimreduction.R")
plot_umap_clusters(seurat_obj, output_dir = "results/umap")
plot_clustering_comparison(
  seurat_obj, resolutions = c(0.4, 0.6, 0.8, 1.0),
  output_dir = "results/umap"
)
plot_feature_umap(
  seurat_obj,
  features = c("nFeature_RNA", "percent.mt", "CD3D", "CD8A", "CD14", "MS4A1"),
  output_dir = "results/umap"
)
```

## Phase 5: Markers and annotation

### Step 5.1 — Find markers

```r
source("scripts/find_markers.R")
markers <- find_all_cluster_markers(
  seurat_obj, resolution = 0.8, method = "wilcox",
  min_pct = 0.25, logfc_threshold = 0.25
)
export_marker_tables(markers, output_dir = "results/markers")
plot_top_markers_heatmap(seurat_obj, markers, n_top = 10, output_dir = "results/markers")
plot_markers_dotplot(seurat_obj, markers, n_top = 5, output_dir = "results/markers")
```

### Step 5.2 — Annotate

Manual:

```r
source("scripts/annotate_celltypes.R")
annotations <- c(
  "0" = "CD4 T cells", "1" = "CD14+ Monocytes",
  "2" = "CD8 T cells", "3" = "B cells",
  "4" = "NK cells",    "5" = "CD16+ Monocytes",
  "6" = "Dendritic cells", "7" = "Megakaryocytes"
)
seurat_obj <- annotate_clusters_manual(
  seurat_obj, annotations, resolution = 0.8
)
plot_annotated_umap(seurat_obj, output_dir = "results/annotation")
```

SingleR:

```r
seurat_obj <- annotate_with_singler(seurat_obj, reference = "HPCA")
```

Azimuth:

```r
seurat_obj <- annotate_with_azimuth(seurat_obj, reference = "pbmc")
```

Hybrid (recommended): SingleR/Azimuth for the first pass → manual refinement against markers.

## Phase 6: Pseudobulk DE

```r
source("scripts/pseudobulk_de.R")

pb <- aggregate_to_pseudobulk(
  seurat_obj, sample_col = "sample_id",
  celltype_col = "cell_type", min_cells = 10
)

de_results <- run_pseudobulk_deseq2(
  pb, formula = "~ batch + condition",
  contrast = c("condition", "treated", "control"),
  alpha = 0.05, output_dir = "results/pseudobulk_de"
)

export_de_results(de_results, output_dir = "results/pseudobulk_de")

for (ct in names(de_results)) {
  plot_volcano(de_results[[ct]], celltype = ct,
                padj_threshold = 0.05, log2fc_threshold = 0.5,
                output_dir = "results/pseudobulk_de")
  plot_ma(de_results[[ct]], celltype = ct, output_dir = "results/pseudobulk_de")
}
```

The DE script refuses to run with N=1 in any condition. Use `FindMarkers` as exploratory only in that case.

## Phase 7: Export

```r
source("scripts/export_results.R")
export_all(
  seurat_obj, output_dir = "results",
  all_markers = markers
)
```

Outputs:

- `seurat_processed.rds` — annotated object.
- `normalized_counts.csv` / `.rds`.
- `cell_metadata.csv`.
- `pca_coordinates.csv`, `umap_coordinates.csv`.
- `cluster_markers_all.csv`.
- `analysis_summary.txt` — dataset shape, retention, integration scores.

## End-to-end example

See `common-patterns.md` for compact recipes covering PBMC, multi-batch, condition contrasts, SingleR + manual hybrid annotation.
