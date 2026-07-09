# Common Patterns (Seurat)

Compact, copy-paste R recipes for common scenarios.

## Pattern 1: Standard PBMC, single batch

```r
source("scripts/setup_and_import.R")
source("scripts/qc_metrics.R")
source("scripts/filter_cells.R")
source("scripts/normalize_data.R")
source("scripts/scale_and_pca.R")
source("scripts/cluster_cells.R")
source("scripts/run_umap.R")
source("scripts/find_markers.R")
source("scripts/annotate_celltypes.R")
source("scripts/plot_dimreduction.R")
source("scripts/export_results.R")

setup_seurat_libraries()
seurat_obj <- import_10x_data("filtered_feature_bc_matrix/")
seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
seurat_obj <- run_doubletfinder(seurat_obj)
seurat_obj <- filter_by_mad_outliers(seurat_obj, remove_doublets = TRUE)
seurat_obj <- run_sctransform(seurat_obj, vars_to_regress = "percent.mt")
seurat_obj <- run_pca_analysis(seurat_obj, n_pcs = 50)
seurat_obj <- cluster_multiple_resolutions(
  seurat_obj, dims = 1:30, resolutions = c(0.4, 0.6, 0.8, 1.0)
)
seurat_obj <- run_umap_reduction(seurat_obj, dims = 1:30)

markers <- find_all_cluster_markers(seurat_obj, resolution = 0.8)

annotations <- c(
  "0" = "CD4 T cells", "1" = "CD14+ Monocytes",
  "2" = "CD8 T cells", "3" = "B cells",
  "4" = "NK cells",    "5" = "CD16+ Monocytes",
  "6" = "Dendritic cells", "7" = "Megakaryocytes"
)
seurat_obj <- annotate_clusters_manual(seurat_obj, annotations, resolution = 0.8)
plot_umap_clusters(seurat_obj, output_dir = "results/umap")
export_all(seurat_obj, output_dir = "results", all_markers = markers)
```

## Pattern 2: Multi-batch with Harmony

```r
source("scripts/integrate_batches.R")
source("scripts/integration_diagnostics.R")

# After loading multiple samples and annotating "batch" column.
seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
seurat_obj <- batch_mad_outlier_detection(seurat_obj, batch_col = "batch", nmads = 5)
seurat_obj <- run_doubletfinder(seurat_obj, batch_col = "batch")
seurat_obj <- filter_by_mad_outliers(seurat_obj, remove_doublets = TRUE)
seurat_obj <- run_sctransform(seurat_obj, vars_to_regress = "percent.mt")
seurat_obj <- run_pca_analysis(seurat_obj, n_pcs = 50)
seurat_obj <- run_harmony_integration(
  seurat_obj, batch_var = "batch", dims_use = 1:30, theta = 2
)
lisi <- compute_lisi_scores(
  seurat_obj, batch_var = "batch", reduction = "harmony"
)
seurat_obj <- cluster_multiple_resolutions(
  seurat_obj, dims = 1:30, reduction = "harmony",
  resolutions = c(0.4, 0.6, 0.8, 1.0)
)
seurat_obj <- run_umap_reduction(
  seurat_obj, dims = 1:30, reduction = "harmony"
)
markers <- find_all_cluster_markers(seurat_obj, resolution = 0.8)
```

## Pattern 3: Multi-batch with Seurat CCA

```r
source("scripts/integrate_batches.R")
seurat_obj <- run_seurat_cca_integration(
  seurat_obj, split_by = "batch", dims = 1:30
)
DefaultAssay(seurat_obj) <- "integrated"
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 50)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
DefaultAssay(seurat_obj) <- "RNA"   # switch back for plotting + DE
```

## Pattern 4: Condition comparison with pseudobulk

```r
source("scripts/pseudobulk_de.R")

# After clustering and annotation on a multi-sample object.
pb <- aggregate_to_pseudobulk(
  seurat_obj, sample_col = "sample_id",
  celltype_col = "cell_type", min_cells = 10
)
de <- run_pseudobulk_deseq2(
  pb, formula = "~ batch + condition",
  contrast = c("condition", "treated", "control"),
  alpha = 0.05, output_dir = "results/pseudobulk_de"
)

for (ct in names(de)) {
  plot_volcano(de[[ct]], celltype = ct, output_dir = "results/pseudobulk_de")
  plot_ma(de[[ct]], celltype = ct, output_dir = "results/pseudobulk_de")
}
```

Requires N ≥ 2 per condition. The script refuses otherwise.

## Pattern 5: SoupX before QC

```r
source("scripts/remove_ambient_rna.R")
source("scripts/setup_and_import.R")
source("scripts/qc_metrics.R")
source("scripts/filter_cells.R")

seurat_obj <- run_soupx_correction(
  raw_matrix_dir = "raw_feature_bc_matrix/",
  filtered_matrix_dir = "filtered_feature_bc_matrix/",
  output_dir = "results/ambient"
)
seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
seurat_obj <- run_doubletfinder(seurat_obj)
seurat_obj <- filter_by_mad_outliers(seurat_obj, remove_doublets = TRUE)
```

## Pattern 6: SingleR annotation

```r
source("scripts/annotate_celltypes.R")
source("scripts/plot_dimreduction.R")

seurat_obj <- annotate_with_singler(seurat_obj, reference = "HPCA")
plot_umap_clusters(
  seurat_obj, cluster_var = "singler_labels",
  output_dir = "results/annotation"
)
```

References to try (via `celldex`):

- `HPCA`: Human Primary Cell Atlas.
- `BlueprintEncode`: human tissue.
- `ImmGen`: mouse immune.
- `MonacoImmuneData`: high-resolution human immune.

## Pattern 7: Hybrid annotation

```r
source("scripts/annotate_celltypes.R")
source("scripts/find_markers.R")
source("scripts/plot_dimreduction.R")

seurat_obj <- annotate_with_singler(seurat_obj, reference = "HPCA")
markers <- find_all_cluster_markers(seurat_obj, resolution = 0.8)

# Inspect each cluster's top markers and SingleR labels, then refine
refined <- c(
  "0" = "CD4 T cells", "1" = "CD14+ Monocytes",
  "2" = "CD8 T cells", "3" = "B cells",
  "4" = "NK cells",    "5" = "CD16+ Monocytes"
)
seurat_obj <- annotate_clusters_manual(
  seurat_obj, refined, resolution = 0.8,
  celltype_col = "cell_type"
)
```

## Pattern 8: Azimuth label transfer

```r
source("scripts/annotate_celltypes.R")
seurat_obj <- annotate_with_azimuth(seurat_obj, reference = "pbmc")
```

Pre-trained references include `pbmc`, `lung`, `heart`, `kidney`, `bone marrow`, etc.

## Pattern 9: Custom tissue thresholds

```r
source("scripts/filter_cells.R")

seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
seurat_obj <- filter_cells_by_tissue(seurat_obj, tissue = "brain")
# Or fully custom:
seurat_obj <- filter_cells_by_qc(
  seurat_obj,
  min_features = 500, max_features = 6000, max_mt_percent = 10
)
seurat_obj <- run_doubletfinder(seurat_obj)
```

## Pattern 10: LogNormalize for very large data

```r
source("scripts/normalize_data.R")
source("scripts/find_variable_features.R")
source("scripts/scale_and_pca.R")

seurat_obj <- run_lognormalize(seurat_obj)
seurat_obj <- find_hvgs(seurat_obj, n_features = 2000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
seurat_obj <- run_pca_analysis(seurat_obj, n_pcs = 50)
```

Use this when SCTransform runs out of memory (typically > 50k cells without big RAM).
