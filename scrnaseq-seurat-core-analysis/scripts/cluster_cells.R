# cluster_cells.R
# ----------------------------------------------------------------------
# Build a neighbor graph and cluster at one or more resolutions.
#
# Inputs:  Seurat object with PCA (or harmony) reduction.
# Outputs: Seurat object with cluster columns in meta.data.
# Example:
#   source("scripts/cluster_cells.R")
#   seurat_obj <- cluster_multiple_resolutions(
#     seurat_obj, dims = 1:30, resolutions = c(0.4, 0.6, 0.8, 1.0)
#   )

#' Cluster at a single resolution.
cluster_seurat <- function(seurat_obj,
                            dims = 1:30,
                            reduction = "pca",
                            resolution = 0.8,
                            algorithm = 1,
                            random_seed = 42) {
  seurat_obj <- Seurat::FindNeighbors(
    seurat_obj, reduction = reduction, dims = dims, verbose = FALSE
  )
  seurat_obj <- Seurat::FindClusters(
    seurat_obj, resolution = resolution, algorithm = algorithm,
    random.seed = random_seed, verbose = FALSE
  )
  n <- length(levels(Seurat::Idents(seurat_obj)))
  message(sprintf("✓ Clustering at resolution %g: %d clusters", resolution, n))
  seurat_obj
}

#' Sweep multiple resolutions; stores each in a `RNA_snn_res.<res>` column.
cluster_multiple_resolutions <- function(seurat_obj,
                                          dims = 1:30,
                                          reduction = "pca",
                                          resolutions = c(0.4, 0.6, 0.8, 1.0),
                                          algorithm = 1,
                                          random_seed = 42) {
  seurat_obj <- Seurat::FindNeighbors(
    seurat_obj, reduction = reduction, dims = dims, verbose = FALSE
  )
  for (r in resolutions) {
    seurat_obj <- Seurat::FindClusters(
      seurat_obj, resolution = r, algorithm = algorithm,
      random.seed = random_seed, verbose = FALSE
    )
  }
  message(sprintf("✓ Clustering completed at %d resolutions",
                  length(resolutions)))
  seurat_obj
}

#' Louvain alias for `cluster_seurat` with algorithm=1.
cluster_louvain <- function(...) cluster_seurat(..., algorithm = 1)

#' Per-cluster QC summary (n_cells, median QC metrics).
calculate_cluster_qc_stats <- function(seurat_obj,
                                         cluster_col = "seurat_clusters") {
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("`", cluster_col, "` not found in meta.data")
  }
  meta <- seurat_obj@meta.data
  metric_cols <- intersect(c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            colnames(meta))
  stats <- stats::aggregate(
    meta[, metric_cols, drop = FALSE],
    by = list(cluster = meta[[cluster_col]]),
    FUN = median
  )
  stats$n_cells <- as.integer(table(meta[[cluster_col]]))
  message("Cluster QC stats:")
  print(stats)
  invisible(stats)
}
