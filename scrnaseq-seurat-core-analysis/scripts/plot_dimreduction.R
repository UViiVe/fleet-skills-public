# plot_dimreduction.R
# ----------------------------------------------------------------------
# UMAP / PCA plots: cluster overlay, feature plots, batch overlay.
#
# Inputs:  Seurat object with reductions$umap and clustering metadata.
# Outputs: PNG (and optional SVG) figures.
# Example:
#   source("scripts/plot_dimreduction.R")
#   plot_umap_clusters(seurat_obj, output_dir = "results/umap")

.save_gg2 <- function(p, base_path, w = 6, h = 5) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(base_path, p, width = w, height = h, dpi = 200)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(sub("\\.png$", ".svg", base_path), p, width = w, height = h)
  }
  message(sprintf("✓ Saved %s", base_path))
}

#' UMAP coloured by cluster (or any metadata column).
plot_umap_clusters <- function(seurat_obj,
                                 cluster_var = NULL,
                                 output_dir = "results/umap") {
  group <- cluster_var %||% "seurat_clusters"
  p <- Seurat::DimPlot(seurat_obj, group.by = group, label = TRUE,
                       repel = TRUE, label.size = 3) +
    ggplot2::ggtitle(paste0("UMAP — ", group))
  .save_gg2(p, file.path(output_dir, sprintf("umap_%s.png", group)),
             w = 7, h = 6)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Side-by-side UMAP at multiple resolutions.
plot_clustering_comparison <- function(seurat_obj,
                                          resolutions = c(0.4, 0.6, 0.8, 1.0),
                                          output_dir = "results/umap") {
  cols <- intersect(
    paste0("RNA_snn_res.", resolutions),
    colnames(seurat_obj@meta.data)
  )
  if (length(cols) == 0) {
    cols <- intersect(
      paste0("SCT_snn_res.", resolutions),
      colnames(seurat_obj@meta.data)
    )
  }
  if (length(cols) == 0) {
    message("  no resolution columns found")
    return(invisible(NULL))
  }
  plots <- lapply(cols, function(c) {
    Seurat::DimPlot(seurat_obj, group.by = c, label = TRUE,
                    repel = TRUE, label.size = 3) +
      ggplot2::ggtitle(c)
  })
  p <- patchwork::wrap_plots(plots, ncol = 2)
  .save_gg2(p, file.path(output_dir, "umap_resolution_comparison.png"),
             w = 12, h = 4 * ceiling(length(plots) / 2))
}

#' UMAP coloured by one or more features (genes or obs columns).
plot_feature_umap <- function(seurat_obj,
                                features,
                                output_dir = "results/umap") {
  valid <- intersect(features, c(rownames(seurat_obj), colnames(seurat_obj@meta.data)))
  if (length(valid) == 0) {
    message("  none of the requested features were found")
    return(invisible(NULL))
  }
  p <- Seurat::FeaturePlot(seurat_obj, features = valid, ncol = min(3, length(valid)),
                            cols = c("lightgrey", "firebrick"))
  .save_gg2(p, file.path(output_dir, "umap_features.png"),
             w = 4 * min(3, length(valid)),
             h = 3 * ceiling(length(valid) / 3))
}

#' PCA scatter coloured by cluster/metadata.
plot_pca_clusters <- function(seurat_obj,
                                 group_by = NULL,
                                 output_dir = "results/pca") {
  group <- group_by %||% "seurat_clusters"
  p <- Seurat::DimPlot(seurat_obj, reduction = "pca", group.by = group, label = TRUE)
  .save_gg2(p, file.path(output_dir, sprintf("pca_%s.png", group)),
             w = 7, h = 6)
}

#' UMAP split by a categorical metadata column.
plot_split_umap <- function(seurat_obj,
                              split_by = "batch",
                              group_by = "seurat_clusters",
                              output_dir = "results/umap") {
  if (!split_by %in% colnames(seurat_obj@meta.data)) {
    message("  ", split_by, " not in meta.data")
    return(invisible(NULL))
  }
  p <- Seurat::DimPlot(seurat_obj, group.by = group_by, split.by = split_by,
                       ncol = 2, label = TRUE, repel = TRUE, label.size = 3)
  n_levels <- length(unique(seurat_obj@meta.data[[split_by]]))
  .save_gg2(p, file.path(output_dir, sprintf("umap_split_%s.png", split_by)),
             w = 6 * min(2, n_levels), h = 5 * ceiling(n_levels / 2))
}
