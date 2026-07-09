# scale_and_pca.R
# ----------------------------------------------------------------------
# Run PCA and choose the number of PCs.
#
# Inputs:  Seurat object that has been normalized (SCT or LogNormalize) and
#          has VariableFeatures populated.
# Outputs: Seurat object with reductions$pca; elbow + loadings plots saved.
# Example:
#   source("scripts/scale_and_pca.R")
#   seurat_obj <- run_pca_analysis(seurat_obj, n_pcs = 50)
#   plot_elbow(seurat_obj, output_dir = "results/pca")

#' Run PCA on the active assay using HVGs.
run_pca_analysis <- function(seurat_obj, n_pcs = 50, seed = 42) {
  set.seed(seed)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
  stdev <- seurat_obj@reductions$pca@stdev
  var_ratio <- stdev^2 / sum(stdev^2)
  cum <- cumsum(var_ratio)
  message(sprintf("PCA cumulative variance: PC1-10=%.2f, PC1-20=%.2f, PC1-30=%.2f",
                  cum[min(10, length(cum))],
                  cum[min(20, length(cum))],
                  cum[min(30, length(cum))]))
  message(sprintf("✓ Ran PCA with %d components", n_pcs))
  seurat_obj
}

#' Suggest the number of PCs from variance knee (default 30 if no PCA).
suggest_n_pcs <- function(seurat_obj,
                           min_pcs = 15, max_pcs = 50,
                           variance_threshold = 0.85) {
  if (is.null(seurat_obj@reductions$pca)) {
    message("  No PCA found; defaulting to 30 PCs.")
    return(30)
  }
  stdev <- seurat_obj@reductions$pca@stdev
  var_ratio <- stdev^2 / sum(stdev^2)
  cum <- cumsum(var_ratio)
  n <- which(cum >= variance_threshold)[1]
  if (is.na(n)) n <- length(cum)
  n <- max(min_pcs, min(n, max_pcs))
  message(sprintf("✓ Suggested n_pcs = %d (cumulative variance ≥ %.0f%%)",
                  n, variance_threshold * 100))
  n
}

#' Elbow plot saved to disk.
plot_elbow <- function(seurat_obj,
                        n_show = 50,
                        output_dir = "results/pca") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  p <- Seurat::ElbowPlot(seurat_obj, ndims = n_show)
  out_png <- file.path(output_dir, "pca_elbow.png")
  ggplot2::ggsave(out_png, p, width = 6, height = 4, dpi = 200)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(sub("\\.png$", ".svg", out_png), p, width = 6, height = 4)
  }
  message(sprintf("✓ Saved %s", out_png))
}

#' Heatmaps of top-loading genes for the first PCs.
plot_pca_heatmaps <- function(seurat_obj,
                                n_pcs = 9,
                                output_dir = "results/pca") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_png <- file.path(output_dir, "pca_heatmaps.png")
  grDevices::png(out_png, width = 1200, height = 800, res = 120)
  on.exit(grDevices::dev.off())
  Seurat::DimHeatmap(seurat_obj, dims = 1:n_pcs, cells = 200, balanced = TRUE)
  message(sprintf("✓ Saved %s", out_png))
}

#' Bar plots of top-loading genes per PC.
plot_pca_loadings <- function(seurat_obj,
                                n_pcs = 5,
                                n_top = 10,
                                output_dir = "results/pca") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  p <- Seurat::VizDimLoadings(seurat_obj, dims = 1:n_pcs, reduction = "pca",
                               nfeatures = n_top)
  out_png <- file.path(output_dir, "pca_loadings.png")
  ggplot2::ggsave(out_png, p, width = 3 * n_pcs, height = 4, dpi = 200)
  message(sprintf("✓ Saved %s", out_png))
}
