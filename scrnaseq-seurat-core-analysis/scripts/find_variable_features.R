# find_variable_features.R
# ----------------------------------------------------------------------
# Identify highly variable genes (HVGs) for downstream PCA and integration.
#
# Inputs:  Seurat object (any normalization). For SCT, HVGs are computed
#          automatically and stored in the SCT assay.
# Outputs: Seurat object with VariableFeatures populated.
# Example:
#   source("scripts/find_variable_features.R")
#   seurat_obj <- find_hvgs(seurat_obj, n_features = 2000)

#' Identify HVGs via the VST method. Works on the RNA assay.
find_hvgs <- function(seurat_obj,
                       n_features = 2000,
                       selection_method = "vst") {
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj,
    selection.method = selection_method,
    nfeatures = n_features,
    verbose = FALSE
  )
  n <- length(Seurat::VariableFeatures(seurat_obj))
  message(sprintf("✓ Selected %d HVGs (method=%s)", n, selection_method))
  seurat_obj
}

#' Plot the mean-variance scatter highlighting HVGs.
plot_variable_features <- function(seurat_obj,
                                     top_n = 20,
                                     output_dir = "results/hvg") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  top <- head(Seurat::VariableFeatures(seurat_obj), top_n)
  p <- Seurat::VariableFeaturePlot(seurat_obj)
  p <- Seurat::LabelPoints(plot = p, points = top, repel = TRUE)
  out_png <- file.path(output_dir, "hvg_selection.png")
  ggplot2::ggsave(out_png, p, width = 6, height = 4, dpi = 200)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(sub("\\.png$", ".svg", out_png), p, width = 6, height = 4)
  }
  message(sprintf("✓ Saved %s", out_png))
}
