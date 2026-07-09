# integrate_batches.R
# ----------------------------------------------------------------------
# Batch integration via Harmony, Seurat CCA, or Seurat RPCA.
#
# Inputs:  Seurat object with PCA computed and a `batch` column in meta.data.
# Outputs: Seurat object with reductions$harmony (or `integrated` assay).
# Example:
#   source("scripts/integrate_batches.R")
#   seurat_obj <- run_harmony_integration(seurat_obj, batch_var = "batch")

#' Pre-flight check before integration.
setup_for_integration <- function(seurat_obj, batch_var) {
  if (!batch_var %in% colnames(seurat_obj@meta.data)) {
    stop("`", batch_var, "` not in meta.data")
  }
  n_batches <- length(unique(seurat_obj@meta.data[[batch_var]]))
  if (n_batches < 2) {
    warning("only 1 batch detected — integration is a no-op")
  }
  message(sprintf("  %d batches", n_batches))
  invisible(seurat_obj)
}

#' Harmony integration on the existing PCA reduction.
run_harmony_integration <- function(seurat_obj,
                                      batch_var = "batch",
                                      dims_use = 1:30,
                                      theta = 2,
                                      lambda = 1,
                                      max_iter = 10) {
  if (!requireNamespace("harmony", quietly = TRUE)) {
    stop("Install with: install.packages('harmony')")
  }
  setup_for_integration(seurat_obj, batch_var)
  message(sprintf("Harmony: theta=%g, dims=%d, max_iter=%d",
                  theta, length(dims_use), max_iter))
  seurat_obj <- harmony::RunHarmony(
    seurat_obj,
    group.by.vars = batch_var,
    dims.use = dims_use,
    theta = theta,
    lambda = lambda,
    max.iter.harmony = max_iter,
    verbose = FALSE
  )
  seurat_obj@misc$integration <- list(
    method = "Harmony", batch_var = batch_var,
    theta = theta, dims = dims_use, lambda = lambda
  )
  message("✓ Harmony integration completed successfully")
  seurat_obj
}

#' Seurat CCA integration (split + IntegrateData).
run_seurat_cca_integration <- function(seurat_obj,
                                          split_by = "batch",
                                          dims = 1:30,
                                          n_features = 2000) {
  setup_for_integration(seurat_obj, split_by)
  message(sprintf("Seurat CCA: dims=%d, n_features=%d",
                  length(dims), n_features))
  obj_list <- Seurat::SplitObject(seurat_obj, split.by = split_by)
  obj_list <- lapply(obj_list, function(x) {
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, nfeatures = n_features, verbose = FALSE)
    x
  })
  features <- Seurat::SelectIntegrationFeatures(obj_list, nfeatures = n_features)
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = obj_list,
    anchor.features = features,
    reduction = "cca",
    dims = dims, verbose = FALSE
  )
  integrated <- Seurat::IntegrateData(anchorset = anchors, dims = dims,
                                       verbose = FALSE)
  Seurat::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
  integrated <- Seurat::RunPCA(integrated, npcs = 50, verbose = FALSE)
  integrated@misc$integration <- list(
    method = "CCA", split_by = split_by, dims = dims
  )
  message("✓ Seurat CCA integration completed successfully")
  integrated
}

#' Seurat RPCA integration (faster than CCA; assumes similar compositions).
run_seurat_rpca_integration <- function(seurat_obj,
                                          split_by = "batch",
                                          dims = 1:30,
                                          n_features = 2000) {
  setup_for_integration(seurat_obj, split_by)
  message(sprintf("Seurat RPCA: dims=%d, n_features=%d",
                  length(dims), n_features))
  obj_list <- Seurat::SplitObject(seurat_obj, split.by = split_by)
  obj_list <- lapply(obj_list, function(x) {
    x <- Seurat::NormalizeData(x, verbose = FALSE)
    x <- Seurat::FindVariableFeatures(x, nfeatures = n_features, verbose = FALSE)
    x <- Seurat::ScaleData(x, verbose = FALSE)
    x <- Seurat::RunPCA(x, npcs = 30, verbose = FALSE)
    x
  })
  features <- Seurat::SelectIntegrationFeatures(obj_list, nfeatures = n_features)
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = obj_list,
    anchor.features = features,
    reduction = "rpca",
    dims = dims, verbose = FALSE
  )
  integrated <- Seurat::IntegrateData(anchorset = anchors, dims = dims,
                                       verbose = FALSE)
  Seurat::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::ScaleData(integrated, verbose = FALSE)
  integrated <- Seurat::RunPCA(integrated, npcs = 50, verbose = FALSE)
  integrated@misc$integration <- list(
    method = "RPCA", split_by = split_by, dims = dims
  )
  message("✓ Seurat RPCA integration completed successfully")
  integrated
}
