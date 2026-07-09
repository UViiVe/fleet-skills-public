# normalize_data.R
# ----------------------------------------------------------------------
# Normalize counts via SCTransform (default) or LogNormalize.
#
# Inputs:  Seurat object with raw counts.
# Outputs: Seurat object with SCT assay (default) or log-normalized RNA assay.
# Example:
#   source("scripts/normalize_data.R")
#   seurat_obj <- run_sctransform(seurat_obj, vars_to_regress = "percent.mt")

#' SCTransform normalization. Default for UMI data.
run_sctransform <- function(seurat_obj,
                             vars_to_regress = NULL,
                             vst_flavor = "v2",
                             verbose = FALSE) {
  message("SCTransform: vst.flavor=", vst_flavor)
  if (!is.null(vars_to_regress)) {
    message("  regressing out: ", paste(vars_to_regress, collapse = ", "))
  }
  seurat_obj <- Seurat::SCTransform(
    seurat_obj,
    vars.to.regress = vars_to_regress,
    vst.flavor = vst_flavor,
    verbose = verbose
  )
  message("✓ SCTransform normalization completed")
  seurat_obj
}

#' Classic LogNormalize + variable features + scale.
run_lognormalize <- function(seurat_obj,
                              scale_factor = 1e4,
                              n_features = 2000,
                              vars_to_regress = c("nCount_RNA", "percent.mt")) {
  seurat_obj <- Seurat::NormalizeData(seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = scale_factor, verbose = FALSE
  )
  seurat_obj <- Seurat::FindVariableFeatures(
    seurat_obj, selection.method = "vst",
    nfeatures = n_features, verbose = FALSE
  )
  vars <- intersect(vars_to_regress, colnames(seurat_obj@meta.data))
  seurat_obj <- Seurat::ScaleData(
    seurat_obj, vars.to.regress = vars, verbose = FALSE
  )
  message("✓ LogNormalize completed")
  seurat_obj
}

#' Side-by-side comparison of LogNormalize vs SCT residuals (if both present).
compare_normalization_methods <- function(seurat_obj,
                                            output_dir = "results/normalization") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(output_dir, "normalization_comparison.csv")
  metrics <- data.frame(
    assay = names(seurat_obj@assays),
    n_features = vapply(names(seurat_obj@assays),
                          function(a) nrow(seurat_obj@assays[[a]]), 0)
  )
  utils::write.csv(metrics, out_path, row.names = FALSE)
  message(sprintf("✓ Saved %s", out_path))
  invisible(metrics)
}
