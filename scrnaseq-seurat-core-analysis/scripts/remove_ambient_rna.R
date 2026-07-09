# remove_ambient_rna.R
# ----------------------------------------------------------------------
# Ambient RNA correction via SoupX.
#
# Inputs:  paths to raw and filtered CellRanger matrix directories.
# Outputs: a Seurat object built from SoupX-corrected counts.
# Example:
#   source("scripts/remove_ambient_rna.R")
#   seurat_obj <- run_soupx_correction("raw_feature_bc_matrix/", "filtered_feature_bc_matrix/")

#' Run SoupX (auto-estimate or manual) and return a Seurat object.
run_soupx_correction <- function(raw_matrix_dir,
                                   filtered_matrix_dir,
                                   contamination_fraction = NULL,
                                   output_dir = "results/ambient") {
  if (!requireNamespace("SoupX", quietly = TRUE)) {
    stop("Install with: install.packages('SoupX')")
  }
  if (!dir.exists(raw_matrix_dir) || !dir.exists(filtered_matrix_dir)) {
    stop("Both `raw_matrix_dir` and `filtered_matrix_dir` must exist.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("SoupX: estimating contamination fraction...")
  sc <- SoupX::load10X(filtered_matrix_dir)
  if (is.null(contamination_fraction)) {
    sc <- SoupX::autoEstCont(sc)
  } else {
    sc <- SoupX::setContaminationFraction(sc, contamination_fraction)
  }
  out <- SoupX::adjustCounts(sc, verbose = 0)

  seurat_obj <- Seurat::CreateSeuratObject(counts = out, project = "soupx")
  saveRDS(seurat_obj, file.path(output_dir, "soupx_corrected.rds"))
  message("✓ SoupX correction completed successfully!")
  seurat_obj
}

#' Quick contamination estimate from before/after total UMIs.
estimate_contamination_fraction <- function(seurat_before, seurat_after) {
  before <- mean(Matrix::colSums(Seurat::GetAssayData(seurat_before, slot = "counts")))
  after  <- mean(Matrix::colSums(Seurat::GetAssayData(seurat_after,  slot = "counts")))
  estimate <- 1 - after / max(before, 1)
  estimate <- max(0, min(0.5, estimate))
  message(sprintf("  contamination estimate: %.1f%%", estimate * 100))
  estimate
}

#' Save a before/after summary CSV plus a histogram.
compare_before_after <- function(seurat_before,
                                   seurat_after,
                                   output_dir = "results/ambient",
                                   markers = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  b <- Matrix::colSums(Seurat::GetAssayData(seurat_before, slot = "counts"))
  a <- Matrix::colSums(Seurat::GetAssayData(seurat_after,  slot = "counts"))
  df <- data.frame(
    metric = c("mean_umi_before", "mean_umi_after", "pct_reduction"),
    value  = c(mean(b), mean(a), 100 * (1 - mean(a) / max(mean(b), 1)))
  )
  utils::write.csv(df, file.path(output_dir, "umi_before_after.csv"),
                    row.names = FALSE)
  message(sprintf("✓ Saved %s", file.path(output_dir, "umi_before_after.csv")))
  invisible(df)
}
