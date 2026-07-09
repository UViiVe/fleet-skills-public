# qc_metrics.R
# ----------------------------------------------------------------------
# Compute QC metrics and flag outliers (fixed thresholds or batch-aware MAD).
#
# Inputs:  Seurat object with raw counts.
# Outputs: Seurat object with QC metrics in meta.data and outlier flags.
# Example:
#   source("scripts/qc_metrics.R")
#   seurat_obj <- calculate_qc_metrics(seurat_obj, species = "human")
#   seurat_obj <- batch_mad_outlier_detection(seurat_obj, batch_col = "batch")

#' Mitochondrial gene pattern for `species`.
get_species_mito_pattern <- function(species) {
  patterns <- list(human = "^MT-", mouse = "^mt-")
  species <- tolower(species)
  if (!species %in% names(patterns)) stop("species must be 'human' or 'mouse'")
  patterns[[species]]
}

#' Ribosomal gene pattern for `species`.
get_species_ribo_pattern <- function(species) {
  patterns <- list(human = "^RP[SL]", mouse = "^Rp[sl]")
  species <- tolower(species)
  if (!species %in% names(patterns)) stop("species must be 'human' or 'mouse'")
  patterns[[species]]
}

#' Add `percent.mt`, `percent.ribo`, and log-transformed counts metrics.
calculate_qc_metrics <- function(seurat_obj,
                                  species = "human",
                                  calculate_ribo = TRUE) {
  mito_pattern <- get_species_mito_pattern(species)
  seurat_obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = mito_pattern)
  message(sprintf("  Mitochondrial genes: %d",
                  sum(grepl(mito_pattern, rownames(seurat_obj)))))

  if (calculate_ribo) {
    ribo_pattern <- get_species_ribo_pattern(species)
    seurat_obj[["percent.ribo"]] <- Seurat::PercentageFeatureSet(seurat_obj, pattern = ribo_pattern)
    message(sprintf("  Ribosomal genes:    %d",
                    sum(grepl(ribo_pattern, rownames(seurat_obj)))))
  }

  seurat_obj$log10_nCount_RNA <- log10(seurat_obj$nCount_RNA + 1)
  seurat_obj$log10_nFeature_RNA <- log10(seurat_obj$nFeature_RNA + 1)

  message(sprintf("✓ QC metrics calculated successfully: %d cells", ncol(seurat_obj)))
  message(sprintf("  median nFeature_RNA: %.0f", median(seurat_obj$nFeature_RNA)))
  message(sprintf("  median nCount_RNA:   %.0f", median(seurat_obj$nCount_RNA)))
  message(sprintf("  median percent.mt:   %.2f%%", median(seurat_obj$percent.mt)))
  seurat_obj
}

#' Tissue-specific recommended thresholds.
get_tissue_qc_thresholds <- function(tissue) {
  thresholds <- list(
    pbmc     = list(min_features = 200, max_features = 2500, max_mt = 5),
    brain    = list(min_features = 200, max_features = 6000, max_mt = 10),
    tumor    = list(min_features = 200, max_features = 5000, max_mt = 20),
    kidney   = list(min_features = 200, max_features = 4000, max_mt = 15),
    liver    = list(min_features = 200, max_features = 4000, max_mt = 15),
    heart    = list(min_features = 200, max_features = 4000, max_mt = 18),
    lung     = list(min_features = 200, max_features = 4000, max_mt = 10),
    pancreas = list(min_features = 200, max_features = 5000, max_mt = 10),
    default  = list(min_features = 200, max_features = 4000, max_mt = 10)
  )
  tissue <- tolower(tissue)
  if (!tissue %in% names(thresholds)) {
    message("Tissue '", tissue, "' not recognized; using default thresholds.")
    tissue <- "default"
  }
  th <- thresholds[[tissue]]
  message(sprintf("Tissue: %s — min_features=%d, max_features=%d, max_mt=%g%%",
                  tissue, th$min_features, th$max_features, th$max_mt))
  th
}

#' Batch-aware MAD outlier detection on QC metrics.
batch_mad_outlier_detection <- function(seurat_obj,
                                         batch_col = "batch",
                                         metrics = c("log10_nCount_RNA",
                                                     "log10_nFeature_RNA",
                                                     "percent.mt"),
                                         nmads = 5) {
  meta <- seurat_obj@meta.data
  if (!batch_col %in% colnames(meta)) {
    message("  No '", batch_col, "' column; treating as single batch.")
    meta[[batch_col]] <- "batch_1"
  }

  meta$outlier <- FALSE
  for (m in metrics) {
    if (!m %in% colnames(meta)) {
      message("  skipping missing metric ", m)
      next
    }
    out_col <- paste0("outlier_", m)
    meta[[out_col]] <- FALSE
    for (b in unique(meta[[batch_col]])) {
      idx <- meta[[batch_col]] == b
      v <- meta[idx, m]
      med <- stats::median(v)
      mad_v <- stats::median(abs(v - med))
      if (mad_v == 0) next
      lo <- med - nmads * mad_v
      hi <- med + nmads * mad_v
      bad <- (v < lo) | (v > hi)
      meta[idx, out_col] <- bad
      meta[idx, "outlier"] <- meta[idx, "outlier"] | bad
      if (any(bad)) {
        message(sprintf("  %s [%s]: %d outliers (%.2f, %.2f)",
                        m, b, sum(bad), lo, hi))
      }
    }
  }
  n_out <- sum(meta$outlier)
  pct <- 100 * n_out / nrow(meta)
  message(sprintf("✓ MAD outlier detection completed: %d outliers (%.1f%%)",
                  n_out, pct))
  seurat_obj@meta.data <- meta
  seurat_obj
}

#' Fixed-threshold outlier flagging using tissue presets (or overrides).
mark_outliers_fixed <- function(seurat_obj,
                                 tissue = "pbmc",
                                 min_features = NULL,
                                 max_features = NULL,
                                 max_mt = NULL) {
  th <- get_tissue_qc_thresholds(tissue)
  if (is.null(min_features)) min_features <- th$min_features
  if (is.null(max_features)) max_features <- th$max_features
  if (is.null(max_mt))       max_mt       <- th$max_mt

  meta <- seurat_obj@meta.data
  meta$outlier_low_features  <- meta$nFeature_RNA < min_features
  meta$outlier_high_features <- meta$nFeature_RNA > max_features
  meta$outlier_high_mt       <- meta$percent.mt   > max_mt
  meta$outlier <- meta$outlier_low_features | meta$outlier_high_features | meta$outlier_high_mt
  seurat_obj@meta.data <- meta

  n_out <- sum(meta$outlier)
  pct <- 100 * n_out / nrow(meta)
  message(sprintf("✓ Fixed-threshold outliers: %d (%.1f%%) — keeping %d",
                  n_out, pct, nrow(meta) - n_out))
  seurat_obj
}
