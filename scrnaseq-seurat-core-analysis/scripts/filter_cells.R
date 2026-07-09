# filter_cells.R
# ----------------------------------------------------------------------
# Filter cells by QC thresholds, run doublet detection, and apply MAD filters.
#
# Inputs:  Seurat object with QC metrics in meta.data.
# Outputs: a smaller Seurat object with predicted_doublet flags set.
# Example:
#   source("scripts/filter_cells.R")
#   seurat_obj <- run_doubletfinder(seurat_obj, batch_col = "batch")
#   seurat_obj <- filter_by_mad_outliers(seurat_obj, remove_doublets = TRUE)

#' Estimate expected doublet rate from cell loading.
#' Approx. 0.8% per 1000 cells loaded (10X-published rate).
.estimate_doublet_rate <- function(n_cells) {
  min(0.10, max(0.005, n_cells / 1000 * 0.008))
}

#' Filter by fixed gene-count and MT% thresholds.
filter_cells_by_qc <- function(seurat_obj,
                                min_features = 200,
                                max_features = NULL,
                                max_mt_percent = NULL) {
  n0 <- ncol(seurat_obj)
  cells_keep <- colnames(seurat_obj)[
    seurat_obj$nFeature_RNA >= min_features &
      (is.null(max_features) | seurat_obj$nFeature_RNA <= max_features) &
      (is.null(max_mt_percent) | seurat_obj$percent.mt <= max_mt_percent)
  ]
  seurat_obj <- subset(seurat_obj, cells = cells_keep)
  message(sprintf("✓ filter_cells_by_qc: kept %d / %d (%.1f%%)",
                  ncol(seurat_obj), n0, 100 * ncol(seurat_obj) / n0))
  seurat_obj
}

#' Apply tissue-appropriate thresholds and optionally remove doublets.
filter_cells_by_tissue <- function(seurat_obj, tissue = "pbmc",
                                    remove_doublets = TRUE) {
  source(file.path(dirname(sys.frame(1)$ofile %||% "."), "qc_metrics.R"))
  th <- get_tissue_qc_thresholds(tissue)
  seurat_obj <- filter_cells_by_qc(
    seurat_obj,
    min_features = th$min_features,
    max_features = th$max_features,
    max_mt_percent = th$max_mt
  )
  if (remove_doublets && "predicted_doublet" %in% colnames(seurat_obj@meta.data)) {
    n0 <- ncol(seurat_obj)
    seurat_obj <- subset(seurat_obj, subset = predicted_doublet == FALSE)
    message(sprintf("✓ removed doublets: kept %d / %d", ncol(seurat_obj), n0))
  }
  seurat_obj
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Run DoubletFinder per batch and store scores/predictions in meta.data.
run_doubletfinder <- function(seurat_obj,
                               batch_col = NULL,
                               expected_doublet_rate = NULL) {
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) {
    message("DoubletFinder not installed; skipping doublet detection.")
    return(seurat_obj)
  }

  do_one <- function(obj, rate) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
    n_pcs <- min(30, ncol(obj@reductions$pca))
    # paramSweep + summarizeSweep + find.pK
    sweep_list <- tryCatch(
      DoubletFinder::paramSweep_v3(obj, PCs = 1:n_pcs, sct = FALSE),
      error = function(e) NULL
    )
    if (is.null(sweep_list)) {
      message("  paramSweep failed; skipping doublet detection.")
      obj$predicted_doublet <- FALSE
      obj$doublet_score <- 0
      return(obj)
    }
    sweep_stats <- DoubletFinder::summarizeSweep(sweep_list, GT = FALSE)
    bcmvn <- DoubletFinder::find.pK(sweep_stats)
    pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
    expected_n <- round(rate * ncol(obj))
    pann <- paste0("pANN_0.25_", pK, "_", expected_n)
    obj <- DoubletFinder::doubletFinder_v3(
      obj, PCs = 1:n_pcs, pN = 0.25, pK = pK,
      nExp = expected_n, reuse.pANN = NULL, sct = FALSE
    )
    cls_col <- grep("^DF.classifications", colnames(obj@meta.data), value = TRUE)[1]
    obj$predicted_doublet <- obj@meta.data[[cls_col]] == "Doublet"
    pann_col <- grep("^pANN", colnames(obj@meta.data), value = TRUE)[1]
    obj$doublet_score <- obj@meta.data[[pann_col]]
    obj
  }

  if (is.null(batch_col) || !batch_col %in% colnames(seurat_obj@meta.data)) {
    rate <- expected_doublet_rate %||% .estimate_doublet_rate(ncol(seurat_obj))
    message(sprintf("DoubletFinder (single run, expected rate = %.3f)", rate))
    seurat_obj <- do_one(seurat_obj, rate)
  } else {
    batches <- unique(seurat_obj@meta.data[[batch_col]])
    seurat_obj$predicted_doublet <- FALSE
    seurat_obj$doublet_score <- 0
    for (b in batches) {
      sub <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[batch_col]] == b])
      rate <- expected_doublet_rate %||% .estimate_doublet_rate(ncol(sub))
      message(sprintf("  DoubletFinder [%s, n=%d, rate=%.3f]", b, ncol(sub), rate))
      sub <- do_one(sub, rate)
      cells <- colnames(sub)
      seurat_obj@meta.data[cells, "predicted_doublet"] <- sub$predicted_doublet
      seurat_obj@meta.data[cells, "doublet_score"] <- sub$doublet_score
    }
  }
  n <- sum(seurat_obj$predicted_doublet)
  message(sprintf("✓ DoubletFinder completed successfully: %d doublets (%.1f%%)",
                  n, 100 * n / ncol(seurat_obj)))
  seurat_obj
}

#' Apply the `outlier` flag and optionally remove doublets.
filter_by_mad_outliers <- function(seurat_obj,
                                    remove_doublets = TRUE,
                                    doublet_score_threshold = NULL) {
  if (!"outlier" %in% colnames(seurat_obj@meta.data)) {
    message("  no 'outlier' column; run batch_mad_outlier_detection first.")
    return(seurat_obj)
  }
  n0 <- ncol(seurat_obj)
  keep <- !seurat_obj@meta.data$outlier
  seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep])
  message(sprintf("  removed outliers: %d", n0 - ncol(seurat_obj)))

  if (remove_doublets && "predicted_doublet" %in% colnames(seurat_obj@meta.data)) {
    m0 <- ncol(seurat_obj)
    if (!is.null(doublet_score_threshold)) {
      keep <- seurat_obj$doublet_score < doublet_score_threshold
    } else {
      keep <- !seurat_obj$predicted_doublet
    }
    seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep])
    message(sprintf("  removed doublets: %d", m0 - ncol(seurat_obj)))
  }
  message(sprintf("✓ Cell filtering completed successfully: kept %d / %d (%.1f%%)",
                  ncol(seurat_obj), n0, 100 * ncol(seurat_obj) / n0))
  seurat_obj
}

#' Filter genes by minimum cells expressing them.
filter_genes <- function(seurat_obj, min_cells = 3) {
  counts <- Seurat::GetAssayData(seurat_obj, slot = "counts")
  keep <- Matrix::rowSums(counts > 0) >= min_cells
  seurat_obj <- subset(seurat_obj, features = rownames(seurat_obj)[keep])
  message(sprintf("✓ filter_genes: kept %d genes", sum(keep)))
  seurat_obj
}

#' Compare per-cell QC distributions before vs after filtering.
compare_before_after_filtering <- function(seurat_before, seurat_after,
                                            output_dir = "results/qc") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  out_path <- file.path(output_dir, "qc_before_after.csv")
  df <- data.frame(
    metric = metrics,
    median_before = vapply(metrics, function(m) median(seurat_before@meta.data[[m]]), 0),
    median_after  = vapply(metrics, function(m) median(seurat_after@meta.data[[m]]), 0)
  )
  utils::write.csv(df, out_path, row.names = FALSE)
  message(sprintf("✓ Saved %s", out_path))
  invisible(df)
}
