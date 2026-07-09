# integration_diagnostics.R
# ----------------------------------------------------------------------
# Compute LISI / ASW metrics for integration quality and produce comparison
# plots between un-integrated PCA and integrated reductions.
#
# Inputs:  Seurat object with one or more reductions (pca, harmony, ...).
# Outputs: a named list of scores plus diagnostic plots.
# Example:
#   source("scripts/integration_diagnostics.R")
#   lisi <- compute_lisi_scores(seurat_obj, batch_var = "batch", reduction = "harmony")

#' Compute per-cell LISI for a categorical column.
.lisi_one <- function(X, labels, perplexity = 30) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Install with: install.packages('RANN')")
  }
  k <- max(perplexity * 3, 30)
  k <- min(k, length(labels) - 1)
  nn <- RANN::nn2(X, X, k = k + 1)
  idx <- nn$nn.idx[, -1, drop = FALSE]
  dists <- nn$nn.dists[, -1, drop = FALSE]

  sigma <- pmax(rowMeans(dists), 1e-8)
  w <- exp(-(dists ^ 2) / (2 * sigma ^ 2))
  w <- w / pmax(rowSums(w), 1e-8)

  unique_labels <- unique(labels)
  lisi <- numeric(nrow(X))
  for (i in seq_len(nrow(X))) {
    neigh <- labels[idx[i, ]]
    p_sq <- 0
    for (lbl in unique_labels) {
      p <- sum(w[i, neigh == lbl])
      p_sq <- p_sq + p ^ 2
    }
    lisi[i] <- 1 / max(p_sq, 1e-8)
  }
  lisi
}

#' Compute iLISI (batch) and cLISI (cell type) on a reduction.
compute_lisi_scores <- function(seurat_obj,
                                  batch_var = "batch",
                                  label_var = NULL,
                                  reduction = "harmony",
                                  perplexity = 30) {
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("reduction `", reduction, "` not found in seurat_obj@reductions")
  }
  emb <- seurat_obj@reductions[[reduction]]@cell.embeddings
  out <- list()

  if (batch_var %in% colnames(seurat_obj@meta.data)) {
    labs <- as.character(seurat_obj@meta.data[[batch_var]])
    val <- mean(.lisi_one(emb, labs, perplexity = perplexity))
    out$batch_LISI <- val
    message(sprintf("  iLISI (%s, %s): %.2f", batch_var, reduction, val))
  }
  if (!is.null(label_var) && label_var %in% colnames(seurat_obj@meta.data)) {
    labs <- as.character(seurat_obj@meta.data[[label_var]])
    val <- mean(.lisi_one(emb, labs, perplexity = perplexity))
    out$celltype_LISI <- val
    message(sprintf("  cLISI (%s, %s): %.2f", label_var, reduction, val))
  }
  out
}

#' Compute batch ASW and cell-type ASW silhouette scores.
compute_asw_scores <- function(seurat_obj,
                                  batch_var = "batch",
                                  label_var = NULL,
                                  reduction = "harmony") {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("Install with: install.packages('cluster')")
  }
  emb <- seurat_obj@reductions[[reduction]]@cell.embeddings
  out <- list()
  dist_mat <- stats::dist(emb)

  if (batch_var %in% colnames(seurat_obj@meta.data)) {
    grp <- as.integer(as.factor(seurat_obj@meta.data[[batch_var]]))
    if (length(unique(grp)) > 1) {
      s <- cluster::silhouette(grp, dist_mat)
      out$batch_ASW <- mean(s[, "sil_width"])
      message(sprintf("  batch ASW: %.3f (target ≈ 0)", out$batch_ASW))
    }
  }
  if (!is.null(label_var) && label_var %in% colnames(seurat_obj@meta.data)) {
    grp <- as.integer(as.factor(seurat_obj@meta.data[[label_var]]))
    if (length(unique(grp)) > 1) {
      s <- cluster::silhouette(grp, dist_mat)
      out$celltype_ASW <- mean(s[, "sil_width"])
      message(sprintf("  cell type ASW: %.3f (target > 0.5)", out$celltype_ASW))
    }
  }
  out
}

#' Save a bar plot of integration metrics.
plot_integration_metrics <- function(seurat_obj,
                                       lisi = list(),
                                       asw = list(),
                                       output_dir = "results/integration") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  df <- data.frame(
    metric = c(names(lisi), names(asw)),
    score  = c(unlist(lisi), unlist(asw)),
    stringsAsFactors = FALSE
  )
  if (nrow(df) == 0) return(invisible(NULL))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = metric, y = score)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)) +
    ggplot2::labs(y = "score", title = "Integration metrics")
  out_png <- file.path(output_dir, "integration_metrics.png")
  ggplot2::ggsave(out_png, p, width = 6, height = 4, dpi = 200)
  message(sprintf("✓ Saved %s", out_png))
}

#' Compute LISI + ASW across multiple reductions for comparison.
compare_integration_methods <- function(seurat_obj,
                                          batch_var = "batch",
                                          label_var = NULL,
                                          reductions = c("pca", "harmony"),
                                          output_dir = "results/integration") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  rows <- list()
  for (red in reductions) {
    if (!red %in% names(seurat_obj@reductions)) {
      message("  skipping missing reduction: ", red)
      next
    }
    message("-- ", red, " --")
    lisi <- compute_lisi_scores(seurat_obj, batch_var = batch_var,
                                 label_var = label_var, reduction = red)
    asw <- compute_asw_scores(seurat_obj, batch_var = batch_var,
                                label_var = label_var, reduction = red)
    rows[[red]] <- c(list(reduction = red), lisi, asw)
  }
  df <- do.call(rbind, lapply(rows, function(r) {
    data.frame(r, stringsAsFactors = FALSE)
  }))
  out_csv <- file.path(output_dir, "integration_comparison.csv")
  utils::write.csv(df, out_csv, row.names = FALSE)
  message(sprintf("✓ Saved %s", out_csv))
  df
}
