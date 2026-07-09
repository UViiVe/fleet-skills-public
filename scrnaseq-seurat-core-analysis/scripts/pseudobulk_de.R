# pseudobulk_de.R
# ----------------------------------------------------------------------
# Pseudobulk aggregation + DESeq2 differential expression per cell type.
#
# Inputs:  Seurat object with raw counts (RNA assay), `sample_col`,
#          `celltype_col`, and a `condition` covariate.
# Outputs: a named list of DE result data frames (one per cell type),
#          plus volcano / MA plots.
# Example:
#   source("scripts/pseudobulk_de.R")
#   pb <- aggregate_to_pseudobulk(seurat_obj, "sample_id", "cell_type", min_cells = 10)
#   de <- run_pseudobulk_deseq2(pb, formula = "~ batch + condition",
#                                contrast = c("condition", "treated", "control"))

#' Sum raw counts per (sample, cell_type), filtering low-cell combinations.
aggregate_to_pseudobulk <- function(seurat_obj,
                                      sample_col = "sample_id",
                                      celltype_col = "cell_type",
                                      min_cells = 10,
                                      min_total_counts = 1) {
  meta <- seurat_obj@meta.data
  if (!all(c(sample_col, celltype_col) %in% colnames(meta))) {
    stop(sprintf("meta.data must contain `%s` and `%s`", sample_col, celltype_col))
  }
  counts <- Seurat::GetAssayData(seurat_obj, slot = "counts", assay = "RNA")
  samples <- unique(as.character(meta[[sample_col]]))
  cts <- unique(as.character(meta[[celltype_col]]))
  out <- list()
  for (ct in cts) {
    counts_rows <- list()
    meta_rows <- list()
    for (s in samples) {
      cells <- rownames(meta)[
        as.character(meta[[sample_col]]) == s &
          as.character(meta[[celltype_col]]) == ct
      ]
      n <- length(cells)
      if (n < min_cells) next
      sub <- counts[, cells, drop = FALSE]
      sums <- Matrix::rowSums(sub)
      if (sum(sums) < min_total_counts) next
      counts_rows[[s]] <- sums
      m <- list(sample = s, cell_type = ct, n_cells = n)
      for (col in c("condition", "batch", "donor", "sex")) {
        if (col %in% colnames(meta)) {
          v <- unique(meta[cells, col])
          m[[col]] <- if (length(v) == 1) as.character(v) else paste(v, collapse = ",")
        }
      }
      meta_rows[[s]] <- m
    }
    if (length(counts_rows) == 0) next
    counts_df <- do.call(cbind, counts_rows)
    meta_df <- do.call(rbind, lapply(meta_rows, as.data.frame, stringsAsFactors = FALSE))
    rownames(meta_df) <- names(counts_rows)
    out[[ct]] <- list(counts = counts_df, meta = meta_df)
    message(sprintf("  %s: %d samples × %d genes",
                    ct, ncol(counts_df), nrow(counts_df)))
  }
  message(sprintf("✓ Pseudobulk aggregation completed: %d cell types", length(out)))
  out
}

#' Confirm at least `require_min` samples per level of the condition factor.
validate_pseudobulk_design <- function(meta, condition, require_min = 2) {
  if (!condition %in% colnames(meta)) {
    message("  WARNING: condition column `", condition, "` missing; skip")
    return(FALSE)
  }
  counts <- table(meta[[condition]])
  bad <- counts[counts < require_min]
  if (length(bad) > 0) {
    message(sprintf("  REFUSING to run DE — levels with N < %d: %s",
                    require_min,
                    paste(sprintf("%s=%d", names(bad), bad), collapse = ", ")))
    return(FALSE)
  }
  TRUE
}

#' Run DESeq2 per cell type. Skips cell types where N<2 in any condition.
run_pseudobulk_deseq2 <- function(pseudobulk,
                                    formula = "~ condition",
                                    contrast = NULL,
                                    alpha = 0.05,
                                    output_dir = NULL) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Install with: BiocManager::install('DESeq2')")
  }
  if (is.null(contrast) || length(contrast) != 3) {
    stop("contrast must be c(factor, level_a, level_b)")
  }
  factor_name <- contrast[1]
  design <- stats::as.formula(formula)

  out <- list()
  for (ct in names(pseudobulk)) {
    counts <- pseudobulk[[ct]]$counts
    meta <- pseudobulk[[ct]]$meta
    if (!validate_pseudobulk_design(meta, factor_name)) next
    message("-- ", ct, " --")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = round(counts),
      colData = meta,
      design = design
    )
    keep <- rowSums(DESeq2::counts(dds) >= 10) >= 3
    dds <- dds[keep, ]
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = contrast, alpha = alpha)
    df <- as.data.frame(res)
    df$gene <- rownames(df)
    df$cell_type <- ct
    out[[ct]] <- df
    message(sprintf("  %s: %d significant at padj<%g",
                    ct, sum(df$padj < alpha, na.rm = TRUE), alpha))
  }
  if (!is.null(output_dir)) {
    export_de_results(out, output_dir)
  }
  message("✓ DESeq2 analysis completed")
  out
}

#' Write per-cell-type CSVs and a combined long table.
export_de_results <- function(de_results, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  combined <- list()
  for (ct in names(de_results)) {
    safe <- gsub("[/ ]", "_", ct)
    utils::write.csv(de_results[[ct]],
                      file.path(output_dir, sprintf("%s_deseq2_results.csv", safe)),
                      row.names = FALSE)
    combined[[ct]] <- de_results[[ct]]
  }
  if (length(combined) > 0) {
    utils::write.csv(do.call(rbind, combined),
                      file.path(output_dir, "all_de_results.csv"), row.names = FALSE)
  }
  message(sprintf("✓ DE results saved to %s", output_dir))
}

#' Volcano plot for one cell type.
plot_volcano <- function(results,
                          celltype = "",
                          padj_threshold = 0.05,
                          log2fc_threshold = 0.5,
                          output_dir = "results/pseudobulk_de") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  df <- results[!is.na(results$padj) & !is.na(results$log2FoldChange), ]
  df$neg_log10_padj <- -log10(pmax(df$padj, 1e-300))
  df$sig <- df$padj < padj_threshold & abs(df$log2FoldChange) > log2fc_threshold
  p <- ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, y = neg_log10_padj, color = sig)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.6) +
    ggplot2::scale_color_manual(values = c("TRUE" = "crimson", "FALSE" = "grey60")) +
    ggplot2::geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Volcano — ", celltype),
                   x = "log2 fold change", y = "-log10 adjusted p")
  safe <- gsub("[/ ]", "_", celltype)
  out <- file.path(output_dir, sprintf("volcano_%s.png", safe))
  ggplot2::ggsave(out, p, width = 6, height = 5, dpi = 200)
  message(sprintf("✓ Saved %s", out))
}

#' MA plot for one cell type.
plot_ma <- function(results,
                     celltype = "",
                     output_dir = "results/pseudobulk_de") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  df <- results[!is.na(results$baseMean) & !is.na(results$log2FoldChange), ]
  df$sig <- if (!is.null(df$padj)) df$padj < 0.05 else abs(df$log2FoldChange) > 1
  p <- ggplot2::ggplot(df, ggplot2::aes(x = baseMean, y = log2FoldChange, color = sig)) +
    ggplot2::geom_point(alpha = 0.5, size = 0.6) +
    ggplot2::scale_color_manual(values = c("TRUE" = "crimson", "FALSE" = "grey60")) +
    ggplot2::scale_x_log10() +
    ggplot2::theme_minimal() +
    ggplot2::geom_hline(yintercept = 0, color = "black", size = 0.3) +
    ggplot2::labs(title = paste0("MA — ", celltype),
                   x = "mean expression", y = "log2 fold change")
  safe <- gsub("[/ ]", "_", celltype)
  out <- file.path(output_dir, sprintf("ma_%s.png", safe))
  ggplot2::ggsave(out, p, width = 6, height = 4.5, dpi = 200)
  message(sprintf("✓ Saved %s", out))
}
