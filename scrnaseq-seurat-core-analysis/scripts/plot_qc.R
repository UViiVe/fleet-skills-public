# plot_qc.R
# ----------------------------------------------------------------------
# QC visualization: violins, scatters, histograms.
#
# Inputs:  Seurat object with QC metrics in meta.data.
# Outputs: PNG figures under `output_dir`.
# Example:
#   source("scripts/plot_qc.R")
#   plot_qc_violin(seurat_obj, batch_col = "batch", output_dir = "results/qc")

.save_gg <- function(p, base_path, w = 6, h = 4) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(base_path, p, width = w, height = h, dpi = 200)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(sub("\\.png$", ".svg", base_path), p, width = w, height = h)
  }
  message(sprintf("✓ Saved %s", base_path))
}

#' Violin plot of QC metrics, optionally split by batch.
plot_qc_violin <- function(seurat_obj,
                            features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                            batch_col = NULL,
                            output_dir = "results/qc") {
  features <- intersect(features, colnames(seurat_obj@meta.data))
  if (length(features) == 0) return(invisible(NULL))
  group_by <- if (!is.null(batch_col) && batch_col %in% colnames(seurat_obj@meta.data)) batch_col else NULL
  p <- Seurat::VlnPlot(seurat_obj, features = features, group.by = group_by,
                       pt.size = 0.1, ncol = length(features))
  .save_gg(p, file.path(output_dir, "qc_violin.png"),
            w = 4 * length(features), h = 4)
}

#' Scatter plot of nCount vs nFeature, optionally with percent.mt.
plot_qc_scatter <- function(seurat_obj, output_dir = "results/qc") {
  p1 <- Seurat::FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA")
  p2 <- Seurat::FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt")
  p <- p1 | p2
  .save_gg(p, file.path(output_dir, "qc_scatter.png"), w = 10, h = 4)
}

#' Histograms of QC metrics.
plot_qc_histogram <- function(seurat_obj,
                                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                output_dir = "results/qc") {
  features <- intersect(features, colnames(seurat_obj@meta.data))
  if (length(features) == 0) return(invisible(NULL))
  meta <- seurat_obj@meta.data
  plots <- lapply(features, function(f) {
    ggplot2::ggplot(meta, ggplot2::aes_string(x = f)) +
      ggplot2::geom_histogram(bins = 80, fill = "steelblue", alpha = 0.7) +
      ggplot2::theme_minimal()
  })
  p <- patchwork::wrap_plots(plots, ncol = length(features))
  .save_gg(p, file.path(output_dir, "qc_histograms.png"),
            w = 4 * length(features), h = 3)
}

#' Top-expressed genes (a quick check for ambient soup contaminants).
plot_highest_expr_genes <- function(seurat_obj,
                                       n_top = 20,
                                       output_dir = "results/qc") {
  counts <- Seurat::GetAssayData(seurat_obj, slot = "counts")
  gene_means <- Matrix::rowSums(counts) / ncol(counts)
  top <- head(sort(gene_means, decreasing = TRUE), n_top)
  df <- data.frame(gene = names(top), mean_counts = as.numeric(top))
  df$gene <- factor(df$gene, levels = rev(df$gene))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = mean_counts, y = gene)) +
    ggplot2::geom_bar(stat = "identity", fill = "steelblue") +
    ggplot2::labs(x = "mean counts per cell", y = "") +
    ggplot2::theme_minimal()
  .save_gg(p, file.path(output_dir, "highest_expr_genes.png"), w = 6, h = 5)
}
