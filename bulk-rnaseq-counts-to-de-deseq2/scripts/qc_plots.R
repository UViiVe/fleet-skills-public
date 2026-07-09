# qc_plots.R
#
# Purpose:
#   Publication-quality QC plots for a fitted DESeqDataSet:
#   dispersion, PCA, MA, volcano. PNG always, SVG when svglite is
#   available (falls back to base svg() device, then skips with a
#   warning if neither works).
#
# Inputs:
#   dds : DESeqDataSet (after DESeq())
#   res : DESeqResults (optional, needed for MA + volcano)
#
# Outputs:
#   PNG/SVG files in output_dir/: dispersion_plot, pca_plot,
#   ma_plot, volcano_plot.
#
# Example:
#   source("scripts/qc_plots.R")
#   run_all_qc(dds, res, output_dir = "qc")

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2",  quietly = TRUE)) stop("Install with: BiocManager::install('DESeq2')")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install with: install.packages('ggplot2')")
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Install with: install.packages('ggrepel')")
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})
.have_ggprism <- requireNamespace("ggprism", quietly = TRUE)
if (.have_ggprism) suppressPackageStartupMessages(library(ggprism))
.have_svglite <- requireNamespace("svglite", quietly = TRUE)

.qc_theme <- function() {
  if (.have_ggprism) {
    ggprism::theme_prism(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.title = element_text(face = "bold"))
  } else {
    theme_bw(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.title = element_text(face = "bold"))
  }
}

.save_plot <- function(p, base_path, width = 8, height = 6, dpi = 300) {
  base_no_ext <- sub("\\.(png|svg)$", "", base_path)
  png_path <- paste0(base_no_ext, ".png")
  ggsave(png_path, plot = p, width = width, height = height,
         dpi = dpi, device = "png")
  cat(sprintf("    wrote %s\n", png_path))

  svg_path <- paste0(base_no_ext, ".svg")
  ok <- FALSE
  if (.have_svglite) {
    ok <- tryCatch({
      ggsave(svg_path, plot = p, width = width, height = height, device = "svg")
      TRUE
    }, error = function(e) FALSE)
  }
  if (!ok) {
    ok <- tryCatch({
      grDevices::svg(svg_path, width = width, height = height)
      print(p); grDevices::dev.off()
      TRUE
    }, error = function(e) FALSE)
  }
  if (ok) cat(sprintf("    wrote %s\n", svg_path)) else
    cat("    SVG skipped (svglite missing and base svg device failed)\n")
}

#' Dispersion estimates plot.
plot_dispersions <- function(dds, output_file) {
  df <- data.frame(
    baseMean    = rowMeans(counts(dds, normalized = TRUE)),
    dispGeneEst = mcols(dds)$dispGeneEst,
    dispFit     = mcols(dds)$dispFit,
    dispersion  = dispersions(dds),
    dispOutlier = mcols(dds)$dispOutlier
  )
  df <- df[!is.na(df$dispGeneEst) & df$baseMean > 0, ]

  p <- ggplot(df, aes(x = baseMean)) +
    geom_point(aes(y = dispGeneEst),  size = 0.5, alpha = 0.3, colour = "black") +
    geom_point(data = df[isTRUE(any(df$dispOutlier)) & !is.na(df$dispOutlier) & df$dispOutlier, ],
               aes(y = dispGeneEst), size = 0.8, colour = "cyan4") +
    geom_line(aes(y = dispFit), colour = "red", linewidth = 1) +
    geom_point(aes(y = dispersion), size = 0.3, alpha = 0.5, colour = "dodgerblue") +
    scale_x_log10() + scale_y_log10() +
    labs(title = "Gene-wise dispersion",
         x = "Mean of normalized counts",
         y = "Dispersion") +
    .qc_theme()

  .save_plot(p, output_file)
  cat("    interpret: black=gene-wise estimate, red=fit, blue=final shrunk, cyan=outlier\n")
}

#' PCA plot with sample labels.
plot_pca <- function(dds, intgroup = "condition", output_file,
                     use_vst = TRUE, label_samples = TRUE) {
  vsd <- if (use_vst) vst(dds, blind = FALSE) else rlog(dds, blind = FALSE)
  pca <- plotPCA(vsd, intgroup = intgroup, returnData = TRUE)
  pv  <- round(100 * attr(pca, "percentVar"))

  p <- ggplot(pca, aes(x = PC1, y = PC2, colour = group, label = name)) +
    geom_point(size = 4, alpha = 0.8) +
    labs(title = "PCA",
         x = sprintf("PC1: %d%%", pv[1]),
         y = sprintf("PC2: %d%%", pv[2]),
         colour = intgroup) +
    .qc_theme()
  if (label_samples) {
    p <- p + geom_text_repel(size = 3, max.overlaps = 20,
                             box.padding = 0.5, point.padding = 0.3,
                             show.legend = FALSE)
  }
  .save_plot(p, output_file)
  cat("    interpret: samples should cluster by", intgroup, "(not by batch)\n")
}

#' MA plot (log2FC vs mean expression) with top genes labeled.
plot_ma <- function(res, output_file, ylim = c(-5, 5),
                    alpha = 0.05, label_top = 10) {
  df <- as.data.frame(res); df$gene <- rownames(df)
  df$significant <- !is.na(df$padj) & df$padj < alpha
  df$y <- pmin(pmax(df$log2FoldChange, ylim[1]), ylim[2])

  top <- df[df$significant, ]
  if (nrow(top) > 0) {
    top <- top[order(top$padj), ]
    df$label <- ifelse(df$gene %in% head(top$gene, label_top), df$gene, "")
  } else df$label <- ""

  p <- ggplot(df, aes(x = baseMean, y = y)) +
    geom_point(aes(colour = significant), size = 0.5, alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "gray40") +
    scale_x_log10() +
    scale_colour_manual(values = c(`FALSE` = "gray60", `TRUE` = "dodgerblue"),
                        labels = c(`FALSE` = "ns",
                                   `TRUE`  = sprintf("padj < %s", alpha)),
                        name = "") +
    labs(title = "MA plot",
         x = "Mean of normalized counts",
         y = "Log2 fold change") +
    coord_cartesian(ylim = ylim) +
    .qc_theme()
  if (any(df$label != "")) {
    p <- p + geom_text_repel(aes(label = label), size = 2.5,
                             max.overlaps = 20, box.padding = 0.5,
                             point.padding = 0.3, segment.size = 0.2,
                             show.legend = FALSE)
  }
  .save_plot(p, output_file)
  cat("    interpret: cloud should be symmetric about y=0\n")
}

#' Volcano plot with thresholds drawn.
plot_volcano <- function(res, output_file, xlim = c(-6, 6),
                         alpha = 0.05, lfc_threshold = 1, label_top = 10) {
  df <- as.data.frame(res); df$gene <- rownames(df)
  df$neg_log10_padj <- -log10(df$padj)
  df$category <- "ns"
  df$category[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange >  lfc_threshold] <- "up"
  df$category[!is.na(df$padj) & df$padj < alpha & df$log2FoldChange < -lfc_threshold] <- "down"
  df$category <- factor(df$category, levels = c("up", "down", "ns"))
  df$x <- pmin(pmax(df$log2FoldChange, xlim[1]), xlim[2])

  hits <- df[df$category != "ns" & !is.na(df$padj), ]
  if (nrow(hits) > 0) {
    hits <- hits[order(hits$padj), ]
    df$label <- ifelse(df$gene %in% head(hits$gene, label_top), df$gene, "")
  } else df$label <- ""

  df_plot <- df[!is.na(df$neg_log10_padj), ]

  p <- ggplot(df_plot, aes(x = x, y = neg_log10_padj)) +
    geom_point(aes(colour = category), size = 0.5, alpha = 0.5) +
    geom_vline(xintercept = c(-lfc_threshold, lfc_threshold),
               linetype = "dashed", colour = "gray40", linewidth = 0.5) +
    geom_hline(yintercept = -log10(alpha),
               linetype = "dashed", colour = "gray40", linewidth = 0.5) +
    scale_colour_manual(values = c(up = "red3", down = "dodgerblue", ns = "gray60"),
                        name = "") +
    labs(title = "Volcano",
         x = "Log2 fold change",
         y = expression(-log[10]~adjusted~p)) +
    coord_cartesian(xlim = xlim) +
    .qc_theme()
  if (any(df_plot$label != "")) {
    p <- p + geom_text_repel(data = df_plot[df_plot$label != "", ],
                             aes(label = label), size = 2.5,
                             max.overlaps = 20, box.padding = 0.5,
                             point.padding = 0.3, segment.size = 0.2,
                             show.legend = FALSE)
  }
  .save_plot(p, output_file)
  n_up   <- sum(df$category == "up",   na.rm = TRUE)
  n_down <- sum(df$category == "down", na.rm = TRUE)
  cat(sprintf("    %d up, %d down at padj < %s, |LFC| > %s\n",
              n_up, n_down, alpha, lfc_threshold))
}

#' Run all four QC plots and write to output_dir.
run_all_qc <- function(dds, res = NULL, output_dir = "qc") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat("\n=== QC plots ===\n")
  cat(sprintf("genes=%d samples=%d design=%s\n",
              nrow(dds), ncol(dds), deparse(design(dds))))

  cat("- dispersion\n")
  plot_dispersions(dds, file.path(output_dir, "dispersion_plot"))

  cat("- PCA\n")
  intgrp <- if ("condition" %in% colnames(colData(dds))) "condition" else colnames(colData(dds))[1]
  plot_pca(dds, intgroup = intgrp,
           output_file = file.path(output_dir, "pca_plot"),
           use_vst = ncol(dds) > 30)

  if (!is.null(res)) {
    cat("- MA\n")
    plot_ma(res, output_file = file.path(output_dir, "ma_plot"))
    cat("- volcano\n")
    plot_volcano(res, output_file = file.path(output_dir, "volcano_plot"))
  }

  cat(sprintf("[OK] QC plots written to %s/\n", output_dir))
}
