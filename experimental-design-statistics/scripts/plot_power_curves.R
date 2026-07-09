# plot_power_curves.R
# Publication-quality power visualizations using ggplot2 + ggprism.
# Always writes PNG; tries SVG (graceful fallback when svglite isn't installed).
#
# Inputs:
#   - cv, fold_changes, sample sizes, etc.
# Outputs:
#   - PNG and (if possible) SVG plots.

.save_plot <- function(plot, base_path, width = 8, height = 6, dpi = 300) {
  png_path <- sub("\\.(svg|png)$", ".png", base_path)
  if (!grepl("\\.png$", png_path)) png_path <- paste0(png_path, ".png")
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  dpi = dpi, device = "png")
  cat("   saved:", png_path, "\n")

  svg_path <- sub("\\.(svg|png)$", ".svg", base_path)
  if (!grepl("\\.svg$", svg_path)) svg_path <- paste0(svg_path, ".svg")
  ok <- tryCatch({
    ggplot2::ggsave(svg_path, plot = plot, width = width, height = height, device = "svg")
    TRUE
  }, error = function(e) FALSE)
  if (!ok) {
    ok2 <- tryCatch({
      grDevices::svg(svg_path, width = width, height = height)
      print(plot); grDevices::dev.off(); TRUE
    }, error = function(e) FALSE)
    if (!ok2) cat("   (SVG export failed; PNG written)\n")
    else      cat("   saved:", svg_path, " (via base svg)\n")
  } else {
    cat("   saved:", svg_path, "\n")
  }
}

.require_plot_pkgs <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggprism", quietly = TRUE))
    stop("Install ggplot2 and ggprism for plotting.")
  if (!requireNamespace("RNASeqPower", quietly = TRUE))
    stop("Install RNASeqPower with BiocManager::install('RNASeqPower')")
}

.palette <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A")

#' Power vs sample size for multiple fold-changes
#'
#' @param cv Coefficient of variation
#' @param fold_changes Vector of fold-changes to plot
#' @param depth Sequencing depth (M reads)
#' @param alpha Significance threshold (default 0.05)
#' @param n_range Range of n per group to evaluate (default 2:15)
#' @param output_file Base output path (PNG + SVG written)
plot_power_vs_samplesize <- function(cv = 0.4,
                                     fold_changes = c(1.5, 2, 3, 4),
                                     depth = 20, alpha = 0.05,
                                     n_range = 2:15,
                                     output_file = "power_vs_samplesize") {
  .require_plot_pkgs()

  df <- do.call(rbind, lapply(fold_changes, function(fc) {
    pw <- vapply(n_range, function(n) {
      RNASeqPower::rnapower(depth = depth, n = n, cv = cv,
                            effect = fc, alpha = alpha)
    }, numeric(1))
    data.frame(n_per_group = n_range, power = pw,
               fold_change = paste0(fc, "x"),
               fc_numeric = fc)
  }))
  df$fold_change <- factor(df$fold_change,
                           levels = paste0(sort(fold_changes), "x"))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = n_per_group, y = power,
                                        color = fold_change)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40") +
    ggplot2::annotate("text", x = max(n_range) - 1, y = 0.82,
                      label = "0.80 target", color = "gray40", size = 3) +
    ggprism::theme_prism() +
    ggplot2::scale_color_manual(values = .palette[seq_along(fold_changes)]) +
    ggplot2::labs(title = "Statistical power vs sample size",
                  subtitle = sprintf("CV = %s, depth = %sM reads, alpha = %s",
                                     cv, depth, alpha),
                  x = "Replicates per group", y = "Power",
                  color = "Fold change") +
    ggplot2::scale_x_continuous(breaks = seq(min(n_range), max(n_range), 2)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   legend.position = "right")

  cat("Saving power-vs-n plot:\n"); .save_plot(p, output_file)
  cat(sprintf("✓ Power-vs-n plot completed: %d combinations\n", nrow(df)))
  invisible(p)
}

#' Power vs fold-change for multiple sample sizes
plot_power_vs_foldchange <- function(n_per_group = c(3, 5, 8, 10),
                                     cv = 0.4, depth = 20, alpha = 0.05,
                                     output_file = "power_vs_foldchange") {
  .require_plot_pkgs()

  fc_range <- seq(1.25, 4, by = 0.25)
  df <- do.call(rbind, lapply(n_per_group, function(n) {
    pw <- vapply(fc_range, function(fc) {
      RNASeqPower::rnapower(depth = depth, n = n, cv = cv,
                            effect = fc, alpha = alpha)
    }, numeric(1))
    data.frame(fold_change = fc_range, power = pw,
               n_label = paste0("n=", n), n_numeric = n)
  }))
  df$n_label <- factor(df$n_label, levels = paste0("n=", sort(n_per_group)))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = fold_change, y = power, color = n_label)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40") +
    ggprism::theme_prism() +
    ggplot2::scale_color_manual(values = .palette[seq_along(n_per_group)]) +
    ggplot2::labs(title = "Statistical power vs fold-change",
                  subtitle = sprintf("CV = %s, depth = %sM reads, alpha = %s",
                                     cv, depth, alpha),
                  x = "Fold change", y = "Power", color = "Sample size") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   legend.position = "right")

  cat("Saving power-vs-FC plot:\n"); .save_plot(p, output_file)
  cat(sprintf("✓ Power-vs-FC plot completed: %d combinations\n", nrow(df)))
  invisible(p)
}

#' Plot CV by tissue from the bundled database (or default categories)
plot_cv_by_tissue <- function(cv_database_path = "references/cv_tissue_database.csv",
                              output_file = "cv_by_tissue") {
  .require_plot_pkgs()

  if (!file.exists(cv_database_path)) {
    df <- data.frame(
      Tissue = c("Cell lines","Inbred mice","Primary cells","Human samples"),
      CV_min = c(0.10, 0.20, 0.30, 0.30),
      CV_max = c(0.20, 0.30, 0.40, 0.50),
      CV_typical = c(0.15, 0.25, 0.35, 0.40)
    )
  } else {
    raw <- utils::read.csv(cv_database_path, stringsAsFactors = FALSE)
    # Aggregate to tissue-level summary
    grouped <- by(raw$cv_median, raw$tissue, function(x) {
      c(min = min(x), max = max(x), typical = median(x))
    })
    rows <- do.call(rbind, lapply(grouped, as.numeric))
    df <- data.frame(Tissue = names(grouped),
                     CV_min = rows[, 1L],
                     CV_max = rows[, 2L],
                     CV_typical = rows[, 3L],
                     stringsAsFactors = FALSE)
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Tissue, y = CV_typical)) +
    ggplot2::geom_point(size = 4, color = "#1B9E77") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = CV_min, ymax = CV_max),
                           width = 0.2, linewidth = 1, color = "#1B9E77") +
    ggprism::theme_prism() +
    ggplot2::labs(title = "Coefficient of variation by tissue",
                  subtitle = "Typical CV ranges for RNA-seq",
                  x = "Tissue", y = "CV",
                  caption = "Point: median CV; bars: min-max") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  cat("Saving CV-by-tissue plot:\n"); .save_plot(p, output_file)
  cat(sprintf("✓ CV plot completed: %d tissues\n", nrow(df)))
  invisible(p)
}

#' Heatmap of required sample size by (CV, fold-change)
plot_sample_size_heatmap <- function(cvs = c(0.2, 0.3, 0.4, 0.5),
                                     fcs = c(1.5, 2, 3, 4),
                                     depth = 20, alpha = 0.05,
                                     target_power = 0.8,
                                     output_file = "sample_size_heatmap") {
  .require_plot_pkgs()

  grid <- expand.grid(CV = cvs, FoldChange = fcs)
  grid$SampleSize <- apply(grid, 1, function(row) {
    ceiling(RNASeqPower::rnapower(depth = depth, cv = row[["CV"]],
                                  effect = row[["FoldChange"]],
                                  alpha = alpha, power = target_power))
  })

  p <- ggplot2::ggplot(grid, ggplot2::aes(x = factor(FoldChange),
                                          y = factor(CV),
                                          fill = SampleSize)) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = SampleSize),
                       size = 5, fontface = "bold") +
    ggplot2::scale_fill_gradient(low = "#1B9E77", high = "#D95F02",
                                 name = "n/group") +
    ggprism::theme_prism() +
    ggplot2::labs(title = sprintf("Required n for %.0f%% power", 100 * target_power),
                  subtitle = sprintf("Depth = %sM reads, alpha = %s", depth, alpha),
                  x = "Fold change to detect", y = "CV") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   plot.subtitle = ggplot2::element_text(hjust = 0.5))

  cat("Saving sample-size heatmap:\n"); .save_plot(p, output_file)
  cat(sprintf("✓ Sample-size heatmap completed: %d cells\n", nrow(grid)))
  invisible(p)
}
