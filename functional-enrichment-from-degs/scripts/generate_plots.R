# generate_plots.R
#
# Purpose:
#   Standard enrichment visualization set:
#     - GSEA dotplot (NES and p.adjust across top pathways)
#     - GSEA running-enrichment plot for the top hit
#     - ORA barplot (up- and downregulated)
#   PNG always; SVG via svglite if installed, else base svg(), else skipped.
#
# Inputs:
#   gsea       : clusterProfiler gseaResult (optional)
#   ora_up     : enrichResult for upregulated genes (optional)
#   ora_down   : enrichResult for downregulated genes (optional)
#   output_dir : directory for output files (default "enrich_out/figures")
#   top_n      : number of pathways to display (default 15)
#
# Outputs:
#   PNG / SVG files in output_dir/.
#
# Example:
#   source("scripts/generate_plots.R")
#   generate_all_plots(gsea, ora_up = ora_up, ora_down = ora_down,
#                      output_dir = "figs")

.have_svglite <- requireNamespace("svglite", quietly = TRUE)

.save_plot <- function(p, base_path, width = 8, height = 6, dpi = 300) {
  if (is.null(p)) return(invisible(NULL))
  base_no_ext <- sub("\\.(png|svg)$", "", base_path)
  png_path <- paste0(base_no_ext, ".png")
  ggplot2::ggsave(png_path, plot = p,
                  width = width, height = height, dpi = dpi, device = "png")
  cat(sprintf("   wrote %s\n", png_path))

  svg_path <- paste0(base_no_ext, ".svg")
  ok <- FALSE
  if (.have_svglite) {
    ok <- tryCatch({
      ggplot2::ggsave(svg_path, plot = p,
                      width = width, height = height, device = "svg")
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
  if (ok) cat(sprintf("   wrote %s\n", svg_path))
  else    cat("   SVG skipped\n")
}

#' GSEA dotplot of top pathways by p.adjust.
plot_gsea_dotplot <- function(gsea, output_file, top_n = 15) {
  if (!requireNamespace("enrichplot", quietly = TRUE))
    stop("Install with: BiocManager::install('enrichplot')")
  if (is.null(gsea) || nrow(as.data.frame(gsea)) == 0) {
    cat("   no GSEA results to plot\n"); return(invisible(NULL))
  }
  suppressPackageStartupMessages(library(enrichplot))
  p <- enrichplot::dotplot(gsea, showCategory = top_n, split = ".sign") +
    ggplot2::facet_grid(. ~ .sign) +
    ggplot2::ggtitle("GSEA: top pathways")
  .save_plot(p, output_file, width = 10, height = 8)
}

#' GSEA running-enrichment plot for the top pathway.
plot_gsea_running_score <- function(gsea, output_file) {
  if (!requireNamespace("enrichplot", quietly = TRUE)) return(invisible(NULL))
  if (is.null(gsea) || nrow(as.data.frame(gsea)) == 0) return(invisible(NULL))
  suppressPackageStartupMessages(library(enrichplot))

  df  <- as.data.frame(gsea)
  top <- df$ID[order(df$p.adjust)][1]
  p   <- enrichplot::gseaplot2(gsea, geneSetID = top, title = top)
  .save_plot(p, output_file, width = 8, height = 6)
}

#' ORA barplot for one direction.
plot_ora_barplot <- function(ora, output_file, direction = "", top_n = 15) {
  if (is.null(ora) || nrow(as.data.frame(ora)) == 0) {
    cat(sprintf("   no ORA-%s results to plot\n", direction))
    return(invisible(NULL))
  }
  if (!requireNamespace("enrichplot", quietly = TRUE))
    stop("Install with: BiocManager::install('enrichplot')")
  suppressPackageStartupMessages(library(enrichplot))
  p <- enrichplot::barplot(ora, showCategory = top_n) +
    ggplot2::ggtitle(sprintf("ORA: %s", direction))
  .save_plot(p, output_file, width = 9, height = 7)
}

#' Run all standard plots.
generate_all_plots <- function(gsea       = NULL,
                               ora_up     = NULL,
                               ora_down   = NULL,
                               output_dir = "enrich_out/figures",
                               top_n      = 15L) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat(sprintf("=== generating plots in %s/ ===\n", output_dir))
  if (!is.null(gsea)) {
    plot_gsea_dotplot      (gsea, file.path(output_dir, "gsea_dotplot"),       top_n)
    plot_gsea_running_score(gsea, file.path(output_dir, "gsea_running_score"))
  }
  if (!is.null(ora_up))
    plot_ora_barplot(ora_up,   file.path(output_dir, "ora_up_barplot"),
                     direction = "upregulated",   top_n = top_n)
  if (!is.null(ora_down))
    plot_ora_barplot(ora_down, file.path(output_dir, "ora_down_barplot"),
                     direction = "downregulated", top_n = top_n)
  cat(sprintf("[OK] plots written to %s/\n", output_dir))
}
