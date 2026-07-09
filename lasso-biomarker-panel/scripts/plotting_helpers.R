# plotting_helpers.R
# Shared helpers for any ad-hoc plotting code that lives outside
# biomarker_plots.R. Used by the example notebooks and one-off variants.

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggprism)
})


#' Apply the workflow's standard ggprism theme to an existing plot.
#'
#' @param plot ggplot object.
#' @param base_size Base font size.
#' @return ggplot object with theme applied.
panel_theme <- function(plot, base_size = 12) {
    plot + ggprism::theme_prism(base_size = base_size) +
        theme(plot.title = element_text(hjust = 0.5,
                                         face = "bold", size = base_size + 2))
}


#' Save a ggplot in both PNG and SVG (with graceful SVG fallback).
#'
#' @param plot ggplot object.
#' @param base_path Path without extension.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param dpi Resolution for PNG.
save_paired <- function(plot, base_path,
                        width = 8, height = 6, dpi = 300) {
    png_path <- paste0(base_path, ".png")
    svg_path <- paste0(base_path, ".svg")
    ggsave(png_path, plot = plot, width = width, height = height,
           dpi = dpi, device = "png")
    tryCatch(ggsave(svg_path, plot = plot, width = width,
                     height = height, device = "svg"),
             error = function(e) {
                 tryCatch({
                     grDevices::svg(svg_path, width = width, height = height)
                     print(plot); grDevices::dev.off()
                 }, error = function(e2) NULL)
             })
    invisible(NULL)
}


#' Centered factor palette generator.
#'
#' Returns a named vector of colours for the levels of `x`. Hand-picked first
#' eight; falls back to viridis for longer factors.
#'
#' @param x A factor or vector of unique values.
#' @return Named character vector mapping each level to a hex colour.
palette_for <- function(x) {
    levels_ <- if (is.factor(x)) levels(x) else sort(unique(x))
    pool <- c("#E74C3C", "#3498DB", "#27AE60", "#F39C12",
              "#9B59B6", "#1ABC9C", "#E67E22", "#34495E")
    if (length(levels_) <= length(pool)) {
        setNames(pool[seq_along(levels_)], levels_)
    } else {
        viridis_pool <- viridisLite::viridis(length(levels_))
        setNames(viridis_pool, levels_)
    }
}
