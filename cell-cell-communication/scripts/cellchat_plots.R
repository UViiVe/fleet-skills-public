# =============================================================================
# Cell-cell communication: visualizations
# -----------------------------------------------------------------------------
# Emits six core figures from a trained CellChat object. Each figure is
# written as PNG and SVG; SVG falls back from `svglite` (when available) to
# the base R `svg()` device automatically.
#
# Inputs : CellChat object from scripts/run_cellchat.R; output directory
# Outputs: PNG + SVG files in `output_dir`
# Example: source("scripts/cellchat_plots.R")
#          generate_all_plots(cellchat, output_dir = "results")
# =============================================================================

suppressPackageStartupMessages({
    library(CellChat)
    library(ggplot2)
    library(ggprism)
})

.has_svglite <- requireNamespace("svglite", quietly = TRUE)
if (.has_svglite) suppressPackageStartupMessages(library(svglite))


# ---- save helpers ------------------------------------------------------------

.save_ggplot <- function(plot_obj, base_path,
                         width = 8, height = 6, dpi = 300) {
    png_path <- sub("\\.(svg|png)$", ".png", base_path)
    ggsave(png_path, plot = plot_obj, width = width, height = height,
           dpi = dpi, device = "png")
    cat("   Saved:", png_path, "\n")

    svg_path <- sub("\\.(svg|png)$", ".svg", base_path)
    tryCatch({
        ggsave(svg_path, plot = plot_obj, width = width, height = height,
               device = "svg")
        cat("   Saved:", svg_path, "\n")
    }, error = function(e) {
        tryCatch({
            grDevices::svg(svg_path, width = width, height = height)
            print(plot_obj)
            grDevices::dev.off()
            cat("   Saved:", svg_path, "\n")
        }, error = function(e2) {
            cat("   (SVG export failed for ggplot)\n")
        })
    })
}

.save_base_plot <- function(draw_fn, base_path,
                            width = 8, height = 6, dpi = 300) {
    png_path <- sub("\\.(svg|png)$", ".png", base_path)
    grDevices::png(png_path, width = width, height = height,
                   units = "in", res = dpi)
    tryCatch(draw_fn(),
             error = function(e) cat("   PNG plot error:",
                                     conditionMessage(e), "\n"))
    grDevices::dev.off()
    cat("   Saved:", png_path, "\n")

    svg_path <- sub("\\.(svg|png)$", ".svg", base_path)
    tryCatch({
        grDevices::svg(svg_path, width = width, height = height)
        draw_fn()
        grDevices::dev.off()
        cat("   Saved:", svg_path, "\n")
    }, error = function(e) {
        tryCatch(grDevices::dev.off(), error = function(e2) NULL)
        cat("   (SVG export failed for base plot)\n")
    })
}


# ---- individual plots --------------------------------------------------------

.plot_interaction_count <- function(cellchat, output_dir) {
    cat("\n   [1/6] interaction count network\n")
    sizes <- as.numeric(table(cellchat@idents))
    mat   <- cellchat@net$count
    .save_base_plot(
        function() {
            par(mar = c(1, 1, 2, 1))
            netVisual_circle(mat,
                             vertex.weight = sizes,
                             weight.scale  = TRUE,
                             label.edge    = FALSE,
                             title.name    = "Number of Interactions")
        },
        file.path(output_dir, "interaction_count_network.png"),
        width = 8, height = 8
    )
}

.plot_interaction_strength <- function(cellchat, output_dir) {
    cat("   [2/6] interaction strength network\n")
    sizes <- as.numeric(table(cellchat@idents))
    mat   <- cellchat@net$weight
    .save_base_plot(
        function() {
            par(mar = c(1, 1, 2, 1))
            netVisual_circle(mat,
                             vertex.weight = sizes,
                             weight.scale  = TRUE,
                             label.edge    = FALSE,
                             title.name    = "Interaction Strength")
        },
        file.path(output_dir, "interaction_strength_network.png"),
        width = 8, height = 8
    )
}

.plot_chord_aggregated <- function(cellchat, output_dir) {
    cat("   [3/6] chord diagram (aggregated)\n")
    weight_net <- cellchat@net$weight
    .save_base_plot(
        function() {
            par(mar = c(1, 1, 2, 1))
            netVisual_chord_cell(cellchat,
                                 net        = weight_net,
                                 title.name = "Cell-Cell Communication Network")
        },
        file.path(output_dir, "chord_aggregated.png"),
        width = 10, height = 10
    )
}

.plot_bubble_LR <- function(cellchat, output_dir) {
    cat("   [4/6] bubble plot of ligand-receptor pairs\n")
    p <- tryCatch(
        netVisual_bubble(cellchat, remove.isolate = TRUE),
        error = function(e) {
            cat("   Full bubble plot failed (", conditionMessage(e),
                "); retrying with top 5 source cell types\n", sep = "")
            top_sources <- levels(cellchat@idents)[seq_len(
                min(5L, length(levels(cellchat@idents))))]
            netVisual_bubble(cellchat,
                             sources.use    = top_sources,
                             remove.isolate = TRUE)
        }
    )

    if (!is.null(p) && inherits(p, "gg")) {
        p <- p + theme_prism(base_size = 10) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
                  axis.text.y = element_text(size = 7),
                  plot.title  = element_text(hjust = 0.5, face = "bold",
                                             size  = 12),
                  legend.position = "right")
        .save_ggplot(p, file.path(output_dir, "bubble_ligand_receptor.png"),
                     width = 12, height = 10)
    }
}

.plot_signaling_heatmaps <- function(cellchat, output_dir) {
    cat("   [5/6] signaling role heatmaps (outgoing + incoming)\n")
    for (pattern in c("outgoing", "incoming")) {
        title_str <- sprintf("%s Signaling Patterns",
                             tools::toTitleCase(pattern))
        png_name  <- sprintf("signaling_%s_heatmap.png", pattern)
        svg_name  <- sub("\\.png$", ".svg", png_name)

        tryCatch({
            ht <- netAnalysis_signalingRole_heatmap(cellchat,
                                                    pattern   = pattern,
                                                    title     = title_str,
                                                    font.size = 8)
            png_path <- file.path(output_dir, png_name)
            grDevices::png(png_path, width = 10, height = 8,
                           units = "in", res = 300)
            ComplexHeatmap::draw(ht)
            grDevices::dev.off()
            cat("   Saved:", png_path, "\n")

            svg_path <- file.path(output_dir, svg_name)
            tryCatch({
                grDevices::svg(svg_path, width = 10, height = 8)
                ComplexHeatmap::draw(ht)
                grDevices::dev.off()
                cat("   Saved:", svg_path, "\n")
            }, error = function(e2) {
                tryCatch(grDevices::dev.off(), error = function(e3) NULL)
                cat("   (SVG export failed for ", pattern, " heatmap)\n",
                    sep = "")
            })
        }, error = function(e) {
            cat("   ", pattern, " heatmap failed: ",
                conditionMessage(e), "\n", sep = "")
        })
    }
}

.plot_signaling_scatter <- function(cellchat, output_dir) {
    cat("   [6/6] signaling role scatter\n")
    tryCatch({
        p <- netAnalysis_signalingRole_scatter(cellchat)
        if (inherits(p, "gg")) {
            p <- p + theme_prism(base_size = 12) +
                theme(plot.title = element_text(hjust = 0.5, face = "bold",
                                                size = 14))
            .save_ggplot(p, file.path(output_dir, "signaling_role_scatter.png"),
                         width = 8, height = 7)
        }
    }, error = function(e) {
        cat("   Signaling scatter failed:", conditionMessage(e), "\n")
    })
}


# ---- main entry point --------------------------------------------------------

#' Generate the standard six-plot CellChat visualization bundle
#'
#' @param cellchat Trained CellChat object.
#' @param output_dir Directory for PNG + SVG output. Created if missing.
generate_all_plots <- function(cellchat, output_dir = "results") {
    cat("\n=== Generating CellChat visualizations ===\n")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    steps <- list(
        list(label = "interaction count network",      fn = .plot_interaction_count),
        list(label = "interaction strength network",   fn = .plot_interaction_strength),
        list(label = "chord diagram",                   fn = .plot_chord_aggregated),
        list(label = "bubble plot",                     fn = .plot_bubble_LR),
        list(label = "signaling role heatmaps",         fn = .plot_signaling_heatmaps),
        list(label = "signaling role scatter",          fn = .plot_signaling_scatter)
    )

    n_ok <- 0L
    for (s in steps) {
        ok <- tryCatch({
            s$fn(cellchat, output_dir)
            TRUE
        }, error = function(e) {
            cat("   ✗ ", s$label, " failed: ",
                conditionMessage(e), "\n", sep = "")
            FALSE
        })
        if (ok) n_ok <- n_ok + 1L
    }

    cat(sprintf("\n✓ All plots generated successfully: %d visualizations saved\n\n",
                n_ok))
}
