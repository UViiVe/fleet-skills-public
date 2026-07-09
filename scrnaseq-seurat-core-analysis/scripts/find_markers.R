# find_markers.R
# ----------------------------------------------------------------------
# Identify cluster markers (Wilcoxon by default). Exploratory — not for
# condition inference (use pseudobulk for that).
#
# Inputs:  Seurat object with clustering metadata.
# Outputs: a data.frame of cluster markers, plus optional plots.
# Example:
#   source("scripts/find_markers.R")
#   markers <- find_all_cluster_markers(seurat_obj, resolution = 0.8)

#' Find markers across all clusters.
find_all_cluster_markers <- function(seurat_obj,
                                       resolution = NULL,
                                       method = "wilcox",
                                       min_pct = 0.25,
                                       logfc_threshold = 0.25,
                                       features_to_exclude = NULL,
                                       only_pos = TRUE) {
  if (!is.null(resolution)) {
    res_col <- grep(paste0("^(RNA|SCT|integrated)_snn_res\\.", resolution, "$"),
                     colnames(seurat_obj@meta.data), value = TRUE)
    if (length(res_col) == 1) {
      Seurat::Idents(seurat_obj) <- res_col
    }
  }

  obj <- seurat_obj
  if (!is.null(features_to_exclude)) {
    keep <- setdiff(rownames(obj), features_to_exclude)
    obj <- subset(obj, features = keep)
  }

  markers <- Seurat::FindAllMarkers(
    obj,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    test.use = method,
    verbose = FALSE
  )
  message(sprintf("✓ Found markers for %d clusters (%d total rows)",
                  length(unique(markers$cluster)), nrow(markers)))
  markers
}

#' Find markers for one cluster vs the rest.
find_markers_for_cluster <- function(seurat_obj,
                                       cluster,
                                       resolution = NULL,
                                       method = "wilcox") {
  if (!is.null(resolution)) {
    res_col <- grep(paste0("^(RNA|SCT|integrated)_snn_res\\.", resolution, "$"),
                     colnames(seurat_obj@meta.data), value = TRUE)
    if (length(res_col) == 1) Seurat::Idents(seurat_obj) <- res_col
  }
  Seurat::FindMarkers(
    seurat_obj, ident.1 = cluster, test.use = method, verbose = FALSE
  )
}

#' Export the markers table (all + per-cluster).
export_marker_tables <- function(all_markers,
                                   output_dir = "results/markers",
                                   top_n = 50) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(all_markers,
                    file.path(output_dir, "cluster_markers_all.csv"),
                    row.names = FALSE)
  for (cluster in unique(all_markers$cluster)) {
    sub <- all_markers[all_markers$cluster == cluster, ]
    sub <- sub[order(-sub$avg_log2FC), ]
    if (!is.null(top_n)) sub <- utils::head(sub, top_n)
    utils::write.csv(sub,
                      file.path(output_dir, sprintf("markers_cluster_%s.csv", cluster)),
                      row.names = FALSE)
  }
  message(sprintf("✓ Exported markers to %s", output_dir))
}

#' Heatmap of top-N markers per cluster.
plot_top_markers_heatmap <- function(seurat_obj,
                                       all_markers,
                                       n_top = 10,
                                       output_dir = "results/markers") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  top <- all_markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = n_top) |>
    dplyr::pull(gene) |>
    unique()
  p <- Seurat::DoHeatmap(seurat_obj, features = top, size = 3) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
  out_png <- file.path(output_dir, "marker_heatmap.png")
  ggplot2::ggsave(out_png, p, width = 10, height = 8, dpi = 200)
  message(sprintf("✓ Saved %s", out_png))
}

#' Dot plot of top-N markers per cluster.
plot_markers_dotplot <- function(seurat_obj,
                                   all_markers,
                                   n_top = 5,
                                   output_dir = "results/markers") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  top <- all_markers |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(avg_log2FC, n = n_top) |>
    dplyr::pull(gene) |>
    unique()
  p <- Seurat::DotPlot(seurat_obj, features = top) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  out_png <- file.path(output_dir, "marker_dotplot.png")
  ggplot2::ggsave(out_png, p, width = 12, height = 5, dpi = 200)
  message(sprintf("✓ Saved %s", out_png))
}

#' Violin plots of selected markers.
plot_markers_violin <- function(seurat_obj,
                                  genes,
                                  output_dir = "results/markers") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  available <- intersect(genes, rownames(seurat_obj))
  if (length(available) == 0) {
    message("  none of the requested genes were found")
    return(invisible(NULL))
  }
  p <- Seurat::VlnPlot(seurat_obj, features = available, pt.size = 0)
  out_png <- file.path(output_dir, "marker_violin.png")
  ggplot2::ggsave(out_png, p, width = 3 * length(available), height = 4, dpi = 200)
  message(sprintf("✓ Saved %s", out_png))
}
