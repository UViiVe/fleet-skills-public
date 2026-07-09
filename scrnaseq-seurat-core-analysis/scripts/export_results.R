# export_results.R
# ----------------------------------------------------------------------
# Export analysis artifacts: RDS, CSV tables, embeddings, markers, summary.
#
# Inputs:  Seurat object (annotated).
# Outputs: results/seurat_processed.rds plus tables and a summary text file.
# Example:
#   source("scripts/export_results.R")
#   export_all(seurat_obj, output_dir = "results", all_markers = markers)

#' One-shot export.
export_all <- function(seurat_obj,
                         output_dir = "results",
                         all_markers = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  export_seurat_object(seurat_obj, output_dir)
  export_count_matrices(seurat_obj, output_dir)
  export_metadata(seurat_obj, output_dir)
  export_dimreductions(seurat_obj, output_dir)
  if (!is.null(all_markers)) export_markers(all_markers, output_dir)
  export_summary_stats(seurat_obj, output_dir)
  message("=== Export Complete ===")
}

#' Alias for `export_all` (matches the original script API).
export_seurat_results <- function(...) export_all(...)

#' Save the Seurat object.
export_seurat_object <- function(seurat_obj,
                                    output_dir,
                                    file_name = "seurat_processed.rds") {
  out <- file.path(output_dir, file_name)
  saveRDS(seurat_obj, out)
  message(sprintf("✓ Wrote %s", out))
}

#' Save normalized counts as CSV + RDS.
export_count_matrices <- function(seurat_obj, output_dir) {
  data <- Seurat::GetAssayData(seurat_obj, slot = "data")
  saveRDS(data, file.path(output_dir, "normalized_counts.rds"))
  if (ncol(data) <= 20000 && nrow(data) <= 50000) {
    df <- as.data.frame(as.matrix(data))
    utils::write.csv(df, file.path(output_dir, "normalized_counts.csv"))
    message(sprintf("✓ Wrote normalized_counts.csv (%d × %d)", nrow(df), ncol(df)))
  } else {
    message("  skipped CSV export (too large); RDS only")
  }
}

#' Save meta.data.
export_metadata <- function(seurat_obj, output_dir, resolution = NULL) {
  meta <- seurat_obj@meta.data
  out <- file.path(output_dir, "cell_metadata.csv")
  utils::write.csv(meta, out)
  message(sprintf("✓ Wrote %s (%d × %d)", out, nrow(meta), ncol(meta)))
}

#' Save PCA / UMAP / Harmony coordinates.
export_dimreductions <- function(seurat_obj, output_dir) {
  for (red in names(seurat_obj@reductions)) {
    coords <- seurat_obj@reductions[[red]]@cell.embeddings
    out <- file.path(output_dir, sprintf("%s_coordinates.csv", red))
    utils::write.csv(coords, out)
    message(sprintf("✓ Wrote %s", out))
  }
}

#' Save the markers data frame.
export_markers <- function(all_markers, output_dir) {
  out <- file.path(output_dir, "cluster_markers_all.csv")
  utils::write.csv(all_markers, out, row.names = FALSE)
  message(sprintf("✓ Wrote %s (%d rows)", out, nrow(all_markers)))
}

#' Write `analysis_summary.txt`.
export_summary_stats <- function(seurat_obj, output_dir) {
  lines <- c(
    "Single-cell RNA-seq Analysis Summary (Seurat)",
    paste(rep("=", 40), collapse = ""),
    "",
    sprintf("Cells: %d", ncol(seurat_obj)),
    sprintf("Genes: %d", nrow(seurat_obj))
  )
  for (m in c("nFeature_RNA", "nCount_RNA", "percent.mt")) {
    if (m %in% colnames(seurat_obj@meta.data)) {
      lines <- c(lines, sprintf("  median %s: %.2f", m,
                                  median(seurat_obj@meta.data[[m]])))
    }
  }
  if (!is.null(seurat_obj@misc$integration)) {
    lines <- c(lines, "", "Integration:")
    for (k in names(seurat_obj@misc$integration)) {
      lines <- c(lines, sprintf("  %s: %s", k,
                                  paste(seurat_obj@misc$integration[[k]], collapse = ",")))
    }
  }
  if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
    lines <- c(lines, "", "Cell type counts:")
    counts <- sort(table(seurat_obj@meta.data$cell_type), decreasing = TRUE)
    for (k in names(counts)) {
      lines <- c(lines, sprintf("  %s: %d", k, counts[[k]]))
    }
  }
  out <- file.path(output_dir, "analysis_summary.txt")
  writeLines(lines, out)
  message(sprintf("✓ Wrote %s", out))
}
