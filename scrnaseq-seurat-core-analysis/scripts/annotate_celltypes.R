# annotate_celltypes.R
# ----------------------------------------------------------------------
# Assign cell type labels (manual, SingleR, or Azimuth).
#
# Inputs:  Seurat object with clustering metadata.
# Outputs: Seurat object with `cell_type` (or `singler_labels`, `azimuth_pred`).
# Example:
#   source("scripts/annotate_celltypes.R")
#   annotations <- c("0" = "CD4 T", "1" = "CD14 mono", ...)
#   seurat_obj <- annotate_clusters_manual(seurat_obj, annotations, resolution = 0.8)

#' Map cluster IDs to cell-type labels.
annotate_clusters_manual <- function(seurat_obj,
                                       annotations,
                                       resolution = NULL,
                                       celltype_col = "cell_type") {
  if (!is.null(resolution)) {
    res_col <- grep(paste0("^(RNA|SCT|integrated)_snn_res\\.", resolution, "$"),
                     colnames(seurat_obj@meta.data), value = TRUE)
    if (length(res_col) == 1) Seurat::Idents(seurat_obj) <- res_col
  }
  clusters <- as.character(Seurat::Idents(seurat_obj))
  labels <- annotations[clusters]
  if (any(is.na(labels))) {
    n_missing <- sum(is.na(labels))
    message(sprintf("  warning: %d cells not in mapping → 'Unassigned'", n_missing))
    labels[is.na(labels)] <- "Unassigned"
  }
  seurat_obj[[celltype_col]] <- labels
  Seurat::Idents(seurat_obj) <- celltype_col
  message(sprintf("✓ Cell type annotation completed: %d unique labels",
                  length(unique(labels))))
  seurat_obj
}

#' Reference-based annotation with SingleR + celldex.
annotate_with_singler <- function(seurat_obj,
                                    reference = "HPCA",
                                    label_col = "singler_labels") {
  if (!requireNamespace("SingleR", quietly = TRUE) ||
      !requireNamespace("celldex", quietly = TRUE) ||
      !requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Install with: BiocManager::install(c('SingleR', 'celldex'))")
  }

  ref <- switch(reference,
    HPCA              = celldex::HumanPrimaryCellAtlasData(),
    BlueprintEncode   = celldex::BlueprintEncodeData(),
    ImmGen            = celldex::ImmGenData(),
    MonacoImmuneData  = celldex::MonacoImmuneData(),
    DICE              = celldex::DatabaseImmuneCellExpressionData(),
    stop("Unknown reference: ", reference)
  )

  test_data <- Seurat::GetAssayData(seurat_obj, slot = "data")
  pred <- SingleR::SingleR(
    test = test_data,
    ref = ref,
    labels = ref$label.main
  )
  seurat_obj[[label_col]] <- pred$pruned.labels
  Seurat::Idents(seurat_obj) <- label_col
  message(sprintf("✓ SingleR annotation completed (reference=%s): %d unique labels",
                  reference, length(unique(pred$pruned.labels))))
  seurat_obj
}

#' Azimuth label transfer from a pre-trained reference.
annotate_with_azimuth <- function(seurat_obj, reference = "pbmc") {
  if (!requireNamespace("Azimuth", quietly = TRUE)) {
    stop("Install with: remotes::install_github('satijalab/azimuth')")
  }
  seurat_obj <- Azimuth::RunAzimuth(seurat_obj, reference = reference)
  message(sprintf("✓ Azimuth annotation completed (reference=%s)", reference))
  seurat_obj
}

#' UMAP coloured by cell-type column.
plot_annotated_umap <- function(seurat_obj,
                                  group_by = "cell_type",
                                  output_dir = "results/annotation") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  p <- Seurat::DimPlot(seurat_obj, group.by = group_by, label = TRUE,
                       repel = TRUE, label.size = 3)
  out_png <- file.path(output_dir, sprintf("umap_%s.png", group_by))
  ggplot2::ggsave(out_png, p, width = 8, height = 6, dpi = 200)
  message(sprintf("✓ Saved %s", out_png))
}

#' Stacked bar plot of cell-type proportions per group.
plot_celltype_proportions <- function(seurat_obj,
                                        group_by = "sample_id",
                                        celltype_col = "cell_type",
                                        output_dir = "results/annotation") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- seurat_obj@meta.data
  if (!all(c(group_by, celltype_col) %in% colnames(meta))) {
    message("  missing ", paste(setdiff(c(group_by, celltype_col), colnames(meta)),
                                  collapse = ", "))
    return(invisible(NULL))
  }
  tab <- as.data.frame(table(meta[[group_by]], meta[[celltype_col]]))
  colnames(tab) <- c("group", "celltype", "n")
  p <- ggplot2::ggplot(tab, ggplot2::aes(x = group, y = n, fill = celltype)) +
    ggplot2::geom_bar(stat = "identity", position = "fill") +
    ggplot2::labs(y = "fraction", x = group_by) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  out_png <- file.path(output_dir, "celltype_proportions.png")
  ggplot2::ggsave(out_png, p, width = 8, height = 5, dpi = 200)
  message(sprintf("✓ Saved %s", out_png))
}
