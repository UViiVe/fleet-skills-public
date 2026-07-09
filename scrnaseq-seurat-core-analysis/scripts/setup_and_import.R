# setup_and_import.R
# ----------------------------------------------------------------------
# Set up Seurat libraries and import scRNA-seq data from common formats.
#
# Inputs:  paths to 10X / H5 / RDS / CSV / TSV.
# Outputs: a Seurat object in memory.
# Example:
#   source("scripts/setup_and_import.R")
#   setup_seurat_libraries()
#   seurat_obj <- import_10x_data("filtered_feature_bc_matrix/")

#' Load core Seurat packages and silence noisy warnings.
setup_seurat_libraries <- function(verbose = TRUE) {
  required <- c("Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "dplyr")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Install missing packages: ", paste(missing, collapse = ", "))
  }
  suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(Matrix)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
  })
  if (verbose) {
    message("Seurat   ", packageVersion("Seurat"))
    message("R        ", R.version.string)
    message("Environment ready.")
  }
  invisible(NULL)
}

#' Import a 10X Genomics CellRanger directory.
#'
#' @param data_dir Path to filtered_feature_bc_matrix/.
#' @param project Project name for the Seurat object.
#' @param min_cells Minimum cells per gene.
#' @param min_features Minimum features per cell.
import_10x_data <- function(data_dir,
                            project = "scRNA",
                            min_cells = 3,
                            min_features = 200) {
  if (!dir.exists(data_dir)) stop("No such directory: ", data_dir)
  counts <- Seurat::Read10X(data.dir = data_dir)
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts, project = project,
    min.cells = min_cells, min.features = min_features
  )
  message(sprintf("✓ Data loaded successfully: %d genes × %d cells",
                  nrow(seurat_obj), ncol(seurat_obj)))
  seurat_obj
}

#' Import a 10X `.h5` file.
import_h5_data <- function(h5_file,
                            project = "scRNA",
                            min_cells = 3,
                            min_features = 200) {
  if (!file.exists(h5_file)) stop("No such file: ", h5_file)
  counts <- Seurat::Read10X_h5(h5_file)
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts, project = project,
    min.cells = min_cells, min.features = min_features
  )
  message(sprintf("✓ Data loaded successfully: %d genes × %d cells",
                  nrow(seurat_obj), ncol(seurat_obj)))
  seurat_obj
}

#' Import a CSV/TSV count matrix. Default orientation is genes × cells.
import_count_matrix <- function(file_path,
                                 transpose = FALSE,
                                 sep = NULL,
                                 project = "scRNA",
                                 min_cells = 3,
                                 min_features = 200) {
  if (!file.exists(file_path)) stop("No such file: ", file_path)
  if (is.null(sep)) {
    sep <- if (grepl("\\.tsv$", file_path, ignore.case = TRUE)) "\t" else ","
  }
  counts <- as.matrix(utils::read.table(
    file_path, sep = sep, header = TRUE, row.names = 1, check.names = FALSE
  ))
  if (transpose) counts <- t(counts)
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = counts, project = project,
    min.cells = min_cells, min.features = min_features
  )
  message(sprintf("✓ Data loaded successfully: %d genes × %d cells",
                  nrow(seurat_obj), ncol(seurat_obj)))
  seurat_obj
}

#' Add per-cell sample metadata from a CSV onto seurat@meta.data.
add_metadata <- function(seurat_obj, metadata_file, merge_on = "row.names") {
  if (!file.exists(metadata_file)) stop("No such file: ", metadata_file)
  meta <- utils::read.csv(metadata_file, stringsAsFactors = FALSE)
  if (merge_on == "row.names" && is.null(rownames(meta))) {
    rownames(meta) <- meta[, 1]; meta <- meta[, -1, drop = FALSE]
  }
  common <- intersect(colnames(seurat_obj), rownames(meta))
  if (length(common) == 0) {
    stop("No matching barcodes between Seurat object and metadata file.")
  }
  for (col in colnames(meta)) {
    seurat_obj@meta.data[common, col] <- meta[common, col]
  }
  message(sprintf("✓ Added %d metadata columns to %d cells",
                  ncol(meta), length(common)))
  seurat_obj
}
