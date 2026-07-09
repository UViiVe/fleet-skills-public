# load_example_data.R
# ----------------------------------------------------------------------
# Load example datasets from the SeuratData package for testing.
#
# Inputs:  dataset name (e.g. "pbmc3k", "ifnb").
# Outputs: a Seurat object with raw counts.
# Example:
#   source("scripts/load_example_data.R")
#   seurat_obj <- load_seurat_data("pbmc3k")

#' Load a named dataset from SeuratData. Installs the package if absent.
#'
#' Available: "pbmc3k" (default), "ifnb", and other SeuratData datasets.
load_seurat_data <- function(dataset_name = "pbmc3k", type = NULL) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))

  if (!requireNamespace("SeuratData", quietly = TRUE)) {
    message("Installing SeuratData (first time only)...")
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
    remotes::install_github("satijalab/seurat-data")
  }
  suppressPackageStartupMessages(library(SeuratData))

  message("Loading ", dataset_name, " from SeuratData...")

  if (!dataset_name %in% SeuratData::InstalledData()$Dataset) {
    message("Installing ", dataset_name, "...")
    SeuratData::InstallData(dataset_name)
  }

  if (is.null(type)) {
    type <- if (dataset_name == "pbmc3k") "default" else "filtered"
  }

  SeuratData::LoadData(dataset_name, type = type)
  seurat_obj <- get(dataset_name)

  if (utils::packageVersion("Seurat") >= "5.0.0") {
    seurat_obj <- Seurat::UpdateSeuratObject(seurat_obj)
  }

  message(sprintf("✓ Data loaded successfully: %d genes × %d cells",
                  nrow(seurat_obj), ncol(seurat_obj)))
  seurat_obj
}
