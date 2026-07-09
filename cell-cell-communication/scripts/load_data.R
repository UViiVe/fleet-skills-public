# =============================================================================
# Cell-cell communication: data loading
# -----------------------------------------------------------------------------
# Provides two entry points:
#   load_cellchat_data(seurat_path = NULL, group.by = "celltype")
#       -> dispatches to the example loader or to the user-file loader
#   load_example_pbmc()
#       -> loads the SeuratData PBMC 3k dataset and returns it pre-validated
#
# Inputs : optional .rds path or NULL for the demo
# Outputs: a Seurat object validated for CellChat (normalized data present,
#          cell-type metadata column present, per-type cell counts logged)
# Example: source("scripts/load_data.R"); obj <- load_cellchat_data()
# =============================================================================

suppressPackageStartupMessages({
    library(Seurat)
})

#' Load the PBMC 3k example dataset (pre-annotated; 2,638 cells, 8 cell types)
#'
#' Installs SeuratData (and devtools) on demand, materializes the
#' `pbmc3k.final` dataset, upgrades it to Seurat v5 format, ensures a
#' `celltype` metadata column, and renormalizes if the data layer is empty.
#'
#' @return A Seurat object ready for CellChat with `celltype` in `meta.data`.
load_example_pbmc <- function() {
    cat("\n=== Loading PBMC 3k example dataset ===\n\n")

    options(repos = c(CRAN = "https://cloud.r-project.org"))

    if (!requireNamespace("SeuratData", quietly = TRUE)) {
        cat("   SeuratData not found; installing...\n")
        if (!requireNamespace("devtools", quietly = TRUE)) {
            install.packages("devtools")
        }
        devtools::install_github("satijalab/seurat-data", upgrade = "never")
    }

    library(SeuratData)

    available <- AvailableData()
    pbmc_row <- grep("^pbmc3k", rownames(available), value = TRUE)[1]
    if (is.na(pbmc_row)) {
        stop("pbmc3k dataset is not registered in SeuratData. ",
             "Reinstall with: remotes::install_github('satijalab/seurat-data')")
    }
    if (!available[pbmc_row, "Installed"]) {
        cat("   Downloading pbmc3k (~6 MB) ...\n")
        InstallData("pbmc3k")
    }

    holder <- new.env()
    utils::data("pbmc3k.final", package = "pbmc3k.SeuratData", envir = holder)
    seurat_obj <- holder$pbmc3k.final

    cat("   Upgrading Seurat object to v5 ...\n")
    seurat_obj <- UpdateSeuratObject(seurat_obj)

    if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj$celltype <- Idents(seurat_obj)
    }

    # Re-normalize defensively if the `data` layer looks unpopulated.
    needs_norm <- tryCatch({
        d <- GetAssayData(seurat_obj, layer = "data")
        nr <- min(5L, nrow(d))
        nc <- min(5L, ncol(d))
        all(d[seq_len(nr), seq_len(nc)] == 0)
    }, error = function(e) TRUE)
    if (needs_norm) {
        cat("   Running NormalizeData() ...\n")
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    }

    .report_summary(seurat_obj, source_label = "PBMC 3k (Satija Lab, 10X Genomics)",
                    species_label = "Human")

    cat(sprintf("\n✓ Data loaded successfully: %d cells, %d cell types\n\n",
                ncol(seurat_obj),
                length(unique(seurat_obj$celltype))))

    invisible(seurat_obj)
}


#' Load a Seurat object for CellChat from a file (or fall back to the demo)
#'
#' @param seurat_path Path to a `.rds` file produced by `saveRDS()`. If NULL,
#'   delegates to `load_example_pbmc()`.
#' @param group.by Metadata column carrying cell-type labels. Default `"celltype"`.
#' @return A validated Seurat object.
load_cellchat_data <- function(seurat_path = NULL, group.by = "celltype") {
    if (is.null(seurat_path)) {
        return(load_example_pbmc())
    }

    cat("\n=== Loading Seurat object from disk ===\n\n")

    if (!file.exists(seurat_path)) {
        stop("File not found: ", seurat_path)
    }
    cat("   Reading:", seurat_path, "\n")
    seurat_obj <- readRDS(seurat_path)

    seurat_obj <- validate_seurat_for_cellchat(seurat_obj, group.by = group.by)

    n_types <- length(unique(seurat_obj[[group.by, drop = TRUE]]))
    cat(sprintf("\n✓ Data loaded successfully: %d cells, %d cell types\n\n",
                ncol(seurat_obj), n_types))

    invisible(seurat_obj)
}


#' Validate that a Seurat object can be handed to CellChat
#'
#' Checks class, the presence of `group.by`, populates normalized data if
#' missing, and warns about cell types with fewer than 10 cells.
#'
#' @param seurat_obj A Seurat object.
#' @param group.by Metadata column with cell-type labels.
#' @return The validated (possibly re-normalized) Seurat object.
validate_seurat_for_cellchat <- function(seurat_obj, group.by = "celltype") {
    if (!inherits(seurat_obj, "Seurat")) {
        stop("Input must be a Seurat object; got class: ", class(seurat_obj)[1])
    }

    if (!group.by %in% colnames(seurat_obj@meta.data)) {
        cols <- colnames(seurat_obj@meta.data)
        likely <- grep("cell.?type|ident|label|annotation|cluster",
                       cols, ignore.case = TRUE, value = TRUE)
        msg <- paste0("Column '", group.by, "' not found in metadata.\n",
                      "  Available columns: ", paste(cols, collapse = ", "))
        if (length(likely)) {
            msg <- paste0(msg, "\n  Candidate cell-type columns: ",
                          paste(likely, collapse = ", "))
        }
        stop(msg)
    }

    norm_data <- tryCatch(GetAssayData(seurat_obj, layer = "data"),
                          error = function(e) NULL)
    if (is.null(norm_data)) {
        cat("   No normalized data found; running NormalizeData() ...\n")
        seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
    }

    ct_counts <- table(seurat_obj[[group.by, drop = TRUE]])
    small_types <- ct_counts[ct_counts < 10L]
    if (length(small_types)) {
        cat("   Warning: cell types with <10 cells (CellChat may drop these):\n")
        for (ct in names(small_types)) {
            cat(sprintf("     %s: %d cells\n", ct, small_types[ct]))
        }
    }

    .report_summary(seurat_obj, source_label = NULL, species_label = NULL,
                    group.by = group.by)
    invisible(seurat_obj)
}


# ---- internal helpers --------------------------------------------------------

.report_summary <- function(seurat_obj,
                            source_label = NULL,
                            species_label = NULL,
                            group.by = "celltype") {
    if (!is.null(source_label)) cat("   Dataset:", source_label, "\n")
    if (!is.null(species_label)) cat("   Species:", species_label, "\n")
    cat("   Cells :", ncol(seurat_obj), "\n")
    cat("   Genes :", nrow(seurat_obj), "\n")
    cat("   Cell types (", group.by, "):\n", sep = "")
    counts <- sort(table(seurat_obj[[group.by, drop = TRUE]]), decreasing = TRUE)
    for (nm in names(counts)) {
        cat(sprintf("     %-25s %d cells\n", nm, counts[nm]))
    }
}
