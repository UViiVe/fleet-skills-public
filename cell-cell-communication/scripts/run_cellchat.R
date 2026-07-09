# =============================================================================
# Cell-cell communication: CellChat v2 inference pipeline
# -----------------------------------------------------------------------------
# Wraps the six core CellChat steps into one entry point:
#   1. createCellChat()
#   2. set CellChatDB (species + optional category filter)
#   3. identifyOverExpressedGenes()/Interactions()
#   4. computeCommunProb() + filterCommunication()
#   5. computeCommunProbPathway() + aggregateNet()
#   6. netAnalysis_computeCentrality()
#
# Inputs : Seurat object validated by scripts/load_data.R
# Outputs: a fully analysed CellChat object
# Example: source("scripts/run_cellchat.R")
#          cellchat <- run_cellchat_analysis(seurat_obj, species = "human")
# =============================================================================

suppressPackageStartupMessages({
    library(CellChat)
    library(Seurat)
})

#' Run the complete CellChat v2 pipeline on a Seurat object
#'
#' @param seurat_obj Annotated Seurat object (see scripts/load_data.R).
#' @param species `"human"` or `"mouse"` — selects CellChatDB.
#' @param group.by Metadata column with cell-type labels (default `"celltype"`).
#' @param db_category One of `"all"`, `"Secreted Signaling"`,
#'   `"ECM-Receptor"`, or `"Cell-Cell Contact"`.
#' @param min.cells Cell-group filter for `filterCommunication()`.
#' @return A CellChat object with all standard slots populated.
run_cellchat_analysis <- function(seurat_obj,
                                  species     = "human",
                                  group.by    = "celltype",
                                  db_category = "all",
                                  min.cells   = 10L) {

    cat("\n=== Running CellChat v2 analysis ===\n\n")

    cellchat <- .step_create(seurat_obj, group.by)
    cellchat <- .step_set_db(cellchat, species, db_category)
    cellchat <- .step_overexpressed(cellchat)
    cellchat <- .step_communication(cellchat, min.cells)
    pathway_info <- .step_pathway_aggregation(cellchat)
    cellchat <- pathway_info$cellchat
    cellchat <- .step_centrality(cellchat)

    # Final inventory.
    df_net <- subsetCommunication(cellchat)
    n_sig <- nrow(df_net)
    n_paths <- length(cellchat@netP$pathways)
    n_LR <- nrow(cellchat@LR$LRsig)

    cat("\n--- Analysis summary ---\n")
    cat("   Species              :", species, "\n")
    cat("   Cell types           :", length(levels(cellchat@idents)), "\n")
    cat("   L-R pairs tested     :", n_LR, "\n")
    cat("   Significant L-R pairs:", n_sig, "\n")
    cat("   Active pathways      :", n_paths, "\n")

    cat(sprintf("\n✓ CellChat analysis completed: %d significant interactions across %d pathways\n\n",
                n_sig, n_paths))

    invisible(cellchat)
}


# ---- internal pipeline steps -------------------------------------------------

.step_create <- function(seurat_obj, group.by) {
    cat("Step 1/6: createCellChat()\n")
    cc <- createCellChat(object = seurat_obj, group.by = group.by)
    cat("   Cell groups:", length(levels(cc@idents)), "\n")
    cc
}

.step_set_db <- function(cellchat, species, db_category) {
    cat("Step 2/6: setting CellChatDB\n")
    db <- switch(tolower(species),
                 human = CellChatDB.human,
                 mouse = CellChatDB.mouse,
                 stop("species must be 'human' or 'mouse'; got: ", species))
    cat(sprintf("   Using CellChatDB.%s\n", tolower(species)))

    if (!identical(db_category, "all")) {
        valid <- unique(db$interaction$annotation)
        if (!db_category %in% valid) {
            stop("Invalid db_category '", db_category, "'. Valid options: ",
                 paste(valid, collapse = ", "))
        }
        db <- subsetDB(db, search = db_category)
        cat(sprintf("   Restricted to: %s\n", db_category))
    } else {
        cat("   Using all signaling categories\n")
    }

    cellchat@DB <- db
    cat("   Database L-R interactions:", nrow(db$interaction), "\n")
    cellchat
}

.step_overexpressed <- function(cellchat) {
    cat("Step 3/6: overexpressed signaling genes\n")
    cellchat <- subsetData(cellchat)

    has_presto <- requireNamespace("presto", quietly = TRUE)
    if (!has_presto) {
        cat("   (presto not installed; using standard Wilcoxon test, slower)\n")
    }
    cellchat <- identifyOverExpressedGenes(cellchat, do.fast = has_presto)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    cat("   Overexpressed L-R pairs:", nrow(cellchat@LR$LRsig), "\n")
    cellchat
}

.step_communication <- function(cellchat, min.cells) {
    cat("Step 4/6: communication probabilities\n")
    cat("   (this step can take 1-3 minutes)\n")
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = min.cells)

    df_net <- subsetCommunication(cellchat)
    cat("   Significant cell-cell interactions:", nrow(df_net), "\n")
    cellchat
}

.step_pathway_aggregation <- function(cellchat) {
    cat("Step 5/6: pathway-level aggregation\n")
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)

    pathways <- cellchat@netP$pathways
    n_paths <- length(pathways)
    cat("   Aggregated into pathways:", n_paths, "\n")

    if (n_paths > 0L) {
        prob <- cellchat@netP$prob
        if (!is.null(prob) && length(dim(prob)) == 3L) {
            strength <- sort(apply(prob, 3, sum), decreasing = TRUE)
            top_n <- min(10L, length(strength))
            cat("\n   Top signaling pathways:\n")
            for (i in seq_len(top_n)) {
                cat(sprintf("     %2d. %-20s strength = %.4f\n",
                            i, names(strength)[i], strength[i]))
            }
        }
    }

    list(cellchat = cellchat)
}

.step_centrality <- function(cellchat) {
    cat("\nStep 6/6: network centrality\n")
    cellchat <- netAnalysis_computeCentrality(cellchat)
    cat("   Computed sender / receiver / mediator / influencer roles\n")
    cellchat
}
