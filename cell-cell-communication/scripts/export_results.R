# =============================================================================
# Cell-cell communication: results export
# -----------------------------------------------------------------------------
# Writes the canonical CellChat result bundle:
#   * CSV tables for interactions, pathways, count/strength matrices, roles
#   * The full trained CellChat object as `.rds`
#   * A markdown analysis report (always)
#   * A PDF analysis report (when rmarkdown + LaTeX/tinytex are available)
#
# Inputs : trained CellChat object; (optional) original Seurat object;
#          output directory
# Outputs: files listed in SKILL.md "Outputs" section
# Example: source("scripts/export_results.R")
#          export_all(cellchat, seurat_obj = seurat_obj,
#                     output_dir = "results")
# =============================================================================

suppressPackageStartupMessages({
    library(CellChat)
})

#' Export the full CellChat results bundle
#'
#' @param cellchat Trained CellChat object.
#' @param seurat_obj Optional Seurat object used for the analysis. Currently
#'   accepted for API symmetry; reserved for future per-sample metadata dumps.
#' @param output_dir Directory for output files. Created if missing.
export_all <- function(cellchat, seurat_obj = NULL, output_dir = "results") {

    cat("\n=== Exporting CellChat results ===\n\n")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    df_net     <- .export_significant_interactions(cellchat, output_dir)
    df_pathway <- .export_pathway_summary(cellchat, output_dir)
    .export_count_matrix(cellchat, output_dir)
    .export_strength_matrix(cellchat, output_dir)
    .export_signaling_roles(cellchat, output_dir)
    .export_top_interactions(df_net, output_dir)
    .export_cellchat_rds(cellchat, output_dir)
    .write_markdown_report(cellchat, df_net, df_pathway, output_dir)
    .render_pdf_report(cellchat, df_net, df_pathway, output_dir)

    cat("\n=== Export Complete ===\n")
    cat("\nFiles written to '", output_dir, "/':\n", sep = "")
    files <- list.files(output_dir, full.names = FALSE)
    for (f in sort(files)) {
        cat(sprintf("  %-45s %s\n", f,
                    .humanize_size(file.size(file.path(output_dir, f)))))
    }
    cat("\n")

    invisible(NULL)
}


# ---- individual writers ------------------------------------------------------

.export_significant_interactions <- function(cellchat, output_dir) {
    cat("1. significant interactions\n")
    df_net <- subsetCommunication(cellchat)
    path <- file.path(output_dir, "significant_interactions.csv")
    write.csv(df_net, path, row.names = FALSE)
    cat(sprintf("   Saved: %s (%d interactions)\n", path, nrow(df_net)))
    df_net
}

.export_pathway_summary <- function(cellchat, output_dir) {
    cat("2. pathway-level summary\n")
    df_pw <- subsetCommunication(cellchat, slot.name = "netP")
    path <- file.path(output_dir, "pathway_summary.csv")
    write.csv(df_pw, path, row.names = FALSE)
    cat(sprintf("   Saved: %s (%d pathway-level interactions)\n",
                path, nrow(df_pw)))
    df_pw
}

.export_count_matrix <- function(cellchat, output_dir) {
    cat("3. interaction count matrix\n")
    path <- file.path(output_dir, "interaction_count_matrix.csv")
    write.csv(cellchat@net$count, path)
    cat(sprintf("   Saved: %s\n", path))
}

.export_strength_matrix <- function(cellchat, output_dir) {
    cat("4. interaction strength matrix\n")
    path <- file.path(output_dir, "interaction_strength_matrix.csv")
    write.csv(cellchat@net$weight, path)
    cat(sprintf("   Saved: %s\n", path))
}

.export_signaling_roles <- function(cellchat, output_dir) {
    cat("5. signaling role centrality\n")
    tryCatch({
        pathways <- cellchat@netP$pathways
        rows <- lapply(pathways, function(pw) {
            cen <- tryCatch(cellchat@netP$centr[[pw]], error = function(e) NULL)
            if (is.null(cen)) return(NULL)
            data.frame(
                pathway          = pw,
                cell_type        = names(cen$outdeg),
                outdeg_sender    = cen$outdeg,
                indeg_receiver   = cen$indeg,
                flowbet_mediator = if (!is.null(cen$flowbet)) cen$flowbet else NA_real_,
                info_influencer  = if (!is.null(cen$info))    cen$info    else NA_real_,
                stringsAsFactors = FALSE
            )
        })
        rows <- Filter(Negate(is.null), rows)
        if (length(rows) > 0L) {
            df_centr <- do.call(rbind, rows)
            rownames(df_centr) <- NULL
            path <- file.path(output_dir, "signaling_roles.csv")
            write.csv(df_centr, path, row.names = FALSE)
            cat(sprintf("   Saved: %s (%d pathways)\n", path, length(pathways)))
        }
    }, error = function(e) {
        cat("   Signaling roles export failed:", conditionMessage(e), "\n")
    })
}

.export_top_interactions <- function(df_net, output_dir) {
    cat("6. top 20 interactions\n")
    tryCatch({
        if (nrow(df_net) > 0L) {
            df_top <- head(df_net[order(df_net$prob, decreasing = TRUE), ], 20L)
            path <- file.path(output_dir, "top_interactions.csv")
            write.csv(df_top, path, row.names = FALSE)
            cat(sprintf("   Saved: %s\n", path))
        }
    }, error = function(e) {
        cat("   Top interactions export failed:", conditionMessage(e), "\n")
    })
}

.export_cellchat_rds <- function(cellchat, output_dir) {
    cat("7. CellChat object (RDS)\n")
    path <- file.path(output_dir, "cellchat_object.rds")
    saveRDS(cellchat, path)
    cat(sprintf("   Saved: %s\n", path))
    cat("   (load with: cellchat <- readRDS('cellchat_object.rds'))\n")
}


# ---- markdown report ---------------------------------------------------------

.write_markdown_report <- function(cellchat, df_net, df_pathway, output_dir) {
    cat("8. markdown report\n")
    pathways   <- cellchat@netP$pathways
    cell_types <- levels(cellchat@idents)

    lines <- c(
        "# Cell-Cell Communication Analysis Report",
        "",
        paste("Date:", Sys.Date()),
        "Tool: CellChat v2",
        "",
        "## Summary",
        "",
        sprintf("- Cell types analyzed: %d", length(cell_types)),
        sprintf("- Significant interactions: %d", nrow(df_net)),
        sprintf("- Active signaling pathways: %d", length(pathways)),
        "",
        "## Introduction",
        "",
        "Cell-cell communication was inferred from single-cell RNA-seq data",
        "using CellChat v2. The analysis identifies intercellular signaling by",
        "modeling ligand-receptor interactions and quantifying communication",
        "probabilities between cell populations.",
        "",
        "### Cell types",
        "",
        paste("-", cell_types),
        "",
        "## Results",
        "",
        "### Top signaling pathways",
        ""
    )

    if (length(pathways) > 0L) {
        prob <- cellchat@netP$prob
        if (!is.null(prob) && length(dim(prob)) == 3L) {
            strength <- sort(apply(prob, 3, sum), decreasing = TRUE)
            top_n <- min(15L, length(strength))
            lines <- c(lines,
                       "| Rank | Pathway | Strength |",
                       "|------|---------|----------|")
            for (i in seq_len(top_n)) {
                lines <- c(lines, sprintf("| %d | %s | %.4f |",
                                          i, names(strength)[i], strength[i]))
            }
        }
    }

    lines <- c(lines, "", "### Top ligand-receptor interactions", "")
    if (nrow(df_net) > 0L) {
        df_top <- head(df_net[order(df_net$prob, decreasing = TRUE), ], 15L)
        lines <- c(lines,
                   "| Source | Target | Interaction | Pathway | Probability |",
                   "|--------|--------|-------------|---------|-------------|")
        for (i in seq_len(nrow(df_top))) {
            lines <- c(lines, sprintf("| %s | %s | %s | %s | %.4f |",
                                      df_top$source[i],
                                      df_top$target[i],
                                      df_top$interaction_name_2[i],
                                      df_top$pathway_name[i],
                                      df_top$prob[i]))
        }
    }

    lines <- c(lines, "",
        "## Methods",
        "",
        "Communication probabilities were estimated with CellChat v2 (Jin et al.,",
        "Nature Communications 2021; Jin et al., Nature Protocols 2024).",
        "Ligand-receptor interactions were drawn from CellChatDB. The triMean",
        "estimator was used for `computeCommunProb()`. Network centrality",
        "(sender, receiver, mediator, influencer) was computed per pathway.",
        "",
        "## Conclusions",
        "",
        sprintf("- %d significant ligand-receptor interactions across %d signaling pathways",
                nrow(df_net), length(pathways)),
        sprintf("- Communication network analysis surfaces dominant senders and receivers among %d cell types",
                length(cell_types)),
        "- Signaling-role analysis identifies key mediator and influencer cell populations",
        "- Results can be used to prioritize targets for downstream functional follow-up",
        "",
        "## References",
        "",
        "- Jin S, et al. Inference and analysis of cell-cell communication using CellChat. Nature Communications. 2021;12:1088.",
        "- Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell communication from single-cell and spatially resolved transcriptomics. Nature Protocols. 2024."
    )

    path <- file.path(output_dir, "analysis_report.md")
    writeLines(lines, path)
    cat(sprintf("   Saved: %s\n", path))
}


# ---- PDF report (optional) ---------------------------------------------------

.render_pdf_report <- function(cellchat, df_net, df_pathway, output_dir) {
    cat("9. PDF report (optional)\n")
    tryCatch({
        candidate_paths <- c(
            "scripts/generate_report.R",
            file.path(dirname(output_dir), "scripts", "generate_report.R")
        )
        for (p in candidate_paths) {
            if (file.exists(p)) {
                source(p)
                break
            }
        }
        if (exists("generate_report", mode = "function")) {
            generate_report(cellchat, df_net, df_pathway,
                            output_dir = output_dir)
        } else {
            cat("   generate_report() not loaded; markdown report is sufficient\n")
        }
    }, error = function(e) {
        cat("   PDF generation skipped:", conditionMessage(e), "\n")
        cat("   (markdown report is still available)\n")
    })
}


# ---- utilities ---------------------------------------------------------------

.humanize_size <- function(bytes) {
    if (is.na(bytes)) return("")
    if (bytes < 1024)        return(paste0(bytes, " B"))
    if (bytes < 1024L^2L)    return(sprintf("%.1f KB", bytes / 1024))
    return(sprintf("%.1f MB", bytes / 1024L^2L))
}
