# =============================================================================
# Cell-cell communication: PDF report (optional)
# -----------------------------------------------------------------------------
# Builds a publication-style PDF report from the analysis using rmarkdown.
# Skips gracefully when rmarkdown or a LaTeX engine is unavailable; the
# markdown report written by export_results.R is always present.
#
# Inputs : CellChat object; significant-interactions data frame;
#          pathway-summary data frame; output directory
# Outputs: <output_dir>/analysis_report.pdf when LaTeX is available
# Example: source("scripts/generate_report.R")
#          generate_report(cellchat, df_net, df_pathway,
#                          output_dir = "results")
# =============================================================================

#' Render a PDF analysis report from CellChat results
#'
#' @param cellchat Trained CellChat object.
#' @param df_net Data frame from `subsetCommunication(cellchat)`.
#' @param df_pathway Data frame from `subsetCommunication(cellchat, slot.name = "netP")`.
#' @param output_dir Directory containing plot PNGs; PDF written here too.
#' @return The path to the PDF (invisibly), or NULL when skipped.
generate_report <- function(cellchat, df_net, df_pathway,
                            output_dir = "results") {

    if (!requireNamespace("rmarkdown", quietly = TRUE)) {
        cat("   rmarkdown not installed; skipping PDF report\n")
        cat("   Install with: install.packages('rmarkdown')\n")
        return(invisible(NULL))
    }

    have_tinytex <- requireNamespace("tinytex", quietly = TRUE) &&
                    tinytex::is_tinytex()
    have_xelatex  <- nchar(Sys.which("xelatex"))  > 0
    have_pdflatex <- nchar(Sys.which("pdflatex")) > 0
    if (!any(c(have_tinytex, have_xelatex, have_pdflatex))) {
        cat("   No LaTeX engine found; skipping PDF report\n")
        cat("   Install with: tinytex::install_tinytex()\n")
        return(invisible(NULL))
    }

    rmd_path <- file.path(output_dir, "_report.Rmd")
    pdf_path <- file.path(output_dir, "analysis_report.pdf")

    pathways   <- cellchat@netP$pathways
    cell_types <- levels(cellchat@idents)
    n_int      <- nrow(df_net)
    n_path     <- length(pathways)
    n_types    <- length(cell_types)

    top_pathways_md     <- .top_pathways_markdown(cellchat, pathways)
    top_interactions_md <- .top_interactions_markdown(df_net)
    plot_includes_md    <- .plot_includes_markdown(output_dir)

    rmd_body <- paste0(
'---
title: "Cell-Cell Communication Analysis Report"
date: "', Sys.Date(), '"
output:
  pdf_document:
    toc: true
    toc_depth: 2
    latex_engine: xelatex
header-includes:
  - \\usepackage{booktabs}
  - \\usepackage{float}
  - \\floatplacement{figure}{H}
---

# Summary

| Metric | Value |
|--------|-------|
| Cell types | ', n_types, ' |
| Significant interactions | ', n_int, ' |
| Active signaling pathways | ', n_path, ' |
| Analysis tool | CellChat v2 |

**Cell types analyzed:** ', paste(cell_types, collapse = ", "), '

# Introduction

This report presents cell-cell communication analysis results inferred from
single-cell RNA-seq data using CellChat v2. CellChat identifies intercellular
signaling by modeling ligand-receptor interactions and quantifying
communication probabilities between cell populations.

# Methods

- Ligand-receptor database: CellChatDB v2
- Communication inference: triMean estimator
- Minimum cells per group: 10 (default `filterCommunication`)
- Network centrality computed for sender, receiver, mediator, influencer roles

# Results

## Top signaling pathways

', top_pathways_md, '

## Top ligand-receptor interactions

', top_interactions_md, '

## Visualizations

', plot_includes_md, '

# Conclusions

- **', n_int, '** significant ligand-receptor interactions identified across **',
n_path, '** signaling pathways
- Communication network analysis reveals dominant senders and receivers among the ',
n_types, ' cell types
- Signaling-role analysis identifies key mediator and influencer populations
- These results can guide target prioritization for follow-up assays

# References

1. Jin S, et al. Inference and analysis of cell-cell communication using CellChat.
   *Nature Communications*. 2021;12:1088.
2. Jin S, Plikus MV, Nie Q. CellChat for systematic analysis of cell-cell
   communication from single-cell and spatially resolved transcriptomics.
   *Nature Protocols*. 2024.
')

    writeLines(rmd_body, rmd_path)

    tryCatch({
        rmarkdown::render(
            rmd_path,
            output_file = basename(pdf_path),
            output_dir  = output_dir,
            quiet       = TRUE
        )
        cat("   Saved:", pdf_path, "\n")
    }, error = function(e) {
        cat("   PDF rendering failed:", conditionMessage(e), "\n")
        cat("   (markdown report is still available)\n")
    })

    unlink(rmd_path)
    invisible(pdf_path)
}


# ---- helpers -----------------------------------------------------------------

.top_pathways_markdown <- function(cellchat, pathways) {
    if (length(pathways) == 0L) return("(no active pathways)")
    prob <- cellchat@netP$prob
    if (is.null(prob) || length(dim(prob)) != 3L) return("(pathway probabilities unavailable)")
    strength <- sort(apply(prob, 3, sum), decreasing = TRUE)
    top_n <- min(10L, length(strength))
    df <- data.frame(
        Pathway  = names(strength)[seq_len(top_n)],
        Strength = round(strength[seq_len(top_n)], 4)
    )
    paste(knitr::kable(df, format = "pipe"), collapse = "\n")
}

.top_interactions_markdown <- function(df_net) {
    if (nrow(df_net) == 0L) return("(no significant interactions)")
    df_top <- head(df_net[order(df_net$prob, decreasing = TRUE), ], 10L)
    df <- data.frame(
        Source       = df_top$source,
        Target       = df_top$target,
        Interaction  = df_top$interaction_name_2,
        Pathway      = df_top$pathway_name,
        Prob         = round(df_top$prob, 4)
    )
    paste(knitr::kable(df, format = "pipe"), collapse = "\n")
}

.plot_includes_markdown <- function(output_dir) {
    pngs <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
    pngs <- pngs[!grepl("_report", pngs)]
    if (length(pngs) == 0L) return("(no plot files found)")
    parts <- character(length(pngs))
    for (i in seq_along(pngs)) {
        nm <- tools::toTitleCase(gsub("_", " ",
                                      tools::file_path_sans_ext(basename(pngs[i]))))
        parts[i] <- paste0("\n### ", nm, "\n\n![", nm, "](",
                           basename(pngs[i]), ")\n")
    }
    paste(parts, collapse = "\n")
}
