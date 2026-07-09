# run_ora.R
#
# Purpose:
#   Run Over-Representation Analysis (ORA) via clusterProfiler::enricher.
#   Always pass an explicit background (universe) — the tested-gene set
#   from your DE analysis, not the whole genome.
#
# Inputs:
#   genes        : character vector of significant gene symbols (one direction at a time)
#   term2gene    : TERM2GENE data.frame (columns term, gene)
#   background   : character vector of all tested genes (the universe)
#   direction    : informational tag for logging ("upregulated"/"downregulated")
#   min_size     : default 10
#   max_size     : default 500
#   pvalue_cutoff: default 0.05
#
# Outputs:
#   clusterProfiler enrichResult object.
#
# Example:
#   source("scripts/run_ora.R")
#   ora_up <- run_ora(sig$up, term2gene, background = sig$background,
#                     direction = "upregulated")

run_ora <- function(genes,
                    term2gene,
                    background,
                    direction     = "",
                    min_size      = 10L,
                    max_size      = 500L,
                    pvalue_cutoff = 0.05) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Install with: BiocManager::install('clusterProfiler')")
  suppressPackageStartupMessages(library(clusterProfiler))

  if (length(genes) == 0L) {
    cat(sprintf("[ORA-%s] no input genes; returning NULL\n", direction))
    return(invisible(NULL))
  }
  if (missing(background) || is.null(background) || length(background) == 0L)
    stop("background (universe) is required — pass the set of all tested genes")

  cat(sprintf("[ORA-%s] %d genes; %d gene sets; background=%d; size=%d..%d\n",
              direction, length(genes), length(unique(term2gene$term)),
              length(background), min_size, max_size))

  res <- enricher(
    gene          = genes,
    TERM2GENE     = term2gene,
    universe      = background,
    minGSSize     = min_size,
    maxGSSize     = max_size,
    pvalueCutoff  = pvalue_cutoff,
    pAdjustMethod = "BH"
  )

  n_sig <- if (is.null(res) || nrow(res@result) == 0) 0 else
           sum(res@result$p.adjust < pvalue_cutoff)
  cat(sprintf("[OK] ORA-%s: %d significant pathways (p.adjust<%.3g)\n",
              direction, n_sig, pvalue_cutoff))
  invisible(res)
}
