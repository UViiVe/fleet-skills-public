# run_gsea.R
#
# Purpose:
#   Run GSEA via clusterProfiler::GSEA on a pre-built ranked gene list and
#   TERM2GENE table. Sets sensible defaults: 10000 permutations, BH adjustment,
#   gene-set sizes 15-500, seed for reproducibility.
#
# Inputs:
#   ranked_genes : named numeric vector, descending; from create_ranked_list()
#   term2gene    : data.frame with columns (term, gene); from get_msigdb_genesets()
#   min_size     : default 15
#   max_size     : default 500
#   pvalue_cutoff: default 0.05
#   n_perm       : default 10000 (use 1000 for quick testing)
#
# Outputs:
#   clusterProfiler gseaResult object.
#
# Example:
#   source("scripts/run_gsea.R")
#   gsea <- run_gsea(ranked, term2gene)

run_gsea <- function(ranked_genes,
                     term2gene,
                     min_size      = 15L,
                     max_size      = 500L,
                     pvalue_cutoff = 0.05,
                     n_perm        = 10000L) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE))
    stop("Install with: BiocManager::install('clusterProfiler')")
  suppressPackageStartupMessages(library(clusterProfiler))

  if (is.null(names(ranked_genes)))
    stop("ranked_genes must be a *named* numeric vector")
  if (is.unsorted(-ranked_genes))
    ranked_genes <- sort(ranked_genes, decreasing = TRUE)

  cat(sprintf("[GSEA] %d genes ranked; %d gene sets; n_perm=%d; size=%d..%d\n",
              length(ranked_genes), length(unique(term2gene$term)),
              n_perm, min_size, max_size))

  res <- GSEA(
    geneList      = ranked_genes,
    TERM2GENE     = term2gene,
    minGSSize     = min_size,
    maxGSSize     = max_size,
    pvalueCutoff  = pvalue_cutoff,
    pAdjustMethod = "BH",
    nPermSimple   = n_perm,
    seed          = TRUE,
    verbose       = FALSE
  )

  n_sig    <- sum(res@result$p.adjust < pvalue_cutoff)
  n_up     <- sum(res@result$p.adjust < pvalue_cutoff & res@result$NES > 0)
  n_down   <- sum(res@result$p.adjust < pvalue_cutoff & res@result$NES < 0)
  cat(sprintf("[OK] GSEA complete: %d significant (p.adjust<%.3g); activated=%d suppressed=%d\n",
              n_sig, pvalue_cutoff, n_up, n_down))
  invisible(res)
}
