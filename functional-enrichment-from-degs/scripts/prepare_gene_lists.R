# prepare_gene_lists.R
#
# Purpose:
#   Build the two inputs used by downstream enrichment functions:
#   - create_ranked_list(de): named numeric vector for GSEA, sorted descending
#   - filter_significant_genes(de, ...): up / down / background gene-symbol lists for ORA
#
# Inputs:
#   de : data.frame with columns gene, log2FC, padj, and optionally stat/pvalue.
#
# Outputs:
#   See per-function descriptions below.
#
# Example:
#   source("scripts/prepare_gene_lists.R")
#   ranked <- create_ranked_list(de)
#   sig    <- filter_significant_genes(de, padj_thresh = 0.05, log2fc_thresh = 1)

#' Build a ranked gene list (descending) for GSEA.
#'
#' Ranking metric preference: stat > sign(log2FC)*-log10(pvalue) > sign(log2FC)*-log10(padj).
#' Drops duplicate gene symbols (keeps best padj) and NA / inf values.
create_ranked_list <- function(de) {
  required <- c("gene", "log2FC", "padj")
  missing <- setdiff(required, colnames(de))
  if (length(missing) > 0) stop("Required columns missing: ", paste(missing, collapse = ", "))

  df <- de
  if (any(duplicated(df$gene))) {
    n_dup <- sum(duplicated(df$gene))
    cat(sprintf("  collapsing %d duplicate symbols (keeping smallest padj)\n", n_dup))
    df <- df[order(df$padj), ]
    df <- df[!duplicated(df$gene), ]
  }

  if ("stat" %in% colnames(df)) {
    df$rank_metric <- df$stat
    cat("  ranking by test statistic (best)\n")
  } else if ("pvalue" %in% colnames(df)) {
    df$rank_metric <- sign(df$log2FC) * -log10(pmax(df$pvalue, 1e-300))
    cat("  ranking by sign(log2FC) * -log10(pvalue)\n")
  } else {
    df$rank_metric <- sign(df$log2FC) * -log10(pmax(df$padj, 1e-300))
    cat("  ranking by sign(log2FC) * -log10(padj)  (last-resort fallback)\n")
  }

  df <- df[is.finite(df$rank_metric), ]
  ranked <- df$rank_metric
  names(ranked) <- df$gene
  ranked <- sort(ranked, decreasing = TRUE)

  cat(sprintf("[OK] ranked list: %d genes  top3: %s  bottom3: %s\n",
              length(ranked),
              paste(head(names(ranked), 3), collapse = ", "),
              paste(tail(names(ranked), 3), collapse = ", ")))
  invisible(ranked)
}

#' Filter DE results into up- and downregulated lists for ORA.
#'
#' Returns list with: up, down, background (all genes in the input table).
filter_significant_genes <- function(de,
                                     padj_thresh   = 0.05,
                                     log2fc_thresh = 1.0) {
  required <- c("gene", "log2FC", "padj")
  missing <- setdiff(required, colnames(de))
  if (length(missing) > 0) stop("Required columns missing: ", paste(missing, collapse = ", "))

  df <- de
  if (any(duplicated(df$gene))) {
    df <- df[order(df$padj), ]
    df <- df[!duplicated(df$gene), ]
  }

  sig  <- df[!is.na(df$padj) & df$padj <= padj_thresh &
             !is.na(df$log2FC) & abs(df$log2FC) >= log2fc_thresh, ]
  up   <- sig$gene[sig$log2FC > 0]
  down <- sig$gene[sig$log2FC < 0]
  bg   <- df$gene

  cat(sprintf("[OK] padj<=%g, |log2FC|>=%g: up=%d  down=%d  background=%d\n",
              padj_thresh, log2fc_thresh, length(up), length(down), length(bg)))
  invisible(list(up = up, down = down, background = bg))
}
