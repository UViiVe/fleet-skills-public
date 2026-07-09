# extract_results.R
#
# Purpose:
#   Helpers for extracting DESeq2 results by coefficient or contrast, applying
#   LFC shrinkage, and filtering by significance.
#
# Inputs:
#   dds : fitted DESeqDataSet
#
# Outputs:
#   DESeqResults objects (returned, not written).
#
# Example:
#   source("scripts/extract_results.R")
#   show_coefficients(dds)
#   res <- extract_by_coefficient(dds, "condition_treated_vs_control")
#   res_shrunk <- apply_lfc_shrinkage(dds, "condition_treated_vs_control")
#   sig <- filter_significant(res, padj_cutoff = 0.05, lfc_cutoff = 1)

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Install with: BiocManager::install('DESeq2')")
  library(DESeq2)
})

#' List coefficients available in a fitted dds.
show_coefficients <- function(dds) {
  coefs <- resultsNames(dds)
  cat("Available coefficients:\n")
  for (i in seq_along(coefs)) cat(sprintf("  %d. %s\n", i, coefs[i]))
  invisible(coefs)
}

#' Pull results by coefficient name (works for any factor whose levels are
#' encoded as factor_level_vs_reference).
extract_by_coefficient <- function(dds, coef_name, alpha = 0.05) {
  cat(sprintf("[extract] coef = %s   alpha = %s\n", coef_name, alpha))
  res <- results(dds, name = coef_name, alpha = alpha)
  print(summary(res))
  invisible(res)
}

#' Pull results by contrast (works for any pair of factor levels, not just
#' those involving the reference).
extract_by_contrast <- function(dds, factor_name, numerator, denominator,
                                alpha = 0.05) {
  cat(sprintf("[extract] contrast = %s vs %s on factor '%s'   alpha = %s\n",
              numerator, denominator, factor_name, alpha))
  res <- results(dds, contrast = c(factor_name, numerator, denominator),
                 alpha = alpha)
  print(summary(res))
  invisible(res)
}

#' Apply LFC shrinkage. apeglm is the default; falls back to ashr if the
#' caller passes contrast = instead of coef =.
apply_lfc_shrinkage <- function(dds, coef_name = NULL,
                                contrast = NULL, shrink_type = "apeglm") {
  if (is.null(coef_name) && is.null(contrast))
    stop("Provide coef_name (preferred) or contrast")

  if (!is.null(contrast) && shrink_type == "apeglm") {
    message("apeglm only supports coef=. Switching to ashr.")
    shrink_type <- "ashr"
  }

  if (shrink_type == "apeglm" && !requireNamespace("apeglm", quietly = TRUE))
    stop("Install with: BiocManager::install('apeglm')")
  if (shrink_type == "ashr"   && !requireNamespace("ashr",   quietly = TRUE))
    stop("Install with: install.packages('ashr')")

  cat(sprintf("[shrink] method = %s\n", shrink_type))
  if (!is.null(coef_name)) {
    res <- lfcShrink(dds, coef = coef_name, type = shrink_type)
  } else {
    res <- lfcShrink(dds, contrast = contrast, type = shrink_type)
  }
  cat("    use shrunk LFC for plots/ranking; use unshrunken for hypothesis tests\n")
  invisible(res)
}

#' Filter a DESeqResults by significance and absolute log2FC.
filter_significant <- function(res, padj_cutoff = 0.05, lfc_cutoff = 1) {
  cat(sprintf("[filter] padj < %s and |log2FC| > %s\n",
              padj_cutoff, lfc_cutoff))
  hits <- subset(res, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
  cat(sprintf("    %d hits  (up=%d, down=%d)\n",
              nrow(hits),
              sum(hits$log2FoldChange > 0, na.rm = TRUE),
              sum(hits$log2FoldChange < 0, na.rm = TRUE)))
  invisible(hits)
}

#' Top N by adjusted p-value.
get_top_genes <- function(res, n = 100) {
  ordered <- res[order(res$padj), ]
  head(ordered, n)
}

#' Convenience: run the standard extract -> shrink -> filter sequence and
#' return all three as a named list.
standard_extraction <- function(dds, coef_name = NULL,
                                apply_shrinkage = TRUE,
                                padj_cutoff = 0.05, lfc_cutoff = 1) {
  if (is.null(coef_name)) {
    show_coefficients(dds)
    stop("Provide coef_name from the list above")
  }
  res <- extract_by_coefficient(dds, coef_name, alpha = padj_cutoff)
  out <- list(results = res)
  if (apply_shrinkage) {
    out$shrunk <- apply_lfc_shrinkage(dds, coef_name = coef_name,
                                      shrink_type = "apeglm")
  }
  out$significant <- filter_significant(res,
                                        padj_cutoff = padj_cutoff,
                                        lfc_cutoff  = lfc_cutoff)
  cat("[OK] extraction complete\n")
  invisible(out)
}
