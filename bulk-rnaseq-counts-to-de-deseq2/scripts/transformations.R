# transformations.R
#
# Purpose:
#   Variance-stabilizing transformations for downstream visualization
#   and clustering. Never feed these back into DESeq() for testing.
#
# Inputs:
#   dds : fitted DESeqDataSet
#
# Outputs:
#   DESeqTransform objects (returned, not written).
#
# Example:
#   source("scripts/transformations.R")
#   vsd <- transform_counts(dds, method = "auto")
#   mat <- get_transformed_matrix(vsd)

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Install with: BiocManager::install('DESeq2')")
  library(DESeq2)
})

#' Size-factor-normalized counts (no log).
get_normalized_counts <- function(dds) {
  if (is.null(sizeFactors(dds)))
    stop("dds has no size factors; run DESeq() first")
  cat("normalized counts shape:",
      nrow(counts(dds)), "genes x", ncol(counts(dds)), "samples\n")
  invisible(counts(dds, normalized = TRUE))
}

#' VST: fast, recommended for n > 30.
apply_vst <- function(dds, blind = FALSE) {
  cat(sprintf("vst(blind = %s)  n=%d samples\n", blind, ncol(dds)))
  invisible(vst(dds, blind = blind))
}

#' rlog: slower, better at very low counts and very small n.
apply_rlog <- function(dds, blind = FALSE) {
  if (ncol(dds) > 100)
    warning("rlog can be slow for n > 100; consider vst() instead")
  cat(sprintf("rlog(blind = %s)  n=%d samples\n", blind, ncol(dds)))
  invisible(rlog(dds, blind = blind))
}

#' Auto-pick vst or rlog by sample size.
transform_counts <- function(dds,
                             method = c("auto", "vst", "rlog"),
                             blind  = FALSE) {
  method <- match.arg(method)
  if (method == "auto") method <- if (ncol(dds) > 30) "vst" else "rlog"
  if (method == "vst")  return(apply_vst (dds, blind = blind))
  if (method == "rlog") return(apply_rlog(dds, blind = blind))
  stop("Unreachable")
}

#' Pull the numeric matrix out of a DESeqTransform.
get_transformed_matrix <- function(transformed) {
  assay(transformed)
}

#' Compare vst and rlog side by side on the same dds. Returns a list of both
#' matrices and prints the correlation between the first two samples.
compare_transformations <- function(dds) {
  cat("computing vst...\n");  vsd <- vst(dds, blind = FALSE)
  cat("computing rlog...\n"); rld <- rlog(dds, blind = FALSE)
  m_v <- assay(vsd);  m_r <- assay(rld)
  if (ncol(m_v) >= 2L) {
    cat(sprintf("Pearson correlation between vst and rlog (sample 1): %.3f\n",
                cor(m_v[, 1], m_r[, 1])))
  }
  invisible(list(vsd = vsd, rld = rld))
}
