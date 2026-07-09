# batch_correction.R
#
# Purpose:
#   Demonstrate a batch-corrected DESeq2 pipeline on simulated data:
#   include the batch term in the design, fit, and extract the
#   condition effect while controlling for batch.
#
# Inputs:
#   None (uses simulated data). Adapt the data-loading section to use
#   your own counts and metadata; the modeling part is the part to keep.
#
# Outputs:
#   Prints fit summary; returns the fitted dds and shrunken results.
#
# Example:
#   source("scripts/batch_correction.R")

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 first")
  if (!requireNamespace("apeglm", quietly = TRUE)) stop("Install apeglm first")
  library(DESeq2); library(apeglm)
})

run_demo <- function(n_genes = 1000, n_per_group = 4, seed = 7L) {
  set.seed(seed)
  n_samples <- n_per_group * 2 * 3   # 3 batches, 2 conditions per batch

  counts <- matrix(
    rnbinom(n_genes * n_samples, mu = 100, size = 10),
    nrow = n_genes,
    dimnames = list(
      paste0("gene", seq_len(n_genes)),
      paste0("sample", seq_len(n_samples))
    )
  )

  coldata <- data.frame(
    condition = factor(rep(c("control", "treated"), times = n_samples / 2)),
    batch     = factor(rep(c("b1", "b2", "b3"),    each  = n_samples / 3)),
    row.names = colnames(counts)
  )

  # Confounding check before building dds
  tab <- table(coldata$batch, coldata$condition)
  if (any(tab == 0)) stop("Simulated coldata is confounded; raise n_per_group")
  cat("batch x condition cross-tab:\n"); print(tab)

  dds <- DESeqDataSetFromMatrix(
    countData = counts, colData = coldata,
    design    = ~ batch + condition
  )
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds$condition <- relevel(dds$condition, ref = "control")
  dds <- DESeq(dds)

  cat("\nresultsNames(dds):\n"); print(resultsNames(dds))

  res        <- results  (dds, name = "condition_treated_vs_control")
  res_shrunk <- lfcShrink(dds, coef = "condition_treated_vs_control",
                          type = "apeglm")
  print(summary(res))
  cat(sprintf("[OK] batch-corrected fit complete: %d genes tested\n", nrow(res)))
  invisible(list(dds = dds, res = res, res_shrunk = res_shrunk))
}

if (sys.nframe() == 0L) {
  run_demo()
}
