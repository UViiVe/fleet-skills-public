# multi_condition.R
#
# Purpose:
#   Demonstrate DESeq2 with a 3-level condition factor: extract two
#   coefficient-form comparisons (each treatment vs control) and one
#   contrast-form comparison (treatmentA vs treatmentB, neither of which
#   is the reference).
#
# Inputs:
#   None (simulated). Replace the data-generating block with your own
#   counts + coldata; keep the modeling pattern.
#
# Outputs:
#   Prints summaries; returns a named list of all three result tables.
#
# Example:
#   source("scripts/multi_condition.R")

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) stop("Install DESeq2 first")
  if (!requireNamespace("apeglm", quietly = TRUE)) stop("Install apeglm first")
  if (!requireNamespace("ashr",   quietly = TRUE))
    stop("Install with: install.packages('ashr')")
  library(DESeq2); library(apeglm); library(ashr)
})

run_demo <- function(n_genes = 1000, n_per_group = 3L, seed = 11L) {
  set.seed(seed)
  n_samples <- n_per_group * 3L
  counts <- matrix(
    rnbinom(n_genes * n_samples, mu = 100, size = 10),
    nrow     = n_genes,
    dimnames = list(
      paste0("gene", seq_len(n_genes)),
      paste0("sample", seq_len(n_samples))
    )
  )
  coldata <- data.frame(
    condition = factor(rep(c("control", "treatmentA", "treatmentB"),
                           each = n_per_group)),
    row.names = colnames(counts)
  )

  dds <- DESeqDataSetFromMatrix(counts, coldata, design = ~ condition)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds$condition <- relevel(dds$condition, ref = "control")
  dds <- DESeq(dds)

  cat("resultsNames(dds):\n"); print(resultsNames(dds))

  # Coefficient-form: each treatment vs the reference.
  res_A <- lfcShrink(dds, coef = "condition_treatmentA_vs_control", type = "apeglm")
  res_B <- lfcShrink(dds, coef = "condition_treatmentB_vs_control", type = "apeglm")

  # Contrast-form: A vs B (neither is the reference) -> ashr.
  res_AB <- lfcShrink(dds,
                      contrast = c("condition", "treatmentA", "treatmentB"),
                      type     = "ashr")

  cat("\n--- treatmentA vs control ---\n"); print(summary(res_A))
  cat("\n--- treatmentB vs control ---\n"); print(summary(res_B))
  cat("\n--- treatmentA vs treatmentB ---\n"); print(summary(res_AB))

  cat("[OK] 3-condition workflow demo complete\n")
  invisible(list(dds = dds, A_vs_control = res_A,
                 B_vs_control = res_B, A_vs_B = res_AB))
}

if (sys.nframe() == 0L) {
  run_demo()
}
