# export_results.R
#
# Purpose:
#   Save DESeq2 outputs to disk: full results table, shrunk results,
#   significant-only table, normalized counts, transformed counts,
#   top-N gene list, and the dds RDS for later re-extraction.
#
# Inputs:
#   dds        : fitted DESeqDataSet
#   res        : DESeqResults (unshrunken)
#   res_shrunk : DESeqResults with shrinkage (optional)
#
# Outputs (default output_dir = "results"):
#   deseq2_results.csv
#   deseq2_results_shrunk.csv          (if res_shrunk provided)
#   deseq2_significant.csv             (padj + |LFC| filter)
#   normalized_counts.csv
#   vst_transformed.csv or rlog_transformed.csv
#   top100_genes.csv
#   dds_object.rds
#
# Example:
#   source("scripts/export_results.R")
#   export_all(dds, res, res_shrunk, output_dir = "results")

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Install with: BiocManager::install('DESeq2')")
  library(DESeq2)
})

.write_results_csv <- function(res, path) {
  df <- as.data.frame(res)
  df$gene_id <- rownames(df)
  df <- df[, c("gene_id", setdiff(colnames(df), "gene_id"))]
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("    wrote %s  (%d rows)\n", path, nrow(df)))
}

#' Write the full DE results table (one row per gene).
export_results <- function(res, output_file = "deseq2_results.csv") {
  .write_results_csv(res, output_file)
}

#' Write a filtered table containing only significant genes.
export_significant <- function(res, padj_threshold = 0.05,
                               lfc_threshold = 0,
                               output_file = "deseq2_significant.csv") {
  sig <- subset(res, padj < padj_threshold & abs(log2FoldChange) > lfc_threshold)
  .write_results_csv(sig, output_file)
  cat(sprintf("    padj<%s and |LFC|>%s   up=%d down=%d\n",
              padj_threshold, lfc_threshold,
              sum(sig$log2FoldChange > 0, na.rm = TRUE),
              sum(sig$log2FoldChange < 0, na.rm = TRUE)))
}

#' Write size-factor-normalized counts.
export_normalized_counts <- function(dds, output_file = "normalized_counts.csv") {
  mat <- counts(dds, normalized = TRUE)
  df  <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE)
  write.csv(df, output_file, row.names = FALSE)
  cat(sprintf("    wrote %s\n", output_file))
}

#' Write transformed counts (vst or rlog).
export_transformed_counts <- function(dds,
                                      output_file = NULL,
                                      method = c("auto", "vst", "rlog")) {
  method <- match.arg(method)
  if (method == "auto") method <- if (ncol(dds) > 30) "vst" else "rlog"
  if (is.null(output_file)) {
    output_file <- if (method == "vst") "vst_transformed.csv" else "rlog_transformed.csv"
  }
  transformed <- if (method == "vst") vst(dds, blind = FALSE) else rlog(dds, blind = FALSE)
  mat <- assay(transformed)
  df  <- data.frame(gene_id = rownames(mat), mat, check.names = FALSE)
  write.csv(df, output_file, row.names = FALSE)
  cat(sprintf("    wrote %s (%s, n=%d)\n", output_file, method, ncol(dds)))
}

#' Top N rows ordered by adjusted p-value.
export_top_genes <- function(res, n = 100, output_file = NULL) {
  if (is.null(output_file)) output_file <- sprintf("top%d_genes.csv", n)
  ord <- res[order(res$padj), ]
  .write_results_csv(head(ord, n), output_file)
}

#' Write everything you typically need from a DESeq2 run in one call.
export_all <- function(dds, res, res_shrunk = NULL,
                       output_dir = "results",
                       padj_threshold = 0.05,
                       lfc_threshold  = 1) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat(sprintf("=== exporting to %s/ ===\n", output_dir))

  cat("[1] full results\n")
  export_results(res, file.path(output_dir, "deseq2_results.csv"))

  if (!is.null(res_shrunk)) {
    cat("[2] shrunk results\n")
    export_results(res_shrunk, file.path(output_dir, "deseq2_results_shrunk.csv"))
  }

  cat("[3] significant\n")
  export_significant(res,
                     padj_threshold = padj_threshold,
                     lfc_threshold  = lfc_threshold,
                     output_file    = file.path(output_dir, "deseq2_significant.csv"))

  cat("[4] normalized counts\n")
  export_normalized_counts(dds, file.path(output_dir, "normalized_counts.csv"))

  cat("[5] dds object\n")
  saveRDS(dds, file.path(output_dir, "dds_object.rds"))
  cat(sprintf("    wrote %s/dds_object.rds  (load with readRDS())\n", output_dir))

  cat("[6] transformed counts\n")
  method   <- if (ncol(dds) > 30) "vst" else "rlog"
  out_path <- file.path(output_dir,
                        if (method == "vst") "vst_transformed.csv" else "rlog_transformed.csv")
  export_transformed_counts(dds, output_file = out_path, method = method)

  cat("[7] top 100\n")
  export_top_genes(res, n = 100,
                   output_file = file.path(output_dir, "top100_genes.csv"))

  cat(sprintf("[OK] export complete: see %s/\n", output_dir))
}
