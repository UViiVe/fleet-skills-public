# load_de_results.R
#
# Purpose:
#   Read a DE results CSV/TSV from DESeq2, edgeR, or limma, auto-detect
#   the source by its column names, and standardize to a common schema:
#       gene, log2FC, pvalue, padj, stat (optional), baseMean (optional)
#
# Inputs:
#   path : path to a CSV or TSV with a header row.
#
# Outputs:
#   data.frame with the standardized columns; prints which source tool
#   was detected.
#
# Example:
#   source("scripts/load_de_results.R")
#   de <- load_de_results("deseq2_results.csv")

load_de_results <- function(path) {
  sep <- if (grepl("\\.tsv$|\\.tab$|\\.txt$", path, ignore.case = TRUE)) "\t" else ","
  df  <- read.table(path, sep = sep, header = TRUE,
                    check.names = FALSE, stringsAsFactors = FALSE)

  # Identify the gene-id column
  gene_col <- intersect(c("gene", "gene_id", "Gene", "GeneID", "ID", "Geneid", "symbol"),
                        colnames(df))
  if (length(gene_col) > 0) {
    df$gene <- df[[gene_col[1]]]
  } else if (!is.null(rownames(df))) {
    df$gene <- rownames(df)
  } else {
    stop("Could not identify a gene-id column in the input")
  }

  cols <- colnames(df)
  # Detect tool and remap columns
  tool <- if (all(c("log2FoldChange", "padj") %in% cols)) {
    "DESeq2"
  } else if (all(c("logFC", "FDR") %in% cols)) {
    "edgeR"
  } else if (all(c("logFC", "adj.P.Val") %in% cols)) {
    "limma"
  } else if (all(c("log2FC", "padj") %in% cols)) {
    "standardized"
  } else {
    stop("Could not identify the source tool from columns: ",
         paste(cols, collapse = ", "))
  }

  out <- data.frame(gene = df$gene, stringsAsFactors = FALSE)
  out$log2FC <- switch(tool,
    "DESeq2"       = df$log2FoldChange,
    "edgeR"        = df$logFC,
    "limma"        = df$logFC,
    "standardized" = df$log2FC
  )
  out$pvalue <- switch(tool,
    "DESeq2"       = df$pvalue,
    "edgeR"        = df$PValue,
    "limma"        = df$P.Value,
    "standardized" = df$pvalue
  )
  out$padj <- switch(tool,
    "DESeq2"       = df$padj,
    "edgeR"        = df$FDR,
    "limma"        = df$adj.P.Val,
    "standardized" = df$padj
  )

  # Optional: test statistic
  if ("stat"  %in% cols) out$stat <- df$stat
  else if ("t" %in% cols) out$stat <- df$t
  else if ("F" %in% cols) out$stat <- df$F

  # Optional: base mean
  if ("baseMean" %in% cols) out$baseMean <- df$baseMean
  else if ("AveExpr" %in% cols) out$baseMean <- df$AveExpr
  else if ("logCPM"  %in% cols) out$baseMean <- df$logCPM

  cat(sprintf("[OK] loaded %d rows from %s (source: %s)\n",
              nrow(out), path, tool))
  invisible(out)
}
