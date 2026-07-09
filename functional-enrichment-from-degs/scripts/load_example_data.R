# load_example_data.R
#
# Purpose:
#   Generate example DE results using the airway dataset and DESeq2.
#   Bundled for testing the enrichment pipeline end-to-end.
#
# Inputs:
#   None.
#
# Outputs:
#   data.frame with columns gene, log2FC, pvalue, padj, stat, baseMean.
#
# Example:
#   source("scripts/load_example_data.R")
#   de <- load_example_de_results()

.ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (length(getOption("repos")) == 0 || identical(unname(getOption("repos")["CRAN"]), "@CRAN@"))
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  invisible(NULL)
}

load_example_de_results <- function() {
  for (p in c("DESeq2", "airway", "AnnotationDbi", "org.Hs.eg.db")) .ensure_pkg(p)
  suppressPackageStartupMessages({
    library(DESeq2); library(airway)
    library(AnnotationDbi); library(org.Hs.eg.db)
  })

  airway <- NULL
  data("airway", package = "airway", envir = environment())
  dds <- DESeqDataSet(airway, design = ~ cell + dex)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds$dex <- relevel(dds$dex, ref = "untrt")
  dds <- DESeq(dds)
  res <- results(dds, name = "dex_trt_vs_untrt")

  # Map Ensembl ids -> symbols
  ens <- sub("\\..*$", "", rownames(res))
  sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL",
                keytype = "ENSEMBL", multiVals = "first")

  out <- data.frame(
    gene     = ifelse(is.na(sym) | sym == "", ens, sym),
    log2FC   = res$log2FoldChange,
    pvalue   = res$pvalue,
    padj     = res$padj,
    stat     = res$stat,
    baseMean = res$baseMean,
    stringsAsFactors = FALSE
  )
  # Drop rows with NA padj (independent filtering removed them) and
  # duplicate symbols (keep best padj).
  out <- out[!is.na(out$padj), ]
  out <- out[order(out$padj), ]
  out <- out[!duplicated(out$gene), ]

  cat(sprintf("[OK] airway DE results loaded: %d genes, %d with padj < 0.05\n",
              nrow(out), sum(out$padj < 0.05)))
  invisible(out)
}
