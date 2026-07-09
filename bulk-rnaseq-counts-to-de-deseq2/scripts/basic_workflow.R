# basic_workflow.R
#
# Purpose:
#   End-to-end DESeq2 pipeline: validate -> construct -> filter -> fit ->
#   extract -> shrink -> summarize. Creates dds, res, res_shrunk in the
#   calling environment.
#
# Inputs (from calling environment):
#   counts  : integer matrix or data.frame (genes x samples)
#   coldata : data.frame with row.names matching colnames(counts), and a
#             column named `condition` (factor with >= 2 levels).
#   design_formula (optional, default ~ condition)
#   reference_level (optional, default first alphabetical level)
#
# Outputs (created in calling environment):
#   dds        : fitted DESeqDataSet
#   res        : DESeqResults (unshrunken)
#   res_shrunk : DESeqResults with apeglm shrinkage
#
# Example:
#   source("scripts/load_example_data.R"); ds <- load_pasilla_data()
#   counts <- ds$counts; coldata <- ds$coldata[, "condition", drop = FALSE]
#   source("scripts/basic_workflow.R")

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE))
    stop("Install with: BiocManager::install('DESeq2')")
  if (!requireNamespace("apeglm", quietly = TRUE))
    stop("Install with: BiocManager::install('apeglm')")
  library(DESeq2)
  library(apeglm)
})

if (!exists("counts")  || is.null(counts))  stop("`counts` must be defined before sourcing this script")
if (!exists("coldata") || is.null(coldata)) stop("`coldata` must be defined before sourcing this script")
if (!exists("design_formula")) design_formula <- ~ condition

# Align metadata to counts and validate
if (!all(colnames(counts) %in% rownames(coldata))) {
  missing_in_meta <- setdiff(colnames(counts), rownames(coldata))
  stop("Sample ID mismatch. Missing in coldata: ",
       paste(head(missing_in_meta, 5), collapse = ", "))
}
coldata <- coldata[colnames(counts), , drop = FALSE]
if (!is.matrix(counts)) counts <- as.matrix(counts)

cat(sprintf("[1/7] data: %d genes x %d samples\n", nrow(counts), ncol(counts)))
if ("condition" %in% colnames(coldata)) {
  cat("      group sizes:",
      paste(sprintf("%s=%d", names(table(coldata$condition)), table(coldata$condition)),
            collapse = ", "), "\n")
}

# Construct DESeqDataSet
cat("[2/7] constructing DESeqDataSet with design ",
    deparse(design_formula), "\n", sep = "")
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design    = design_formula)

# Pre-filter low-count genes
pre_n <- nrow(dds)
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep, ]
cat(sprintf("[3/7] pre-filter: kept %d of %d genes (rowSums >= 10)\n",
            nrow(dds), pre_n))

# Set reference level
if ("condition" %in% colnames(colData(dds))) {
  if (!exists("reference_level")) reference_level <- levels(dds$condition)[1]
  dds$condition <- relevel(dds$condition, ref = reference_level)
  cat(sprintf("[4/7] reference level: '%s'\n", reference_level))
  cat(sprintf("      positive log2FC means higher in '%s'\n",
              setdiff(levels(dds$condition), reference_level)[1]))
}

# Fit
cat("[5/7] running DESeq() (size factors + dispersions + Wald)\n")
dds <- DESeq(dds)

# Extract unshrunken results
coef_name <- resultsNames(dds)[2]
cat(sprintf("[6/7] extracting results for coefficient: %s\n", coef_name))
res <- results(dds, name = coef_name)

# Shrink (apeglm)
cat("[7/7] applying apeglm LFC shrinkage\n")
res_shrunk <- lfcShrink(dds, coef = coef_name, type = "apeglm")

# Summarize
cat("\n--- summary ---\n")
print(summary(res, alpha = 0.05))

n_sig <- sum(res$padj < 0.05 & abs(res$log2FoldChange) >= 1, na.rm = TRUE)
n_up  <- sum(res$padj < 0.05 & res$log2FoldChange >=  1, na.rm = TRUE)
n_dn  <- sum(res$padj < 0.05 & res$log2FoldChange <= -1, na.rm = TRUE)

cat(sprintf("\n[OK] DESeq2 workflow completed: %d significant (padj<0.05, |LFC|>=1) | up=%d down=%d\n",
            n_sig, n_up, n_dn))
cat("     Next: run scripts/qc_plots.R then scripts/export_results.R\n")

invisible(list(dds = dds, res = res, res_shrunk = res_shrunk))
