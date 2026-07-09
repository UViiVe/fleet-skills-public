# load_example_data.R
#
# Purpose:
#   Provides loaders for two well-known example RNA-seq count datasets
#   (pasilla, airway) and a validator for user-supplied count + metadata.
#
# Inputs:
#   - load_pasilla_data() and load_airway_data() take no arguments;
#     they install the relevant Bioconductor data packages if missing.
#   - validate_input_data(counts, coldata, design_col) checks user-supplied
#     count matrices against DESeq2's requirements.
#
# Outputs:
#   - load_*_data() return list(counts = matrix, coldata = data.frame,
#     description = character).
#   - validate_input_data() returns the same list with counts re-aligned
#     to coldata row order; stops with an informative message on failure.
#
# Example:
#   source("scripts/load_example_data.R")
#   ds <- load_pasilla_data()
#   counts  <- ds$counts
#   coldata <- ds$coldata

.ensure_cran_mirror <- function() {
  current <- getOption("repos")
  if (length(current) == 0L || identical(unname(current["CRAN"]), "@CRAN@")) {
    options(repos = c(CRAN = "https://cloud.r-project.org"))
  }
  invisible(NULL)
}

.ensure_bioc_pkg <- function(pkg) {
  .ensure_cran_mirror()
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager...")
    install.packages("BiocManager")
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing %s ...", pkg))
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
  invisible(NULL)
}

#' Validate a user-supplied count matrix + sample metadata.
#'
#' @param counts integer matrix or data.frame, genes in rows, samples in columns.
#' @param coldata data.frame with row.names matching colnames(counts).
#' @param design_col character; name of the factor column in coldata used
#'        as the primary comparison.
#' @return list(counts = matrix aligned to coldata, coldata = aligned data.frame).
validate_input_data <- function(counts, coldata, design_col = "condition") {
  problems <- character(0)

  if (missing(counts)  || is.null(counts))  stop("`counts` is required")
  if (missing(coldata) || is.null(coldata)) stop("`coldata` is required")
  if (!is.matrix(counts) && !is.data.frame(counts))
    problems <- c(problems, "`counts` must be a matrix or data.frame")
  if (!is.data.frame(coldata))
    problems <- c(problems, "`coldata` must be a data.frame")

  if (is.null(colnames(counts)) || is.null(rownames(coldata)))
    problems <- c(problems, "`counts` needs colnames; `coldata` needs rownames")

  if (length(problems) == 0L) {
    missing_in_meta <- setdiff(colnames(counts), rownames(coldata))
    if (length(missing_in_meta) > 0L) {
      problems <- c(problems, paste0(
        "These count columns have no matching row in coldata: ",
        paste(head(missing_in_meta, 5), collapse = ", "),
        if (length(missing_in_meta) > 5) " ..." else ""
      ))
    }
    if (!design_col %in% colnames(coldata)) {
      problems <- c(problems, sprintf(
        "Design column '%s' not in coldata. Available: %s",
        design_col, paste(colnames(coldata), collapse = ", ")
      ))
    }
    if (any(is.na(counts)))
      problems <- c(problems, "`counts` contains NA (not allowed by DESeq2)")
    if (any(as.matrix(counts) < 0, na.rm = TRUE))
      problems <- c(problems, "`counts` contains negative values")
  }

  if (length(problems) > 0L) {
    stop("Input validation failed:\n  - ",
         paste(problems, collapse = "\n  - "))
  }

  coldata <- coldata[colnames(counts), , drop = FALSE]
  if (is.data.frame(counts)) counts <- as.matrix(counts)

  cat(sprintf("[OK] %d genes x %d samples; design column: %s\n",
              nrow(counts), ncol(counts), design_col))
  if (is.factor(coldata[[design_col]]) || is.character(coldata[[design_col]])) {
    cat("    Groups:", paste(unique(coldata[[design_col]]), collapse = ", "), "\n")
  }
  invisible(list(counts = counts, coldata = coldata))
}

#' Load the pasilla example dataset.
#'
#' Brooks et al. (2011), Drosophila pasilla RNAi knockdown vs untreated.
#' Strips the "fb" suffix from metadata rownames so they match count columns.
load_pasilla_data <- function() {
  .ensure_bioc_pkg("pasilla")
  pasilla_dir <- system.file("extdata", package = "pasilla", mustWork = TRUE)

  counts <- as.matrix(read.table(
    file.path(pasilla_dir, "pasilla_gene_counts.tsv"),
    header = TRUE, row.names = 1, check.names = FALSE
  ))
  coldata <- read.csv(
    file.path(pasilla_dir, "pasilla_sample_annotation.csv"),
    row.names = 1
  )
  rownames(coldata) <- sub("fb$", "", rownames(coldata))
  coldata$condition <- factor(coldata$condition)
  coldata$type      <- factor(coldata$type)
  coldata <- coldata[colnames(counts), , drop = FALSE]

  cat(sprintf("[OK] pasilla loaded: %d genes x %d samples\n",
              nrow(counts), ncol(counts)))
  cat("    condition: treated / untreated;  type: single-read / paired-end\n")
  invisible(list(
    counts      = counts,
    coldata     = coldata,
    description = "Drosophila pasilla RNAi knockdown (Brooks et al. 2011)"
  ))
}

#' Load the airway example dataset.
#'
#' Himes et al. (2014), human airway smooth muscle cells, dexamethasone vs untreated.
load_airway_data <- function() {
  .ensure_bioc_pkg("airway")
  # SummarizedExperiment is a dependency of airway
  airway <- NULL
  data("airway", package = "airway", envir = environment())
  if (is.null(airway)) stop("airway dataset failed to load")

  counts  <- SummarizedExperiment::assay(airway)
  coldata <- as.data.frame(SummarizedExperiment::colData(airway))
  coldata <- coldata[, c("cell", "dex"), drop = FALSE]
  coldata$cell <- factor(coldata$cell)
  coldata$dex  <- factor(coldata$dex)

  cat(sprintf("[OK] airway loaded: %d genes x %d samples\n",
              nrow(counts), ncol(counts)))
  cat("    dex: trt / untrt;  cell: 4 cell lines  (consider ~ cell + dex)\n")
  invisible(list(
    counts      = counts,
    coldata     = coldata,
    description = "Human airway smooth muscle cells, dexamethasone (Himes et al. 2014)"
  ))
}
