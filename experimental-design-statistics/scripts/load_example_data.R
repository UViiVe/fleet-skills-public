# load_example_data.R
# Provide example pilot data for testing the experimental-design-statistics workflow,
# and a loader for the bundled tissue-CV database.
#
# Inputs:
#   - None for load_example_data() (generates synthetic RNA-seq counts).
#   - For load_cv_database(): path to references/cv_tissue_database.csv.
# Outputs:
#   - load_example_data(): list of counts, metadata, DESeqDataSet, dispersions, CV summary.
#   - load_cv_database(): data frame.
# Example:
#   source("scripts/load_example_data.R")
#   pilot <- load_example_data()
#   cvs   <- load_cv_database("references/cv_tissue_database.csv")

options(repos = c(CRAN = "https://cloud.r-project.org"))

.require_or_install <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg)
    }
  }
}

#' Generate synthetic RNA-seq pilot data with realistic count statistics
#'
#' Builds an 8-sample, 2000-gene synthetic dataset suitable for testing
#' the power, sample-size, and batch-assignment workflows.
#'
#' @param n_genes Total simulated genes (default 2000)
#' @param n_per_group Replicates per condition (default 4)
#' @param de_fraction Fraction of genes to make differentially expressed (default 0.20)
#' @param seed Random seed for reproducibility
#' @return list with counts, metadata, dds (DESeqDataSet), dispersions, cv
load_example_data <- function(n_genes = 2000,
                              n_per_group = 4,
                              de_fraction = 0.20,
                              seed = 42L) {

  .require_or_install("DESeq2", bioc = TRUE)
  suppressPackageStartupMessages({
    library(DESeq2)
  })

  message("=== Generating synthetic pilot RNA-seq data ===")

  set.seed(seed)

  n_samples <- 2L * n_per_group

  metadata <- data.frame(
    sample_id = sprintf("sample_%02d", seq_len(n_samples)),
    condition = rep(c("control", "treated"), each = n_per_group),
    batch     = rep(c("A", "B"), times = n_per_group),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$sample_id

  # Negative-binomial counts with gene-specific mean and dispersion
  mu_per_gene   <- runif(n_genes, 50, 5000)
  size_per_gene <- runif(n_genes, 0.1, 2)

  counts <- matrix(
    rnbinom(n_genes * n_samples,
            mu = rep(mu_per_gene, n_samples),
            size = rep(size_per_gene, n_samples)),
    nrow = n_genes, ncol = n_samples
  )

  # Inject DE for de_fraction of genes
  de_genes <- sample.int(n_genes, size = round(de_fraction * n_genes))
  treated_cols <- seq.int(n_per_group + 1L, n_samples)
  for (g in de_genes) {
    fc <- sample(c(0.5, 2.0), 1)
    counts[g, treated_cols] <- round(counts[g, treated_cols] * fc)
  }

  rownames(counts) <- sprintf("gene_%05d", seq_len(n_genes))
  colnames(counts) <- metadata$sample_id

  message(sprintf("  samples: %d  |  genes: %d  |  conditions: %s",
                  n_samples, n_genes,
                  paste(unique(metadata$condition), collapse = ", ")))

  dds <- DESeqDataSetFromMatrix(countData = counts, colData = metadata,
                                design = ~ condition)
  message("  estimating dispersions...")
  dds <- DESeq(dds, quiet = TRUE)
  disp <- dispersions(dds)

  norm_counts <- counts(dds, normalized = TRUE)
  gene_means <- rowMeans(norm_counts)
  gene_sds   <- apply(norm_counts, 1, sd)
  gene_cvs   <- gene_sds / gene_means
  gene_cvs   <- gene_cvs[is.finite(gene_cvs) & gene_cvs > 0]

  cv_summary <- list(
    median   = median(gene_cvs, na.rm = TRUE),
    mean     = mean(gene_cvs, na.rm = TRUE),
    per_gene = gene_cvs
  )

  message(sprintf("  median CV = %.3f  |  median dispersion = %.4f",
                  cv_summary$median, median(disp, na.rm = TRUE)))
  cat("\n✓ Example pilot data loaded successfully!\n")

  invisible(list(
    counts      = counts,
    metadata    = metadata,
    dds         = dds,
    dispersions = disp,
    cv          = cv_summary
  ))
}


#' Load the bundled tissue-CV database
#'
#' @param path Path to references/cv_tissue_database.csv (default looks for it).
#' @param tissue Optional filter on tissue column.
#' @return data frame.
load_cv_database <- function(path = "references/cv_tissue_database.csv",
                             tissue = NULL) {
  if (!file.exists(path)) {
    stop(sprintf("CV database not found at: %s", path))
  }
  cv_db <- read.csv(path, stringsAsFactors = FALSE)

  if (!is.null(tissue)) {
    cv_db <- cv_db[tolower(cv_db$tissue) == tolower(tissue), , drop = FALSE]
    if (nrow(cv_db) == 0) {
      stop(sprintf("No rows for tissue: %s", tissue))
    }
  }

  cat(sprintf("✓ Loaded CV database (%d rows)\n", nrow(cv_db)))
  invisible(cv_db)
}


#' Smoke-test the workflow with example data
test_workflow <- function() {
  message("=== Testing experimental-design-statistics workflow ===")
  pilot <- load_example_data()
  message("\nNext steps:")
  message("  source('scripts/power_rnaseq.R')")
  message("  source('scripts/sample_size_de.R')")
  message("  source('scripts/batch_assignment.R')")
  invisible(pilot)
}
