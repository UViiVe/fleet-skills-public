# sample_size_scrna.R
# Sample-size guidance for single-cell RNA-seq experiments.
# Two-level variation: between-sample (biological) and within-sample (cell-to-cell).
# The cardinal rule: cells per sample don't substitute for samples.
#
# Inputs:
#   - Optional SingleCellExperiment pilot, fold-change, power, cell-type proportion.
# Outputs:
#   - Recommended samples per group, cells per sample, and total cells to sequence.

.try_pkg <- function(pkg) requireNamespace(pkg, quietly = TRUE)

#' Estimate pilot scRNA-seq parameters
#'
#' @param pilot_sce A SingleCellExperiment with a `sample` column in colData
#' @return list with dropout, mean expression, median CV, sample/cell counts
estimate_scrna_params <- function(pilot_sce) {
  if (!.try_pkg("SingleCellExperiment"))
    stop("Install SingleCellExperiment with BiocManager::install('SingleCellExperiment')")
  if (!methods::is(pilot_sce, "SingleCellExperiment"))
    stop("pilot_sce must be a SingleCellExperiment")

  cm <- SingleCellExperiment::counts(pilot_sce)

  dropout <- sum(cm == 0) / length(cm)
  gene_means <- rowMeans(cm)
  mean_expr <- mean(gene_means[gene_means > 0])

  cvs <- apply(cm, 1, function(x) {
    mu <- mean(x); if (mu > 0) sd(x) / mu else NA_real_
  })
  med_cv <- median(cvs, na.rm = TRUE)

  n_cells <- ncol(pilot_sce)
  n_samples <- if ("sample" %in% colnames(SingleCellExperiment::colData(pilot_sce)))
    length(unique(pilot_sce$sample)) else NA_integer_

  list(
    dropout_rate     = dropout,
    mean_expression  = mean_expr,
    median_cv        = med_cv,
    n_cells          = n_cells,
    n_samples        = n_samples,
    n_genes          = nrow(cm),
    interpretation   = .scrna_param_msg(dropout, med_cv, n_cells)
  )
}

#' Sample-size recommendation for scRNA-seq
#'
#' If `powsimR` is installed it may be wrapped externally; this function returns
#' a heuristic recommendation based on Vieth et al. (2017) and assay defaults.
#'
#' @param pilot_sce Optional SingleCellExperiment for parameter extraction
#' @param fold_change Minimum fold-change to detect (default 2)
#' @param power Target power (default 0.8)
#' @param min_cells_per_sample Minimum cells to sequence per sample (default 500)
#' @param target_cell_type_prop Target cell-type proportion for rare-type studies
#' @param dropout_rate Default if no pilot (0.6)
#' @param mean_expression Default if no pilot (10)
samplesize_scrna <- function(pilot_sce = NULL,
                             fold_change = 2,
                             power = 0.8,
                             min_cells_per_sample = 500,
                             target_cell_type_prop = NULL,
                             dropout_rate = 0.6,
                             mean_expression = 10) {

  if (!is.null(pilot_sce)) {
    params <- estimate_scrna_params(pilot_sce)
    dropout_rate    <- params$dropout_rate
    mean_expression <- params$mean_expression
  }

  base_n <- if (fold_change >= 3) 3 else if (fold_change >= 2) 6 else 10
  if (power >= 0.90) base_n <- base_n + 2

  if (!is.null(target_cell_type_prop)) {
    if (target_cell_type_prop < 0.05) {
      base_n <- base_n + 4
      min_cells_per_sample <- max(min_cells_per_sample, 2000)
    } else if (target_cell_type_prop < 0.15) {
      base_n <- base_n + 2
      min_cells_per_sample <- max(min_cells_per_sample, 1000)
    }
  }

  if (dropout_rate > 0.7) {
    base_n <- base_n + 2
    min_cells_per_sample <- min_cells_per_sample + 200
  }

  out <- list(
    recommended_n_per_group     = base_n,
    recommended_cells_per_sample = min_cells_per_sample,
    total_samples               = 2L * base_n,
    total_cells                 = 2L * base_n * min_cells_per_sample,
    parameters = list(fold_change = fold_change, power = power,
                      dropout_rate = dropout_rate,
                      mean_expression = mean_expression,
                      target_cell_type_prop = target_cell_type_prop),
    recommendation = .scrna_design_msg(base_n, min_cells_per_sample, target_cell_type_prop)
  )
  cat(sprintf("✓ scRNA-seq design: %d samples/group, ~%d cells/sample, %s total cells\n",
              base_n, min_cells_per_sample, format(out$total_cells, big.mark = ",")))
  out
}

#' Quick reference table for scRNA-seq sample sizes by study type
quick_reference_scrna <- function() {
  data.frame(
    study_type = c("Cell-type discovery",
                   "DE between conditions (common cell types)",
                   "DE between conditions (rare cell types <5%)",
                   "Trajectory inference",
                   "Spatial integration"),
    min_samples_per_group = c(3, 6, 10, 5, 4),
    recommended_samples_per_group = c(5, 8, 12, 8, 6),
    cells_per_sample = c("1000-2000", "500-1000", "2000-5000", "1000-2000", "1000-2000"),
    note = c("Prioritize diversity", "Balance samples vs cells",
             "More samples critical", "Continuous sampling",
             "Pilot for integration"),
    stringsAsFactors = FALSE
  )
}

#' Compute cells per sample to capture a rare cell type
#'
#' @param cell_type_proportion Expected fraction of target cell type (e.g., 0.01)
#' @param min_cells_needed Minimum target cells required (default 50)
#' @param samples_per_group Samples per condition (default 6)
cells_for_rare_type <- function(cell_type_proportion,
                                min_cells_needed = 50,
                                samples_per_group = 6) {
  cells_per_sample <- ceiling(min_cells_needed /
                                (cell_type_proportion * samples_per_group))
  cells_per_sample <- ceiling(cells_per_sample * 1.2)  # 20% buffer

  list(
    recommended_cells_per_sample = cells_per_sample,
    expected_target_per_sample   = ceiling(cells_per_sample * cell_type_proportion),
    total_cells_to_sequence      = cells_per_sample * samples_per_group * 2,
    expected_total_target_cells  = ceiling(cells_per_sample * cell_type_proportion *
                                             samples_per_group * 2),
    note = sprintf("To obtain ~%d cells from a %.2f%% population",
                   min_cells_needed, cell_type_proportion * 100)
  )
}

.scrna_param_msg <- function(dropout, cv, n_cells) {
  parts <- sprintf("dropout=%.1f%%, CV=%.2f, cells=%d. ",
                   dropout * 100, cv, n_cells)
  if (dropout > 0.8)  parts <- paste0(parts, "High dropout — deeper sequencing or more cells. ")
  if (cv > 1.5)        parts <- paste0(parts, "High cell-to-cell variability — more samples needed. ")
  if (n_cells < 500)   parts <- paste0(parts, "Pilot is small — estimates uncertain.")
  parts
}

.scrna_design_msg <- function(n_samples, n_cells, rare_prop) {
  base <- sprintf("scRNA-seq design: %d samples/group, ~%d cells/sample (total %s cells). ",
                  n_samples, n_cells,
                  format(n_samples * 2 * n_cells, big.mark = ","))
  if (!is.null(rare_prop) && rare_prop < 0.10)
    base <- paste0(base, "Rare cell type — prioritize more samples over more cells. ")
  paste0(base, "Key principle: BIOLOGICAL SAMPLES > CELLS/SAMPLE.")
}
