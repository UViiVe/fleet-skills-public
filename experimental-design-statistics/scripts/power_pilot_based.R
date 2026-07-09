# power_pilot_based.R
# Pilot-driven power and sample size. Extracts dispersion from a DESeq2 object
# and translates it into a CV that RNASeqPower can consume.
#
# Inputs:
#   - pilot_dds: a DESeqDataSet that has been through DESeq() (dispersions estimated).
# Outputs:
#   - estimated power or required n based on the pilot's variability structure.

.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Package required: %s. Install via BiocManager::install('%s').",
                 pkg, pkg), call. = FALSE)
}

#' Summarize dispersion from a DESeq2 object
#'
#' @param pilot_dds A DESeqDataSet (DESeq2 v1.40+)
#' @return list with summary statistics
extract_dispersion_deseq2 <- function(pilot_dds) {
  .need("DESeq2")
  if (!methods::is(pilot_dds, "DESeqDataSet"))
    stop("pilot_dds must be a DESeqDataSet object")

  disp <- DESeq2::dispersions(pilot_dds)
  if (is.null(disp) || all(is.na(disp))) {
    message("Estimating dispersions on supplied pilot_dds...")
    pilot_dds <- DESeq2::estimateSizeFactors(pilot_dds)
    pilot_dds <- DESeq2::estimateDispersions(pilot_dds)
    disp <- DESeq2::dispersions(pilot_dds)
  }

  med <- median(disp, na.rm = TRUE)
  list(
    median_dispersion = med,
    mean_dispersion   = mean(disp, na.rm = TRUE),
    iqr_dispersion    = stats::quantile(disp, probs = c(0.25, 0.75), na.rm = TRUE),
    n_genes           = length(disp),
    fitted            = disp,
    interpretation    = .disp_msg(med)
  )
}

#' Compute power from a pilot DESeqDataSet
#'
#' Converts dispersion to a depth-aware CV: CV ≈ sqrt(disp + 1/depth).
#'
#' @param pilot_dds DESeqDataSet
#' @param n_per_group Replicates per group for the proposed main study
#' @param fold_change Minimum fold-change to detect
#' @param alpha Significance threshold (default 0.05)
#' @param depth Sequencing depth (M reads) — defaults to pilot mean library size
calc_power_from_pilot <- function(pilot_dds, n_per_group, fold_change,
                                  alpha = 0.05, depth = NULL) {
  .need("DESeq2"); .need("RNASeqPower")

  disp_info <- extract_dispersion_deseq2(pilot_dds)
  med_disp <- disp_info$median_dispersion

  if (is.null(depth)) {
    lib_sizes <- colSums(DESeq2::counts(pilot_dds))
    depth <- mean(lib_sizes) / 1e6
  }

  cv <- sqrt(med_disp + 1 / depth)
  power <- RNASeqPower::rnapower(depth = depth, n = n_per_group,
                                 cv = cv, effect = fold_change, alpha = alpha)

  cat(sprintf("✓ Pilot-based power completed: power=%.3f (CV=%.3f, depth=%.1fM)\n",
              power, cv, depth))

  list(
    power = power,
    parameters = list(n_per_group = n_per_group, fold_change = fold_change,
                      alpha = alpha, depth = depth, estimated_cv = cv,
                      median_dispersion = med_disp),
    pilot_info = list(n_samples = ncol(pilot_dds), n_genes = nrow(pilot_dds)),
    interpretation = .power_msg(power, cv)
  )
}

#' Solve for required n from a pilot DESeqDataSet
#'
#' @inheritParams calc_power_from_pilot
#' @param target_power Target power (default 0.8)
calc_samplesize_from_pilot <- function(pilot_dds, fold_change,
                                       target_power = 0.8, alpha = 0.05,
                                       depth = NULL) {
  .need("DESeq2"); .need("RNASeqPower")

  disp_info <- extract_dispersion_deseq2(pilot_dds)
  med_disp <- disp_info$median_dispersion

  if (is.null(depth)) {
    lib_sizes <- colSums(DESeq2::counts(pilot_dds))
    depth <- mean(lib_sizes) / 1e6
  }
  cv <- sqrt(med_disp + 1 / depth)

  n_raw <- RNASeqPower::rnapower(depth = depth, cv = cv, effect = fold_change,
                                 alpha = alpha, power = target_power)
  n <- ceiling(n_raw)

  cat(sprintf("✓ Pilot-based sample-size completed: required n=%d per group (CV=%.3f)\n",
              n, cv))

  list(
    required_n_per_group = n,
    total_samples = 2L * n,
    parameters = list(fold_change = fold_change, target_power = target_power,
                      alpha = alpha, depth = depth, estimated_cv = cv,
                      median_dispersion = med_disp),
    pilot_info = list(n_samples = ncol(pilot_dds), n_genes = nrow(pilot_dds)),
    recommendation = .from_pilot_msg(n, cv, fold_change)
  )
}

#' Thin wrapper noting that PROPER simulation needs explicit configuration
#'
#' Returns guidance rather than running a partial PROPER pipeline.
calc_power_proper <- function(pilot_dds, n_per_group_range = 3:8,
                              fold_change_range = c(1.5, 2, 3), nsims = 10) {
  if (!requireNamespace("PROPER", quietly = TRUE)) {
    message("PROPER not installed; install via BiocManager::install('PROPER').")
  }
  message("PROPER simulation requires more configuration than fits a wrapper.")
  message("See the PROPER vignette for the proper estParam + simRun + comparePower workflow.")
  list(
    note = "PROPER simulation deferred to user — see vignette.",
    pilot_info = list(n_samples = ncol(pilot_dds),
                      n_genes   = nrow(pilot_dds))
  )
}

.disp_msg <- function(d) {
  if (d < 0.1) "Low dispersion: cell lines or controlled conditions"
  else if (d < 0.3) "Moderate dispersion: inbred model organisms"
  else if (d < 0.5) "High dispersion: outbred animals or human samples"
  else "Very high dispersion: check for batch effects or outliers"
}

.power_msg <- function(p, cv) {
  base <- sprintf("Power=%.3f based on pilot (CV=%.3f). ", p, cv)
  if (p >= 0.80) paste0(base, "Adequate.")
  else paste0(base, "Low — pilot variability is higher than expected.")
}

.from_pilot_msg <- function(n, cv, fc) {
  base <- sprintf("Required n=%d (CV=%.3f, FC=%.2f). ", n, cv, fc)
  if (cv > 0.4) base <- paste0(base, "Note: high pilot CV; consider variance reduction. ")
  if (n > 15)   base <- paste0(base, "Large n required — consider larger effects or paired design.")
  base
}
