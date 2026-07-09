# power_rnaseq.R
# Bulk RNA-seq power and sample size calculations using RNASeqPower.
# Reference: Hart SN et al. (2013), J Comput Biol 20:970-978.
#
# Inputs:
#   - depth (M reads), n_per_group, cv, fold_change, alpha
# Outputs:
#   - power estimate or required sample size, with interpretation.

.has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

.require_rnaseqpower <- function() {
  if (!.has("RNASeqPower")) {
    stop("Install with: BiocManager::install('RNASeqPower')", call. = FALSE)
  }
}

#' Calculate power for a bulk RNA-seq design
#'
#' @param depth Sequencing depth in millions of reads per sample
#' @param n_per_group Replicates per condition (must be >= 2)
#' @param cv Coefficient of variation (biological variability)
#' @param fold_change Minimum fold-change to detect (> 1)
#' @param alpha Significance threshold (default 0.05)
#' @return list with power, parameters, interpretation
calc_power_rnaseq <- function(depth, n_per_group, cv, fold_change, alpha = 0.05) {
  .require_rnaseqpower()
  if (depth <= 0 || n_per_group < 2 || cv <= 0 || fold_change <= 1)
    stop("depth, n_per_group, cv must be positive; n >= 2; fold_change > 1")
  if (alpha <= 0 || alpha >= 1) stop("alpha must lie in (0, 1)")

  power <- RNASeqPower::rnapower(depth = depth, n = n_per_group,
                                 cv = cv, effect = fold_change, alpha = alpha)
  out <- list(
    power = power,
    parameters = list(depth = depth, n_per_group = n_per_group, cv = cv,
                      fold_change = fold_change, alpha = alpha),
    interpretation = .interpret_power(power)
  )
  cat(sprintf("✓ Power analysis completed: %d rows / power=%.3f for n=%d, FC=%.1f, CV=%.2f\n",
              1, power, n_per_group, fold_change, cv))
  out
}

#' Solve for the required sample size
#'
#' @param depth Sequencing depth (M reads)
#' @param cv Coefficient of variation
#' @param fold_change Minimum detectable fold-change
#' @param target_power Target statistical power (default 0.8)
#' @param alpha Significance threshold (default 0.05)
#' @return list with required_n_per_group, total_samples, recommendation
calc_sample_size_rnaseq <- function(depth, cv, fold_change,
                                    target_power = 0.8, alpha = 0.05) {
  .require_rnaseqpower()
  if (depth <= 0 || cv <= 0 || fold_change <= 1)
    stop("depth, cv positive; fold_change > 1")
  if (target_power <= 0 || target_power >= 1) stop("power must lie in (0, 1)")

  n_raw <- RNASeqPower::rnapower(depth = depth, cv = cv, effect = fold_change,
                                 alpha = alpha, power = target_power)
  n <- ceiling(n_raw)

  out <- list(
    required_n_per_group = n,
    total_samples = 2L * n,
    parameters = list(depth = depth, cv = cv, fold_change = fold_change,
                      target_power = target_power, alpha = alpha),
    recommendation = .recommend_n(n, cv, fold_change)
  )
  cat(sprintf("✓ Sample-size calculation completed: required n=%d per group (%d total) for power=%.2f\n",
              n, 2L * n, target_power))
  out
}

#' Estimate CV from a raw pilot count matrix
#'
#' Filters low-count genes (mean < 10) then takes the median per-gene CV.
#'
#' @param pilot_counts Matrix or data frame of raw counts (genes x samples)
#' @return list with estimated_cv, cv_range (25-75%), n_genes_used
estimate_cv_from_pilot <- function(pilot_counts) {
  if (!is.matrix(pilot_counts) && !is.data.frame(pilot_counts))
    stop("pilot_counts must be a matrix or data frame")
  m <- as.matrix(pilot_counts)
  m <- m[rowMeans(m) >= 10, , drop = FALSE]
  if (nrow(m) == 0) stop("All genes below mean count threshold of 10")

  gene_means <- rowMeans(m)
  gene_sds   <- apply(m, 1, sd)
  cvs <- gene_sds / gene_means
  cvs <- cvs[is.finite(cvs)]

  list(
    estimated_cv  = median(cvs, na.rm = TRUE),
    cv_range      = stats::quantile(cvs, probs = c(0.25, 0.75), na.rm = TRUE),
    n_genes_used  = nrow(m),
    recommendation = .interpret_cv(median(cvs, na.rm = TRUE))
  )
}

#' Compute power across a range of sample sizes (for plotting)
#'
#' @param depth Sequencing depth (M reads)
#' @param n_range Integer vector of sample sizes
#' @param cv Coefficient of variation
#' @param fold_change Fold-change threshold
#' @param alpha Significance threshold
#' @return data frame with n_per_group, power
calc_power_range <- function(depth, n_range, cv, fold_change, alpha = 0.05) {
  .require_rnaseqpower()
  powers <- vapply(n_range, function(n) {
    RNASeqPower::rnapower(depth = depth, n = n, cv = cv,
                          effect = fold_change, alpha = alpha)
  }, numeric(1))
  data.frame(n_per_group = n_range, power = powers,
             depth = depth, cv = cv, fold_change = fold_change)
}

.interpret_power <- function(p) {
  if (p >= 0.90) "Excellent: high power"
  else if (p >= 0.80) "Good: adequate power for most applications"
  else if (p >= 0.70) "Moderate: consider increasing n"
  else "Low: increase n or relax effect size"
}

.recommend_n <- function(n, cv, fc) {
  msg <- sprintf("Recommended: n=%d per group (%d total) for CV=%.2f, FC=%.2f\n",
                 n, 2L * n, cv, fc)
  if (n < 3) msg <- paste0(msg, "Warning: n<3 is too low for valid inference\n")
  if (n > 20) msg <- paste0(msg, "Note: large n required — verify CV assumption\n")
  msg
}

.interpret_cv <- function(cv) {
  if (cv < 0.2) "Low CV: cell lines or controlled conditions"
  else if (cv < 0.3) "Moderate CV: inbred animals, primary cells"
  else if (cv < 0.5) "High CV: outbred animals, human cohorts"
  else "Very high CV: consider reducing variability or paired design"
}
