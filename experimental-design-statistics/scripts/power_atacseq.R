# power_atacseq.R
# ATAC-seq (and ChIP-seq / peak-based assays) power and sample-size via ssizeRNA.
# Reference: Liu et al. (2014) Bioinformatics 30:301-304; ssizeRNA package vignette.
#
# Inputs:
#   - m (total peaks), m1 (expected differential peaks),
#     fold_change, FDR, mean count per peak (mu), dispersion (disp),
#     either n_per_group (-> power) or power (-> required n).
# Outputs:
#   - power or required_n_per_group with interpretation.

.req <- function(pkg, bioc = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Install with: %s::install('%s')",
                 if (bioc) "BiocManager" else "install.packages", pkg), call. = FALSE)
}

#' Power or required n for ATAC-seq (peak-based) differential analysis
#'
#' @param m Total peaks tested
#' @param m1 Expected number of differentially accessible peaks
#' @param fold_change Minimum fold-change to detect (> 1)
#' @param fdr FDR threshold (default 0.05)
#' @param n_per_group If provided, computes power; if NULL, requires `power`
#' @param mu Mean read count per peak per sample (default 10)
#' @param disp Dispersion parameter (default 0.1; typical ATAC range 0.1-0.5)
#' @param power Target power if computing required n
#' @return list with either `power` or `required_n_per_group`, plus parameters
calc_power_atac <- function(m, m1, fold_change,
                            fdr = 0.05,
                            n_per_group = NULL,
                            mu = 10, disp = 0.1,
                            power = NULL) {
  .req("ssizeRNA")

  if (m <= 0 || m1 <= 0 || m1 > m)
    stop("m, m1 positive; m1 <= m")
  if (fold_change <= 1 || mu <= 0 || disp < 0)
    stop("fold_change > 1; mu positive; disp non-negative")

  pi0 <- 1 - (m1 / m)

  if (!is.null(n_per_group)) {
    res <- ssizeRNA::ssizeRNA_vary(
      nGenes = m, pi0 = pi0, m = n_per_group,
      mu = mu, disp = disp, fc = fold_change, fdr = fdr, power = NULL
    )
    out <- list(
      power = res$power,
      n_per_group = n_per_group,
      parameters = list(total_peaks = m, differential_peaks = m1,
                        fold_change = fold_change, fdr = fdr,
                        mean_count = mu, dispersion = disp),
      interpretation = .atac_power_msg(res$power, n_per_group)
    )
    cat(sprintf("✓ ATAC power completed: power=%.3f at n=%d\n",
                res$power, n_per_group))
    return(out)
  }

  if (!is.null(power)) {
    res <- ssizeRNA::ssizeRNA_vary(
      nGenes = m, pi0 = pi0,
      mu = mu, disp = disp, fc = fold_change, fdr = fdr,
      power = power, m = NULL
    )
    n <- ceiling(res$m)
    out <- list(
      required_n_per_group = n,
      total_samples = 2L * n,
      target_power = power,
      parameters = list(total_peaks = m, differential_peaks = m1,
                        fold_change = fold_change, fdr = fdr,
                        mean_count = mu, dispersion = disp),
      recommendation = .atac_n_msg(n, fold_change)
    )
    cat(sprintf("✓ ATAC sample-size completed: required n=%d per group for power=%.2f\n",
                n, power))
    return(out)
  }

  stop("Provide either n_per_group (-> power) or power (-> required n)")
}

#' Estimate ATAC mean and dispersion from pilot peak-count matrix
#'
#' @param pilot_counts Matrix or data frame of raw counts (peaks x samples)
#' @return list with estimated_mu, estimated_dispersion, n_peaks_used
estimate_atac_params <- function(pilot_counts) {
  .req("DESeq2")

  m <- as.matrix(pilot_counts)
  m <- m[rowMeans(m) >= 5, , drop = FALSE]
  if (nrow(m) == 0) stop("All peaks below mean count threshold of 5")

  est_mu <- mean(rowMeans(m))

  condition <- factor(rep(c("A", "B"), length.out = ncol(m)))
  coldata <- data.frame(condition = condition, row.names = colnames(m))
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = m, colData = coldata,
                                        design = ~ condition)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::estimateDispersions(dds)
  est_disp <- median(DESeq2::dispersions(dds), na.rm = TRUE)

  list(
    estimated_mu          = est_mu,
    estimated_dispersion  = est_disp,
    n_peaks_used          = nrow(m),
    recommendation        = .atac_params_msg(est_mu, est_disp)
  )
}

.atac_power_msg <- function(p, n) {
  base <- sprintf("Power = %.3f at n=%d per group. ", p, n)
  if (p >= 0.80) paste0(base, "Adequate for detecting differential peaks.")
  else paste0(base, "Low — consider n=", n + 2, " or larger effects.")
}

.atac_n_msg <- function(n, fc) {
  base <- sprintf("Recommended: n=%d per group for fold-change %.2f. ", n, fc)
  if (n < 3) paste0(base, "Warning: ATAC typically needs >= 3 replicates.")
  else if (n <= 5) paste0(base, "Reasonable ATAC design.")
  else base
}

.atac_params_msg <- function(mu, disp) {
  parts <- sprintf("mu=%.1f, disp=%.3f. ", mu, disp)
  if (mu < 5) parts <- paste0(parts, "Low coverage — consider deeper sequencing. ")
  else if (mu > 50) parts <- paste0(parts, "Good coverage. ")
  if (disp < 0.2) parts <- paste0(parts, "Low variability (favorable).")
  else if (disp > 0.5) parts <- paste0(parts, "High variability — more replicates needed.")
  parts
}
