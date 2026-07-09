# sample_size_de.R
# Sample size estimation for differential-expression studies via ssizeRNA,
# accounting for the multiple-testing burden.
#
# Inputs:
#   - pi0 (proportion non-DE), m1 (expected DE genes), fold_change, fdr, power,
#     mean_count, dispersion.
# Outputs:
#   - Required sample size per group with rationale.

.req <- function(pkg, bioc = TRUE) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop(sprintf("Install %s with %s.",
                 pkg, if (bioc) "BiocManager::install" else "install.packages"),
         call. = FALSE)
}

#' Required sample size for DE with FDR control
#'
#' @param pi0 Proportion of non-DE genes (default 0.9)
#' @param m1 Expected number of truly DE genes (default 2000)
#' @param fold_change Minimum fold-change to detect (> 1)
#' @param fdr Target FDR (default 0.05)
#' @param power Target power (default 0.8)
#' @param mean_count Mean read count per gene (default 20)
#' @param dispersion Dispersion parameter (default 0.1)
calc_samplesize_de <- function(pi0 = 0.9, m1 = 2000,
                               fold_change = 1.5, fdr = 0.05, power = 0.8,
                               mean_count = 20, dispersion = 0.1) {
  .req("ssizeRNA")
  if (pi0 <= 0 || pi0 >= 1) stop("pi0 must lie in (0, 1)")
  if (m1 <= 0 || fold_change <= 1) stop("m1 > 0; fold_change > 1")
  if (power <= 0 || power >= 1) stop("power must lie in (0, 1)")

  m <- ceiling(m1 / (1 - pi0))

  res <- ssizeRNA::ssizeRNA_vary(
    nGenes = m, pi0 = pi0, m = NULL,
    mu = mean_count, disp = dispersion,
    fc = fold_change, fdr = fdr, power = power
  )
  n <- ceiling(res$m)

  cat(sprintf("✓ Sample size estimated: %d per group (%d total) for %d DE genes at FC=%.2f\n",
              n, 2L * n, m1, fold_change))

  list(
    required_n_per_group = n,
    total_samples = 2L * n,
    parameters = list(total_genes = m, de_genes = m1, pi0 = pi0,
                      fold_change = fold_change, fdr = fdr,
                      target_power = power, mean_count = mean_count,
                      dispersion = dispersion),
    recommendation = .rec_de(n, fold_change)
  )
}

#' Sample size driven by pilot DESeq2 data
#'
#' @param pilot_dds DESeqDataSet from pilot experiment (DESeq run)
#' @param fold_change Minimum fold-change to detect
#' @param power Target power (default 0.8)
#' @param fdr Target FDR (default 0.05)
#' @param expected_de_prop Expected proportion of DE genes (default 0.10)
samplesize_from_pilot <- function(pilot_dds, fold_change,
                                  power = 0.8, fdr = 0.05,
                                  expected_de_prop = 0.10) {
  .req("DESeq2"); .req("ssizeRNA")

  if (is.null(DESeq2::dispersions(pilot_dds))) {
    pilot_dds <- DESeq2::estimateSizeFactors(pilot_dds)
    pilot_dds <- DESeq2::estimateDispersions(pilot_dds)
  }
  med_disp <- median(DESeq2::dispersions(pilot_dds), na.rm = TRUE)

  norm_counts <- DESeq2::counts(pilot_dds, normalized = TRUE)
  mean_count <- mean(rowMeans(norm_counts))

  n_genes <- nrow(pilot_dds)
  pi0 <- 1 - expected_de_prop
  m1  <- ceiling(n_genes * expected_de_prop)

  res <- ssizeRNA::ssizeRNA_vary(
    nGenes = n_genes, pi0 = pi0, m = NULL,
    mu = mean_count, disp = med_disp,
    fc = fold_change, fdr = fdr, power = power
  )
  n <- ceiling(res$m)

  cat(sprintf("✓ Pilot-based sample size: %d per group (median disp=%.3f)\n",
              n, med_disp))

  list(
    required_n_per_group = n,
    total_samples = 2L * n,
    pilot_parameters = list(median_dispersion = med_disp,
                            mean_count = mean_count,
                            n_genes = n_genes,
                            expected_de_genes = m1),
    parameters = list(fold_change = fold_change, power = power, fdr = fdr),
    recommendation = .rec_from_pilot_de(n, med_disp, fold_change)
  )
}

#' Sample size for paired designs (within-subject correlation reduces variance)
#'
#' @param pilot_dds DESeqDataSet from pilot paired data
#' @param fold_change Minimum fold-change to detect
#' @param power Target power (default 0.8)
#' @param fdr Target FDR (default 0.05)
#' @param disp_reduction Variance reduction factor for pairing (default 0.7)
paired_samplesize <- function(pilot_dds, fold_change,
                              power = 0.8, fdr = 0.05,
                              disp_reduction = 0.7) {
  .req("DESeq2"); .req("ssizeRNA")

  if (is.null(DESeq2::dispersions(pilot_dds))) {
    pilot_dds <- DESeq2::estimateSizeFactors(pilot_dds)
    pilot_dds <- DESeq2::estimateDispersions(pilot_dds)
  }
  med_disp <- median(DESeq2::dispersions(pilot_dds), na.rm = TRUE)
  adj_disp <- med_disp * disp_reduction

  norm_counts <- DESeq2::counts(pilot_dds, normalized = TRUE)
  mean_count <- mean(rowMeans(norm_counts))
  n_genes <- nrow(pilot_dds)

  res <- ssizeRNA::ssizeRNA_vary(
    nGenes = n_genes, pi0 = 0.9, m = NULL,
    mu = mean_count, disp = adj_disp,
    fc = fold_change, fdr = fdr, power = power
  )
  pairs <- ceiling(res$m)

  cat(sprintf("✓ Paired sample size: %d pairs (%d total samples)\n", pairs, 2L * pairs))

  list(
    required_pairs = pairs,
    total_samples = 2L * pairs,
    pilot_parameters = list(median_dispersion = med_disp,
                            adjusted_dispersion = adj_disp,
                            mean_count = mean_count, n_genes = n_genes),
    parameters = list(fold_change = fold_change, power = power, fdr = fdr,
                      disp_reduction = disp_reduction),
    note = "Pairs are the unit of replication, not independent samples.",
    recommendation = .rec_paired(pairs, fold_change)
  )
}

#' Quick reference table of sample sizes at various effect sizes
#'
#' @param cv Coefficient of variation (default 0.4)
#' @param depth Sequencing depth in M reads (default 20)
quick_reference_samplesize <- function(cv = 0.4, depth = 20) {
  if (!requireNamespace("RNASeqPower", quietly = TRUE))
    stop("Install RNASeqPower with BiocManager::install('RNASeqPower')")

  fcs <- c(1.25, 1.5, 2, 3, 4)
  ns <- vapply(fcs, function(fc) {
    ceiling(RNASeqPower::rnapower(depth = depth, cv = cv,
                                  effect = fc, alpha = 0.05, power = 0.8))
  }, numeric(1))

  difficulty <- c("Very challenging", "Challenging", "Moderate", "Easy", "Very easy")
  data.frame(fold_change = fcs, required_n_per_group = ns,
             total_samples = 2L * ns, difficulty = difficulty)
}

.rec_de <- function(n, fc) {
  base <- sprintf("Recommended: %d per group for FC=%.2f. ", n, fc)
  if (n < 3) paste0(base, "Warning: < 3 replicates is not recommended for RNA-seq.")
  else if (n <= 6) paste0(base, "Standard RNA-seq design.")
  else if (n > 12) paste0(base, "Large n required — consider larger effects or paired design.")
  else base
}

.rec_from_pilot_de <- function(n, disp, fc) {
  base <- sprintf("Pilot-derived: n=%d per group; observed dispersion %.3f, target FC %.2f. ",
                  n, disp, fc)
  if (disp > 0.5) base <- paste0(base, "High dispersion — investigate batch effects.")
  base
}

.rec_paired <- function(n_pairs, fc) {
  base <- sprintf("Paired design: %d pairs (%d total) for FC=%.2f. ",
                  n_pairs, 2L * n_pairs, fc)
  if (n_pairs < 5) base <- paste0(base, "Warning: <5 pairs may be underpowered even with pairing.")
  base
}
