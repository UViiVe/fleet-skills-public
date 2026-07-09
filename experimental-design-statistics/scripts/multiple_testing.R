# multiple_testing.R
# Multiple-testing correction methods and a method recommender.
#
# Inputs:
#   - p-values (vector); optional covariates (for IHW).
# Outputs:
#   - adjusted p / q values, significance counts, interpretation strings.

.has <- function(pkg) requireNamespace(pkg, quietly = TRUE)

#' Benjamini-Hochberg FDR
#'
#' @param pvalues Numeric vector of raw p-values (NA allowed)
#' @param q_threshold FDR threshold (default 0.05)
#' @return list with adjusted p-values, significance flags, counts
apply_bh_fdr <- function(pvalues, q_threshold = 0.05) {
  ok <- !is.na(pvalues)
  adj <- p.adjust(pvalues[ok], method = "BH")
  sig <- adj < q_threshold

  full_adj <- rep(NA_real_, length(pvalues)); full_adj[ok] <- adj
  full_sig <- rep(NA, length(pvalues));      full_sig[ok] <- sig

  out <- list(
    method = "Benjamini-Hochberg FDR",
    adjusted_pvalues = full_adj,
    is_significant = full_sig,
    n_significant = sum(sig, na.rm = TRUE),
    n_tested = sum(ok),
    threshold = q_threshold,
    interpretation = sprintf(
      "BH-FDR: %d of %d significant (%.2f%%); expected FP ≈ %d",
      sum(sig, na.rm = TRUE), sum(ok),
      100 * sum(sig, na.rm = TRUE) / sum(ok),
      round(sum(sig, na.rm = TRUE) * q_threshold)
    )
  )
  cat(sprintf("✓ BH-FDR completed: %d significant of %d\n",
              out$n_significant, out$n_tested))
  out
}

#' Bonferroni FWER correction
apply_bonferroni <- function(pvalues, alpha = 0.05) {
  ok <- !is.na(pvalues)
  adj <- p.adjust(pvalues[ok], method = "bonferroni")
  sig <- adj < alpha

  full_adj <- rep(NA_real_, length(pvalues)); full_adj[ok] <- adj
  full_sig <- rep(NA, length(pvalues));      full_sig[ok] <- sig

  out <- list(
    method = "Bonferroni",
    adjusted_pvalues = full_adj,
    is_significant = full_sig,
    n_significant = sum(sig, na.rm = TRUE),
    n_tested = sum(ok),
    threshold = alpha,
    bonferroni_threshold = alpha / sum(ok),
    interpretation = sprintf(
      "Bonferroni: %d significant at p < %.2e (FWER controlled at %.2f)",
      sum(sig, na.rm = TRUE), alpha / sum(ok), alpha
    )
  )
  cat(sprintf("✓ Bonferroni completed: %d significant\n", out$n_significant))
  out
}

#' Storey q-value
apply_qvalue <- function(pvalues, fdr_threshold = 0.05) {
  if (!.has("qvalue"))
    stop("Install qvalue with BiocManager::install('qvalue')")

  ok <- !is.na(pvalues)
  qobj <- qvalue::qvalue(pvalues[ok])
  qv <- qobj$qvalues
  sig <- qv < fdr_threshold

  full_qv <- rep(NA_real_, length(pvalues)); full_qv[ok] <- qv
  full_sig <- rep(NA, length(pvalues));      full_sig[ok] <- sig

  out <- list(
    method = "q-value (Storey)",
    qvalues = full_qv,
    is_significant = full_sig,
    n_significant = sum(sig, na.rm = TRUE),
    n_tested = sum(ok),
    threshold = fdr_threshold,
    pi0 = qobj$pi0,
    estimated_true_positives = round((1 - qobj$pi0) * sum(ok)),
    interpretation = sprintf(
      "q-value: %d significant; estimated pi0=%.2f (≈%d true positives among tests)",
      sum(sig, na.rm = TRUE), qobj$pi0,
      round((1 - qobj$pi0) * sum(ok))
    )
  )
  cat(sprintf("✓ q-value completed: %d significant\n", out$n_significant))
  out
}

#' Independent Hypothesis Weighting (IHW)
#'
#' @param pvalues Vector of p-values
#' @param covariates Vector of covariates (must be independent of p-value under null)
#' @param alpha Significance threshold (default 0.05)
apply_ihw <- function(pvalues, covariates, alpha = 0.05) {
  if (!.has("IHW"))
    stop("Install IHW with BiocManager::install('IHW')")
  ok <- !is.na(pvalues) & !is.na(covariates)
  res <- IHW::ihw(pvalues[ok], covariates[ok], alpha = alpha)

  adj <- IHW::adj_pvalues(res)
  rej <- IHW::rejections(res)
  full_adj <- rep(NA_real_, length(pvalues)); full_adj[ok] <- adj
  full_rej <- rep(NA,        length(pvalues)); full_rej[ok] <- rej

  bh_n <- sum(p.adjust(pvalues[ok], method = "BH") < alpha)
  out <- list(
    method = "Independent Hypothesis Weighting (IHW)",
    adjusted_pvalues = full_adj,
    is_significant = full_rej,
    n_significant = sum(rej),
    n_tested = sum(ok),
    threshold = alpha,
    bh_rejections = bh_n,
    additional_discoveries = sum(rej) - bh_n,
    interpretation = if (sum(rej) > bh_n)
      sprintf("IHW: %d significant; BH would find %d (+%.1f%% with informative covariate)",
              sum(rej), bh_n, 100 * (sum(rej) - bh_n) / max(bh_n, 1))
    else
      sprintf("IHW: %d significant; similar to BH (%d) — covariate not informative",
              sum(rej), bh_n)
  )
  cat(sprintf("✓ IHW completed: %d significant (vs BH %d)\n",
              out$n_significant, bh_n))
  out
}

#' Recommend a correction method given test characteristics
#'
#' @param n_tests Number of tests
#' @param expected_true_positives Optional estimate of true positives
#' @param effect_heterogeneity One of "low", "moderate", "high"
#' @param sample_size Sample size per group (optional)
#' @return list with recommended_method and justification
recommend_method <- function(n_tests,
                             expected_true_positives = NULL,
                             effect_heterogeneity = "moderate",
                             sample_size = NULL) {

  method <- NA_character_; justification <- NA_character_

  if (n_tests < 100) {
    method <- "Bonferroni"
    justification <- sprintf(
      "Small number of tests (n=%d) — Bonferroni gives strong FWER control. Use p < %.2e.",
      n_tests, 0.05 / n_tests
    )
  } else if (n_tests > 500000L) {
    method <- "Bonferroni"
    justification <- sprintf(
      "GWAS-scale tests (n=%s) — use the standard 5e-8 threshold.",
      format(n_tests, big.mark = ",")
    )
  } else if (effect_heterogeneity == "high" && n_tests > 5000L) {
    method <- "IHW"
    justification <- "Heterogeneous effects + informative covariate (e.g., mean expression) — IHW typically finds 10-20% more discoveries than BH. Fall back to BH if no covariate is available."
  } else if (!is.null(expected_true_positives) &&
             expected_true_positives > n_tests * 0.05) {
    method <- "q-value"
    justification <- sprintf(
      "Many true positives expected (%d of %d) — q-value is more powerful than BH when pi0 < 0.9.",
      expected_true_positives, n_tests
    )
  } else if (!is.null(sample_size) && sample_size < 6L) {
    method <- "Permutation-based FDR"
    justification <- sprintf(
      "Small sample size (n=%d) — parametric assumptions may fail; consider permutation FDR.",
      sample_size
    )
  } else {
    method <- "Benjamini-Hochberg FDR"
    justification <- "Default for genomics; well understood, broadly accepted. Use FDR < 0.05."
  }

  cat("\n=== Multiple-testing recommendation ===\n")
  cat(sprintf("Method: %s\n", method))
  cat(justification, "\n")

  invisible(list(recommended_method = method,
                 justification = justification,
                 parameters = list(n_tests = n_tests,
                                   expected_true_positives = expected_true_positives,
                                   effect_heterogeneity = effect_heterogeneity,
                                   sample_size = sample_size)))
}

#' Compare multiple correction methods on the same vector of p-values
#'
#' @param pvalues Numeric vector
#' @param covariates Optional covariates for IHW
#' @param threshold Significance threshold (default 0.05)
compare_methods <- function(pvalues, covariates = NULL, threshold = 0.05) {
  bh   <- apply_bh_fdr(pvalues, threshold)
  bonf <- apply_bonferroni(pvalues, threshold)
  qv   <- tryCatch(apply_qvalue(pvalues, threshold), error = function(e) NULL)
  ihw  <- if (!is.null(covariates))
            tryCatch(apply_ihw(pvalues, covariates, threshold), error = function(e) NULL)
          else NULL

  rows <- list(
    c(bh$method, bh$n_significant, bh$n_tested),
    c(bonf$method, bonf$n_significant, bonf$n_tested)
  )
  if (!is.null(qv))  rows <- c(rows, list(c(qv$method, qv$n_significant, qv$n_tested)))
  if (!is.null(ihw)) rows <- c(rows, list(c(ihw$method, ihw$n_significant, ihw$n_tested)))

  out <- do.call(rbind, lapply(rows, function(r) {
    data.frame(Method = r[[1L]],
               Significant = as.integer(r[[2L]]),
               Tested = as.integer(r[[3L]]),
               stringsAsFactors = FALSE)
  }))
  out$Percent <- sprintf("%.2f%%", 100 * out$Significant / out$Tested)
  cat("\n=== Method comparison ===\n"); print(out)
  invisible(out)
}
