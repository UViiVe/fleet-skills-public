# lasso_workflow.R
# Nested-CV LASSO / elastic net panel selection with stability analysis.
# Reference: Ali et al. Nature Medicine 2025
# (github.com/NeuroGenomicsAndInformatics/NatMed_2025_GNPC)
#
# Exposes:
#   run_lasso_panel(X, y, ...)
#   predict_biomarker_panel(model_result, new_X)

suppressPackageStartupMessages({
    library(glmnet)
    library(pROC)
})


#' Run repeated nested-CV LASSO with stability selection.
#'
#' Per Ali et al. (2025): n_repeats class-balanced 70/30 splits, each with
#' an inner cv.glmnet 10-fold loop for lambda; the panel is the set of
#' features selected in >= stability_threshold of the outer folds, refit
#' on the full discovery cohort.
#'
#' @param X Samples x features numeric matrix (already standardized).
#' @param y Binary outcome vector (0/1) aligned to rows of X.
#' @param n_repeats Number of outer-fold repeats (default 50).
#' @param train_fraction Fraction held in the training partition (default 0.7).
#' @param alpha Elastic-net mixing (1 = LASSO, 0.5 = elastic net, 0 = ridge).
#' @param n_inner_folds Folds for inner cv.glmnet (default 10).
#' @param stability_threshold Minimum selection frequency for a stable feature.
#' @param top_n_features Optional cap on selected features per iteration.
#' @param seed Random seed for reproducibility.
#' @return list with final_model, stable_features, feature_importance,
#'         fold_aucs, mean_auc, auc_ci, cv_predictions, parameters.
run_lasso_panel <- function(X, y,
                            n_repeats = 50,
                            train_fraction = 0.7,
                            alpha = 1.0,
                            n_inner_folds = 10,
                            stability_threshold = 0.8,
                            top_n_features = NULL,
                            seed = 42) {

    cat("\n=== LASSO biomarker panel selection ===\n\n")
    cat(sprintf("  samples:     %d\n", nrow(X)))
    cat(sprintf("  features:    %d\n", ncol(X)))
    cat(sprintf("  outcome:     %d positives, %d negatives\n",
                sum(y == 1L), sum(y == 0L)))
    cat(sprintf("  n_repeats:   %d\n", n_repeats))
    cat(sprintf("  alpha:       %.2f (%s)\n", alpha,
                if (alpha == 1) "pure LASSO" else if (alpha == 0) "ridge"
                else "elastic net"))
    cat(sprintf("  stability:   %.2f\n\n", stability_threshold))

    set.seed(seed)
    feature_names <- colnames(X)
    selection_count <- setNames(integer(ncol(X)), feature_names)
    coef_matrix <- matrix(0, nrow = n_repeats, ncol = ncol(X),
                          dimnames = list(NULL, feature_names))
    fold_aucs <- numeric(n_repeats)
    fold_sensitivities <- numeric(n_repeats)
    fold_specificities <- numeric(n_repeats)
    cv_predictions <- vector("list", n_repeats)

    cat("Running outer folds...\n")
    pb_step <- max(1L, n_repeats %/% 10L)

    idx_pos <- which(y == 1L)
    idx_neg <- which(y == 0L)

    for (i in seq_len(n_repeats)) {
        n_train_pos <- round(length(idx_pos) * train_fraction)
        n_train_neg <- round(length(idx_neg) * train_fraction)
        train_idx <- c(sample(idx_pos, n_train_pos),
                        sample(idx_neg, n_train_neg))
        test_idx <- setdiff(seq_len(nrow(X)), train_idx)

        X_train <- X[train_idx, , drop = FALSE]
        y_train <- y[train_idx]
        X_test <- X[test_idx, , drop = FALSE]
        y_test <- y[test_idx]

        cv_fit <- tryCatch(
            cv.glmnet(x = X_train, y = y_train, family = "binomial",
                       alpha = alpha,
                       nfolds = min(n_inner_folds, length(y_train)),
                       type.measure = "auc",
                       standardize = FALSE),
            error = function(e) NULL
        )
        if (is.null(cv_fit)) {
            fold_aucs[i] <- NA_real_
            next
        }

        coefs <- as.vector(coef(cv_fit, s = "lambda.min"))[-1]
        names(coefs) <- feature_names
        selected <- which(coefs != 0)
        if (length(selected) > 0L) {
            selection_count[selected] <- selection_count[selected] + 1L
            coef_matrix[i, selected] <- coefs[selected]
            if (!is.null(top_n_features) &&
                length(selected) > top_n_features) {
                top <- order(abs(coefs[selected]), decreasing = TRUE)[
                    seq_len(top_n_features)]
                selected <- selected[top]
            }
        }

        pred_probs <- as.vector(predict(cv_fit, newx = X_test,
                                          s = "lambda.min", type = "response"))
        if (length(unique(y_test)) == 2L) {
            roc_obj <- tryCatch(pROC::roc(y_test, pred_probs, quiet = TRUE),
                                error = function(e) NULL)
            if (!is.null(roc_obj)) {
                fold_aucs[i] <- pROC::auc(roc_obj)
                co <- pROC::coords(roc_obj, "best",
                                   ret = c("sensitivity", "specificity"))
                fold_sensitivities[i] <- co$sensitivity[1]
                fold_specificities[i] <- co$specificity[1]
            } else {
                fold_aucs[i] <- NA_real_
            }
        }

        cv_predictions[[i]] <- data.frame(
            sample_id      = rownames(X_test),
            true_label     = y_test,
            predicted_prob = pred_probs,
            fold           = i,
            stringsAsFactors = FALSE
        )

        if (i %% pb_step == 0L || i == 1L) {
            cat(sprintf("  fold %d/%d  AUC=%.3f  selected=%d\n",
                        i, n_repeats, fold_aucs[i],
                        length(selected)))
        }
    }

    valid <- sum(!is.na(fold_aucs))
    selection_freq <- selection_count / valid
    stable_features <- names(selection_freq[selection_freq >= stability_threshold])

    cat(sprintf("\nValid outer folds: %d / %d\n", valid, n_repeats))
    cat(sprintf("Features ever selected: %d\n", sum(selection_count > 0L)))
    cat(sprintf("Features at >= %.0f%% stability: %d\n",
                100 * stability_threshold, length(stable_features)))

    if (length(stable_features) < 2L) {
        cat("WARNING: <2 stable features — relaxing to top 10 by frequency.\n")
        top10 <- names(sort(selection_freq, decreasing = TRUE))[
            seq_len(min(10, length(selection_freq)))]
        stable_features <- top10[selection_freq[top10] > 0]
        if (length(stable_features) < 2L) {
            stable_features <- names(sort(selection_freq, decreasing = TRUE))[1:2]
        }
    }

    cat("Stable feature panel:\n")
    for (feat in stable_features) {
        cat(sprintf("  - %s  (selected in %.1f%% of folds)\n",
                    feat, 100 * selection_freq[feat]))
    }

    feature_importance <- data.frame(
        feature              = feature_names,
        selection_frequency  = as.numeric(selection_freq),
        mean_coefficient     = colMeans(coef_matrix, na.rm = TRUE),
        sd_coefficient       = apply(coef_matrix, 2, sd, na.rm = TRUE),
        is_stable            = feature_names %in% stable_features,
        stringsAsFactors     = FALSE
    )
    feature_importance <- feature_importance[
        order(-feature_importance$selection_frequency), ]

    final_model <- NULL
    if (length(stable_features) >= 2L) {
        X_stable <- X[, stable_features, drop = FALSE]
        final_model <- cv.glmnet(x = X_stable, y = y,
                                  family = "binomial", alpha = alpha,
                                  nfolds = min(n_inner_folds, length(y)),
                                  type.measure = "auc",
                                  standardize = FALSE)
        final_preds <- as.vector(predict(final_model, newx = X_stable,
                                          s = "lambda.min", type = "response"))
        final_auc <- pROC::auc(pROC::roc(y, final_preds, quiet = TRUE))
        cat(sprintf("\nFinal-model in-sample AUC (reference only): %.3f\n",
                    final_auc))
        cat("Cite the CV AUC below for the unbiased performance estimate.\n")
    } else {
        cat("\nWARNING: no stable features — final model not fit.\n")
    }

    valid_aucs <- fold_aucs[!is.na(fold_aucs)]
    mean_auc <- mean(valid_aucs)
    auc_ci <- quantile(valid_aucs, probs = c(0.025, 0.975))

    cat("\nCross-validation performance:\n")
    cat(sprintf("  mean AUC:       %.3f  (95%% CI %.3f - %.3f)\n",
                mean_auc, auc_ci[1], auc_ci[2]))
    cat(sprintf("  mean sens:      %.3f\n",
                mean(fold_sensitivities, na.rm = TRUE)))
    cat(sprintf("  mean spec:      %.3f\n",
                mean(fold_specificities, na.rm = TRUE)))

    result <- list(
        final_model            = final_model,
        stable_features        = stable_features,
        feature_importance     = feature_importance,
        fold_aucs              = fold_aucs,
        mean_auc               = mean_auc,
        auc_ci                 = auc_ci,
        fold_sensitivities     = fold_sensitivities,
        fold_specificities     = fold_specificities,
        cv_predictions         = do.call(rbind, cv_predictions),
        selection_frequencies  = selection_freq,
        coef_matrix            = coef_matrix,
        parameters             = list(
            n_repeats           = n_repeats,
            train_fraction      = train_fraction,
            alpha               = alpha,
            n_inner_folds       = n_inner_folds,
            stability_threshold = stability_threshold,
            seed                = seed,
            n_samples           = nrow(X),
            n_features          = ncol(X),
            n_stable            = length(stable_features)
        )
    )

    cat("\nokay LASSO panel selection completed successfully!\n")
    cat(sprintf("  panel size: %d features\n", length(stable_features)))
    cat(sprintf("  CV AUC:     %.3f\n", mean_auc))
    invisible(result)
}


#' Apply the locked panel to a new feature matrix.
#'
#' @param model_result Output of run_lasso_panel.
#' @param new_X Samples x features matrix using the discovery feature names.
#' @return Named numeric vector of predicted probabilities.
predict_biomarker_panel <- function(model_result, new_X) {
    if (is.null(model_result$final_model)) {
        stop("Final model is missing — panel selection failed.")
    }
    stable <- model_result$stable_features
    available <- intersect(stable, colnames(new_X))
    if (length(available) == 0L) {
        stop("None of the panel features are present in new_X.")
    }
    if (length(available) < length(stable)) {
        cat(sprintf("WARNING: %d / %d panel features missing from new_X\n",
                    length(stable) - length(available), length(stable)))
    }
    preds <- as.vector(predict(model_result$final_model,
                                newx = new_X[, available, drop = FALSE],
                                s = "lambda.min", type = "response"))
    names(preds) <- rownames(new_X)
    preds
}
