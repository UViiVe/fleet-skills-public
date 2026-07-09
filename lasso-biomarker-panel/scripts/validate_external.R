# validate_external.R
# Apply a locked biomarker panel to an external cohort.

validate_panel <- function(model,
                           X,
                           y,
                           cohort_name = "validation") {
    if (!exists("predict_biomarker_panel")) {
        stop("source('scripts/lasso_workflow.R') before validate_external.R")
    }

    y_factor <- factor(y)
    if (nlevels(y_factor) != 2L) {
        stop("External validation outcome must be binary.")
    }
    y_binary <- ifelse(y_factor == levels(y_factor)[2], 1L, 0L)
    pred <- predict_biomarker_panel(model, X)

    roc_obj <- NULL
    auc_value <- NA_real_
    if (requireNamespace("pROC", quietly = TRUE)) {
        roc_obj <- pROC::roc(y_binary, pred, quiet = TRUE)
        auc_value <- as.numeric(pROC::auc(roc_obj))
    }

    result <- list(
        cohort_name = cohort_name,
        n = length(y_binary),
        predictions = data.frame(
            sample_id = names(pred),
            true_label = y_binary,
            predicted_prob = as.numeric(pred),
            stringsAsFactors = FALSE
        ),
        auc = auc_value,
        roc = roc_obj
    )

    cat(sprintf("okay external validation completed: %s, n=%d, AUC=%s\n",
                cohort_name, result$n,
                ifelse(is.na(auc_value), "NA", sprintf("%.3f", auc_value))))
    invisible(result)
}
