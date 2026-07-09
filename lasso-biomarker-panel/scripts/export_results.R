# export_results.R
# Persist LASSO biomarker-panel outputs.

export_all <- function(model,
                       output_dir = "results",
                       data = NULL,
                       features = NULL,
                       validation_result = NULL) {
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    saveRDS(model, file.path(output_dir, "lasso_model.rds"))

    if (!is.null(model$feature_importance)) {
        write.csv(model$feature_importance,
                  file.path(output_dir, "feature_importance.csv"),
                  row.names = FALSE)
    }

    if (!is.null(model$stable_features)) {
        write.csv(data.frame(feature = model$stable_features),
                  file.path(output_dir, "stable_features.csv"),
                  row.names = FALSE)
    }

    if (!is.null(model$cv_predictions)) {
        write.csv(model$cv_predictions,
                  file.path(output_dir, "cv_predictions.csv"),
                  row.names = FALSE)
    }

    summary <- data.frame(
        metric = c("mean_auc", "auc_ci_low", "auc_ci_high", "panel_size"),
        value = c(
            model$mean_auc,
            model$auc_ci[[1]],
            model$auc_ci[[2]],
            length(model$stable_features)
        )
    )
    write.csv(summary, file.path(output_dir, "performance_summary.csv"),
              row.names = FALSE)

    if (!is.null(validation_result)) {
        saveRDS(validation_result, file.path(output_dir, "validation_result.rds"))
        write.csv(
            data.frame(
                cohort = validation_result$cohort_name,
                auc = validation_result$auc,
                n = validation_result$n
            ),
            file.path(output_dir, "validation_summary.csv"),
            row.names = FALSE
        )
    }

    cat("=== Export Complete ===\n")
    invisible(output_dir)
}
