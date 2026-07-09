# biomarker_plots.R
# ROC, stability, coefficient forest, calibration, AUC distribution, panel
# heatmap. Each plot saves as PNG and SVG with a graceful svglite fallback.

suppressPackageStartupMessages({
    library(ggplot2)
    library(ggprism)
    library(ggrepel)
    library(pROC)
    library(ComplexHeatmap)
    library(circlize)
})

.has_svglite <- requireNamespace("svglite", quietly = TRUE)
if (.has_svglite) suppressPackageStartupMessages(library(svglite))


.save_ggplot <- function(plot, base, width = 8, height = 6, dpi = 300) {
    png_path <- paste0(base, ".png")
    svg_path <- paste0(base, ".svg")
    ggsave(png_path, plot = plot, width = width, height = height,
           dpi = dpi, device = "png")
    cat(sprintf("   saved %s\n", basename(png_path)))
    saved <- tryCatch({
        ggsave(svg_path, plot = plot, width = width, height = height,
               device = "svg")
        TRUE
    }, error = function(e) FALSE)
    if (!saved) {
        saved <- tryCatch({
            grDevices::svg(svg_path, width = width, height = height)
            print(plot); grDevices::dev.off()
            TRUE
        }, error = function(e) FALSE)
    }
    if (!saved) cat("   (SVG export skipped)\n")
}


.save_heatmap <- function(ht, base, width = 9, height = 7, dpi = 300) {
    png_path <- paste0(base, ".png")
    grDevices::png(png_path, width = width, height = height,
                   units = "in", res = dpi)
    ComplexHeatmap::draw(ht); grDevices::dev.off()
    cat(sprintf("   saved %s\n", basename(png_path)))
    tryCatch({
        grDevices::svg(paste0(base, ".svg"), width = width, height = height)
        ComplexHeatmap::draw(ht); grDevices::dev.off()
    }, error = function(e) cat("   (SVG export skipped)\n"))
}


#' ROC curve(s) — discovery + optional validation overlay.
plot_roc_curve <- function(model_result, validation_result = NULL,
                           output_dir = "results") {
    cv_pred <- model_result$cv_predictions
    roc_disc <- pROC::roc(cv_pred$true_label, cv_pred$predicted_prob,
                          quiet = TRUE)
    disc_df <- data.frame(
        fpr = 1 - roc_disc$specificities,
        tpr = roc_disc$sensitivities,
        cohort = sprintf("Discovery (CV AUC = %.3f)", model_result$mean_auc)
    )
    roc_df <- disc_df
    if (!is.null(validation_result) && !is.null(validation_result$roc)) {
        roc_val <- validation_result$roc
        val_df <- data.frame(
            fpr = 1 - roc_val$specificities,
            tpr = roc_val$sensitivities,
            cohort = sprintf("Validation (AUC = %.3f)",
                              validation_result$auc)
        )
        roc_df <- rbind(disc_df, val_df)
    }
    p <- ggplot(roc_df, aes(fpr, tpr, color = cohort)) +
        geom_line(size = 1.2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    color = "grey60") +
        scale_color_manual(values = c("#E74C3C", "#3498DB")) +
        theme_prism(base_size = 12) +
        coord_equal() +
        labs(title = "ROC curve", x = "False positive rate",
             y = "True positive rate", color = NULL)
    .save_ggplot(p, file.path(output_dir, "roc_curve"))
}


#' Feature selection-frequency bar plot.
plot_stability_barplot <- function(model_result, output_dir = "results") {
    fi <- model_result$feature_importance
    top <- head(fi, 25)
    top$feature <- factor(top$feature, levels = rev(top$feature))
    p <- ggplot(top, aes(selection_frequency, feature,
                          fill = is_stable)) +
        geom_col() +
        scale_fill_manual(values = c(`TRUE` = "#27AE60", `FALSE` = "#7F8C8D")) +
        geom_vline(xintercept = model_result$parameters$stability_threshold,
                   linetype = "dashed", color = "grey40") +
        theme_prism(base_size = 12) +
        labs(title = "Feature selection frequency", x = "fraction of folds",
             y = NULL, fill = "stable")
    .save_ggplot(p, file.path(output_dir, "stability_barplot"),
                 width = 8, height = 7)
}


#' LASSO coefficient forest with ±1 SD bars.
plot_coefficient_forest <- function(model_result, output_dir = "results") {
    fi <- model_result$feature_importance
    sub <- fi[fi$is_stable, ]
    if (nrow(sub) == 0L) {
        cat("   no stable features — skipping forest plot\n"); return()
    }
    sub <- sub[order(sub$mean_coefficient), ]
    sub$feature <- factor(sub$feature, levels = sub$feature)
    p <- ggplot(sub, aes(mean_coefficient, feature)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        geom_errorbarh(aes(xmin = mean_coefficient - sd_coefficient,
                            xmax = mean_coefficient + sd_coefficient),
                        height = 0.25, color = "grey30") +
        geom_point(size = 3, color = "#2E86C1") +
        theme_prism(base_size = 12) +
        labs(title = "LASSO coefficients (panel features)",
             x = "mean coefficient (+/- SD)", y = NULL)
    .save_ggplot(p, file.path(output_dir, "coefficient_forest"),
                 width = 8, height = max(4, 0.4 * nrow(sub) + 3))
}


#' Calibration curve via 10 probability bins.
plot_calibration_curve <- function(model_result, output_dir = "results") {
    cv_pred <- model_result$cv_predictions
    bins <- cut(cv_pred$predicted_prob,
                breaks = seq(0, 1, by = 0.1), include.lowest = TRUE)
    cal <- aggregate(cv_pred$true_label, list(bin = bins),
                     function(x) c(mean = mean(x), n = length(x)))
    cal_df <- data.frame(
        midpoint = seq(0.05, 0.95, by = 0.1),
        observed = as.numeric(cal$x[, "mean"]),
        n        = as.integer(cal$x[, "n"])
    )
    p <- ggplot(cal_df, aes(midpoint, observed, size = n)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    color = "grey60") +
        geom_point(color = "#E74C3C", alpha = 0.85) +
        scale_size_continuous(range = c(2, 8)) +
        theme_prism(base_size = 12) +
        labs(title = "Calibration curve",
             x = "Predicted probability bin midpoint",
             y = "Observed event rate", size = "n")
    .save_ggplot(p, file.path(output_dir, "calibration_curve"))
}


#' Distribution of CV AUC across outer folds.
plot_auc_distribution <- function(model_result, output_dir = "results") {
    aucs <- model_result$fold_aucs
    df <- data.frame(auc = aucs[!is.na(aucs)])
    p <- ggplot(df, aes(auc)) +
        geom_histogram(binwidth = 0.02, fill = "#3498DB", alpha = 0.85,
                       color = "white") +
        geom_vline(xintercept = model_result$mean_auc,
                   linetype = "dashed", color = "#E74C3C") +
        theme_prism(base_size = 12) +
        labs(title = sprintf("CV AUC distribution (mean %.3f, 95%% CI %.3f-%.3f)",
                              model_result$mean_auc,
                              model_result$auc_ci[1],
                              model_result$auc_ci[2]),
             x = "AUC", y = "outer folds")
    .save_ggplot(p, file.path(output_dir, "auc_distribution"))
}


#' Clustered heatmap of panel-feature expression by class.
plot_feature_heatmap <- function(X, y, stable_features,
                                  output_dir = "results") {
    avail <- intersect(stable_features, colnames(X))
    if (length(avail) == 0L) {
        cat("   no panel features in X — skipping heatmap\n"); return()
    }
    mat <- t(X[, avail, drop = FALSE])
    z <- (mat - rowMeans(mat)) / apply(mat, 1, sd)
    order_y <- order(y)
    z <- z[, order_y]
    y_ord <- y[order_y]
    col_anno <- ComplexHeatmap::HeatmapAnnotation(
        outcome = factor(y_ord, levels = c(0, 1)),
        col = list(outcome = c(`0` = "#3498DB", `1` = "#E74C3C")))
    col_fun <- circlize::colorRamp2(c(-2, 0, 2),
                                     c("#2166AC", "#F7F7F7", "#B2182B"))
    ht <- ComplexHeatmap::Heatmap(z, col = col_fun, name = "z",
                                  top_annotation = col_anno,
                                  show_column_names = FALSE,
                                  column_title = "Panel features")
    .save_heatmap(ht, file.path(output_dir, "feature_heatmap"),
                  width = 9, height = max(4, 0.3 * length(avail) + 3))
}


#' Orchestrator that produces every plot in the suite.
generate_all_plots <- function(model_result,
                               validation_result = NULL,
                               X = NULL, y = NULL,
                               output_dir = "results") {
    cat("\n=== Generating biomarker plots ===\n")
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    plot_roc_curve(model_result, validation_result, output_dir)
    plot_stability_barplot(model_result, output_dir)
    plot_coefficient_forest(model_result, output_dir)
    plot_calibration_curve(model_result, output_dir)
    plot_auc_distribution(model_result, output_dir)
    if (!is.null(X) && !is.null(y)) {
        plot_feature_heatmap(X, y, model_result$stable_features, output_dir)
    }
    cat("\nokay All biomarker plots generated successfully!\n")
}
