# prepare_features.R
# Build the (X, y) pair used by run_lasso_panel from an expression matrix and
# metadata table. Optionally filters by differential expression first.
#
# Exposes:
#   prepare_feature_matrix(expression, metadata, outcome_col, top_n_variable = 2000)
#   filter_by_de(expression, metadata, outcome_col, padj_threshold = 0.05)
#   prepare_validation_features(validation_expression, discovery_features)


#' Build the standardized feature matrix and outcome vector for LASSO.
#'
#' Steps:
#'   1. Align expression columns with metadata rows.
#'   2. Coerce the outcome to 0/1.
#'   3. Drop genes with NA or zero variance.
#'   4. Optionally keep the top N most variable features.
#'   5. Standardize features to mean 0, SD 1 (LASSO is scale-sensitive).
#'
#' @param expression Matrix or data.frame: features (rows) x samples (columns).
#' @param metadata data.frame; row.names must match colnames(expression).
#' @param outcome_col Name of the binary outcome column in metadata.
#' @param top_n_variable Cap on the number of features kept (NULL = no cap).
#' @return list(X, y, feature_names, n_samples, n_features, outcome_levels)
prepare_feature_matrix <- function(expression,
                                   metadata,
                                   outcome_col,
                                   top_n_variable = 2000) {

    if (!outcome_col %in% colnames(metadata)) {
        stop(sprintf("Outcome column '%s' not found in metadata.", outcome_col))
    }

    expression <- as.matrix(expression)
    storage.mode(expression) <- "double"

    samples <- intersect(colnames(expression), rownames(metadata))
    if (length(samples) == 0L) {
        stop("No shared samples between expression columns and metadata rows.")
    }
    expression <- expression[, samples, drop = FALSE]
    metadata <- metadata[samples, , drop = FALSE]

    y_raw <- metadata[[outcome_col]]
    levels_seen <- sort(unique(na.omit(y_raw)))
    if (length(levels_seen) != 2L) {
        stop(sprintf("Outcome must be binary; found levels: %s",
                     paste(levels_seen, collapse = ", ")))
    }
    y <- ifelse(y_raw == levels_seen[2], 1L, 0L)
    keep <- !is.na(y)
    if (sum(keep) < length(y)) {
        cat(sprintf("  Dropping %d samples with missing outcome\n",
                    length(y) - sum(keep)))
    }
    expression <- expression[, keep, drop = FALSE]
    y <- y[keep]

    variances <- apply(expression, 1, var, na.rm = TRUE)
    keep_features <- !is.na(variances) & variances > 0
    expression <- expression[keep_features, , drop = FALSE]
    variances <- variances[keep_features]

    if (!is.null(top_n_variable) && nrow(expression) > top_n_variable) {
        top <- order(variances, decreasing = TRUE)[seq_len(top_n_variable)]
        expression <- expression[top, , drop = FALSE]
    }

    X <- t(expression)
    X <- scale(X, center = TRUE, scale = TRUE)
    X[is.nan(X)] <- 0
    rownames(X) <- colnames(expression)

    cat(sprintf("okay prepared features: %d samples x %d features\n",
                nrow(X), ncol(X)))
    cat(sprintf("  outcome levels: %s -> 0, %s -> 1\n",
                levels_seen[1], levels_seen[2]))
    cat(sprintf("  class balance:  %d positives, %d negatives\n",
                sum(y == 1L), sum(y == 0L)))

    invisible(list(
        X              = X,
        y              = y,
        feature_names  = colnames(X),
        n_samples      = nrow(X),
        n_features     = ncol(X),
        outcome_levels = levels_seen
    ))
}


#' Pre-filter features by univariate t-test / Wilcoxon p-value.
#'
#' Useful to reduce a giant matrix before LASSO. Caller is responsible for
#' avoiding circularity (filter on an independent slice, not the LASSO data).
#'
#' @param expression Feature matrix (features x samples).
#' @param metadata data.frame aligned to expression columns.
#' @param outcome_col Binary outcome column name.
#' @param padj_threshold BH-adjusted p threshold (q) for keeping a feature.
#' @return character vector of feature names that passed the q-value threshold.
filter_by_de <- function(expression, metadata, outcome_col,
                         padj_threshold = 0.05) {
    samples <- intersect(colnames(expression), rownames(metadata))
    expression <- expression[, samples, drop = FALSE]
    metadata <- metadata[samples, , drop = FALSE]

    y <- metadata[[outcome_col]]
    pos <- which(y == sort(unique(y))[2])
    neg <- which(y == sort(unique(y))[1])

    pvals <- apply(expression, 1, function(row) {
        tryCatch(t.test(row[pos], row[neg])$p.value,
                 error = function(e) NA_real_)
    })
    q <- p.adjust(pvals, method = "BH")
    keep <- names(q)[q < padj_threshold & !is.na(q)]

    cat(sprintf("okay filter_by_de retained %d / %d features (q < %g)\n",
                length(keep), length(q), padj_threshold))
    keep
}


#' Prepare the validation X matrix for an already-locked panel.
#'
#' Subsets the validation expression matrix to features present in the
#' discovery feature list and applies the same per-feature scaling as
#' the discovery pipeline.
#'
#' @param validation_expression features x samples matrix.
#' @param discovery_features Character vector from features$feature_names.
#' @return list(X, n_matched, feature_names)
prepare_validation_features <- function(validation_expression,
                                        discovery_features) {
    validation_expression <- as.matrix(validation_expression)
    storage.mode(validation_expression) <- "double"

    available <- intersect(discovery_features, rownames(validation_expression))
    if (length(available) == 0L) {
        stop("None of the discovery features are present in the validation matrix.")
    }
    sub <- validation_expression[available, , drop = FALSE]
    X <- t(sub)
    X <- scale(X, center = TRUE, scale = TRUE)
    X[is.nan(X)] <- 0

    cat(sprintf("okay prepared validation features: %d / %d discovery features matched\n",
                length(available), length(discovery_features)))
    invisible(list(
        X             = X,
        n_matched     = length(available),
        feature_names = available
    ))
}
