# load_example_data.R
# Small synthetic demos for the public quickstart.

make_synthetic_biomarker_data <- function(n_samples = 120,
                                          n_features = 500,
                                          n_signal = 8,
                                          outcome = "endotype",
                                          seed = 42) {
    set.seed(seed)
    sample_ids <- sprintf("sample_%03d", seq_len(n_samples))
    feature_ids <- sprintf("gene_%04d", seq_len(n_features))

    y <- rbinom(n_samples, 1, 0.5)
    expression <- matrix(
        rnorm(n_features * n_samples),
        nrow = n_features,
        ncol = n_samples,
        dimnames = list(feature_ids, sample_ids)
    )
    signal_features <- seq_len(min(n_signal, n_features))
    expression[signal_features, y == 1] <- expression[signal_features, y == 1] + 1.2

    metadata <- data.frame(
        sample_id = sample_ids,
        outcome = factor(ifelse(y == 1, "case", "control")),
        endotype = factor(ifelse(y == 1, "endotype_high", "endotype_low")),
        mortality = factor(ifelse(y == 1, "event", "no_event")),
        subtype = factor(ifelse(y == 1, "subtype_a", "subtype_b")),
        stringsAsFactors = FALSE,
        row.names = sample_ids
    )

    list(
        expression = expression,
        metadata = metadata,
        outcome_col = outcome,
        source = "synthetic public demo"
    )
}

load_sepsis_data <- function(outcome = "endotype") {
    make_synthetic_biomarker_data(n_samples = 120, n_features = 500,
                                  outcome = outcome, seed = 101)
}

load_breast_cancer_pcr_data <- function(outcome = "subtype") {
    make_synthetic_biomarker_data(n_samples = 100, n_features = 400,
                                  outcome = outcome, seed = 202)
}

load_imvigor210_data <- function() {
    make_synthetic_biomarker_data(n_samples = 110, n_features = 450,
                                  outcome = "outcome", seed = 303)
}

load_unifi_data <- function() {
    make_synthetic_biomarker_data(n_samples = 130, n_features = 450,
                                  outcome = "outcome", seed = 404)
}

load_breast_cancer_validation_data <- function(cohort_name = "synthetic-validation") {
    data <- make_synthetic_biomarker_data(n_samples = 80, n_features = 400,
                                          outcome = "subtype", seed = 505)
    data$cohort_name <- cohort_name
    data
}
