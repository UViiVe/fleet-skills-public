#!/usr/bin/env Rscript

## glmGamPoi batch-corrected differential expression for pooled CRISPR screens.
##
## Reads a subset h5ad written by the python wrapper, fits an NB GLM per
## perturbation with a donor/batch fixed effect, and writes one CSV per
## perturbation.
##
## Usage:
##   Rscript run_glmgampoi.R <h5ad> <perturbations> <control> <gene_col> \
##                           <donor_col> <min_expr> <output_dir>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop("Usage: Rscript run_glmgampoi.R <h5ad> <perts_csv> <control> <gene_col> <donor_col> <min_expr> <out_dir>", call. = FALSE)
}

h5ad_file       <- args[1]
perturbations_s <- args[2]
control_group   <- args[3]
gene_col        <- args[4]
donor_col       <- args[5]
min_expression  <- as.numeric(args[6])
output_dir      <- args[7]

suppressPackageStartupMessages({
    library(glmGamPoi)
    library(SingleCellExperiment)
    library(zellkonverter)
})

cat("glmGamPoi DE — pooled CRISPR screen\n")
cat(sprintf("  input          : %s\n", h5ad_file))
cat(sprintf("  control group  : %s\n", control_group))
cat(sprintf("  gene column    : %s\n", gene_col))
cat(sprintf("  donor column   : %s\n", donor_col))
cat(sprintf("  min expression : %d\n", min_expression))
cat(sprintf("  output dir     : %s\n\n", output_dir))

if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
}

sce <- tryCatch(readH5AD(h5ad_file),
                error = function(e) stop(sprintf("Failed to read %s: %s", h5ad_file, conditionMessage(e)), call. = FALSE))
cat(sprintf("Loaded %d cells × %d genes\n", ncol(sce), nrow(sce)))

for (req in c(gene_col, donor_col)) {
    if (!(req %in% colnames(colData(sce)))) {
        stop(sprintf("Required column '%s' missing from cell metadata", req), call. = FALSE)
    }
}

perturbations <- strsplit(perturbations_s, ",")[[1]]
cat(sprintf("Testing %d perturbations\n\n", length(perturbations)))

for (idx in seq_along(perturbations)) {
    gene <- perturbations[idx]
    cat(sprintf("[%d/%d] %s\n", idx, length(perturbations), gene))

    tryCatch({
        keep <- colData(sce)[[gene_col]] %in% c(gene, control_group)
        sub  <- sce[, keep]
        n_p  <- sum(colData(sub)[[gene_col]] == gene)
        n_c  <- sum(colData(sub)[[gene_col]] == control_group)
        cat(sprintf("  cells: %d perturbed, %d control\n", n_p, n_c))
        if (n_p < 20 || n_c < 20) { cat("  skipping (insufficient cells)\n"); next }

        gtotals <- rowSums(counts(sub))
        sub     <- sub[gtotals > min_expression, ]
        cat(sprintf("  genes after expression filter: %d\n", nrow(sub)))
        if (nrow(sub) < 100) { cat("  skipping (too few genes)\n"); next }

        colData(sub)[[gene_col]] <- factor(
            colData(sub)[[gene_col]],
            levels = c(control_group, gene)
        )
        design  <- as.formula(sprintf("~ %s + %s - 1", gene_col, donor_col))
        cat("  fitting glm_gp …\n")
        fit <- glm_gp(sub, design = design, on_disk = FALSE, verbose = FALSE)

        gene_lvl <- sprintf("%s%s", gene_col, gene)
        ctrl_lvl <- sprintf("%s%s", gene_col, control_group)
        contrast <- as.formula(sprintf("~ `%s` - `%s`", gene_lvl, ctrl_lvl))
        cat("  testing differential expression …\n")
        res <- test_de(fit, contrast = contrast, sort_by = "pval", verbose = FALSE)

        out <- data.frame(
            gene             = res$name,
            log2fc           = res$lfc,
            pval             = res$pval,
            adj_pval         = res$adj_pval,
            mean_expression  = res$mean,
            stringsAsFactors = FALSE
        )
        n_sig <- sum(out$adj_pval < 0.05, na.rm = TRUE)
        cat(sprintf("  significant (FDR<0.05): %d\n", n_sig))
        out_file <- file.path(output_dir, sprintf("%s_glmgampoi.csv", gene))
        write.csv(out, out_file, row.names = FALSE, quote = FALSE)
        cat(sprintf("  wrote %s\n", out_file))
    }, error = function(e) {
        cat(sprintf("  error: %s\n", conditionMessage(e)))
    })
    cat("\n")
}

cat("glmGamPoi tier complete\n")
