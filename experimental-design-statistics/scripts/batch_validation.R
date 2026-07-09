# batch_validation.R
# Validate batch designs before any sample is processed.
# Confounding test (chi-square) is the cardinal check.
#
# Inputs:
#   - batch_assignment: data frame with batch + condition + (optional) covariates.
# Outputs:
#   - PASS/FAIL on confounding; balance summary; PNG+SVG visualization.

#' Save a plot in both PNG and SVG formats (with graceful svglite fallback).
#'
#' Tries ggsave for both formats; if SVG via ggsave fails (often missing svglite),
#' falls back to base R svg(). PNG is always written.
.save_plot <- function(plot, base_path, width = 8, height = 6, dpi = 300) {
  png_path <- sub("\\.(svg|png)$", ".png", base_path)
  if (!grepl("\\.png$", png_path)) png_path <- paste0(png_path, ".png")
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height,
                  dpi = dpi, device = "png")
  cat("   saved:", png_path, "\n")

  svg_path <- sub("\\.(svg|png)$", ".svg", base_path)
  if (!grepl("\\.svg$", svg_path)) svg_path <- paste0(svg_path, ".svg")
  ok <- tryCatch({
    ggplot2::ggsave(svg_path, plot = plot, width = width, height = height, device = "svg")
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    ok2 <- tryCatch({
      grDevices::svg(svg_path, width = width, height = height)
      print(plot)
      grDevices::dev.off()
      TRUE
    }, error = function(e) FALSE)
    if (!ok2) cat("   (SVG export failed; PNG written)\n")
    else      cat("   saved:", svg_path, " (via base svg)\n")
  } else {
    cat("   saved:", svg_path, "\n")
  }
}

#' Test whether batch is confounded with condition (the cardinal check)
#'
#' @param batch_assignment Data frame with batch + condition columns
#' @param condition_var Name of the condition column
#' @param batch_var Name of the batch column (default "batch")
#' @return list with status PASS/FAIL, p-value, contingency table, message
check_confounding <- function(batch_assignment, condition_var,
                              batch_var = "batch") {
  if (!condition_var %in% colnames(batch_assignment))
    stop(sprintf("Condition '%s' not found", condition_var))
  if (!batch_var %in% colnames(batch_assignment))
    stop(sprintf("Batch '%s' not found", batch_var))

  ct <- table(batch_assignment[[batch_var]],
              batch_assignment[[condition_var]])
  chi <- suppressWarnings(chisq.test(ct))
  confounded <- chi$p.value < 0.05

  if (confounded) {
    msg <- sprintf(
      "WARNING: Batch is CONFOUNDED with %s\nChi-square p = %s < 0.05\nRegenerate the batch assignment before processing samples.",
      condition_var, format(chi$p.value, digits = 4)
    )
    status <- "FAIL"
  } else {
    msg <- sprintf(
      "PASS: No confounding detected between batch and %s\nChi-square p = %s > 0.05",
      condition_var, format(chi$p.value, digits = 4)
    )
    status <- "PASS"
  }
  cat(msg, "\n")
  invisible(list(status = status, is_confounded = confounded,
                 p_value = chi$p.value, contingency_table = ct,
                 message = msg, test = chi))
}

#' Test covariate balance across batches
#'
#' @param batch_assignment Data frame with batch + covariate columns
#' @param covariates Character vector of covariate column names
#' @param batch_var Batch column name (default "batch")
#' @return data frame: covariate, chi_square_p, is_balanced, interpretation
check_balance <- function(batch_assignment, covariates, batch_var = "batch") {
  if (!batch_var %in% colnames(batch_assignment))
    stop(sprintf("Batch '%s' not found", batch_var))
  missing <- setdiff(covariates, colnames(batch_assignment))
  if (length(missing) > 0)
    stop(sprintf("Covariates not found: %s", paste(missing, collapse = ", ")))

  rows <- lapply(covariates, function(co) {
    ct <- table(batch_assignment[[batch_var]], batch_assignment[[co]])
    chi <- suppressWarnings(chisq.test(ct))
    balanced <- chi$p.value > 0.05
    interp <- if (balanced) "Well balanced"
              else if (chi$p.value > 0.01) "Moderately balanced"
              else "Poorly balanced — consider regenerating"
    data.frame(covariate = co, chi_square_p = chi$p.value,
               is_balanced = balanced, interpretation = interp,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)

  cat("\nCovariate balance summary:\n")
  for (i in seq_len(nrow(out)))
    cat(sprintf("  %s: %s (p = %.4f)\n",
                out$covariate[i], out$interpretation[i],
                out$chi_square_p[i]))
  invisible(out)
}

#' Visualize the batch design as a bar chart
#'
#' @param batch_assignment Data frame with batch + condition columns
#' @param condition_var Condition column name
#' @param batch_var Batch column name
#' @param output_file Base output path (PNG + SVG written)
visualize_batch_design <- function(batch_assignment, condition_var,
                                   batch_var = "batch", output_file = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("ggprism", quietly = TRUE))
    stop("Install ggplot2 and ggprism")

  df <- as.data.frame(table(
    Batch = batch_assignment[[batch_var]],
    Condition = batch_assignment[[condition_var]]
  ))

  p <- ggplot2::ggplot(df, ggplot2::aes(x = Batch, y = Freq, fill = Condition)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggprism::theme_prism() +
    ggplot2::labs(title = "Batch design: samples per batch by condition",
                  x = "Batch", y = "Number of samples", fill = "Condition") +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  if (!is.null(output_file)) {
    cat("Saving batch design plot:\n")
    .save_plot(p, output_file)
  } else {
    print(p)
  }
  invisible(p)
}

#' Comprehensive batch validation report (confounding + balance + summary)
#'
#' @param batch_assignment Data frame with batch + condition + covariate columns
#' @param condition_var Condition column name
#' @param covariates Optional vector of additional covariates to check
#' @param output_file Optional Markdown file to save report
generate_batch_report <- function(batch_assignment, condition_var,
                                  covariates = NULL, output_file = NULL) {
  cat("\n=== Confounding check (critical) ===\n")
  conf <- check_confounding(batch_assignment, condition_var)

  bal <- NULL
  if (!is.null(covariates)) {
    cat("\n=== Covariate balance ===\n")
    bal <- check_balance(batch_assignment, covariates)
  }

  cat("\n=== Sample distribution ===\n")
  sample_summary <- as.data.frame(table(
    Batch = batch_assignment$batch,
    Condition = batch_assignment[[condition_var]]
  ))
  print(sample_summary)

  cat("\n=== Overall assessment ===\n")
  if (conf$status == "PASS") {
    cat("PASS: design has no batch-condition confounding.\n")
    if (!is.null(bal)) {
      if (all(bal$is_balanced)) cat("All covariates well balanced.\n")
      else cat("Some covariates poorly balanced: ",
               paste(bal$covariate[!bal$is_balanced], collapse = ", "),
               " — include as covariates in the analysis model.\n", sep = "")
    }
  } else {
    cat("FAIL: confounding detected. Do not proceed.\n")
  }

  if (!is.null(output_file)) {
    lines <- c(
      "# Batch design validation report",
      sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
      "",
      "## Confounding check",
      conf$message,
      ""
    )
    if (!is.null(bal)) {
      lines <- c(lines, "## Covariate balance",
                 paste(capture.output(print(bal)), collapse = "\n"))
    }
    writeLines(lines, output_file)
    message(sprintf("Report written: %s", output_file))
  }

  invisible(list(confounding = conf, balance = bal,
                 sample_summary = sample_summary))
}

#' Check that recommended metadata fields are present
check_metadata_completeness <- function(batch_assignment) {
  recommended <- c("sample_id","batch","processing_date","operator",
                   "reagent_lot","instrument","plate","well",
                   "library_prep_batch","sequencing_run","flowcell")
  present <- intersect(recommended, colnames(batch_assignment))
  missing <- setdiff(recommended, present)

  cat(sprintf("\nMetadata completeness: %d / %d recommended fields present\n",
              length(present), length(recommended)))
  if (length(missing) > 0)
    cat(sprintf("Missing: %s\n", paste(missing, collapse = ", ")))
  else
    cat("All recommended fields present.\n")

  invisible(list(present = present, missing = missing,
                 completeness = length(present) / length(recommended)))
}
