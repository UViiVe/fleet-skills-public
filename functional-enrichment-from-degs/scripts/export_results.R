# export_results.R
#
# Purpose:
#   Write the enrichment outputs to disk: CSVs of result tables, RDS of
#   the cluster Profiler objects, and a Markdown summary with the top hits.
#
# Inputs:
#   gsea       : gseaResult (optional)
#   ranked_genes : named numeric vector used for GSEA (optional)
#   ora_up     : enrichResult for upregulated genes (optional)
#   ora_down   : enrichResult for downregulated genes (optional)
#   output_dir : directory for all outputs (default "enrich_out")
#
# Outputs:
#   gsea_results.csv             gsea_result.rds
#   ora_up_results.csv           ora_up_result.rds
#   ora_down_results.csv         ora_down_result.rds
#   ranked_genes.rds             enrichment_summary.md
#
# Example:
#   source("scripts/export_results.R")
#   export_all(gsea, ranked,
#              ora_up = ora_up, ora_down = ora_down,
#              output_dir = "enrich_out")

.write_enrichment_csv <- function(obj, path) {
  if (is.null(obj)) return(invisible(NULL))
  df <- as.data.frame(obj)
  if (nrow(df) == 0) {
    cat(sprintf("   %s: 0 rows (no significant enrichment)\n", path))
    return(invisible(NULL))
  }
  write.csv(df, path, row.names = FALSE)
  cat(sprintf("   wrote %s  (%d rows)\n", path, nrow(df)))
}

.format_top_pathways <- function(obj, top_n = 10, label = "") {
  if (is.null(obj)) return("(no results)")
  df <- as.data.frame(obj)
  if (nrow(df) == 0) return("(no significant pathways)")
  df <- head(df[order(df$p.adjust), ], top_n)
  rows <- character(0)
  if ("NES" %in% colnames(df)) {
    for (i in seq_len(nrow(df))) {
      rows[i] <- sprintf("| %s | %.2f | %.2e |",
                         df$ID[i], df$NES[i], df$p.adjust[i])
    }
    header <- "| Pathway | NES | p.adjust |\n|---|---|---|"
  } else {
    for (i in seq_len(nrow(df))) {
      rows[i] <- sprintf("| %s | %s | %.2e |",
                         df$ID[i], df$GeneRatio[i], df$p.adjust[i])
    }
    header <- "| Pathway | GeneRatio | p.adjust |\n|---|---|---|"
  }
  paste(c(header, rows), collapse = "\n")
}

export_all <- function(gsea         = NULL,
                       ranked_genes = NULL,
                       ora_up       = NULL,
                       ora_down     = NULL,
                       output_dir   = "enrich_out") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  cat(sprintf("=== exporting to %s/ ===\n", output_dir))

  # CSVs
  .write_enrichment_csv(gsea,     file.path(output_dir, "gsea_results.csv"))
  .write_enrichment_csv(ora_up,   file.path(output_dir, "ora_up_results.csv"))
  .write_enrichment_csv(ora_down, file.path(output_dir, "ora_down_results.csv"))

  # RDS
  if (!is.null(gsea))         { saveRDS(gsea,         file.path(output_dir, "gsea_result.rds"));     cat("   wrote gsea_result.rds\n") }
  if (!is.null(ora_up))       { saveRDS(ora_up,       file.path(output_dir, "ora_up_result.rds"));   cat("   wrote ora_up_result.rds\n") }
  if (!is.null(ora_down))     { saveRDS(ora_down,     file.path(output_dir, "ora_down_result.rds")); cat("   wrote ora_down_result.rds\n") }
  if (!is.null(ranked_genes)) { saveRDS(ranked_genes, file.path(output_dir, "ranked_genes.rds"));    cat("   wrote ranked_genes.rds\n") }

  # Summary
  sections <- character(0)
  sections <- c(sections, "# Enrichment Summary\n")
  sections <- c(sections, sprintf("- Date: %s", Sys.Date()))
  if (!is.null(ranked_genes))
    sections <- c(sections, sprintf("- Ranked genes: %d", length(ranked_genes)))
  if (!is.null(gsea)) {
    df <- as.data.frame(gsea)
    sections <- c(sections,
      sprintf("\n## GSEA (%d significant pathways)\n", sum(df$p.adjust < 0.05)),
      .format_top_pathways(gsea, top_n = 10, label = "GSEA"))
  }
  if (!is.null(ora_up)) {
    df <- as.data.frame(ora_up)
    sections <- c(sections,
      sprintf("\n## ORA upregulated (%d significant pathways)\n",
              sum(df$p.adjust < 0.05)),
      .format_top_pathways(ora_up, top_n = 10, label = "ORA-up"))
  }
  if (!is.null(ora_down)) {
    df <- as.data.frame(ora_down)
    sections <- c(sections,
      sprintf("\n## ORA downregulated (%d significant pathways)\n",
              sum(df$p.adjust < 0.05)),
      .format_top_pathways(ora_down, top_n = 10, label = "ORA-down"))
  }
  writeLines(sections, file.path(output_dir, "enrichment_summary.md"))
  cat("   wrote enrichment_summary.md\n")

  cat(sprintf("[OK] export complete: %s/\n", output_dir))
  invisible(NULL)
}
