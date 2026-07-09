# get_msigdb_genesets.R
#
# Purpose:
#   Wrapper around msigdbr that returns a TERM2GENE data.frame ready for
#   clusterProfiler's GSEA() or enricher(). Supports the standard MSigDB
#   collections (H, KEGG, REACTOME, GO:BP/MF/CC, C6, C7) for human or mouse.
#
# Inputs:
#   species    : "human" (default) or "mouse"
#   categories : character vector of MSigDB collection keys. Accepted forms:
#                  "H", "C2:CP:KEGG", "C2:CP:REACTOME", "C5:GO:BP",
#                  "C5:GO:MF", "C5:GO:CC", "C6", "C7".
#                Or shorthand: "KEGG", "REACTOME", "GO:BP", etc.
#
# Outputs:
#   data.frame with columns (term, gene); suitable for clusterProfiler.
#
# Example:
#   source("scripts/get_msigdb_genesets.R")
#   tg <- get_msigdb_genesets(species = "human", categories = c("H", "C2:CP:REACTOME"))

.ensure_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (length(getOption("repos")) == 0 || identical(unname(getOption("repos")["CRAN"]), "@CRAN@"))
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    install.packages(pkg)
  }
}

.normalize_category <- function(cat) {
  cat_upper <- toupper(cat)
  switch(cat_upper,
    "H"               = list(category = "H",  subcategory = NULL),
    "KEGG"            = list(category = "C2", subcategory = "CP:KEGG"),
    "C2:CP:KEGG"      = list(category = "C2", subcategory = "CP:KEGG"),
    "REACTOME"        = list(category = "C2", subcategory = "CP:REACTOME"),
    "C2:CP:REACTOME"  = list(category = "C2", subcategory = "CP:REACTOME"),
    "C2:CP:BIOCARTA"  = list(category = "C2", subcategory = "CP:BIOCARTA"),
    "BIOCARTA"        = list(category = "C2", subcategory = "CP:BIOCARTA"),
    "GO:BP"           = list(category = "C5", subcategory = "GO:BP"),
    "C5:GO:BP"        = list(category = "C5", subcategory = "GO:BP"),
    "GO:MF"           = list(category = "C5", subcategory = "GO:MF"),
    "C5:GO:MF"        = list(category = "C5", subcategory = "GO:MF"),
    "GO:CC"           = list(category = "C5", subcategory = "GO:CC"),
    "C5:GO:CC"        = list(category = "C5", subcategory = "GO:CC"),
    "C6"              = list(category = "C6", subcategory = NULL),
    "C7"              = list(category = "C7", subcategory = NULL),
    stop("Unsupported MSigDB category: ", cat)
  )
}

get_msigdb_genesets <- function(species    = "human",
                                categories = c("H")) {
  .ensure_pkg("msigdbr")
  library(msigdbr)

  species_full <- switch(tolower(species),
    "human" = "Homo sapiens",
    "mouse" = "Mus musculus",
    stop("Unsupported species: ", species)
  )

  parts <- lapply(categories, .normalize_category)
  rows  <- lapply(parts, function(p) {
    args <- list(species = species_full, category = p$category)
    if (!is.null(p$subcategory)) args$subcategory <- p$subcategory
    do.call(msigdbr, args)
  })

  combined <- do.call(rbind, rows)
  out <- unique(combined[, c("gs_name", "gene_symbol")])
  colnames(out) <- c("term", "gene")
  out <- out[!is.na(out$gene) & out$gene != "", ]

  cat(sprintf("[OK] %d unique gene-set/gene rows; %d distinct gene sets; species=%s; categories=%s\n",
              nrow(out), length(unique(out$term)), species,
              paste(categories, collapse = ", ")))
  invisible(out)
}
