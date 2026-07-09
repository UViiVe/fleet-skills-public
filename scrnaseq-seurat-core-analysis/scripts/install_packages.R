# install_packages.R
# ----------------------------------------------------------------------
# Install all packages needed for the scRNA-seq Seurat workflow.
#
# Run once after cloning the skill:
#   source("scripts/install_packages.R")
#   install_seurat_skill()

#' Install CRAN, Bioconductor, and GitHub dependencies in order.
install_seurat_skill <- function(force = FALSE) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))

  cran_pkgs <- c(
    "Seurat", "SeuratObject", "Matrix", "ggplot2", "patchwork", "dplyr",
    "ggrepel", "harmony", "RANN", "cluster", "remotes", "BiocManager"
  )
  bioc_pkgs <- c("DESeq2", "muscat", "SingleR", "celldex", "SummarizedExperiment")
  github_pkgs <- list(
    DoubletFinder = "chris-mcginnis-ucsf/DoubletFinder"
  )

  need_cran <- cran_pkgs[force | !vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need_cran) > 0) {
    message("Installing CRAN: ", paste(need_cran, collapse = ", "))
    install.packages(need_cran)
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  need_bioc <- bioc_pkgs[force | !vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need_bioc) > 0) {
    message("Installing Bioconductor: ", paste(need_bioc, collapse = ", "))
    BiocManager::install(need_bioc, ask = FALSE, update = FALSE)
  }

  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  for (pkg in names(github_pkgs)) {
    if (force || !requireNamespace(pkg, quietly = TRUE)) {
      message("Installing GitHub: ", pkg, " (", github_pkgs[[pkg]], ")")
      remotes::install_github(github_pkgs[[pkg]])
    }
  }

  message("✓ Installation completed")
  invisible(NULL)
}
