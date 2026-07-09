# run_umap.R
# ----------------------------------------------------------------------
# Compute UMAP and/or t-SNE embeddings.
#
# Inputs:  Seurat object with a reduction (pca, harmony, integrated.cca).
# Outputs: Seurat object with reductions$umap (and optionally $tsne).
# Example:
#   source("scripts/run_umap.R")
#   seurat_obj <- run_umap_reduction(seurat_obj, dims = 1:30)

#' Run UMAP on the chosen reduction. Always sets a seed.
run_umap_reduction <- function(seurat_obj,
                                 dims = 1:30,
                                 reduction = "pca",
                                 min_dist = 0.3,
                                 spread = 1,
                                 seed = 42) {
  seurat_obj <- Seurat::RunUMAP(
    seurat_obj, reduction = reduction, dims = dims,
    min.dist = min_dist, spread = spread, seed.use = seed,
    verbose = FALSE
  )
  message(sprintf("✓ UMAP completed successfully (min.dist=%g)", min_dist))
  seurat_obj
}

#' Run t-SNE (slow; UMAP is preferred for production).
run_tsne_reduction <- function(seurat_obj,
                                 dims = 1:30,
                                 reduction = "pca",
                                 perplexity = 30,
                                 seed = 42) {
  seurat_obj <- Seurat::RunTSNE(
    seurat_obj, reduction = reduction, dims = dims,
    perplexity = perplexity, seed.use = seed,
    verbose = FALSE
  )
  message(sprintf("✓ t-SNE completed (perplexity=%g)", perplexity))
  seurat_obj
}

#' Run both UMAP and t-SNE.
run_both_reductions <- function(seurat_obj, dims = 1:30) {
  seurat_obj <- run_umap_reduction(seurat_obj, dims = dims)
  seurat_obj <- run_tsne_reduction(seurat_obj, dims = dims)
  seurat_obj
}
