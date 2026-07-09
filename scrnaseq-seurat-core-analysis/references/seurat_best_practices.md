# Seurat Best Practices

Defaults for the Seurat v5 pipeline. Maps directly to the scripts under `scripts/`.

## Normalization

### SCTransform (recommended)

```r
seurat_obj <- SCTransform(
  seurat_obj,
  vst.flavor = "v2",
  vars.to.regress = c("percent.mt"),
  verbose = FALSE
)
```

`run_sctransform()` wraps this. SCTransform fits a regularized negative-binomial GLM per gene, accounting for sequencing depth. It replaces `NormalizeData()` + `ScaleData()` + `FindVariableFeatures()` in one step and produces Pearson residuals that are good for clustering.

When to use:

- UMI-based data (10X, Drop-seq, inDrop).
- Multi-batch experiments (better variance stabilization).
- Most production workflows.

Caveat: memory-heavy on > 30k cells. Use 32+ GB RAM, or drop to LogNormalize.

### LogNormalize (classic)

```r
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 1e4)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
```

`run_lognormalize()` wraps this. Library-size + log1p, the original Seurat workflow.

When to use:

- Very large datasets (> 50k cells) where SCT runs out of memory.
- Smart-seq2 / full-length data.
- When you need exact reproducibility with older Seurat publications.

## Variable features

SCTransform selects HVGs automatically as part of the model. For LogNormalize:

```r
FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
```

`vst` (variance-stabilizing transform of the mean–variance relationship) is the default and works well. 2000 features is the convention; bump to 3000 for atlases, drop to 1000 for very small data.

Visualize: `VariableFeaturePlot(seurat_obj)` shows the mean–variance scatter. Look for a clean separation of the selected HVGs from the bulk.

## Dimensionality reduction

### Number of PCs

```r
seurat_obj <- RunPCA(seurat_obj, npcs = 50)
ElbowPlot(seurat_obj, ndims = 50)
```

`run_pca_analysis()` runs 50 PCs and prints cumulative variance for the first 10 / 20 / 30. `suggest_n_pcs()` returns a variance-knee–based recommendation.

**Never use < 15 PCs.** Standard is 20–30. The "elbow" heuristic chronically undercounts; trust `suggest_n_pcs()` over your eye.

Cumulative variance should reach 80–90% by your chosen `n_pcs`. If it doesn't, the HVG set or QC is the problem, not the PC count.

### JackStraw (optional)

`JackStraw()` provides a permutation-based test of PC significance. Useful as a cross-check, but slow (15–30 min on standard datasets) and not part of the default workflow.

## Clustering

```r
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, reduction = "pca")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, algorithm = 1)
```

`cluster_multiple_resolutions()` runs `FindNeighbors` once and `FindClusters` across multiple resolutions, storing results in `meta.data$RNA_snn_res.<n>` (or `SCT_snn_res.<n>` after SCTransform).

### Algorithm choice

- `algorithm = 1` — Louvain (default).
- `algorithm = 4` — Leiden (newer, well-connected guarantee). Requires `leiden` Python module via `reticulate`.

Both give similar results in practice. Stick with Louvain unless you have a reason.

### Resolution

There's no single right value. Run a sweep:

```r
for (res in c(0.4, 0.6, 0.8, 1.0)) {
  seurat_obj <- FindClusters(seurat_obj, resolution = res)
}
```

Choose by:

1. **Marker coherence.** Each cluster has ≥ 3 strong markers that make biological sense.
2. **Stability.** Adjacent resolutions give consistent assignments for most cells.
3. **Biology.** What's the downstream question? "How many CD4 substates?" needs finer resolution than "what are the major lineages?"

### Neighbor graph

`k.param = 20` is the Seurat default and works for ~ 95% of datasets. Drop to 10 for small datasets (< 5k cells), raise to 30 for very large ones.

Always use the same `dims` for `FindNeighbors`, `FindClusters` (implicit), and `RunUMAP` — discrepancies cause visualization artifacts.

## UMAP

```r
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, min.dist = 0.3, spread = 1, seed.use = 42)
```

`run_umap_reduction()` wraps this and always sets `seed.use`.

- `min.dist`: 0.1 → tight, 0.3 (default) → standard, 0.5+ → loose.
- `spread`: 1 is fine almost always.
- `seed.use`: always set. UMAP is non-deterministic.

UMAP is a visualization tool, not a quantitative coordinate system. Don't compute pseudo-distances on UMAP coordinates.

## Marker gene identification

```r
markers <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)
```

`find_all_cluster_markers()` wraps this and returns a tidy data frame.

### Test selection

- `wilcox` (default): non-parametric, robust, fast enough.
- `t` / `LR` (logistic regression): faster, can handle covariates.
- `MAST`: zero-inflated model for scRNA-seq, slower but principled.
- `bimod`: Seurat's likelihood-ratio test on bimodal expression.
- `DESeq2`: cell-level DESeq2 — slow and **not the right tool** for inferential condition contrasts (use pseudobulk).

For exploratory cluster characterization, Wilcoxon is the right default. None of these tests are valid for condition contrasts — use pseudobulk.

### Filtering

```r
top_per_cluster <- markers |>
  group_by(cluster) |>
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) |>
  slice_max(avg_log2FC, n = 10)
```

## Data management

### Assays in Seurat v5

- `RNA` — raw counts and (after `NormalizeData`) log-normalized.
- `SCT` — SCTransform output: `counts` (corrected), `data` (log of corrected), `scale.data` (Pearson residuals).
- `integrated` — output of `IntegrateData()` for CCA/RPCA workflows.

Keep RNA around even after SCTransform. Pseudobulk DE needs raw counts.

### Memory

- Use sparse matrices (`Matrix::CsparseMatrix`) — Seurat does this by default.
- For very large data, `BPCells` and `DiskAnn` provide on-disk variants.
- Drop intermediate assays you don't need: `DefaultAssay(seurat_obj) <- "RNA"; seurat_obj[["SCT"]] <- NULL`.

## Visualization style

Seurat's built-in `DimPlot`, `FeaturePlot`, `VlnPlot`, `DotPlot` are publication-quality enough for most figures. For polish:

- Combine `ggplot2` + `ggprism` for clean styling.
- `patchwork` for multi-panel layouts.
- `ggrepel` for non-overlapping labels in volcano plots.

Always save PNG and SVG. The skill's plot scripts include SVG fallback handling (graceful PNG-only when `svglite` isn't installed).

## Performance tips

- Set `options(future.globals.maxSize = 8 * 1024^3)` (8 GB) for big SCTransform / DESeq2 runs.
- Parallelize with `future::plan("multisession", workers = N)`.
- For very large datasets, use `BPCells` to keep counts on disk.

## Workflow checklist

- [ ] Loaded data, RNA assay populated with raw counts.
- [ ] Calculated QC metrics, visualized distributions.
- [ ] Filtered with MAD or tissue-fixed thresholds, retained ≥ 70%.
- [ ] Ran DoubletFinder (per batch for multi-batch).
- [ ] Chose SCTransform (default) or LogNormalize.
- [ ] Ran PCA, chose `n_pcs` from elbow / `suggest_n_pcs()`.
- [ ] Integrated batches (Harmony / CCA / RPCA) when multi-batch.
- [ ] Built neighbor graph on the correct reduction.
- [ ] Clustered at multiple resolutions, chose one.
- [ ] Ran UMAP with `seed.use`.
- [ ] Found markers (Wilcoxon, exploratory).
- [ ] Annotated cell types (manual ± SingleR/Azimuth).
- [ ] Pseudobulk DE if condition contrast required.
- [ ] Exported `.rds` + CSVs + figures via `export_all()`.

## Resources

- Seurat v5 vignettes: <https://satijalab.org/seurat/articles/>
- Single-cell best practices book: <https://www.sc-best-practices.org/>
- Seurat GitHub: <https://github.com/satijalab/seurat>

## Citations

1. Hao Y, Stuart T, Kowalski MH, et al. *Dictionary learning for integrative, multimodal and scalable single-cell analysis.* Nat Biotechnol 42:293–304 (2024).
2. Hafemeister C, Satija R. *Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression.* Genome Biol 20:296 (2019).
3. Choudhary S, Satija R. *Comparison and evaluation of statistical error models for scRNA-seq.* Genome Biol 23:27 (2022).
4. Traag VA, Waltman L, van Eck NJ. *From Louvain to Leiden: guaranteeing well-connected communities.* Sci Rep 9:5233 (2019).
