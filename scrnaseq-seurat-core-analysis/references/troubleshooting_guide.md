# Seurat Troubleshooting

Symptoms, causes, fixes. Grouped by phase.

## Installation

### `there is no package called 'Seurat'`

```r
install.packages("Seurat")
```

Ensure R ≥ 4.2. Older R hits dependency walls.

### `Error: package 'SeuratObject' could not be loaded`

Version mismatch. Reinstall together:

```r
install.packages(c("Seurat", "SeuratObject", "Matrix"))
```

### `installation of package 'irlba' had non-zero exit status`

Compiler issue. On macOS:

```bash
brew install gcc
```

On Linux:

```bash
sudo apt-get install r-base-dev libcurl4-openssl-dev libssl-dev
```

### `BiocManager::install("DESeq2")` fails

Update Bioconductor first:

```r
BiocManager::install(update = TRUE, ask = FALSE)
```

## Loading

### `Read10X` errors with "Barcode file missing"

CellRanger output directory must contain `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`. Verify:

```r
list.files("path/to/filtered_feature_bc_matrix/")
```

If your CellRanger version is old, file names may be `genes.tsv.gz` instead of `features.tsv.gz`. Symlink or rename.

### Out of memory loading a large matrix

Use `BPCells` for on-disk counts:

```r
library(BPCells)
mat <- open_matrix_dir("path/")
```

Or downsample for prototyping:

```r
seurat_obj <- subset(seurat_obj, cells = sample(colnames(seurat_obj), 10000))
```

## QC

### Cell retention < 50%

Visualize before deciding:

```r
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
```

Switch to MAD:

```r
seurat_obj <- batch_mad_outlier_detection(seurat_obj, batch_col = "batch", nmads = 7)
```

Or relax tissue thresholds — see `qc_guidelines.md`.

### High `percent.mt` everywhere

Check the distribution shape:

```r
hist(seurat_obj$percent.mt, breaks = 100)
```

- Bimodal: filter at the trough.
- Unimodal at 10–15%: tissue-specific (cardiomyocytes, hepatocytes). Don't filter blindly.
- Unimodal at > 20%: sample probably degraded.

### `Mitochondrial genes detected: 0`

Wrong species pattern or Ensembl IDs:

```r
head(rownames(seurat_obj))
grep("^MT-|^mt-", rownames(seurat_obj), value = TRUE)
```

Convert Ensembl IDs to symbols with `biomaRt`, `annotables`, or `clusterProfiler::bitr`.

### DoubletFinder predicts ~ 100% doublets

Expected rate misconfigured. Check the auto-rate is sensible (~ 0.8% per 1000 cells). If overriding, pass a fraction not a percent: 0.075 means 7.5%.

## Normalization

### SCTransform crashes with `memory not mapped`

RAM exhausted. Switch to LogNormalize:

```r
seurat_obj <- run_lognormalize(seurat_obj)
```

Or run on a server with > 64 GB RAM. SCTransform is the most memory-heavy step in the pipeline.

### `Error: variances must be positive`

A few cells have zero variance (all the same gene at all UMIs after subsetting). Filter first:

```r
seurat_obj <- subset(seurat_obj, subset = nCount_RNA > 100)
```

### Negative values after `ScaleData`

Expected — `ScaleData` z-scores the data. Use this for PCA, not for plotting expression.

## PCA and dimensionality

### `Cannot find variable features`

```r
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
```

Or for SCT-normalized data, variable features are stored in `SCTAssay`; ensure `DefaultAssay(seurat_obj) <- "SCT"` before `RunPCA`.

### PCA variance very low

Causes: too many lowly expressed genes, batch effects dominating, not using HVGs.

- Verify HVGs computed.
- `RunPCA(seurat_obj, features = VariableFeatures(seurat_obj))`.
- Color UMAP by batch — if batch dominates PC1–2, integrate before clustering.

### UMAP is one giant blob

- Too few PCs → `dims = 1:50`.
- Wrong reduction → if you ran Harmony, use `reduction = "harmony"` everywhere downstream.
- Resolution too low → bump to 1.0.

### UMAP changes every run

```r
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, seed.use = 42)
```

Always pass `seed.use`.

## Integration

### `Error: 'harmony' not found in reductions`

You skipped the integration step. Run `run_harmony_integration()` first, or fall back to `reduction = "pca"` for single-batch.

### Batches still cluster after "integration"

- Harmony: raise `theta` to 3 or 4.
- CCA: increase `anchor.features` to 3000.
- Switch methods — Harmony to CCA, or try scVI via `reticulate`.

### Lost biological variation

Over-corrected. Lower `theta` (Harmony) or reduce `dims` (CCA). Check that batch and biology aren't confounded.

## Clustering

### Too many small clusters (50+)

Lower resolution:

```r
seurat_obj <- FindClusters(seurat_obj, resolution = 0.4)
```

Or fewer neighbors:

```r
seurat_obj <- FindNeighbors(seurat_obj, k.param = 10)
```

### Clusters driven by QC metrics

- Tighten MT filter.
- Re-run SCTransform with `vars.to.regress = c("percent.mt", "nCount_RNA")`.
- Drop the problem cluster post-hoc:
  ```r
  seurat_obj <- subset(seurat_obj, idents = setdiff(levels(seurat_obj), "5"))
  ```

### Known cell types fail to separate

- Increase resolution sweep.
- Use more PCs (`dims = 1:50`).
- Check for batch effects (color UMAP by batch).
- Try integration if multi-batch.

## Marker identification

### No significant markers

- Clusters too similar — check `BuildClusterTree(seurat_obj)`.
- Relax thresholds: `logfc.threshold = 0.1`, `min.pct = 0.05`.
- Try `test.use = "LR"` or `"MAST"`.

### Top markers all `RPS*` / `RPL*` / `MT-*`

Filter before ranking:

```r
exclude <- grep("^(MT-|RPS|RPL)", rownames(seurat_obj), value = TRUE)
markers <- find_all_cluster_markers(
  seurat_obj, resolution = 0.8,
  features_to_exclude = exclude
)
```

## Plotting

### `Error: 'svglite' is not installed`

Optional. Either install `svglite`, or accept PNG-only output — the skill's plot scripts handle this gracefully.

### Blurry figures

```r
ggsave("figure.png", dpi = 300, width = 6, height = 5, units = "in")
ggsave("figure.svg", width = 6, height = 5, units = "in")
```

### `Error: Faceting variables must have at least one value`

`split.by` column doesn't exist or is all NA. Verify:

```r
table(seurat_obj@meta.data$batch, useNA = "ifany")
```

## SingleR / Azimuth

### SingleR returns mostly `NA`

Reference doesn't match the tissue. Pick the right one — `HPCA` for human blood/immune, `BlueprintEncode` for human tissue, `ImmGen` for mouse immune.

### Azimuth fails to download reference

Cache directory may be locked. Clear:

```r
unlink("~/.cache/Azimuth", recursive = TRUE)
```

## Pseudobulk DE

### `Error in DESeqDataSetFromMatrix(...)`

Common causes:

- Mean aggregation instead of sum — use `aggregate_to_pseudobulk()`.
- Normalized counts passed in — use raw counts from `RNA$counts`.
- Sample metadata has rows not matching counts columns — check `all(colnames(counts) == rownames(meta))`.

### Refuses to run with N=1

Correct behavior. Use `FindMarkers` as exploratory only and label it. Report N=1 as a design limitation.

### `Error: every gene contains at least one zero`

Filter low-count genes:

```r
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
```

## Export

### `saveRDS` very slow

Drop intermediates first:

```r
seurat_obj[["SCT"]] <- NULL    # if not needed downstream
seurat_obj[["integrated"]] <- NULL
saveRDS(seurat_obj, "results/seurat_processed.rds")
```

Or save raw and processed assays separately.

### `Error: object 'seurat_obj' not found` after kernel restart

```r
seurat_obj <- readRDS("results/seurat_processed.rds")
```

The export is idempotent — safe to re-run upstream steps if needed.

### Converting Seurat → AnnData

```r
library(SeuratDisk)
SaveH5Seurat(seurat_obj, "seurat_obj.h5Seurat")
Convert("seurat_obj.h5Seurat", dest = "h5ad")
```

Or via `zellkonverter`:

```r
library(zellkonverter)
sce <- as.SingleCellExperiment(seurat_obj)
writeH5AD(sce, "seurat_obj.h5ad")
```

## Environment

### Conflicting package versions

Fresh env is the most reliable fix. With `renv`:

```r
renv::init()
renv::install("Seurat")
```

Or a conda env:

```bash
conda create -n seurat-env r-base=4.3
conda activate seurat-env
R -e 'install.packages(c("Seurat", "BiocManager"))'
```

### `future` package limits

For SCTransform and DESeq2 on large data:

```r
options(future.globals.maxSize = 8 * 1024^3)
library(future)
plan("multisession", workers = 4)
```

## When to ask for help

If stuck:

1. Check the Seurat docs and vignettes: <https://satijalab.org/seurat/articles/>
2. Search GitHub issues: <https://github.com/satijalab/seurat/issues>
3. Ask on Seurat Discourse: <https://github.com/satijalab/seurat/discussions>
4. Bioconductor support forum.

Include:

- `packageVersion("Seurat")`, `R.version.string`, OS.
- Full error traceback.
- A minimal reproducible example.

## Citations

The fixes here are drawn from the Seurat documentation, the satijalab GitHub issues, and the Heumos et al. *Best Practices* review (Nat Rev Genet 2023).
