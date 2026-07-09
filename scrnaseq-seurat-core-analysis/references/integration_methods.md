# Batch Integration (Seurat)

When to integrate, which method to pick, and how to know whether the integration actually worked.

## When to integrate

Always when:

- Multiple sequencing runs or library-prep dates.
- Multi-donor or multi-condition cohorts.
- Initial PCA/UMAP shows clear batch-driven clusters.

Skip when:

- Single batch, single donor.
- Batch is perfectly confounded with biology — integration will erase the signal.
- Un-integrated UMAP shows batches already overlapping cleanly.

Always look at the un-integrated UMAP first. If batch is a minor secondary effect within biological clusters, integration may over-correct.

## Method choice

```
Have cell-type labels (even partial)?
  yes → Harmony with prior labels, or Azimuth/transfer learning
  no  → Harmony (default), CCA, or RPCA

Speed vs accuracy?
  speed   → Harmony
  accuracy → Seurat CCA (slower, more thorough)

> 100k cells, similar compositions?
  RPCA (faster than CCA, scales better)

> 200k cells?
  Harmony or sketch-based methods
```

## Comparison

| Method      | Pros                                                                          | Cons                                                       | Best for                                                |
| ----------- | ----------------------------------------------------------------------------- | ---------------------------------------------------------- | ------------------------------------------------------- |
| **Harmony** | Fast (seconds to minutes); hard to over-correct; no GPU; widely used          | Less powerful than CCA on complex batches                  | Default for most multi-batch studies                    |
| **Seurat CCA** | Well-established; learns shared structure across batches                   | Slow; anchor finding finicky on highly heterogeneous data  | Moderate batch effects, < 100k cells                    |
| **Seurat RPCA** | Faster than CCA; scales to > 100k cells                                   | Requires similar cell-type compositions across batches     | Large datasets with shared cell types                   |
| **scVI**    | State-of-the-art uncertainty quantification (via reticulate)                  | GPU strongly recommended; Python dependency                | Python-savvy labs with GPU access                       |
| **fastMNN** | Fast, principled mutual nearest neighbors                                     | Less popular than Harmony; older                           | Legacy pipelines                                        |

## Harmony

The default choice. Iterative soft k-means + linear correction in PCA space.

```r
library(harmony)
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "batch",
  dims.use = 1:30,
  theta = 2,
  lambda = 1,
  max.iter.harmony = 10
)
```

`run_harmony_integration()` wraps this.

Parameters:

- `theta`: diversity penalty. 0 = no correction, 1 = gentle, 2 = standard (default), 4 = aggressive. Start at 2, raise if under-corrected.
- `lambda`: ridge regression penalty, default 1, leave it.
- `max.iter.harmony`: 10 is usually enough.

Downstream clustering and UMAP must use `reduction = "harmony"`:

```r
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:30)
```

## Seurat CCA

```r
seurat_list <- SplitObject(seurat_obj, split.by = "batch")
seurat_list <- lapply(seurat_list, function(x) {
  x |> NormalizeData() |> FindVariableFeatures(nfeatures = 2000)
})
features <- SelectIntegrationFeatures(seurat_list, nfeatures = 2000)
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  reduction = "cca",
  dims = 1:30
)
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(seurat_integrated) <- "integrated"
seurat_integrated <- ScaleData(seurat_integrated)
seurat_integrated <- RunPCA(seurat_integrated, npcs = 50)
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:30)
```

`run_seurat_cca_integration()` wraps this with sane defaults. Note that downstream analysis runs on the `integrated` assay, which has only the anchor-feature genes; switch back to `RNA` for gene-level visualization and DE.

## Seurat RPCA

Same pipeline as CCA but pass `reduction = "rpca"` to `FindIntegrationAnchors`. Each batch's data is projected onto the others' PCA space — faster, scales better, but assumes similar cell-type compositions across batches.

```r
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30
)
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
```

`run_seurat_rpca_integration()` wraps this.

## Seurat v5 streamlined integration

Seurat v5 introduced a unified `IntegrateLayers()` API that works on a single layered RNA assay (no `SplitObject`):

```r
seurat_obj[["RNA"]] <- split(seurat_obj[["RNA"]], f = seurat_obj$batch)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- IntegrateLayers(
  object = seurat_obj,
  method = HarmonyIntegration,   # or CCAIntegration, RPCAIntegration, scVIIntegration
  orig.reduction = "pca",
  new.reduction = "integrated.harmony"
)
seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.harmony", dims = 1:30)
```

The skill scripts default to the layered v5 API where available and fall back to the older split-list approach if not.

## Integration quality

Visual inspection first: UMAP coloured by batch should overlap, by cell type should separate, marker genes should remain specific.

### LISI

`compute_lisi_scores()` returns:

- **iLISI (batch)**: 1 = no mixing, `n_batches` = perfect mixing. Higher is better.
- **cLISI (cell type)**: 1 = perfect separation, `n_celltypes` = total chaos. Lower is better.

Target: iLISI close to `n_batches`, cLISI close to 1.

### ASW

`compute_asw_scores()` returns silhouette widths:

- **Batch ASW**: should be close to 0 (well-mixed).
- **Cell-type ASW**: should be high, > 0.5 ideal.

### Numbers to aim for

| Metric             | Target           |
| ------------------ | ---------------- |
| Batch iLISI        | ≈ `n_batches`    |
| Cell-type cLISI    | ≈ 1              |
| Batch ASW          | ≈ 0              |
| Cell-type ASW      | > 0.5            |

These aren't absolute. The trade-off between batch mixing and biology preservation is real. Ask: *do my downstream conclusions change?* If markers look right and cell types are recoverable, you're done.

## Choosing by data size

| Dataset size       | Recommendation                            |
| ------------------ | ----------------------------------------- |
| < 10k cells        | Harmony (overkill but free since small)   |
| 10–50k cells       | Harmony or CCA                            |
| 50–200k cells      | Harmony or RPCA                           |
| > 200k cells       | Harmony or RPCA with sketch               |

## Failure modes

### Over-correction

Symptoms: cell types collapse together, markers lose specificity, cLISI rises.

Fixes:

- Harmony: lower `theta` to 1.
- CCA/RPCA: reduce `dims` used in integration.

### Under-correction

Symptoms: UMAP still shows batch-driven clusters, iLISI close to 1.

Fixes:

- Harmony: raise `theta` to 3 or 4.
- CCA/RPCA: increase `anchor.features`, raise `k.anchor`.

### Batch–condition confounding

If one condition has only one sample, no method can perfectly separate batch noise from biology for that sample. Clustering and annotation remain valid; quantitative composition comparisons and pseudobulk DE do not. `run_harmony_integration()` warns when this is detected.

### Memory issues with CCA

CCA is memory-heavy. On > 50k cells, switch to RPCA or Harmony.

### Cell types only in one batch

A tumor-specific cell type in tumor samples only is not a batch effect. Don't integrate it away — either exclude it from integration or accept that it won't mix.

## Best practices

1. **Always check the un-integrated UMAP first.**
2. **Compare methods.** Run Harmony and CCA, compare LISI/ASW, pick the one that preserves your biology.
3. **Validate with markers.** Cell-type markers should remain specific after integration.
4. **Use integrated coordinates for clustering / UMAP, raw counts (RNA assay) for DE.** Never run DESeq2 on integrated data.
5. **Document parameters.** Save them into `seurat_obj@misc$integration`.

## Citations

1. Korsunsky I, Millard N, Fan J, et al. *Fast, sensitive and accurate integration of single-cell data with Harmony.* Nat Methods 16:1289–1296 (2019).
2. Hao Y, Hao S, Andersen-Nissen E, et al. *Integrated analysis of multimodal single-cell data.* Cell 184:3573–3587.e29 (2021).
3. Stuart T, Butler A, Hoffman P, et al. *Comprehensive integration of single-cell data.* Cell 177:1888–1902.e21 (2019).
4. Luecken MD, Büttner M, Chaichoompu K, et al. *Benchmarking atlas-level data integration in single-cell genomics.* Nat Methods 19:41–50 (2022).
5. Tran HTN, Ang KS, Chevrier M, et al. *A benchmark of batch-effect correction methods for single-cell RNA sequencing data.* Genome Biol 21:12 (2020).
