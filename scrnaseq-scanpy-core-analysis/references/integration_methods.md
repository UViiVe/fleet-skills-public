# Batch Integration

Batch effects are systematic technical differences between samples — different chip runs, different library-prep dates, different reagents, different labs. Left uncorrected, cells cluster by batch instead of by biology, cell types fragment across batches, and downstream DE becomes invalid.

This guide covers when to integrate, how to choose a method, and how to assess whether the integration actually worked.

## When to integrate

Always integrate when:

- You have multiple 10X runs or sequencing batches.
- You have multi-donor or multi-condition data.
- The initial PCA-only UMAP shows clear batch-driven clusters.

Skip integration when:

- You have a single batch and a single donor.
- Batch is perfectly confounded with biology (e.g., one batch = one tissue) — integration will erase the biology you care about. Re-design the analysis.
- The initial UMAP shows batches already overlapping cleanly.

Always look at the un-integrated UMAP first. If the structure looks like biology and batch is a minor effect within clusters, integration may over-correct and hurt you.

## Method choice

```
Have cell-type labels (even partial)?
  yes → scANVI (best), or Harmony (faster)
  no  → scVI (best), or Harmony (faster)

Speed vs accuracy?
  speed   → Harmony
  accuracy → scVI / scANVI

> 100k cells?
  Harmony (or sketch-based methods)

Stuck in the R / Seurat ecosystem?
  Seurat CCA / RPCA, or Harmony via SeuratWrappers
```

### Comparison

| Method      | Pros                                                                       | Cons                                                            | Best for                                                  |
| ----------- | -------------------------------------------------------------------------- | --------------------------------------------------------------- | --------------------------------------------------------- |
| **scVI**    | State-of-the-art; gene-specific batch effects; uncertainty estimates       | Needs GPU for practical speed; 10–30 min; black-box             | Most use cases, especially complex batch structures       |
| **scANVI**  | Semi-supervised; leverages labels; preserves rare cell types               | Needs labels; same cost as scVI                                 | When you already have a coarse annotation                 |
| **Harmony** | Very fast (minutes); interpretable; hard to over-correct; no GPU           | Less powerful than scVI; can mis-mix when batch ≠ technical     | Large datasets, quick iteration, simple batch structure   |
| **Seurat CCA** | Well-established in Seurat; works on moderate batch effects             | Slow on large data; anchor finding is finicky                   | R users with moderate batch effects                       |
| **Seurat RPCA** | Faster than CCA; scales well                                          | Requires similar cell-type compositions across batches          | Large Seurat datasets with shared cell types              |

## scVI

scVI fits a deep generative model on raw counts with batch as a categorical covariate. The output is a low-dimensional latent representation (`X_scVI`) that you feed into the neighbor graph in place of PCA.

```python
import scvi

scvi.model.SCVI.setup_anndata(adata, batch_key="batch", layer="counts")
model = scvi.model.SCVI(
    adata,
    n_latent=30,
    n_layers=2,
    n_hidden=128,
    dropout_rate=0.1,
)
model.train(max_epochs=400, early_stopping=True, accelerator="gpu")
adata.obsm["X_scVI"] = model.get_latent_representation()
```

Parameter notes:

- `n_latent`: 30 for complex data; 10 (default) for very simple cases. Don't exceed 50 — you'll overfit.
- `max_epochs`: 400–600 with early stopping. The model saves itself when validation loss stops improving.
- `gene_likelihood`: `"nb"` (negative binomial) for UMI data — the default and what you want. `"zinb"` for Smart-seq2 / full-length.
- `n_hidden`: drop from 128 → 64 if GPU memory is tight.

Always train on raw counts. Don't feed normalized data — scVI assumes the count generative model.

## scANVI

scANVI extends scVI by adding a classifier head on cell-type labels. It's semi-supervised: you can have labels for some cells and `"Unknown"` for the rest, and it will propagate.

The recommended workflow is two-stage: train scVI first to get good unsupervised embeddings, then fine-tune as scANVI:

```python
scvi.model.SCANVI.setup_anndata(
    adata,
    batch_key="batch",
    labels_key="cell_type",
    unlabeled_category="Unknown",
)
scanvi = scvi.model.SCANVI.from_scvi_model(
    scvi_model,
    unlabeled_category="Unknown",
    labels_key="cell_type",
)
scanvi.train(max_epochs=200, accelerator="gpu")
adata.obsm["X_scANVI"] = scanvi.get_latent_representation()
adata.obs["scanvi_pred"] = scanvi.predict()
```

This consistently improves rare cell-type preservation over unsupervised scVI when even coarse labels are available.

## Harmony

Harmony performs iterative soft k-means clustering and linear correction in the PCA space. Fast (seconds to minutes), interpretable, robust to over-correction, no GPU required.

```python
import harmonypy as hm
import scanpy as sc

sc.tl.pca(adata, n_comps=50)
ho = hm.run_harmony(
    adata.obsm["X_pca"],
    adata.obs,
    "batch",
    theta=2,
    max_iter_harmony=10,
)
adata.obsm["X_pca_harmony"] = ho.Z_corr.T
```

Parameter notes:

- `theta`: the diversity penalty. 0 = no correction, 1 = gentle, 2 = standard, 4 = aggressive. Start at 2; raise if under-corrected.
- `lambda`: ridge regression penalty, default 1, leave it.
- `max_iter_harmony`: 10 is usually enough.

R equivalent in Seurat:

```r
seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "batch",
  dims.use = 1:30,
  theta = 2,
)
```

## Seurat CCA / RPCA

For R users staying in the Seurat ecosystem. CCA finds shared correlation structure across batches; RPCA uses reciprocal PCA projections and is faster but assumes similar cell-type compositions.

CCA:

```r
seurat_list <- SplitObject(seurat_obj, split.by = "batch")
seurat_list <- lapply(seurat_list, function(x) {
  x |> NormalizeData() |> FindVariableFeatures(nfeatures = 2000)
})
features <- SelectIntegrationFeatures(seurat_list)
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  anchor.features = features,
  reduction = "cca",
  dims = 1:30,
)
seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
```

RPCA: identical pipeline with `reduction = "rpca"`.

Use RPCA when the dataset exceeds ~100k cells or when batches share most cell types. Use CCA when batches differ in composition or you prioritize accuracy.

## Integration quality

Visual inspection first: UMAP colored by batch should overlap; colored by cell type should still separate cleanly. Marker genes should remain specific to their canonical cell types.

### LISI

Local Inverse Simpson's Index measures local diversity around each cell.

- **iLISI (batch)**: 1 = no mixing, `n_batches` = perfect mixing. Higher is better. Target: close to `n_batches`.
- **cLISI (cell type)**: 1 = perfect separation, `n_celltypes` = total chaos. Lower is better. Target: close to 1.

```python
from scripts.integration_diagnostics import compute_lisi_scores

lisi = compute_lisi_scores(adata, batch_key="batch", label_key="cell_type", use_rep="X_scVI")
print(f"iLISI: {lisi['batch_LISI']:.2f}")
print(f"cLISI: {lisi['celltype_LISI']:.2f}")
```

### ASW

Average Silhouette Width. Treats batch labels like cluster assignments:

- **Batch ASW**: should be close to 0 (batches well-mixed within cell-type clusters).
- **Cell type ASW**: should be high (cell types well-separated). > 0.5 is good.

### What "good" looks like

| Metric             | Target           |
| ------------------ | ---------------- |
| Batch iLISI        | ≈ `n_batches`    |
| Cell-type cLISI    | ≈ 1              |
| Batch ASW          | ≈ 0              |
| Cell-type ASW      | > 0.5            |

These numbers are not absolute. The trade-off between batch mixing and biology preservation is real; integration always loses some signal. The right question is: *do my downstream conclusions change?* If marker genes look right and cell types are recoverable, you're done.

## Choosing the method

| Dataset size       | Recommendation                                |
| ------------------ | --------------------------------------------- |
| < 10k cells        | scVI (overkill but free since it's small)     |
| 10–50k cells       | scVI or Harmony                               |
| 50–200k cells      | Harmony or scANVI                             |
| > 200k cells       | Harmony or Seurat RPCA (sketch)               |

Computational resources:

- GPU available, Python: scVI/scANVI.
- CPU only: Harmony.
- Limited time: Harmony.
- Need uncertainty: scVI (model.get_normalized_expression returns posteriors).

## Common failure modes

### Over-correction

Symptoms: cell types collapse together, markers lose specificity, cLISI rises.

Fixes:

- scVI: reduce `n_latent` to 20 or `max_epochs` to 200.
- Harmony: lower `theta` to 1.
- Seurat: reduce dimensions used in integration.

### Under-correction

Symptoms: UMAP still shows batch-driven clusters, iLISI close to 1.

Fixes:

- scVI: bump `n_latent` to 40, train longer.
- Harmony: raise `theta` to 3 or 4.
- Seurat: increase `k.anchor` and anchor features.

### Batch–condition confounding (N=1 in a condition)

If a biological condition has only one sample, scVI cannot perfectly separate batch noise from biology for that sample. iLISI < `n_batches` is *expected* — perfect mixing would erase the real biological difference.

What stays valid:

- Clustering and cell-type annotation.
- Marker gene identification.
- Within-condition analyses (e.g., T cell substates within COVID samples).

What needs caution:

- Quantitative composition comparisons across conditions (N=1 vs N=5 — you can't compute a meaningful contrast).
- Pseudobulk DE — DESeq2 needs N ≥ 2 per group and will refuse to run.

`run_scvi_integration` accepts a `condition_key` argument and auto-detects N=1 conditions. The warning is informational, not blocking.

### GPU out of memory

- Reduce `n_hidden` 128 → 64.
- Reduce `train_size` to 0.9.
- Subsample temporarily to confirm it's a memory issue, not a parameter issue.
- Drop to CPU: `accelerator="cpu"` (10× slower).

### Rare cell types split by batch

- Switch to scANVI with even coarse labels.
- Increase `n_latent` to 40–50.
- Don't downsample — that makes it worse.

### Cell types only in one batch

Genuinely batch-confined populations (a tumor-specific cell type in tumor samples only) are not a batch effect. Don't integrate them away — either exclude them from integration, or accept the warning that they don't mix.

## Best practices

1. **Always check the un-integrated UMAP first.** If batches mix already, integration is a no-op at best and harmful at worst.
2. **Compare methods.** Run scVI and Harmony, look at LISI/ASW for both, pick the one that preserves your biology.
3. **Validate with markers.** Cell-type markers should remain specific after integration.
4. **Use integrated coordinates for clustering / UMAP, raw counts for DE.** Never run DESeq2 on scVI-corrected counts.
5. **Document parameters.** Save them into `adata.uns["integration"]` so future-you can reproduce.

## Citations

1. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. *Deep generative modeling for single-cell transcriptomics.* Nat Methods 15:1053–1058 (2018).
2. Xu C, Lopez R, Mehlman E, et al. *Probabilistic harmonization and annotation of single-cell transcriptomics data with deep generative models.* Mol Syst Biol 17:e9620 (2021).
3. Korsunsky I, Millard N, Fan J, et al. *Fast, sensitive and accurate integration of single-cell data with Harmony.* Nat Methods 16:1289–1296 (2019).
4. Hao Y, Hao S, Andersen-Nissen E, et al. *Integrated analysis of multimodal single-cell data.* Cell 184:3573–3587.e29 (2021).
5. Luecken MD, Büttner M, Chaichoompu K, et al. *Benchmarking atlas-level data integration in single-cell genomics.* Nat Methods 19:41–50 (2022).
6. Tran HTN, Ang KS, Chevrier M, et al. *A benchmark of batch-effect correction methods for single-cell RNA sequencing data.* Genome Biol 21:12 (2020).
