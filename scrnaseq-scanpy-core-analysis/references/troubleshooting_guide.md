# Scanpy Troubleshooting

Symptoms, causes, fixes. Grouped by phase.

## Loading

### `ModuleNotFoundError: No module named 'scanpy'`

```bash
pip install "scanpy[leiden]" anndata scrublet
```

### `FileNotFoundError: barcodes.tsv.gz not found`

CellRanger output directories must contain three files (any of these naming conventions):

- `barcodes.tsv.gz` (or `.tsv`)
- `features.tsv.gz` (or `genes.tsv.gz` for older CellRanger)
- `matrix.mtx.gz` (or `.mtx`)

Confirm:

```python
import os
print(os.listdir("path/to/filtered_feature_bc_matrix/"))
```

If you have a `.h5` file, use `import_h5_data()` instead.

### Out-of-memory loading a large dataset

```python
adata = sc.read_h5ad("data.h5ad", backed="r")
```

Or subset for prototyping:

```python
adata = sc.read_10x_mtx("path/", make_unique=True)
adata = adata[:10000, :].copy()
```

## QC

### Too few cells retained (< 50%)

Diagnose with violins first:

```python
sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
             jitter=0.4, multi_panel=True)
```

If the distributions look fine and you're losing cells anyway, your thresholds are too strict. Switch to MAD:

```python
median = adata.obs["n_genes_by_counts"].median()
mad = np.median(np.abs(adata.obs["n_genes_by_counts"] - median))
upper = median + 3 * mad
lower = median - 3 * mad
```

Or use tissue-appropriate thresholds — see `qc_guidelines.md`.

### All cells have high mitochondrial percentage

Check the distribution shape:

```python
import matplotlib.pyplot as plt
plt.hist(adata.obs["pct_counts_mt"], bins=100)
plt.xlabel("% mitochondrial")
plt.show()
```

- Bimodal (low peak < 10%, high peak > 20%): mix of good and bad cells. Filter at the trough.
- Unimodal at 10–15%: could be tissue-specific (cardiomyocytes, hepatocytes). Don't filter blindly.
- Unimodal at > 20% across the board: sample is degraded. Consider rerunning the wet lab.

### `Mitochondrial genes identified: 0`

Wrong species pattern or unusual gene symbols:

```python
print(adata.var_names[adata.var_names.str.contains("MT", case=False)][:10])

mt_human = adata.var_names.str.startswith("MT-")
mt_mouse = adata.var_names.str.startswith("mt-")
print(f"Human: {mt_human.sum()}, mouse: {mt_mouse.sum()}")
```

Pass the correct `species` arg to `calculate_qc_metrics`.

## Normalization

### `ValueError: inplace=True requires raw to be None`

`adata.raw` is already populated.

```python
adata.raw = None
sc.pp.normalize_total(adata, inplace=True)
# or
adata = sc.pp.normalize_total(adata, inplace=False)
```

### Negative values after normalization

You ran regression before log transform. Correct order:

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])
sc.pp.scale(adata)
```

### Memory error during `sc.pp.scale`

Scaling densifies. Subset to HVGs first:

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata)
```

## PCA and embedding

### PC1–30 explain only 20% of variance

Causes: too many lowly expressed genes, batch effects dominating, not using HVGs.

```python
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.tl.pca(adata, use_highly_variable=True)
```

Check for batch effects: `sc.pl.pca(adata, color="batch")`. If batch dominates PC1–2, integrate before clustering.

### UMAP is one giant blob

- Not enough PCs: `sc.pp.neighbors(adata, n_pcs=50)`.
- Need more neighbors: `n_neighbors=15`.
- Resolution too low: `sc.tl.leiden(adata, resolution=1.0)`.

### UMAP changes every run

UMAP is stochastic:

```python
sc.tl.umap(adata, random_state=42)
```

## Clustering

### 50+ tiny clusters

Resolution too high:

```python
sc.tl.leiden(adata, resolution=0.4)
```

Or fewer neighbors:

```python
sc.pp.neighbors(adata, n_neighbors=5)
```

### Clusters driven by QC metrics

Symptoms: cluster identity correlates with `pct_counts_mt` or `total_counts`, no biological markers.

Fixes:

- Tighten QC.
- Regress during scaling: `sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])`.
- Remove the problem cluster post-hoc.

### Known cell types fail to separate

- Resolution sweep — go higher.
- Use more PCs.
- Check for batch effects (color UMAP by batch).
- Try integration (scVI / Harmony).

## Marker identification

### No significant markers

- Clusters too similar — check dendrogram: `sc.tl.dendrogram(adata, groupby="leiden")`.
- Thresholds too strict — relax `min_in_group_fraction` to 0.10, `min_fold_change` to 0.5.
- Try `t-test_overestim_var` instead of `wilcoxon`.

### Top markers are all RPS/RPL/MT genes

Filter them out before ranking:

```python
adata_for_markers = adata[:, ~adata.var_names.str.startswith(("MT-", "RPS", "RPL"))].copy()
sc.tl.rank_genes_groups(adata_for_markers, groupby="leiden_0.8")
```

## Plotting

### `ModuleNotFoundError: No module named 'seaborn'`

```bash
pip install seaborn matplotlib
```

If still failing, fall back to scanpy's built-in plots — they're ugly but always work.

### Blurry figures

```python
import matplotlib.pyplot as plt
plt.rcParams["figure.dpi"] = 300
plt.savefig("figure.svg", dpi=300, bbox_inches="tight")
```

### Can't see cluster labels on UMAP

```python
sc.pl.umap(adata, color="leiden",
           legend_loc="on data",
           legend_fontsize="x-small",
           legend_fontoutline=2)
```

## Memory and performance

### `MemoryError` or kernel crash

- Confirm sparsity: `from scipy.sparse import issparse; print(issparse(adata.X))`.
- Subset to HVGs before scaling/PCA.
- Use backed mode for huge files: `sc.read_h5ad("data.h5ad", backed="r+")`.
- Increase swap, or move to a bigger machine.

### Analysis is slow

```python
import os, scanpy as sc
sc.settings.n_jobs = max(1, os.cpu_count() - 2)
sc.tl.rank_genes_groups(adata, groupby="leiden", n_jobs=sc.settings.n_jobs)
```

## Integration

### Batches still cluster after "integration"

- Try a different method (scVI ↔ Harmony).
- For scVI: bump `n_latent` to 40, `max_epochs` to 600.
- For Harmony: raise `theta` to 3 or 4.

### Lost biological variation

Over-corrected. Lower `n_latent` (scVI) or `theta` (Harmony). Check that batch and biology aren't confounded — if they are, integration can't save you.

## Export

### H5AD write hangs

Large object with default gzip is slow. The script defaults to `lzf` which is ~5–10× faster:

```python
save_h5ad(adata, "output.h5ad")              # compression="lzf"
save_h5ad(adata, "output.h5ad", compression=None)   # fastest, biggest
```

Expected write times (very approximate):

| Cells     | lzf    | gzip      | None    |
| --------- | ------ | --------- | ------- |
| 3 000     | 2–5 s  | 10–30 s   | 1–2 s   |
| 10 000    | 5–15 s | 30–90 s   | 3–8 s   |
| 50 000    | 15–45 s| 2–5 min   | 10–30 s |
| 100 000+  | 30–90 s| 5–15 min  | 20–60 s |

If you interrupt, the export is idempotent — just re-run that step.

For very large datasets, drop the counts layer before saving if you don't need it downstream:

```python
if "counts" in adata.layers:
    del adata.layers["counts"]
save_h5ad(adata, "output.h5ad")
```

### CSV files too large

Use H5AD instead, or save sparse:

```python
import scipy.io
scipy.io.mmwrite("normalized_counts.mtx", adata.X)
```

### Can't open H5AD in R

```r
library(Seurat)
library(SeuratDisk)
Convert("data.h5ad", dest = "h5seurat")
obj <- LoadH5Seurat("data.h5seurat")
```

Or use `zellkonverter` in Bioconductor for direct H5AD ↔ SCE conversion.

## Environment

### Conflicting package versions

Fresh env is the most reliable fix:

```bash
conda create -n scanpy-env python=3.11
conda activate scanpy-env
pip install "scanpy[leiden]" python-igraph leidenalg scvi-tools harmonypy celltypist
```

### `ImportError: cannot import name 'UMAP' from 'umap'`

```bash
pip uninstall umap umap-learn
pip install umap-learn
```

## When to ask for help

If you're stuck:

1. Check the scanpy docs: <https://scanpy.readthedocs.io/>
2. Search GitHub issues: <https://github.com/scverse/scanpy/issues>
3. Ask on scverse Discourse: <https://discourse.scverse.org/>
4. Stack Overflow with `scanpy` and `single-cell` tags.

When asking, include:

- `sc.__version__`, Python version, OS.
- Full error traceback.
- A minimal reproducible example (10 cells, 100 genes, simulated counts).

## Citations

The fixes here are drawn from the scanpy documentation, the scverse Discourse, and the Heumos et al. *Best Practices* review (Nat Rev Genet 2023).
