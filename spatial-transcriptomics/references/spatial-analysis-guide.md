# Spatial Transcriptomics: Parameter Reference

In-depth notes on the parameters that govern 10x Visium analysis with Squidpy + Scanpy. The default workflow ships with sensible defaults; this guide explains when to deviate and why.

## Spatial Neighbors Graph

The neighbors graph encodes which spots are considered adjacent in tissue space. Every downstream spatial statistic (Moran's I, neighborhood enrichment, co-occurrence) is computed against this graph, so the choice of `coord_type`, `n_neighs`, and `n_rings` materially affects the results.

### `coord_type`

| Value | Use for | Notes |
|---|---|---|
| `grid` | Visium (hexagonal array) | Uses `array_row` / `array_col` indices; `n_neighs=6` matches the hex geometry. |
| `generic` | Slide-seq, MERFISH, custom | Operates on `.obsm['spatial']` directly; requires `n_neighs` **or** `radius`. |

Visium default:

```python
sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
```

### `n_neighs` and `n_rings`

- `n_neighs` (default 6 for `grid`): number of neighbors per spot. Raise to 12 or 18 to include second-ring neighbors. Useful for noisy data where one-ring smoothing is too local.
- `n_rings` (default 1): number of hexagonal rings around each spot. `n_rings=2` captures broader spatial context at ~2× the runtime.

For non-Visium layouts:

```python
sq.gr.spatial_neighbors(adata, coord_type="generic", radius=150)        # fixed distance
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=10)       # KNN
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)     # Delaunay triangulation
```

## Spatially Variable Genes (Moran's I)

Moran's I quantifies spatial autocorrelation: whether expression is clustered, dispersed, or random across tissue.

| Moran's I | Interpretation |
|---|---|
| Near +1 | Strongly clustered — high near high, low near low. The interesting signal. |
| Near 0 | No spatial structure. |
| Near −1 | Dispersed / checkerboard. Rare in real data; usually a sign of preprocessing artifacts. |

### Significance

- `n_perms` (default 100): permutations for the null distribution. Higher → tighter p-values, longer runtime.
  - Exploration: 100 (about 1–2 min on the demo).
  - Publication: 1000 (about 10–15 min).
- FDR: Benjamini–Hochberg via `pval_norm_fdr_bh`. Standard cutoff 0.05; tighten to 0.01 for stringent claims.
- Alternative statistic: Geary's C (`mode="geary"`) is more sensitive to local differences. Useful as a sanity check; the two often correlate strongly.

### Output columns

`adata.uns["moranI"]` is a DataFrame with:

| Column | Meaning |
|---|---|
| `I` | Moran's I statistic |
| `pval_norm` | p-value from the normal approximation |
| `pval_norm_fdr_bh` | FDR-corrected p-value |
| `var_norm` | variance under the normality assumption |

Top SVGs typically include tissue-specific markers, ECM genes at tissue boundaries, and immune markers in infiltrated regions.

## Neighborhood Enrichment

Quantifies whether two cluster types are spatially co-localized (positive z) or segregated (negative z) more than chance.

| z-score | Reading |
|---|---|
| > 2 | Significantly co-localized — clusters tend to be adjacent. |
| −2 to 2 | No significant preference. |
| < −2 | Significantly segregated — clusters avoid each other. |

Diagonal entries indicate the spatial cohesion of each cluster individually (how often a cluster's spots have same-cluster neighbors).

## Co-occurrence

Computes co-occurrence probability of cluster pairs over a range of distances. Unlike neighborhood enrichment (which is binary adjacent/not), co-occurrence captures distance-dependent relationships.

Tips and known pitfalls:
- Use the default `n_splits` — overriding it can produce NaN in some squidpy versions.
- Very small clusters (< 10 spots) give noisy estimates; consider merging or excluding them.
- Distance bins are tissue-coordinate units; the absolute scale depends on the Visium scale factors stored in `.uns["spatial"]`.

## QC Thresholds

| Parameter | Typical range | Heart tissue | Notes |
|---|---|---|---|
| `min_genes` | 200–500 | 200 | Lower for sparse FFPE; raise for high-quality fresh tissue. |
| `min_cells` | 3–20 | 10 | Filters genes detected in too few spots. |
| `max_pct_mito` | 10–30% | 50% | Mito-rich tissues (heart, muscle, brown adipose) need a permissive cap. |

Cardiac tissue specifically: cardiomyocytes are mitochondria-dense, so the median `pct_counts_mt` sits at ~30–40%. A 50% cap retains > 95% of spots. For brain, liver, immune, or epithelial tissues, drop to 20%.

## Clustering

| Parameter | Range | Default | Effect |
|---|---|---|---|
| `resolution` | 0.2–2.0 | 0.8 | Higher → more clusters. |
| `n_neighbors` | 5–30 | 15 | Higher → smoother UMAP and clusters. |
| `n_pcs` | 10–50 | 30 | Higher → captures more variance. |

For spatial data, lower resolutions (0.5–0.8) tend to produce clusters that map cleanly onto visible anatomical regions. Higher resolutions fragment domains in ways that may not correspond to biology. Always validate against the H&E image.

## Highly Variable Gene Selection

| Parameter | Range | Default | Notes |
|---|---|---|---|
| `n_top_genes` | 1000–3000 | 2000 | Standard for a Visium dataset (~33k genes). |
| `flavor` | `seurat_v3` | `seurat_v3` | Operates on raw counts; recommended for spatial data. |

## Advanced

### Geary's C

```python
sq.gr.spatial_autocorr(adata, mode="geary", n_perms=100)
# Results in adata.uns['gearyC']
```

Run both Moran's I and Geary's C and reconcile the gene lists — Geary's C catches local differences Moran's I can miss.

### Subsetting for performance

For very large sections (> 10k spots) consider analyzing a region of interest:

```python
mask = (
    (adata.obsm["spatial"][:, 0] > x_min) &
    (adata.obsm["spatial"][:, 0] < x_max)
)
adata_sub = adata[mask].copy()
```

This is preferable to downsampling when you have a defined anatomical region in mind.

### Custom spatial graphs

```python
sq.gr.spatial_neighbors(adata, coord_type="generic", radius=150)
sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=10)
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
```

Delaunay triangulation is robust to irregular layouts and is sometimes a better choice for non-Visium imaging-based platforms.

## Downstream Analyses to Layer On

1. **Cell-type deconvolution** — combine SVG output with a matched scRNA-seq reference using cell2location, RCTD, or STdeconvolve to produce per-spot cell-type probabilities.
2. **Spatial ligand–receptor analysis** — use CellChat v2's spatial mode or COMMOT on top of the cluster assignments to find spatially constrained signaling.
3. **Spatial gene–gene correlation** — compute spatially weighted correlation between top SVGs to surface co-regulated programs.
4. **Multi-sample comparison** — integrate sections (e.g. disease vs control) with Harmony / Scanorama / scVI before re-running differential SVG and differential abundance analysis.

## References

- Palla G, et al. **Squidpy: a scalable framework for spatial omics analysis.** *Nature Methods* (2022) 19, 171–178.
- Wolf FA, et al. **SCANPY: large-scale single-cell gene expression data analysis.** *Genome Biology* (2018) 19, 15.
- Moran PAP. **Notes on continuous stochastic phenomena.** *Biometrika* (1950) 37, 17–23.
- Geary RC. **The Contiguity Ratio and Statistical Mapping.** *The Incorporated Statistician* (1954) 5(3), 115–146.
- Anselin L. **Local indicators of spatial association — LISA.** *Geographical Analysis* (1995) 27(2), 93–115.
- Squidpy documentation: https://squidpy.readthedocs.io
- Scanpy documentation: https://scanpy.readthedocs.io
