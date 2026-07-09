# Scanpy Best Practices

Opinionated defaults for the core scanpy pipeline. These map directly to the script defaults in `scripts/`.

## Normalization

### Target-sum + log1p (default)

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
```

This is what `run_standard_normalization()` does. Works for almost everything droplet-based. Interpretable and stable.

Target sum of `1e4` is convention — it sets each cell's total to 10 000 (think "CPMs per cell"). You can pass `target_sum=None` to use the median library size, which avoids the dilution-of-deep-cells effect, but downstream tools mostly assume `1e4` so stay there unless you have a reason.

### Pearson residuals (alternative)

```python
sc.experimental.pp.normalize_pearson_residuals(adata)
```

Replaces the log1p approach with analytic Pearson residuals under a negative-binomial model. Better variance stabilization on highly heteroscedastic data (large dynamic range in mean expression). Output is already on a scaled scale, so you skip `sc.pp.scale`. Still flagged experimental in scanpy.

When to consider it: very deep sequencing (Smart-seq3 + UMI), or when standard normalization leaves PC1 dominated by library size.

## Highly variable genes

```python
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=2000,
    flavor="seurat_v3",
    layer="counts",
    batch_key="batch",   # optional
)
```

- `n_top_genes = 2000` is the sweet spot. Bump to 3000 for atlases with many cell types; drop to 1000 for very small datasets.
- `flavor="seurat_v3"` (operates on raw counts, robust to batch effects) is preferred when you have a `layer="counts"`. Otherwise `flavor="seurat"` works on log-normalized data.
- Pass `batch_key` to compute HVGs per batch and intersect/union — protects against batch-specific dropouts.

Visualize with `sc.pl.highly_variable_genes(adata)` before scaling. If the cloud has a flat shape with no clear "highly variable" tail, your data is too sparse and HVG selection is going to be noisy.

## Dimensionality reduction

### Number of PCs

Compute 50 PCs upfront. Use the elbow or `suggest_n_pcs()` (variance-ratio knee detection) to choose how many to feed into the neighbor graph.

**Hard rule: never use < 15 PCs.** Heuristics like "last PC where the variance ratio exceeds 5% of PC1" routinely undercount — they're calibrated for bulk RNA-seq, not single-cell.

Defaults that work:

- PBMC and most immune tissues: 20–30 PCs.
- Brain, kidney, lung, tumor (more cell types): 30–50 PCs.
- Unknown: 30 PCs is a safe default.

Cumulative variance explained should hit 80–90% by your chosen `n_pcs`. If it doesn't, your HVGs or QC is the problem, not the PC count.

### Regression during scaling

```python
sc.pp.regress_out(adata, keys=["total_counts", "pct_counts_mt"])
sc.pp.scale(adata, max_value=10)
```

Regress out technical variation (library size, MT%), not biological signal. Don't regress out `n_genes_by_counts` and `total_counts` together — they're collinear and the regression becomes unstable.

Some practitioners now skip regression entirely on the theory that scVI/Harmony will catch the structure downstream. That's reasonable for multi-batch integrated workflows. For single-batch PBMC, keep the regression.

## Clustering

### Leiden vs Louvain

Use Leiden. It guarantees well-connected communities (Louvain can produce disconnected ones), gives more reproducible results, and is comparably fast. Louvain is in scanpy only for backward compatibility.

```python
sc.tl.leiden(adata, resolution=0.8, random_state=42)
```

### Resolution

There is no "right" resolution. Run a sweep:

```python
for res in [0.4, 0.6, 0.8, 1.0]:
    sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{res}")
```

Pick by:

1. **Marker coherence.** Each cluster should have ≥ 3 strong, biologically sensible markers.
2. **Stability.** Adjacent resolutions should give consistent assignments for most cells.
3. **Biology.** If your downstream question is "how many subtypes of T cells?" you need fine resolution.

Resolution scale rules of thumb:

- 0.3–0.5: major lineages only (T, B, myeloid).
- 0.6–0.8: standard granularity (CD4 vs CD8, classical vs non-classical mono).
- 1.0–1.5: fine substates (naive vs memory subsets, transitional populations).
- > 2.0: usually overshoots and you start clustering noise.

### Neighbor graph

```python
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
```

- `n_neighbors=10` is the default and works for ~95% of datasets.
- Bigger datasets (> 50k cells): try `n_neighbors=15` or `20`.
- Small datasets (< 5k cells): try `n_neighbors=5` or `8`.

Use the same `n_neighbors` for the neighbor graph and UMAP — discrepancies cause visualization artifacts.

## UMAP

```python
sc.tl.umap(adata, min_dist=0.5, spread=1.0, random_state=42)
```

- `min_dist` controls tightness. 0.1 → tight clusters, gaps; 0.5 → standard; 1.0 → loose, broad spread.
- `spread=1.0` is fine for almost everything.
- Always set `random_state` — UMAP is non-deterministic.

UMAP is a visualization tool, not a quantitative coordinate system. Don't compute pseudo-distances on UMAP coordinates; use the neighbor graph or the integrated representation for that.

## Marker gene identification

```python
sc.tl.rank_genes_groups(adata, groupby="leiden_0.8", method="wilcoxon")
```

- Wilcoxon: the right default — non-parametric, robust, fast enough.
- t-test: fastest, assumes normality (wrong for count data); the `t-test_overestim_var` variant is what scanpy recommends if you want speed.
- Logistic regression: when you need to control for covariates (donor, batch). Slow on large datasets.

Filter the markers afterwards:

```python
sc.tl.filter_rank_genes_groups(
    adata,
    min_in_group_fraction=0.25,
    min_fold_change=1.0,
    max_out_group_fraction=0.5,
)
```

These are **exploratory** statistics — never use the resulting p-values as if they were inferential. Pseudoreplication inflates Type I error massively.

## Data management

### Sparsity

Keep `adata.X` sparse (`scipy.sparse.csr_matrix`) until you have to densify. Most scanpy functions handle sparse fine. The conversion happens automatically in some PCA/scaling paths; subset to HVGs first if memory matters.

### Layers convention

- `adata.X` = log-normalized data (what most plotting and DE functions assume).
- `adata.layers["counts"]` = raw integer UMIs (required for pseudobulk DE and scVI).
- `adata.raw` = optional snapshot of normalized data before HVG subsetting (so you can plot any gene later, not just HVGs).

### Memory-frugal patterns

- Subset to HVGs before scaling/regression: `adata = adata[:, adata.var.highly_variable].copy()`.
- Use backed mode for huge files you only need to slice: `sc.read_h5ad("data.h5ad", backed="r+")`.
- Delete heavy `uns` entries after you're done with them: `del adata.uns["rank_genes_groups"]`.

## Visualization style

Scanpy's built-in plots are fast but not publication-ready. For figures:

- Set `sc.settings.set_figure_params(dpi=300, facecolor="white", frameon=False)`.
- Use seaborn/matplotlib directly for custom panels:

```python
fig, ax = plt.subplots(figsize=(8, 6))
umap_df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
umap_df["cell_type"] = adata.obs["cell_type"].values
for ct, sub in umap_df.groupby("cell_type"):
    ax.scatter(sub["UMAP1"], sub["UMAP2"], s=0.5, alpha=0.7, label=ct)
ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, markerscale=5)
sns.despine(ax=ax)
```

Save both PNG (for slides) and SVG (for posters and final figures).

## Performance tips

- Set `sc.settings.n_jobs = max(1, os.cpu_count() - 2)` once at the top of your notebook.
- For very large datasets, the rapids-based backend (`method="rapids"` in `sc.pp.neighbors`) gives ~10× speedups on GPU.
- Save your `AnnData` periodically. A 30-minute scVI run that you lose to a kernel crash is expensive.

## Workflow checklist

- [ ] Loaded data, stored raw counts in `layers["counts"]`.
- [ ] Calculated QC metrics, visualized distributions.
- [ ] Filtered low-quality cells, kept ≥ 70%.
- [ ] Normalized (target-sum + log1p), log1p-transformed.
- [ ] Identified HVGs (Seurat v3 flavor preferred).
- [ ] Scaled, optionally regressed out technical variation.
- [ ] Ran PCA, chose `n_pcs` from elbow / `suggest_n_pcs()`.
- [ ] Built neighbor graph using PCA or integrated rep.
- [ ] Clustered at multiple resolutions, picked one.
- [ ] Ran UMAP with `random_state=42`.
- [ ] Found markers (Wilcoxon, exploratory).
- [ ] Annotated cell types (manual ± CellTypist).
- [ ] Exported `.h5ad` + CSV tables + figures.

## Resources

- Docs: <https://scanpy.readthedocs.io/>
- Single-Cell Best Practices book (Heumos et al.): <https://www.sc-best-practices.org/>
- Scanpy tutorials: <https://scanpy.readthedocs.io/en/stable/tutorials.html>
- scverse: <https://scverse.org/>

## References

1. Wolf FA, Angerer P, Theis FJ. *SCANPY: large-scale single-cell gene expression data analysis.* Genome Biol 19:15 (2018).
2. Heumos L, Schaar AC, Lance C, et al. *Best practices for single-cell analysis across modalities.* Nat Rev Genet 24:550–572 (2023).
3. Traag VA, Waltman L, van Eck NJ. *From Louvain to Leiden: guaranteeing well-connected communities.* Sci Rep 9:5233 (2019).
4. Lause J, Berens P, Kobak D. *Analytic Pearson residuals for normalization of single-cell RNA-seq UMI data.* Genome Biol 22:258 (2021).
