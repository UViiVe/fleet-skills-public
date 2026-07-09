# QC Thresholds and Filtering Strategy

The three QC metrics that drive cell filtering in scRNA-seq are:

- `n_genes_by_counts` — genes detected per cell. Low values are empty droplets and dying cells; very high values are usually doublets.
- `total_counts` — total UMIs per cell. Tracks library size and correlates with `n_genes_by_counts`.
- `pct_counts_mt` — fraction of UMIs mapping to mitochondrial genes. High values mark stressed or lysed cells, but several tissues legitimately run hot.

There is no single universal threshold for any of these. The defaults below are starting points; always inspect the distribution first.

## Tissue-specific defaults

| Tissue           | `min_genes` | `max_genes` | `max_pct_mt` | Notes                                                                |
| ---------------- | ----------- | ----------- | ------------ | -------------------------------------------------------------------- |
| PBMC             | 200         | 2500        | 5            | Small immune cells, low metabolic demand.                            |
| Brain (cortex)   | 200         | 6000        | 10           | Neurons legitimately express many genes.                             |
| Tumor            | 200         | 5000        | 20           | Dissociation stress + necrotic cells push MT% up.                    |
| Kidney           | 200         | 4000        | 15           | Proximal tubule cells are large; podocytes are small.                |
| Liver            | 200         | 4000        | 15           | Hepatocytes carry many mitochondria by nature.                       |
| Heart            | 200         | 4000        | 15–20        | Cardiomyocytes can have 15–25% MT — bimodal distribution expected.   |
| Skeletal muscle  | 200         | 4000        | 15           | Similar to heart, high energy demand.                                |
| Lung             | 200         | 4000        | 10           | Mixed cell sizes; type II pneumocytes shed surfactant into the soup. |
| Pancreas         | 200         | 5000        | 10           | Acinar cells are large and gene-rich.                                |

A `max_pct_mt = 5` filter on heart tissue will throw away most of your cardiomyocytes. Bias the threshold to the tissue.

## Platform considerations

- **10X Chromium v3:** ~1500–3000 genes/cell. v2 chemistry is half that.
- **Drop-seq:** ~500–2000 genes/cell. Lower capture efficiency — tighter empty-droplet filtering.
- **Smart-seq2:** full-length, ~5000–10 000 genes/cell. Use higher `min_genes` (500–1000); doublets are rare with plate-based capture.

## Adaptive vs fixed thresholds

**Fixed thresholds** (one number per metric) are interpretable and good for single-batch data with strong tissue priors. Use `mark_outliers_fixed(tissue=...)` for this.

**Median Absolute Deviation (MAD) thresholds** adapt to each batch's distribution and avoid penalizing metabolically active cell types in tissues where they're a minority. This is the default for multi-batch experiments:

```python
metric_values = adata.obs["log1p_total_counts"]
median = np.median(metric_values)
mad = np.median(np.abs(metric_values - median))
lower = median - nmads * mad
upper = median + nmads * mad
```

`nmads = 5` is conservative (keeps more cells). `nmads = 3` is aggressive (more like the IQR rule). The `batch_mad_outlier_detection` script applies the rule independently per batch.

Always compute MAD on log-transformed counts/genes — these metrics are heavy-tailed.

## Doublet detection

Always run for droplet methods. Plate-based methods (Smart-seq2) rarely produce doublets so it's optional there.

`run_scrublet_detection` auto-scales the expected doublet rate to the cell loading: ~0.8% per 1000 cells captured, which mirrors 10X's published rates. Override `expected_doublet_rate` if you loaded much more or less than standard.

Inspect `adata.obs["doublet_score"].hist()` after running. The bimodal shape — singlets clustered near 0.05–0.15, doublets in a separate hump near 0.3+ — is the diagnostic. If you see no second peak, the data is clean and the threshold is fine. If the boundary looks wrong, override `doublet_score_threshold` manually in `filter_by_mad_outliers`.

## Workflow

1. Visualize before filtering:
   ```python
   sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"])
   sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")
   sc.pl.scatter(adata, "total_counts", "pct_counts_mt")
   ```
2. Decide between MAD (multi-batch) and fixed thresholds (single batch with strong tissue prior).
3. Run Scrublet per batch.
4. Filter, then re-check the distribution.
5. Confirm cell retention > 70%. If you lost more, ask why — too-strict thresholds or genuinely poor sample.

## Cell type–specific filtering

For tissues with very heterogeneous cell sizes (heart: cardiomyocytes vs endothelial, brain: neurons vs glia), permissive initial filtering followed by post-hoc filtering by cluster works better than a single global threshold:

```python
cm_mask = adata.obs["cell_type"] == "Cardiomyocyte"
adata.obs.loc[cm_mask, "pass_qc"] = adata.obs.loc[cm_mask, "pct_counts_mt"] < 20
adata.obs.loc[~cm_mask, "pass_qc"] = adata.obs.loc[~cm_mask, "pct_counts_mt"] < 10
```

This requires preliminary clustering and is iterative.

## When filtering fails

- **< 50% retained:** thresholds too strict, or sample is genuinely degraded. Look at the QC distributions; if MT% has a single peak at 20% across all cells the sample is probably toast.
- **Clusters driven by QC metrics, not biology:** tighten MT filter, regress out `total_counts` and `pct_counts_mt` during scaling, or drop the cluster post-hoc.
- **Unexpectedly high MT% everywhere:** check the dissociation protocol with the wet-lab side. Cold dissociation (e.g., subtilisin A on ice) gives much lower MT% than 37°C collagenase.

## Checklist

- [ ] Visualized QC distributions before setting thresholds.
- [ ] Used tissue-appropriate values, not PBMC defaults blindly.
- [ ] Ran Scrublet (per batch) for droplet data.
- [ ] Confirmed retention > 70%.
- [ ] Verified no entire population was deleted (compare cluster counts before vs after filtering on a coarse pass).
- [ ] Inspected for batch-specific QC differences.
- [ ] Documented all decisions and rationale in the report.

## Citations

1. Heumos L, Schaar AC, Lance C, et al. *Best practices for single-cell analysis across modalities.* Nat Rev Genet. 24:550–572 (2023).
2. Osorio D, Cai JJ. *Systematic determination of the mitochondrial proportion in human and mouse tissues for single-cell RNA-sequencing data quality control.* Bioinformatics 37:963–967 (2021).
3. Germain P-L, Lun A, Garcia Meixide C, et al. *Doublet identification in single-cell sequencing data using scDblFinder.* F1000Research 10:979 (2021).
4. Wolock SL, Lopez R, Klein AM. *Scrublet: computational identification of cell doublets in single-cell transcriptomic data.* Cell Syst 8:281–291 (2019).
