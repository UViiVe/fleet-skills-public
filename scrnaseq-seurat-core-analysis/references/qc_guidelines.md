# QC Thresholds and Filtering Strategy (Seurat)

The three primary QC metrics in Seurat:

- `nFeature_RNA` — genes detected per cell.
- `nCount_RNA` — total UMIs per cell.
- `percent.mt` — fraction of UMIs mapping to mitochondrial genes (added by `PercentageFeatureSet` or `calculate_qc_metrics()`).

The defaults below are starting points — visualize the distributions before picking numbers.

## Tissue-specific defaults

| Tissue           | `min_features` | `max_features` | `max_mt_percent` | Notes                                                          |
| ---------------- | -------------- | -------------- | ---------------- | -------------------------------------------------------------- |
| PBMC             | 200            | 2500           | 5                | Standard immune cells; baseline.                               |
| Brain (cortex)   | 200            | 6000           | 10               | Neurons legitimately express many genes.                       |
| Tumor            | 200            | 5000           | 20               | Necrosis + dissociation stress push MT% up.                    |
| Kidney           | 200            | 4000           | 15               | Proximal tubule cells are large; podocytes are small.          |
| Liver            | 200            | 4000           | 15               | Hepatocytes carry many mitochondria naturally.                 |
| Heart            | 200            | 4000           | 15–20            | Cardiomyocytes ~ 15–25% MT — bimodal expected.                  |
| Skeletal muscle  | 200            | 4000           | 15               | Similar to heart.                                              |
| Lung             | 200            | 4000           | 10               | Type II pneumocytes shed surfactant into the soup.             |
| Pancreas         | 200            | 5000           | 10               | Acinar cells large and gene-rich.                              |

`max_mt_percent = 5` on heart data discards most cardiomyocytes — bias the threshold to the tissue.

## Platform notes

- **10X Chromium v3:** typical 1500–3000 features per cell. v2 is half.
- **Drop-seq:** 500–2000 features per cell. Less efficient; tighter empty-droplet filtering.
- **Smart-seq2:** 5000–10 000 features per cell. Higher `min_features` (500–1000). Doublets are rare with plate-based capture.

## MAD vs fixed thresholds

**Fixed thresholds** are one number per metric — interpretable, good for single-batch data with strong tissue priors. The skill's `filter_cells_by_qc()` and `filter_cells_by_tissue()` use this approach.

**Median Absolute Deviation (MAD) thresholds** adapt to each batch's own distribution and protect minority cell types in heterogeneous tissues:

```r
metric_values <- seurat_obj$nFeature_RNA
median_val <- median(metric_values)
mad_val <- median(abs(metric_values - median_val))
lower <- median_val - nmads * mad_val
upper <- median_val + nmads * mad_val
```

`batch_mad_outlier_detection()` applies this per batch on log-transformed metrics. `nmads = 5` keeps more cells; `nmads = 3` is aggressive.

Always log-transform counts and feature counts before MAD — both metrics are heavy-tailed.

## Doublet detection

Always run for droplet-based methods. `run_doubletfinder()` auto-scales the expected doublet rate to cell loading (~ 0.8% per 1000 cells), matching the 10X published rates. Override `expected_doublet_rate` if you loaded much more or less than standard.

After running, inspect `seurat_obj$doublet_score` — you should see a bimodal histogram with singlets near 0.05–0.15 and doublets in a hump beyond 0.3. If the boundary looks wrong, override the threshold manually.

DoubletFinder requires preliminary normalization and clustering — the script handles this internally. If it fails with cryptic errors, run `setup_for_integration()` first to ensure data is in the expected shape.

## Recommended workflow

1. Visualize first:
   ```r
   VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
   FeatureScatter(seurat_obj, "nCount_RNA", "nFeature_RNA")
   FeatureScatter(seurat_obj, "nCount_RNA", "percent.mt")
   ```
   Look for bimodal MT distributions (good vs bad cells), upper-tail doublet candidates, batch-specific shifts.

2. Decide between MAD (multi-batch) and fixed (single batch with strong tissue prior).

3. Run DoubletFinder per batch.

4. Filter, then re-check distributions.

5. Confirm retention ≥ 70%. Under that, ask why.

## Cell type–specific filtering

For very heterogeneous tissues (heart, brain), permissive initial filtering plus post-hoc cluster-level filtering beats a single global threshold:

```r
cm_cells <- WhichCells(seurat_obj, idents = "Cardiomyocyte")
non_cm_cells <- setdiff(colnames(seurat_obj), cm_cells)
keep <- c(
  cm_cells[seurat_obj$percent.mt[cm_cells] < 20],
  non_cm_cells[seurat_obj$percent.mt[non_cm_cells] < 10]
)
seurat_obj <- subset(seurat_obj, cells = keep)
```

This requires preliminary clustering and is iterative.

## When filtering fails

- **< 50% retained:** thresholds too strict, or sample genuinely degraded. Compare the violins; if `percent.mt` is unimodal at 20% across the board, the sample is toast.
- **Clusters driven by QC metrics, not biology:** tighten MT, pass `vars.to.regress = c("nCount_RNA", "percent.mt")` to SCTransform, or drop the problem cluster post-hoc.
- **MT% high everywhere:** check the wet-lab dissociation protocol. Cold dissociation (subtilisin A on ice) gives much lower MT% than 37°C collagenase.

## Checklist

- [ ] Visualized QC distributions before picking thresholds.
- [ ] Used tissue-appropriate values rather than PBMC defaults blindly.
- [ ] Ran DoubletFinder per batch (droplet data only).
- [ ] Confirmed retention > 70%.
- [ ] Verified no entire population was deleted (compare cluster counts pre vs post a coarse pass).
- [ ] Inspected for batch-specific QC differences.
- [ ] Documented decisions in the report.

## Citations

1. Hao Y, Stuart T, Kowalski MH, et al. *Dictionary learning for integrative, multimodal and scalable single-cell analysis.* Nat Biotechnol 42:293–304 (2024).
2. Heumos L, Schaar AC, Lance C, et al. *Best practices for single-cell analysis across modalities.* Nat Rev Genet 24:550–572 (2023).
3. McGinnis CS, Murrow LM, Gartner ZJ. *DoubletFinder: doublet detection in single-cell RNA sequencing data using artificial nearest neighbors.* Cell Syst 8:329–337 (2019).
4. Osorio D, Cai JJ. *Systematic determination of the mitochondrial proportion in human and mouse tissues for single-cell RNA-sequencing data quality control.* Bioinformatics 37:963–967 (2021).
