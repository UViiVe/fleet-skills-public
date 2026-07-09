# Decision Guide

Six recurring choices you'll make in any Seurat workflow, and how to make them.

## 1. Ambient RNA correction

| Skip                                          | Run                                               |
| --------------------------------------------- | ------------------------------------------------- |
| Filtered 10X output, low-soup tissue          | Raw 10X output                                    |
| Fresh PBMC                                    | High-soup tissue (brain, lung, tumor, lymph node) |
| Pre-corrected data                            | Long processing time before droplet generation    |
| Nuclei-based methods (snRNA-seq)              | Suspected ambient contamination in marker plots   |

Method:

- **SoupX** (R-native): fast, interpretable, no GPU. Default for Seurat pipelines.
- **CellBender** (Python): better on complex soup, needs GPU.

## 2. QC threshold strategy

| Fixed thresholds                              | MAD                                               |
| --------------------------------------------- | ------------------------------------------------- |
| Single batch with strong tissue prior         | Multi-batch experiment                            |
| You know exactly what to filter               | Cell types vary in QC distribution                |
| `filter_cells_by_qc()` / `filter_cells_by_tissue()` | `batch_mad_outlier_detection()`              |

Both methods support the same downstream pipeline. Default to MAD for safety.

## 3. Normalization

| SCTransform                                   | LogNormalize                                      |
| --------------------------------------------- | ------------------------------------------------- |
| UMI-based data (10X, Drop-seq)                | Smart-seq2, full-length data                      |
| Multi-batch (better variance stabilization)   | Very large datasets (> 50k cells) on tight RAM    |
| Default for most production workflows         | Need exact reproducibility with classic workflows |

If SCTransform crashes from memory, fall back to LogNormalize. The downstream pipeline is identical.

## 4. Batch integration

```
Multi-batch?
  no  → skip integration, use PCA
  yes → continue

> 100k cells with similar compositions?  → RPCA
Complex batch effects, < 100k cells?      → CCA
Speed and simplicity?                     → Harmony (default)
Have partial cell-type labels?            → Azimuth label transfer or scANVI (via reticulate)
```

Always check the un-integrated UMAP first. If batches mix already, skip integration.

## 5. Clustering resolution

| Goal                        | Resolution range |
| --------------------------- | ---------------- |
| Major lineages only         | 0.3–0.5          |
| Standard subtypes           | 0.6–0.8          |
| Fine substates              | 1.0–1.5          |
| Over-clustering (testing)   | > 2.0            |

Run a sweep (`cluster_multiple_resolutions(resolutions = c(0.4, 0.6, 0.8, 1.0))`) and choose by marker coherence and stability.

## 6. Cell type annotation

| Manual                                        | Automated (SingleR / Azimuth)                     |
| --------------------------------------------- | ------------------------------------------------- |
| Accurate when you have markers expertise      | Fast, reasonable first pass                       |
| Customizable to your specific cells           | Limited to reference's cell types                 |
| Time-intensive                                | Known failure modes (ILC→NK, HSC over-call)       |
| Required for novel/rare populations           | Pick the matching reference (HPCA, Azimuth pbmc)  |

Hybrid is best: automated for the first pass → manual refinement against your markers.

## 7. Cell-level DE vs pseudobulk DE

| Cell-level (`FindMarkers`)                    | Pseudobulk (`run_pseudobulk_deseq2`)              |
| --------------------------------------------- | ------------------------------------------------- |
| Cluster characterization                      | Condition contrasts across replicates             |
| Exploratory analysis                          | Inferential statistics for publication            |
| Single dataset, no condition contrast         | Multi-sample, multi-condition design              |
| **Never report p-values for condition effects** | Reports valid p-values with FDR control         |

This is the single most important decision in the pipeline. Get it wrong and your paper doesn't replicate.

## When in doubt

- Run the un-integrated UMAP first.
- Visualize QC distributions before setting thresholds.
- Test multiple clustering resolutions, pick by biology.
- Use the manual + automated hybrid for annotation.
- Use pseudobulk for any condition contrast.
- Always set seeds.
- Document everything in the report.
