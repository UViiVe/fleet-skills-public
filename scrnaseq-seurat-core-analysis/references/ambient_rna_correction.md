# Ambient RNA Correction (Seurat)

Droplet-based scRNA-seq captures cell-free RNA along with cells. The "soup" — lysed-cell content in the suspension — contaminates every droplet with the same background, producing false-positive expression of cell-type-specific markers in the wrong cells.

In Seurat workflows, the standard correction tool is **SoupX**. CellBender (Python) is an alternative when a GPU is available.

## When to correct

Always:

- Raw CellRanger output (`raw_feature_bc_matrix/`).
- High-soup tissues: brain, lung, tumor, lymph nodes, anything heavily dissociated.
- Long processing time between dissociation and droplet generation.
- Poor sample quality observed in the wet lab.

Skip:

- Pre-filtered data already corrected upstream.
- Low-soup tissues processed fresh (PBMC, sorted cells).
- Nuclei-based methods (snRNA-seq) — much less ambient RNA.

## SoupX (R-native)

`run_soupx_correction()` wraps the SoupX workflow. Needs both raw and filtered matrices.

```r
library(SoupX)

sc <- load10X("path/to/cellranger/outs/")
sc <- autoEstCont(sc)               # auto-estimate contamination fraction
corrected_counts <- adjustCounts(sc)

# Make a Seurat object from the corrected counts
seurat_obj <- CreateSeuratObject(counts = corrected_counts)
```

Parameters:

- `autoEstCont`: automatic estimation — default and correct most of the time.
- `setContaminationFraction`: manual override when auto-estimate misses (high-soup tissues sometimes need ~ 0.15–0.20).
- `forceAccept`: TRUE if the auto-estimator throws a "high contamination" warning for known high-soup tissues.

Interpreting contamination level:

| Estimate     | Interpretation        | Action                          |
| ------------ | --------------------- | ------------------------------- |
| < 5%         | Low                   | Continue                        |
| 5–15%        | Moderate              | Typical for many tissues        |
| 15–25%       | High                  | Review sample quality           |
| > 25%        | Very high             | Consider excluding the sample   |

## CellBender (Python alternative)

CellBender's deep generative model often outperforms SoupX on complex contamination (especially brain and tumor). Requires a GPU.

```bash
cellbender remove-background \
  --input raw_feature_bc_matrix.h5 \
  --output cellbender_output.h5 \
  --expected-cells 10000 \
  --total-droplets-included 30000 \
  --epochs 200 \
  --cuda
```

Read the output back into R via `Read10X_h5("cellbender_output_filtered.h5")` and wrap in a Seurat object.

## Validation

After correction, confirm:

1. **Total UMI reduction** of 5–20%:
   ```r
   reduction <- 1 - mean(colSums(GetAssayData(after, slot = "counts"))) /
                    mean(colSums(GetAssayData(before, slot = "counts")))
   message(sprintf("Mean UMI reduction: %.1f%%", reduction * 100))
   ```

2. **Marker specificity** — cell-type markers should be more restricted to the right populations:
   ```r
   FeaturePlot(before, features = c("CD3D", "MS4A1"))
   FeaturePlot(after,  features = c("CD3D", "MS4A1"))
   ```

3. **Hemoglobin removal** — in non-erythroid tissues, `HBA1`/`HBB` should drop sharply in T cells, B cells, etc.

## Tissue-specific notes

### Brain (15–30% contamination expected)

Neurons lyse easily; synaptic vesicle content is rich. Common contaminants: `SYT1`, `SNAP25`, `NRXN1/3`. Manual contamination fraction of 0.15–0.20 often works better than auto-estimation.

### Lung (10–25%)

Type II pneumocytes are abundant; surfactant proteins (`SFTPC`, `SFTPB`, `SFTPA1`) dominate the soup. SoupX usually adequate; CellBender preferred if dissociation was harsh.

### Tumor (10–30%)

Necrotic cells release RNA; dissociation often harsh. Common contaminants: cell-cycle, heat-shock, hemoglobin. CellBender preferred for complex patterns.

### PBMC (2–8%)

Fresh blood is the easiest case. Contamination is mostly hemoglobin from a few contaminating erythrocytes. Correction is optional; SoupX is sufficient if you do correct.

## Troubleshooting

### `autoEstCont` errors with "highly contaminated"

Either set the contamination fraction manually:

```r
sc <- setContaminationFraction(sc, contamination = 0.15)
```

Or pass `forceAccept = TRUE` if you've verified the estimate is sensible.

### Over-correction

- Estimate > 30%? Suspicious — verify by inspecting top contaminating genes.
- Loss of biological variation: validate with markers; re-run with manual lower contamination fraction.

### CellBender doesn't converge

- Increase `--epochs` to 300–400.
- Use `--total-droplets-included` at 2.5× expected instead of 3×.
- Confirm input is the raw matrix, not filtered.

### GPU out of memory

- Reduce `--total-droplets-included`.
- Drop `--cuda` and run on CPU (slow but works).

## Best practices

- Start with raw matrices when correction is needed.
- Document estimated contamination in the report.
- Compare before/after with markers, not just summary statistics.
- Use tissue-appropriate parameters (don't use PBMC settings on brain).
- Run correction per sample / batch.
- Do not double-correct.
- Do not use normalized counts as input — both methods assume raw counts.

## Citations

1. Young MD, Behjati S. *SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data.* GigaScience 9:giaa151 (2020).
2. Fleming SJ, Chaffin MD, Arduini A, et al. *Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender.* Nat Methods 20:1323–1335 (2023).
3. Xi NM, Li JJ. *Benchmarking computational doublet-detection methods for single-cell RNA sequencing data.* Cell Syst 12:176–194.e6 (2021).
