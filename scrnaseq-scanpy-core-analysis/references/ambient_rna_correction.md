# Ambient RNA Correction

Droplet-based scRNA-seq captures cell-free RNA along with cells. Lysed cells in the suspension contaminate every downstream droplet with the same background — the "soup." When the soup is rich (high-expression genes from lysed cells), it shows up as false-positive expression of cell-type-specific markers in wrong cell types: neurons appearing to express immune genes, T cells looking erythroid, etc.

Correcting for this is the first preprocessing step when the data is raw and the tissue is soupy.

## When to correct

Always:

- Raw CellRanger output (`raw_feature_bc_matrix/`).
- High-soup tissues: brain, lung, tumor, lymph nodes, anything heavily dissociated.
- Long processing time between dissociation and droplet generation.
- Poor sample quality (high cell death observed during prep).

Skip:

- Pre-filtered data already corrected upstream.
- Low-soup tissues processed fresh: PBMC from fresh blood, sorted cells.
- Nuclei-based methods (snRNA-seq) — much less ambient RNA.

## CellBender vs SoupX

| Method     | Pros                                                                                | Cons                                                          | Best for                                       |
| ---------- | ----------------------------------------------------------------------------------- | ------------------------------------------------------------- | ---------------------------------------------- |
| CellBender | Deep generative model; auto-estimates contamination; handles complex patterns       | GPU recommended; 30–60 min/sample; less interpretable         | Most cases, especially brain/tumor             |
| SoupX      | Fast (minutes); interpretable; CPU only                                             | Needs both filtered and raw matrices; threshold selection     | Quick analyses, when CellBender unavailable    |

## CellBender

```bash
cellbender remove-background \
    --input raw_feature_bc_matrix.h5 \
    --output cellbender_output.h5 \
    --expected-cells 10000 \
    --total-droplets-included 30000 \
    --epochs 200 \
    --cuda
```

Key parameters:

- `--expected-cells`: pull from the CellRanger web summary, or estimate from the knee plot.
- `--total-droplets-included`: should be 2–3× expected cells. Too low → poor ambient profile estimation. Too high → slower, risks pulling in low-quality cells. For high-soup tissues, use 3× expected.
- `--epochs`: 150 default; 200–300 for complex/brain data.
- `--fpr`: cell-calling false positive rate, default 0.01. Lower → more conservative.
- `--cuda`: needed for practical speed. CPU mode is 10–20× slower.

Outputs:

- `cellbender_output_filtered.h5` — corrected counts, filtered cells.
- `cellbender_output_cell_barcodes.csv` — called barcodes.
- `cellbender_output_metrics.csv` — QC metrics.
- `cellbender_output.pdf` — diagnostic plots.

Python wrapper:

```python
from scripts.remove_ambient_rna import run_cellbender, estimate_contamination

adata = run_cellbender(
    raw_h5="raw_feature_bc_matrix.h5",
    expected_cells=10000,
    total_droplets=30000,
    output_dir="results/cellbender",
    epochs=200,
    use_cuda=True,
)
print(f"Contamination: {estimate_contamination(adata):.1%}")
```

## SoupX

Needs both the raw and filtered matrices. The raw matrix gives the ambient profile (from empty droplets); the filtered matrix gives the cell barcodes to correct.

R workflow:

```r
library(SoupX)
sc <- load10X("data/raw_feature_bc_matrix/")
sc <- autoEstCont(sc)            # auto-estimate contamination fraction
out <- adjustCounts(sc)          # corrected counts
```

Python wrapper:

```python
from scripts.remove_ambient_rna import run_soupx_python

adata = run_soupx_python(
    raw_matrix_dir="raw_feature_bc_matrix/",
    filtered_matrix_dir="filtered_feature_bc_matrix/",
)
```

Parameters:

- `auto_estimate=True` is the default and usually correct.
- `contamination_range=c(0.01, 0.5)` — bound the estimate. Use `c(0.05, 0.5)` for high-soup tissue.

## Interpreting contamination level

| Estimate     | Interpretation              | Action                            |
| ------------ | --------------------------- | --------------------------------- |
| < 5%         | Low                         | Continue                          |
| 5–15%        | Moderate                    | Expected for many tissues         |
| 15–25%       | High                        | Review sample quality             |
| > 25%        | Very high                   | Consider excluding the sample     |

## Validation

After correction, check:

1. **Total UMI reduction:** 5–20% is typical.
   ```python
   before = adata_before.X.sum(axis=1).mean()
   after = adata_after.X.sum(axis=1).mean()
   print(f"Reduction: {(1 - after / before) * 100:.1f}%")
   ```

2. **Marker specificity:** cell-type markers should be more restricted to the right populations after correction. Compare before/after UMAPs colored by `CD3D` (T cells), `MS4A1` (B cells), etc.

3. **Hemoglobin removal:** in non-erythroid samples, `HBA1`/`HBB` should drop dramatically in T cells, B cells, etc.

## Tissue-specific notes

### Brain (high soup, 15–30%)

Neurons lyse easily during dissociation; synaptic vesicle content is rich. Common contaminants: `SYT1`, `SNAP25`, `NRXN1`, `NRXN3`. Always use CellBender; set `--total-droplets-included` to 3× expected; bump epochs to 250–300.

### Lung (high soup, 10–25%)

Type II pneumocytes are abundant and release surfactant proteins. Common contaminants: `SFTPC`, `SFTPB`, `SFTPA1`. CellBender preferred; SoupX works for milder cases.

### Tumor (variable, 10–30%)

Necrotic and dying cells release RNA; dissociation is often harsh. Common contaminants: cell-cycle genes, heat-shock proteins, hemoglobin (from contaminating blood). CellBender needed for the complex contamination patterns; SoupX may under-correct.

### PBMC (low soup, 2–8%)

Fresh blood is the easiest case. Contamination is usually just hemoglobin from a small number of contaminating erythrocytes. Correction is optional; SoupX is sufficient if you do correct.

## Troubleshooting

### CellBender doesn't converge

- Increase `--epochs` to 300–400.
- Try `--total-droplets-included` at 2.5× instead of 3×.
- Confirm input is truly the raw matrix, not filtered.

### Over-correction

- Estimate > 30%? Lower with `--fpr 0.005` (CellBender) or manually set contamination fraction (SoupX).
- Lost biological variation: validate with markers; you may need to re-run with conservative parameters.

### GPU out of memory

- Reduce `--total-droplets-included`.
- Drop `--cuda` and run on CPU (slow but doable).
- Use a smaller batch size if exposed by the version.

## Best practices

- Start with raw matrices when correction is needed.
- Document estimated contamination in your report.
- Compare before/after using marker genes, not just summary statistics.
- Use tissue-appropriate parameters (don't use PBMC settings on brain).
- Run correction per sample/batch.
- Do not double-correct — if you got data that was already CellBender-processed, stop.
- Do not use normalized counts as input — both methods assume raw counts.

## Citations

1. Fleming SJ, Chaffin MD, Arduini A, et al. *Unsupervised removal of systematic background noise from droplet-based single-cell experiments using CellBender.* Nat Methods 20:1323–1335 (2023).
2. Young MD, Behjati S. *SoupX removes ambient RNA contamination from droplet-based single-cell RNA sequencing data.* GigaScience 9:giaa151 (2020).
3. Xi NM, Li JJ. *Benchmarking computational doublet-detection methods for single-cell RNA sequencing data.* Cell Syst 12:176–194.e6 (2021).
