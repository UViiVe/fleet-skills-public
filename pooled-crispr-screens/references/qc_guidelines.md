# Quality Control Guidelines

QC for a pooled CRISPR screen has three layers: standard scRNA-seq QC (per cell), sgRNA-assignment QC (per cell with its mapped guide), and screen-level QC (per perturbation, per batch, per library). This reference covers all three with cell-type-specific recommendations.

---

## Cell-Type-Specific Thresholds

These are starting points. Always inspect distributions and adjust to the inflection points.

### Primary cortical neurons

```python
qc = {"min_genes": 2500, "max_genes": 8000,
      "min_counts": 8000, "max_counts": 50000,
      "max_mito_pct": 0.12, "max_ribo_pct": 0.40}
```

Reasoning: neurons are large, transcriptionally complex (8–12k genes detected at saturation), and metabolically sensitive. Mito% above ~12% almost always means apoptosis.

### iPSC-derived neurons

```python
qc = {"min_genes": 2000, "max_genes": 7000,
      "min_counts": 5000, "max_counts": 40000,
      "max_mito_pct": 0.18, "max_ribo_pct": 0.40}
```

Reasoning: differentiation heterogeneity broadens every distribution. Looser thresholds preserve real biology at the cost of some debris.

### PBMCs

```python
qc = {"min_genes": 500, "max_genes": 5000,
      "min_counts": 1000, "max_counts": 30000,
      "max_mito_pct": 0.20, "max_ribo_pct": 0.50}
```

Reasoning: PBMCs span sizes — small naive T cells to large activated monocytes — so the floor is low.

### Macrophages (THP-1, primary, BMDM)

```python
qc = {"min_genes": 1000, "max_genes": 6000,
      "min_counts": 3000, "max_counts": 40000,
      "max_mito_pct": 0.22, "max_ribo_pct": 0.45}
```

Activated (M1) macrophages run high mito% from oxidative metabolism — do not mistake activation for death.

### Cancer cell lines (K562, HEK293T, HeLa, etc.)

```python
qc = {"min_genes": 1500, "max_genes": 6000,
      "min_counts": 5000, "max_counts": 40000,
      "max_mito_pct": 0.22, "max_ribo_pct": 0.45}
```

---

## Standard scRNA-seq QC Metrics

### `n_genes_by_counts`

- <500 in any cell type: debris, lysed cells, or empty droplets
- 1k–5k: healthy range (cell-type-dependent)
- >8k in most populations: doublet candidate (two cells per droplet)

### `total_counts` (UMI)

- <2k: low-quality cell, poor capture
- 3k–30k: typical, depending on cell type and chemistry
- >50k in most populations: doublet candidate

Saturation depends on read depth:

| Reads / cell | Expected UMI / cell | Saturation |
|---|---|---|
| 20k | 3k–5k | 60–70% |
| 50k | 8k–15k | 70–80% |
| 100k | 15k–30k | 80–90% |

### `pct_counts_mt`

- <10%: clearly healthy
- 10–20%: moderate stress, still usable in most cell types
- >20%: dying or dead in most cell types

Watch out for cell types with intrinsically high mito content — cardiomyocytes, brown adipocytes, activated macrophages.

### Ribosomal fraction

Useful as a secondary signal. The combination **high ribo% AND high mito%** flags dying-but-still-translating cells better than either alone.

```python
ribo_genes = adata.var_names.str.startswith(("RPS", "RPL"))
adata.obs["pct_counts_ribo"] = (
    np.asarray(adata[:, ribo_genes].X.sum(axis=1)).flatten()
    / adata.obs["total_counts"]
)
```

---

## sgRNA Assignment QC

### Mapping rate

How many of your scRNA-seq cells carry an assigned sgRNA?

| Technology | Expected | Investigate if |
|---|---|---|
| 10X Feature Barcoding | 50–70% | <40% |
| Direct capture (Replogle) | 60–80% | <50% |
| Separate sgRNA-seq library | 70–90% | <60% |

Diagnostics for low mapping:
- Run the Feature Barcode library on Bioanalyzer — confirm the right insert size
- Verify the sgRNA reference matches the actual library design (typos kill this)
- Check viral MOI — too low means few cells were ever transduced
- Inspect cell-barcode overlap between mapping file and AnnData

### Singlet criteria

Default: a cell is a singlet if **>70% of its sgRNA UMIs come from one guide**, with **>=2 UMIs total**. Drop cells where the top guide is <70% — these are doublets at the guide level even if not at the cell level.

Relaxed alternative (use only when capture is genuinely poor): 60% / 1 UMI. Expect more false assignments.

### Doublet rate

A function of MOI and ambient RNA. Track it per library:

| MOI | Expected doublets | Acceptable |
|---|---|---|
| 0.3 | 5–8% | <10% |
| 0.5 | 8–12% | <15% |
| 0.8 | 15–20% | <20% |

---

## Per-Perturbation Coverage

Cells per perturbation after QC determines downstream statistical power:

| Cells per perturbation | What you can do | Recommendation |
|---|---|---|
| <10 | Nothing reliable | Exclude or flag |
| 10–20 | Underpowered | Include with low confidence |
| 20–50 | Adequate for t-test, marginal for GLM | Standard |
| 50–100 | Comfortable for everything | Ideal |
| >100 | Detects small effect sizes | Best |

When you see a perturbation with <10 cells, ask:
- Is it an essential gene that dropped out? (Expected for CRISPRko; flag, do not delete.)
- Was the sgRNA poorly represented in the library?
- Is it absent in only some batches? (Then it is a batch-specific assignment artefact.)

```python
cells_per_perturbation = adata.obs.groupby("gene").size()
crosstab = pd.crosstab(adata.obs["gene"], adata.obs["batch"])
print(crosstab[crosstab.min(axis=1) == 0])   # perturbations missing in some batch
```

---

## Control Distribution

Non-targeting controls should be:
- 5–15% of total cells after QC (depends on library design)
- Roughly evenly distributed across batches
- Tightly clustered in PCA (a diffuse control cloud is a red flag for batch effects)

```python
n_controls = (adata.obs["gene"] == "non-targeting").sum()
print(f"Controls: {n_controls} ({n_controls / adata.n_obs:.1%})")
```

---

## Target-Gene Validation

After the initial DE screen, ask: did the perturbation actually move its own target gene? This is your strongest filter against non-functional sgRNAs.

### CRISPRi expected target effect

| log2FC | Fold change | Interpretation |
|---|---|---|
| < -1.0 | >50% reduction | Strong knockdown |
| -1.0 to -0.5 | 30–50% reduction | Moderate, acceptable |
| -0.5 to 0 | <30% reduction | Weak — likely non-functional |
| >= 0 | No change or increase | Inactive sgRNA, remove |

Validation thresholds in practice: `min_log2fc = -0.5`, `max_pval = 0.05`.

### CRISPRa expected target effect

| log2FC | Fold change | Interpretation |
|---|---|---|
| > 2.0 | >4x increase | Strong activation |
| 1.0–2.0 | 2–4x increase | Moderate, acceptable |
| 0.5–1.0 | 1.4–2x | Weak |
| <= 0 | No change or decrease | Inactive sgRNA |

### Library validation rate

| Validation rate | Library quality | Action |
|---|---|---|
| 60–80% | Good | Proceed |
| 50–60% | Acceptable | Proceed but filter to validated |
| <50% | Poor | Troubleshoot the perturbation system before trusting results |

---

## Batch Effects

### Detection

```python
import scanpy as sc, statsmodels.api as sm, pandas as pd
sc.tl.pca(adata)
X = pd.get_dummies(adata.obs["batch"])
r2 = [sm.OLS(adata.obsm["X_pca"][:, k], X).fit().rsquared for k in range(10)]
print(f"Mean variance in PC1–10 explained by batch: {sum(r2)/10:.1%}")
```

- <5% — no correction needed
- 5–15% — consider correction; check that biology is preserved
- >15% — correction recommended

### Correction options

- **Harmony** — fast, preserves biological structure (`scanpy.external.pp.harmony_integrate`)
- **Seurat v3 integration** — strong for small numbers of batches
- **BBKNN** — k-NN-level correction
- **`sc.pp.regress_out(["batch"])`** — linear, blunt; useful only if effects are small

### Per-perturbation batch balance

Check whether each perturbation appears in all batches:

```python
crosstab = pd.crosstab(adata.obs["gene"], adata.obs["batch"])
cv = crosstab.std(axis=1) / crosstab.mean(axis=1)
print(cv.nlargest(10))   # most batch-imbalanced perturbations
```

Imbalanced perturbations (CV > 1) should be analyzed within-batch or excluded from cross-batch comparisons.

---

## Outlier-Detector QC

After running outlier detection, check the control false-positive rate:

```python
control_outlier_rate = (
    results["perturbation_summary"]
    .query("gene == 'non-targeting'")["outlier_fraction"]
    .iloc[0]
)
```

- <5%: well-calibrated detector
- 5–10%: somewhat elevated — try `contamination=0.03`
- >10%: false-positive-driven — investigate batch effects first

Hit-rate expectations:

| Screen type | Expected hit rate |
|---|---|
| Targeted essential-gene screen | 5–10% |
| Targeted pathway screen | 10–30% |
| Candidate-gene screen | 20–50% |
| Unbiased discovery | 1–5% |

---

## QC Checklist (Pre-Analysis)

- [ ] sgRNA mapping rate >= 40%
- [ ] Doublet rate <= 15%
- [ ] Mean cells per perturbation >= 30
- [ ] <10% of perturbations have <20 cells
- [ ] Controls >= 5% of total cells
- [ ] Batch variance on PC1–10 <= 15%
- [ ] QC distributions match expectations for the cell type
- [ ] No obvious technical artefacts in PCA (rings, stripes, isolated clusters of unidentified origin)
