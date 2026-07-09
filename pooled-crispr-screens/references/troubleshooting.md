# Troubleshooting

Symptoms, diagnoses, and fixes for the failure modes that come up in practice.

---

## Data Loading

### `ValueError: Variable names are not unique`

10X feature matrices sometimes carry duplicate gene symbols (alias collisions). Make them unique before anything else:

```python
adata.var_names_make_unique()
```

To inspect what is duplicated:

```python
dups = adata.var_names[adata.var_names.duplicated()].unique()
```

### Low sgRNA mapping rate (<30%)

Diagnose first:

```python
df = pd.read_table("sgrna_mapping.txt", header=None)
print(df.head(), df[1].nunique(), len(df))
overlap = set(adata.obs_names) & set(df.index)
print(f"Barcode overlap: {len(overlap)}/{adata.n_obs} ({len(overlap)/adata.n_obs:.1%})")
```

Likely causes and fixes:
- **Barcode-format mismatch**: re-process CellRanger with the correct feature reference
- **Low viral titer**: re-transduce at MOI 0.3–0.5
- **Library-prep issue**: increase PCR cycles for the Feature Barcode amplification (12 → 15)
- **sgRNA reference wrong**: verify the mapping file's sgRNA IDs exactly match the library design

### High doublet rate (>15%)

```python
sns.histplot(adata.obs["sgrna_umi_count"], bins=50)   # diagnostic
```

Causes: MOI too high; ambient RNA contamination; barcode collision. Fixes: tighten singlet threshold to >70% UMI from one guide; dead-cell removal pre-capture; for future runs, lower MOI.

---

## QC

### Too many cells dropped (>50% loss)

Plot the distributions first:

```python
fig, axes = plt.subplots(1, 3, figsize=(15, 4))
axes[0].hist(adata.obs["n_genes"],     bins=50); axes[0].axvline(2000, c="r", ls="--")
axes[1].hist(adata.obs["n_counts"],    bins=50); axes[1].axvline(8000, c="r", ls="--")
axes[2].hist(adata.obs["percent_mito"], bins=50); axes[2].axvline(0.15, c="r", ls="--")
```

If the inflection point sits below the threshold, the threshold is too strict. Relax to cell-type-appropriate values (see `qc_guidelines.md`).

Two-stage QC for perturbations that cause cellular stress:

```python
control = apply_qc_filters(adata[adata.obs["gene"] == "non-targeting"],
                           min_genes=2500, max_mito_pct=0.12)
perturbed = apply_qc_filters(adata[adata.obs["gene"] != "non-targeting"],
                             min_genes=2000, max_mito_pct=0.20)
adata_filt = ad.concat([control, perturbed])
```

### Specific perturbations with <10 cells

```python
counts = adata.obs.groupby("gene").size().sort_values()
print(counts.head(10))
```

If they are known essential genes, expect depletion in CRISPRko. Otherwise: check the plasmid pool by NGS, check the viral library, check per-batch cross-tabulation for batch-specific dropout.

---

## Normalization

### "Highly-expressed genes are not considered"

This is `sc.pp.normalize_total` doing the right thing. `MALAT1`, ribosomal genes, and a handful of others dominate total counts and bias normalization if left in. To make it explicit:

```python
sc.pp.normalize_total(adata, target_sum=1e6, exclude_highly_expressed=True)
```

To exclude a custom set:

```python
adata = adata[:, ~adata.var_names.isin(["MALAT1", "NEAT1"])].copy()
sc.pp.normalize_total(adata, target_sum=1e6)
```

---

## Outlier Detection

### >50% of cells flagged as outliers

Check control outlier rate first:

```python
print(results["perturbation_summary"]
      .query("gene == 'non-targeting'")["outlier_fraction"]
      .iloc[0])
# should be < 0.05
```

If high: tighten contamination (`LocalOutlierFactor(contamination=0.03)`), check batch effects (run Harmony on `X_pca`), and ensure ≥100 control cells per batch.

### Too few outliers, known positives missed

Check DE-gene counts per perturbation:

```python
de_counts = results["perturbation_summary"]["n_de_genes"]
print(de_counts.describe(), (de_counts < 10).sum())
```

Check the target-gene shift directly:

```python
target = "SOX5"
ctrl = adata[adata.obs["gene"] == "non-targeting", target].X.mean()
prt  = adata[adata.obs["gene"] == target, target].X.mean()
print(f"{target}: {prt/ctrl:.2f}x  (CRISPRi expect < 0.5)")
```

Loosen the detector (`contamination=0.10`) or switch to IsolationForest. If target genes are not moving at all, the perturbation system itself is suspect.

---

## Differential Expression

### Few DE genes for most perturbations

Causes:
- <20 cells per perturbation → underpowered. Aggregate to pseudobulk if you must.
- Low sequencing depth → median UMI per cell <2000 means dropout dominates
- Weak perturbation effects → check target validation; the perturbation may just be that subtle

Adaptations:
- Pseudobulk by perturbation: `adata.to_df().groupby(adata.obs["gene"]).sum()`
- Restrict DE to outlier cells only (improves signal when most perturbed cells are unaffected)
- Loosen p-value to 0.10 for the screening tier

### `MemoryError` during DE

```python
adata = sc.read_h5ad("file.h5ad", backed="r")   # disk-backed
```

Or batch processing 100 perturbations at a time and write incrementally. Or downsample to 1000 cells per perturbation for the screening tier.

---

## Batch Effects

### Strong batch separation on PCA

Quantify:

```python
from sklearn.linear_model import LinearRegression
X = pd.get_dummies(adata.obs["batch"])
y = adata.obsm["X_pca"][:, 0]
print(f"PC1 R² by batch: {LinearRegression().fit(X, y).score(X, y):.1%}")
```

Above ~15% → correct. Options:

- **Harmony** (fast, biology-preserving):
  ```python
  import scanpy.external as sce
  sce.pp.harmony_integrate(adata, "batch", basis="X_pca")
  ```
- **`sc.pp.regress_out(["batch"])`** + rescale (linear, blunt; useful for small per-batch shifts)
- **Per-batch analysis**: if batches are biologically distinct and you cannot share controls, analyze each separately and intersect hits

After correction, verify that positive controls still validate. If they don't, you over-corrected.

---

## Visualization

### Featureless UMAP

Probable causes:
- Insufficient preprocessing — confirm normalize + log + scale ran
- Too few PCs for neighbor graph — use `n_pcs=30–50`
- Use HVG-selected matrix:

  ```python
  sc.pp.highly_variable_genes(adata, n_top_genes=2000)
  adata_hvg = adata[:, adata.var["highly_variable"]].copy()
  sc.tl.pca(adata_hvg); sc.pp.neighbors(adata_hvg); sc.tl.umap(adata_hvg)
  ```

- Adjust UMAP: `min_dist=0.5` for global structure, `n_neighbors=30` for smoother embedding

### KeyError on a target gene

Gene-symbol drift between the sgRNA library and the expression matrix:

```python
from scripts.gene_name_corrections import detect_mismatches, suggest_corrections, correct_gene_names
mismatches = detect_mismatches(adata)
corrections = suggest_corrections(adata, use_synonyms=True)   # needs mygene
adata = correct_gene_names(adata, corrections)
```

Common manual corrections: `TMEM55A → PIP4P2`, `ATP5C1 → ATP5F1C`, `ATP5H → ATP5PD`.

---

## Export

### Excel file too large

Excel sheets max out at 256 in most readers and the file balloons beyond ~100 MB. Split:

```python
for i, batch in enumerate(np.array_split(list(de_results.keys()), 10)):
    with pd.ExcelWriter(f"de_batch_{i+1}.xlsx") as writer:
        for gene in batch:
            de_results[gene].head(50).to_excel(writer, sheet_name=gene[:31])
```

Or just emit CSVs — that is what downstream tools usually want anyway.

### SVG export silently absent

PNG is always saved as fallback. If you actually need SVG, install a matplotlib SVG backend (`pip install matplotlib --upgrade` usually suffices, otherwise check `matplotlib.get_backend()`).

---

## glmGamPoi Subprocess Failures

`Rscript: command not found` → install R 4.0+, ensure `Rscript` is on PATH.

`there is no package called 'glmGamPoi'` → in R:

```R
install.packages("BiocManager")
BiocManager::install(c("glmGamPoi", "zellkonverter"))
```

`Error in readH5AD: file not found` → the python helper writes a temporary h5ad and passes its path. Verify temp directory is writable.

Convergence failures on individual genes are normal for low-count genes — glmGamPoi will skip them and continue.

---

## Getting Help

When the workflow still does not behave:

1. Re-examine CellRanger output logs — sequencing depth, saturation, percent reads in cells
2. Inspect raw Bioanalyzer / TapeStation traces of the Feature Barcode library
3. Confirm dCas9 expression in the cell line (Western, fluorescent reporter) — positive controls failing is usually a perturbation-system failure
4. Compare against published reference datasets (Papalexi 2021, Replogle 2022) using the same code path
5. Validate top hits experimentally before trusting them — qPCR / Western is much cheaper than chasing a computational artefact through follow-up biology
