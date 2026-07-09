# Statistical Methods

The screen splits into two statistical problems: (1) outlier detection — flagging cells whose transcriptome no longer looks like a non-targeting control; and (2) differential expression — quantifying per-gene effects between a perturbation and the control. This document covers both with parameter advice.

---

## Outlier Detection in PCA Space

Strategy: for each perturbation, find genes differentially expressed against control, project both groups into PCA on those genes, fit an outlier model to controls, then classify perturbed cells.

### Why not just cluster?

Clustering requires that perturbation effects be large enough to form distinct clusters. Outlier methods detect *cell-level shifts* and tolerate sparse phenotypes — e.g., a perturbation that affects 30% of cells will be invisible in clustering but obvious in outlier-fraction scoring.

### LocalOutlierFactor (default)

Density-based: for each perturbed cell, compare its local k-NN density to its neighbors' local densities. Outliers sit in sparser regions than their nearest controls.

```python
from sklearn.neighbors import LocalOutlierFactor
clf = LocalOutlierFactor(n_neighbors=20, contamination="auto", novelty=True).fit(control_pca)
preds = clf.predict(perturbed_pca)   # -1 outlier, 1 inlier
```

Strengths: sensitive to localized phenotypes; fast; small-sample-friendly. Weaknesses: k-sensitive; high dimensional collapse (use PCA first, n_pcs=4–10).

Tuning:
- Conservative: `n_neighbors=30, contamination=0.03`
- Default: `n_neighbors=20, contamination="auto"`
- Sensitive: `n_neighbors=10, contamination=0.10`

### IsolationForest

Random splits isolate outliers in fewer splits than inliers. Tree-based, ensemble.

```python
from sklearn.ensemble import IsolationForest
clf = IsolationForest(n_estimators=100, contamination="auto", random_state=42).fit(control_pca)
preds = clf.predict(perturbed_pca)
```

Strengths: scales well; handles high dimensions; detects *global* outliers. Weaknesses: misses subtle localized phenotypes; needs ≥100 controls; stochastic.

### OneClassSVM

Learns a non-linear decision boundary around the control distribution (RBF kernel).

```python
from sklearn.svm import OneClassSVM
clf = OneClassSVM(kernel="rbf", gamma="scale", nu=0.1).fit(control_pca)
preds = clf.predict(perturbed_pca)
```

Strengths: complex non-linear boundaries; small-sample-friendly. Weaknesses: O(n²)–O(n³); hyperparameter-sensitive.

### When to use which

| Situation | Method |
|---|---|
| Default starting point | LocalOutlierFactor |
| Too many control outliers | IsolationForest (or lower LOF contamination) |
| Few hits, weak phenotypes | Lower LOF contamination, or try IsolationForest |
| Small screen, complex manifold | OneClassSVM |

### Ensemble call

For maximum confidence, require agreement:

```python
pred_lof = LocalOutlierFactor(novelty=True).fit(control_pca).predict(perturbed_pca)
pred_iso = IsolationForest(random_state=42).fit(control_pca).predict(perturbed_pca)
consensus = (pred_lof == -1) & (pred_iso == -1)
```

---

## Differential Expression

### t-test (Welch)

```python
import diffxpy.api as de
test = de.test.t_test(data=adata, grouping="gene", is_logged=True)
```

Closed-form, fast (sub-second per perturbation), log-data-friendly. Use for the **screening tier**. Assumes approximate normality (log1p of CPM gets you most of the way there) and unequal variances (Welch).

### Wilcoxon rank-sum

```python
test = de.test.wilcoxon(data=adata, grouping="gene")
```

Non-parametric, robust to outliers, slower than t-test, weaker power when normality holds. Use when cells per group are small (<50) or distributions are obviously non-normal.

### Negative-binomial GLM (Wald)

```python
test = de.test.wald(data=adata, formula_loc="~ 1 + gene", factor_loc_totest="gene")
```

Models raw counts directly with NB dispersion. Slower but principled. Use for the **validation tier** on the top ~50 hits.

### glmGamPoi (recommended for the validation tier)

`glmGamPoi` is a fast NB GLM implementation with shrinkage of dispersion estimates. It handles fixed effects (batch, donor) cleanly. Operates on raw counts.

```R
fit <- glm_gp(sce, design = ~ gene + donor - 1)
res <- test_de(fit, contrast = ~ `geneTARGET` - `genenon-targeting`)
```

The python helper `run_glmgampoi_batch()` shells out to R via subprocess and reads back per-gene CSVs.

### MAST (zero-inflated hurdle model)

Two-part GLM: logistic for the detection rate plus linear for the conditional expression. Was the field standard for early scRNA-seq DE but is overkill for most CRISPR-screen comparisons today and is slow.

### Speed comparison (1k genes, 100 cells/group)

| Method | Wall time | RAM |
|---|---|---|
| t-test | <1 s | low |
| Wilcoxon | ~5 s | low |
| Negative binomial (diffxpy) | ~30 s | medium |
| glmGamPoi (R) | ~1 min | medium |
| MAST | ~5 min | high |

### Recommended pipeline

```
Screening (all perturbations)   → t-test
↓ filter to hits with ≥ 10 DE genes
Target validation                → check target log2FC direction
↓ retain validated only
Refinement (top ~50)             → glmGamPoi with batch correction
```

---

## Multiple-Testing Correction

### Benjamini–Hochberg FDR (default)

```python
from statsmodels.stats.multitest import multipletests
_, qvals, _, _ = multipletests(pvals, method="fdr_bh")
```

Controls the expected false-discovery rate. Appropriate for screening — you can tolerate 5% of "discoveries" being false in exchange for power.

### Bonferroni

```python
alpha = 0.05 / len(pvals)
```

Controls family-wise error rate. Use for small, high-confidence screens (<100 perturbations) where any false positive is costly.

### Permutation-based FDR

When the test's null distribution is suspect (e.g. heavily zero-inflated data with unusual variance structure), permute the perturbation label N times, redo the test, and compute the empirical fraction of null statistics that exceed the observed threshold.

---

## Effect-Size Metrics

### log2 fold-change

```python
log2fc = np.log2(mean_perturbed / mean_control)
```

Thresholds:
- |log2FC| > 0.5: modest (1.4x)
- |log2FC| > 1.0: strong (2x)
- |log2FC| > 2.0: very strong (4x)

Always couple a fold-change cut with a p-value cut.

### Cohen's d

```python
def cohens_d(a, b):
    pooled = ((len(a)-1)*a.var(ddof=1) + (len(b)-1)*b.var(ddof=1)) / (len(a)+len(b)-2)
    return (a.mean() - b.mean()) / pooled**0.5
```

- |d| < 0.2: small
- |d| = 0.5: medium
- |d| > 0.8: large

Useful for cross-perturbation comparison when raw fold-changes are confounded by mean expression.

---

## Practical Recommendations

Outlier detection:
1. Start LOF with `n_neighbors=20, contamination="auto"`
2. Compute control outlier rate — target < 5%
3. If too high, lower contamination to 0.03 or switch to IsolationForest
4. If too few hits, raise contamination to 0.10 or lower `min_de_genes` from 10 to 5

Differential expression:
1. Screen with t-test
2. Validate target genes (direction + magnitude)
3. Refine top hits with glmGamPoi if batch correction matters
4. Always pair p-values with effect-size thresholds
5. FDR-correct genome-wide; Bonferroni only for small focused screens

### Reporting checklist

- Statistical test used at each stage
- Multiple-testing-correction method
- Effect-size metric and threshold
- Number of comparisons performed
- Sample sizes per group
- Batch / donor structure of the data and how it was modelled

Example: *"Differential expression was assessed with Welch's t-test on log1p(CPM) and then validated for the top 50 hits with `glmGamPoi` (NB GLM with batch as a fixed effect). P-values were adjusted with Benjamini–Hochberg; genes with FDR < 0.05 and |log2FC| > 0.5 were called significant."*
