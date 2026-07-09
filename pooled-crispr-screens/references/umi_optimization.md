# UMI-Threshold Optimization for sgRNA Assignment

The UMI cutoff for declaring "this cell carries this sgRNA" trades sensitivity (more cells assigned) against specificity (fewer false assignments from ambient RNA). When capture efficiency is the bottleneck — large genome-wide screens or library-prep quality issues — sweeping the threshold and picking the elbow is worth the compute time.

Method based on Wang et al. (2025) "Integrative analysis of pooled CRISPR screens provide functional insights into AD GWAS risk genes".

---

## When This Is Worth Doing

Use it when:
- The screen is large (>500 perturbations) — the payoff scales with the number of perturbations whose call quality shifts
- Mapping rate is <40% under standard thresholds
- You have compute headroom for ~5× the normal runtime
- The publication / decision turns on hit-list confidence

Skip when:
- The screen is small (<100 perturbations) — the noise floor of "is this an optimal cutoff?" exceeds the gain
- Mapping is already >60% under default thresholds
- It's a pilot

---

## Sweep

Test 3–5 thresholds, run the full pipeline at each, and tabulate hit counts.

```python
umi_thresholds = [3, 5, 10, 15, 20]
results = {}

for t in umi_thresholds:
    adata_t = assign_guides_with_threshold(adata, mapping_file, min_umi=t)
    cells_per_gene = adata_t.obs.groupby("gene").size()
    results[t] = {
        "n_cells": adata_t.n_obs,
        "mean_cells_per_gene": cells_per_gene.mean(),
        "median_cells_per_gene": cells_per_gene.median(),
        "genes_with_50plus_cells": (cells_per_gene >= 50).sum(),
        "genes_with_low_cells": (cells_per_gene < 20).sum(),
    }
```

Then run QC, normalization, and the screening DE at each threshold:

```python
for t, adata_t in adata_by_threshold.items():
    adata_f = apply_qc_filters(adata_t)
    adata_n = normalize_only(adata_f)
    de = screen_all_perturbations(adata_n, control_group="non-targeting")
    n_hits = sum(1 for _, df in de.items() if (df["qval"] < 0.05).sum() >= 10)
    results[t]["n_hits"] = n_hits
```

Summarize:

```python
summary = pd.DataFrame(results).T.rename_axis("UMI_threshold")
print(summary)
```

A typical outcome:

| UMI threshold | Cells | Cells/gene mean | Genes ≥50 cells | Hits |
|---|---|---|---|---|
| 3 | 45,000 | 65 | 850 | 125 |
| **5** | **42,000** | **58** | **820** | **148** |
| 10 | 35,000 | 48 | 720 | 142 |
| 15 | 28,000 | 38 | 580 | 118 |
| 20 | 22,000 | 30 | 450 | 95 |

The optimum is usually the elbow where hits stop increasing while cells-per-gene starts collapsing — here, 5.

---

## Interpretation Guide

Lower thresholds (3–5 UMI):
- More cells assigned → better statistical power
- Higher risk of ambient-RNA misassignment
- Use when capture is poor

Higher thresholds (10–20 UMI):
- Fewer false assignments
- Less ambient-RNA contamination
- Use when capture is high and ambient RNA is the concern

Balanced (5–8 UMI):
- The standard regime for clean 10X Feature Barcoding screens

---

## Quality Monitors at Each Threshold

```python
for t, ad_t in adata_by_threshold.items():
    mapping = ad_t.n_obs / total_cells_pre_assignment
    doublets = doublets_detected[t] / cells_with_sgrna[t]
    cv = ad_t.obs.groupby("gene").size().std() / ad_t.obs.groupby("gene").size().mean()
    print(f"UMI={t}: map={mapping:.1%}  doublets={doublets:.1%}  cv={cv:.2f}")
```

Control-cell DE counts also catch bad thresholds:

```python
for t, de in de_by_threshold.items():
    ctrl_de = (de.get("non-targeting", pd.DataFrame())["qval"] < 0.05).sum()
    print(f"UMI={t}: control DE genes = {ctrl_de}  (expect < 500)")
```

A threshold that yields >500 control "DE" genes is contaminated with ambient signal.

---

## Computational Cost

Five thresholds ≈ 5× the standard runtime.

| Screen size | Default runtime | Sweep runtime |
|---|---|---|
| 100 perturbations | ~5 min | +25 min |
| 500 perturbations | ~30 min | +2 h |
| 1000 perturbations | ~1 h | +4 h |

RAM: each threshold makes its own AnnData (~2–5 GB per 50k cells), so process sequentially if you are RAM-constrained. With cores to spare:

```python
from concurrent.futures import ProcessPoolExecutor
with ProcessPoolExecutor(max_workers=5) as pool:
    results = dict(zip(umi_thresholds, pool.map(analyze_threshold, [adata]*5, umi_thresholds)))
```

---

## After Picking a Threshold

Validate the choice:
1. Target-gene validation rate at the chosen threshold should be >60%
2. Hit overlap with the adjacent thresholds should be >70% — if not, the optimum is unstable
3. Cell-count distributions should be balanced across genes
4. Control outlier rate should be <5%

---

## Alternatives

### Adaptive per-gene thresholding

Different thresholds for different expression levels — high-expression genes carry more ambient RNA, so they justify a higher cutoff:

```python
def adaptive(expression):
    if expression > 100: return 10
    if expression > 10:  return 5
    return 3
```

Better calibrated, harder to validate. Reserve for screens where you've already done the global sweep and want to push further.

### Probabilistic assignment

Treat assignment as a probability rather than a hard cut:

```python
p_assign = umi_count / (umi_count + background_umi)
# include cells with p > 0.7
```

Retains more information but requires downstream tools that propagate uncertainty — most don't, so the benefit is rarely realized.

---

## References

- Wang et al. (2025) "Integrative analysis of pooled CRISPR screens provide functional insights into AD GWAS risk genes" *bioRxiv*. https://github.com/Alector-BIO/CRISPR_PerturbSeq_pHM_publication
- Replogle JM, et al. (2020) "Combinatorial single-cell CRISPR screens by direct guide RNA capture and targeted sequencing" *Nature Biotechnology* 38(8):954–961
- Storey JD, Tibshirani R (2003) "Statistical significance for genomewide studies" *PNAS* 100(16):9440–9445 — multiple-testing under variable sample sizes
