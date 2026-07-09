# Multiple Testing — Pre-Specified Methods

When you compute many p-values, an uncorrected 0.05 threshold gives mostly noise. This document is the short-form decision guide; the deeper treatment lives in [experimental-design-statistics](../../experimental-design-statistics/references/multiple_testing_guide.md).

---

## The numbers

| Number of tests at α = 0.05 | Expected false positives if no effects |
|---|---|
| 20 | 1 |
| 100 | 5 |
| 1 000 | 50 |
| 20 000 (RNA-seq) | 1 000 |
| 850 000 (EPIC methylation) | 42 500 |
| 10 000 000 (GWAS) | 500 000 |

Any analysis with > 20 tests needs correction.

---

## Methods at a glance

| Method | Controls | Power | Best for |
|---|---|---|---|
| Bonferroni | FWER | Low | < 100 tests, validation |
| BH-FDR | FDR | Moderate | Default for genomics |
| q-value (Storey) | FDR | Moderate-high | Many true positives expected |
| IHW | FDR | High | Have informative covariate |
| Local FDR | local | Moderate-high | Per-test prioritization |

---

## Decision tree

```
m < 100 tests          → Bonferroni
m ≥ 100, GWAS-scale    → Bonferroni at 5e-8
m ≥ 100, otherwise:
    Have informative covariate? → IHW
    Expect many true positives? → q-value
    Otherwise                   → BH-FDR  (default)
```

When in doubt, BH-FDR at q < 0.05.

---

## Pre-specification rule

Choose the method and threshold **before** seeing the result count. Switching methods after the fact ("we tried BH but got 7 hits, so we used q < 0.10") is method-shopping and inflates the false-discovery rate. Reviewers spot it.

---

## Implementation

### Python

```python
from statsmodels.stats.multitest import multipletests
reject, qvals, _, _ = multipletests(pvalues, alpha=0.05, method="fdr_bh")
```

Methods supported by `multipletests`:

- `bonferroni`
- `fdr_bh` (Benjamini-Hochberg)
- `fdr_by` (Benjamini-Yekutieli, for dependent tests)
- `holm`, `holm-sidak`, `hommel`

### R

```r
adjusted <- p.adjust(pvalues, method = "BH")
adjusted <- p.adjust(pvalues, method = "bonferroni")

library(qvalue)
qobj <- qvalue(pvalues)
qobj$qvalues

library(IHW)
ihw_res <- ihw(pvalues, covariates = mean_expression, alpha = 0.05)
adj_pvalues(ihw_res)
```

### DESeq2 / edgeR / limma

Their `padj` / `FDR` / `adj.P.Val` columns are already corrected. Don't re-correct.

---

## Stratified vs pooled correction

If tests come from multiple distinct strata (e.g., DE in 10 cell types), two choices:

- **Within-stratum correction:** correct each stratum separately. Treats each as an independent question. Good for per-cell-type biology.
- **Pooled correction:** combine all tests, one correction. Controls overall FDR. More conservative.

Choose based on the inferential framing. Within-stratum is usually preferred when downstream interpretation is per-stratum.

---

## Common mistakes

- **No correction at all.** Almost all hits in any genomic-scale analysis are noise without it.
- **Correcting twice.** DESeq2's `padj` is already corrected. Don't apply BH on top.
- **Mixing scales inconsistently.** RNA-seq discovery uses FDR; downstream qPCR validation of 5 genes uses Bonferroni.
- **Treating q-value like a p-value.** "q < 0.05" doesn't mean "this gene has < 5% chance of being a false positive" — that's local FDR.

---

## Reporting

For methods sections:

> P-values were adjusted using the Benjamini-Hochberg procedure; genes with adjusted p (q-value) < 0.05 were considered significant.

For IHW:

> Multiple testing was controlled with Independent Hypothesis Weighting using mean expression as the covariate, alpha = 0.05.

For GWAS:

> Genome-wide significance was defined as p < 5e-8 (Bonferroni).

---

## References

- Benjamini Y, Hochberg Y. *Controlling the false discovery rate.* J R Stat Soc B 1995;57:289–300.
- Storey JD. *A direct approach to false discovery rates.* J R Stat Soc B 2002;64:479–498.
- Ignatiadis N et al. *Data-driven hypothesis weighting increases detection power.* Nat Methods 2016;13:577–580.
- Efron B. *Local false discovery rates.* Stanford Tech Report, 2005.
