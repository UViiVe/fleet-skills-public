# Multiple-Testing Correction Guide

How to choose and apply a correction method when testing thousands of genes simultaneously.

---

## Why correction is mandatory

Test 20 000 genes at α = 0.05 with no correction: even if no genes are truly differentially expressed, you expect 1 000 false positives. The result is almost entirely noise.

Scale at which this matters:

| Assay | Tests |
|---|---|
| Bulk RNA-seq | ~15 K – 20 K genes |
| ATAC-seq | 50 K – 150 K peaks |
| Methylation EPIC array | ~850 K CpGs |
| WGBS | millions of CpGs |
| GWAS | 0.5 M – 10 M SNPs |

Family-wise error rate (FWER) at uncorrected α = 0.05 with 100 independent tests is already ~99%. Correction is not optional.

---

## Methods

### Bonferroni

Use α′ = α / m where m is the number of tests. Controls FWER.

Pros: dead simple, no assumptions about test dependence, well understood.
Cons: extremely conservative; massive false-negative rate at genomics scale.

Use when: m < 100 (candidate panels), GWAS (m × 5 × 10⁻⁸ threshold), or in validation experiments where false positives are very costly.

```r
adjusted <- p.adjust(pvalues, method = "bonferroni")
```

### Benjamini-Hochberg FDR (BH)

Controls the expected proportion of false positives among rejections (FDR), not the probability of any false positive.

Algorithm: sort p-values, reject any p_(k) ≤ (k/m) × α.

Pros: dramatically more powerful than Bonferroni at genomics scale; standard in the literature.
Cons: rejections are a *set* — you cannot point at any one and say "this is 95% likely real".

Use as default for RNA-seq, ATAC-seq, ChIP-seq, methylation, anywhere with > 100 tests.

```r
adjusted <- p.adjust(pvalues, method = "BH")
```

### Storey's q-value

BH augmented with an estimate of π₀ — the proportion of truly null hypotheses.

Pros: more powerful than BH when many true effects exist (π₀ < 0.9).
Cons: π₀ estimation can be unstable with small n; needs well-calibrated p-values.

Use when: many true positives expected (treatments with broad effects, big tissue comparisons).

```r
library(qvalue)
qobj <- qvalue(pvalues)
qobj$qvalues
```

### Independent Hypothesis Weighting (IHW)

Data-driven weighting of hypotheses by an independent covariate (e.g., mean expression). Can gain 10–50% more discoveries than BH while preserving FDR control.

Pros: strictly dominates BH if you have an informative covariate.
Cons: needs a covariate that's independent of the test statistic under the null. Mean expression is good; fold-change is not (correlated with the test statistic).

Use when: you have an informative covariate and want maximum power.

```r
library(IHW)
ihw_res <- ihw(pvalues, covariates = mean_expression, alpha = 0.05)
adj_pvalues(ihw_res)
```

### Local FDR

Estimates the probability that *this specific test* is a false positive.

Pros: gene-specific error estimate; better for prioritizing follow-up.
Cons: harder to interpret; needs good p-value density estimation.

Use when: prioritizing individual hits for validation, comfortable with Bayesian framing.

---

## Method comparison

| Method | Controls | Power | Complexity | Default for |
|---|---|---|---|---|
| Bonferroni | FWER | Very low at scale | Trivial | < 100 tests, validation |
| BH-FDR | FDR | Moderate | Trivial | Genomics standard |
| q-value | FDR | Moderate-high | Low | Many true positives expected |
| IHW | FDR | High | Low | Want max power with covariate |
| Local FDR | Local | Moderate-high | Medium | Per-gene prioritization |

---

## Decision tree

```
m < 100 → Bonferroni
m ≥ 100, GWAS → Bonferroni at 5×10⁻⁸
m ≥ 100, otherwise:
    Have informative covariate (e.g., mean expression)? → IHW
    Expect > 10% of tests to be true positives? → q-value
    Otherwise → BH-FDR (the default)
```

When in doubt, use BH-FDR at q < 0.05. It's the field default, well-understood by reviewers, and rarely the wrong choice.

---

## Significance thresholds

| q threshold | Use when |
|---|---|
| q < 0.05 | Standard discovery — recommended default |
| q < 0.10 | Pilot, exploratory, hypothesis-generating |
| q < 0.01 | Confirmatory, expensive follow-up, clinical |

GWAS uses Bonferroni at p < 5 × 10⁻⁸ as a community convention.

---

## Tool-specific implementations

### DESeq2

Default uses BH-FDR; the `padj` column is already adjusted. Don't re-correct.

```r
res <- results(dds, alpha = 0.05)          # BH-FDR
res_ihw <- results(dds, filterFun = ihw)   # IHW
res_shrunk <- lfcShrink(dds, coef = "condition_treated_vs_control", type = "apeglm")
```

`lfcShrink` shrinks log fold-changes toward zero, improving power for genes with high estimation variance.

### edgeR

```r
fit <- glmQLFit(dge, design)
qlf <- glmQLFTest(fit, coef = 2)
results <- topTags(qlf, n = Inf, adjust.method = "BH")
```

### limma

```r
fit <- lmFit(eset, design)
fit2 <- eBayes(fit)
results <- topTable(fit2, coef = 2, adjust.method = "BH", n = Inf)
```

### Manual

```r
p.adjust(pvalues, method = "BH")
p.adjust(pvalues, method = "bonferroni")
```

---

## Interpreting q < 0.05

Correct reading: "Among genes I'm calling significant, ~5% are expected to be false positives."

Common misinterpretations:

- ❌ "Gene X has a 5% chance of being a false positive." (That's local FDR, not FDR.)
- ❌ "95% confidence the effect is real." (Not a confidence statement.)
- ❌ "p < 0.05 after correction." (q-values aren't p-values.)

When reporting:

> Differential expression was assessed with DESeq2 (v1.40); p-values were adjusted with the Benjamini-Hochberg procedure and genes with adjusted p < 0.05 were considered significant.

Or with IHW:

> Multiple testing was controlled using Independent Hypothesis Weighting (IHW v1.28) with mean expression as the covariate, yielding 1 512 significant genes at FDR < 0.05 — a 21% increase over standard BH-FDR (1 247 genes).

---

## Special situations

### Nothing significant

Common causes, in order:

1. Underpowered study (most likely).
2. No real effect.
3. High variability or batch effects.
4. Wrong contrast or wrong timepoint.

What not to do: relax the threshold post-hoc, switch to uncorrected p-values, report "trends".

What to do: examine top-ranked genes for biological coherence, check QC and batch effects, consider whether the effect size assumption was realistic, plan a better follow-up.

### Too many significant (> 50%)

Possible causes:

- Real broad effect (e.g., gross tissue type comparison, strong drug).
- Unresolved batch effect inflating test statistics.
- Library size normalization issue.

Diagnose with PCA, sample-sample correlation heatmaps, before/after batch correction comparison.

### Stratified testing (multiple cell types or tissues)

Two options:

- **Within-stratum correction:** apply BH-FDR within each cell type separately. Treats each as an independent question.
- **Across-stratum correction:** pool all tests, apply one BH-FDR. More conservative; controls overall FDR.

Choose based on the inferential question. Per-cell-type biology is often best handled within-stratum, with a global-FDR-controlled summary table on top.

### Pre-filtering to reduce m

Filtering tests based on a covariate independent of the test statistic under the null is valid and can boost power. DESeq2 does this automatically (independent filtering on mean expression). The filter must be independent — never filter on fold-change or p-value.

---

## Discovery → validation

Validation experiments test fewer hypotheses → use FWER (Bonferroni), not FDR.

Typical workflow:

1. Discovery RNA-seq: 20 K genes, BH-FDR at 0.05 → 1 000 hits.
2. Select top 20 for qPCR validation.
3. qPCR validation: 20 tests, Bonferroni at α = 0.05 / 20 = 0.0025.

The validation set asks a different question (is *this specific gene* real?), which is what FWER controls.

---

## Common mistakes

- **Double-correcting:** running BH inside DESeq2's `results()` and then re-running BH on the output. Use `padj` directly.
- **Choosing threshold after seeing results:** "q < 0.05 doesn't work, let's try 0.10." This is p-hacking.
- **Reporting "trends":** non-significant means non-significant; "trend toward significance" is not statistical language.
- **Mixing methods inconsistently:** FDR on the discovery RNA-seq but raw p < 0.05 on validation qPCR. Use the right correction for each scale.

---

## Citations

- Benjamini Y, Hochberg Y. *Controlling the false discovery rate: a practical and powerful approach to multiple testing.* J R Stat Soc B 1995;57:289–300.
- Storey JD, Tibshirani R. *Statistical significance for genomewide studies.* PNAS 2003;100:9440–9445.
- Ignatiadis N, Klaus B, Zaugg JB, Huber W. *Data-driven hypothesis weighting increases detection power in genome-scale multiple testing.* Nat Methods 2016;13:577–580.
- Efron B. *Local false discovery rates.* Stanford Tech Report, 2005.

---

## Cross-references

- Power calculations including correction burden — [power_analysis_guidelines.md](power_analysis_guidelines.md).
- Design overview — [experimental_design_best_practices.md](experimental_design_best_practices.md).
- Script — `scripts/multiple_testing.R`.
