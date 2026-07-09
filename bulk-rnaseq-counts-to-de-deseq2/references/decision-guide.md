# DESeq2 Decision Guide

Five decision points that have to be made before you can call yourself done. None of them have a universally correct default; all of them are explainable in three sentences.

---

## Decision 1: Transformation (vst vs rlog)

**You need this for** PCA, heatmaps, clustering, WGCNA, anything that wants homoscedastic values.

| Method | Use when | Pros | Cons |
|--------|----------|------|------|
| `vst(dds, blind = FALSE)` | n > 30, large studies | Fast on 1000+ samples | Slightly less accurate at very low counts and very small n |
| `rlog(dds, blind = FALSE)` | n < 30, especially small replicate counts | Better behaved at low expression | Slow (minutes to hours) once n > 100 |

```
n samples?
├── > 30 ────────────→ vst()
├── 10–30 ──────────→ either; vst() if you want speed, rlog() for accuracy
└── < 10 ──────────→ rlog()
```

`blind` controls whether the dispersion estimate used by the transformation knows about the design. Use `blind = TRUE` only for fully exploratory PCA where you don't want the labels to influence the variance estimate; for downstream clustering use `blind = FALSE`.

---

## Decision 2: LFC shrinkage method

**You need this for** plotting, ranking, GSEA — anywhere you'd be embarrassed to show a |log2FC| > 8 from a gene with two reads in one sample.

### apeglm (default)

```r
res_shrnk <- lfcShrink(dds, coef = resultsNames(dds)[2], type = "apeglm")
```

Heavy-tailed prior. Best preservation of true large effects, cleanest top-N lists.

Constraint: needs `coef =`. You cannot do `lfcShrink(..., contrast = c(factor, num, denom), type = "apeglm")`. If you must use a contrast (e.g., treatmentA vs treatmentB when neither is the reference), either re-level and refit, or use ashr.

### ashr

```r
res_shrnk <- lfcShrink(dds, contrast = c("condition", "A", "B"), type = "ashr")
```

Adaptive mixture-of-normals prior. Slightly more shrinkage of large effects than apeglm, but works with `contrast =`, is fast, and is the right choice when you have many pairwise comparisons in a multi-level factor.

### normal (legacy)

The original DESeq2 prior. Inferior to both above; only here for reproducing older analyses.

```
Need contrast= form?
├── YES → ashr
└── NO ──→ apeglm
```

### Shrunk vs unshrunk in your output table

Report both. Use unshrunk LFC + SE for "the treated samples had a log2FC of 1.4 (SE 0.3, padj 1.2e-5)". Use shrunk LFC for plots and ranked lists. They are not redundant — they serve different purposes.

---

## Decision 3: Design formula

**You need this for** every single analysis. The simplest possible design that captures the structure you care about beats a complicated design that uses degrees of freedom on noise.

### Workflow

1. Build `dds` with `design = ~ condition`.
2. Run `vsd <- vst(dds, blind = TRUE); plotPCA(vsd, intgroup = c("condition","batch"))`.
3. Look at the plot.

What you see decides what to do next:

| Observation | Design |
|-------------|--------|
| Samples cluster by condition cleanly | `~ condition` — keep it simple |
| Samples cluster by batch *and* condition | `~ batch + condition` — adjust |
| Samples cluster by batch, condition mixed | Check confounding with `table(batch, condition)`. If confounded you cannot adjust |
| Paired structure (same subject pre/post) | `~ subject + condition` |
| Testing if treatment effect varies by genotype | `~ genotype + treatment + genotype:treatment` |
| Multi-level factor and you want average | `~ condition` with LRT vs `~ 1` |

### Confounding check

```r
table(coldata$batch, coldata$condition)
#         control treated
# batch1        4       0
# batch2        0       4
```
A zero in any cell means you cannot adjust for batch — it perfectly predicts condition. Options: collect more samples to balance, or drop the batch term and report the limitation honestly.

### Decision tree

```
Known batch / covariate?
├── YES → balanced (every cell of batch×condition > 0)?
│   ├── YES → ~ batch + condition
│   └── NO  → confounded. Either drop, or merge into single composite factor ~ 0 + group
└── NO ───→ paired samples?
    ├── YES → ~ subject + condition
    └── NO ─→ testing interaction?
        ├── YES → ~ A + B + A:B
        └── NO ─→ ~ condition
```

Always sanity-check: `Matrix::rankMatrix(model.matrix(design(dds), colData(dds))) == ncol(model.matrix(design(dds), colData(dds)))` should be `TRUE`. If false you're confounded.

---

## Decision 4: Pre-filtering rule

**You need this for** speed and memory; it has minor impact on FDR.

| Rule | Reason |
|------|--------|
| `rowSums(counts(dds)) >= 10` | Default. Removes genes with essentially no signal. |
| `rowSums(counts(dds) >= 10) >= n_per_group` | Stricter — gene must be expressed in at least n samples. Useful when you suspect technical zeros dominate. |
| `rowMeans(counts(dds)) >= 5` | Filter by average expression. Roughly equivalent to default on balanced designs. |

For most analyses the default `>= 10` total is fine. DESeq2's internal *independent filtering* downstream prunes low-mean genes from BH adjustment regardless of what you pre-filter, so the visible impact on padj is small.

---

## Decision 5: Significance thresholds

| padj cutoff | LFC cutoff | When |
|-------------|-----------|------|
| < 0.05 | \|log2FC\| ≥ 1 | **Standard.** Solid headline list. |
| < 0.1 | none | Discovery screen, follow-up assay planning |
| < 0.01 | \|log2FC\| ≥ 2 | High-stringency — small validation panel, downstream qPCR/Western |
| < 0.05 | \|log2FC\| ≥ 0.5 | In-vivo tissue where biological effects are real but modest |

LFC threshold is study-dependent: cell line knockouts produce dramatic LFCs (2× to 10×), patient tissue often shows real biology at 1.3-1.5× because of heterogeneity dilution. Don't carry over a cell-line LFC threshold to a tissue analysis without thinking.

Report both the padj and LFC cutoff used in your methods section. A "significant" gene with no LFC threshold means something very different from one with |LFC| ≥ 1.

See also: [bulk-rnaseq-differential-expression](../../bulk-rnaseq-differential-expression/SKILL.md) for the methodology-level discussion of FDR control, why padj ≠ raw p, and common interpretation errors.

---

## Quick reference card

| Decision | Default | Override when |
|----------|---------|---------------|
| Transformation | `vst()` if n > 30 else `rlog()` | n ~ 30 — either is fine |
| Shrinkage | `apeglm` | Need contrast=, use `ashr` |
| Design | `~ condition` | PCA shows batch separation |
| Pre-filter | `rowSums >= 10` | Want a stricter expression-in-N-samples rule |
| padj cutoff | 0.05 | Screening → 0.1; small validation → 0.01 |
| LFC cutoff | 1 (2×) | Tissue → 0.5; cell line + KO → 2 |

---

## Things that look like decisions but aren't

- **"Should I normalize before DESeq2?"** No. DESeq2 does median-of-ratios internally. Passing in CPM/TPM is wrong.
- **"Should I remove the batch effect from counts before DESeq2?"** No. Include batch in the design. Removing it in the count matrix corrupts the dispersion estimate.
- **"Should I log-transform first?"** No, raw integer counts are required.
- **"Should I round my Salmon outputs?"** No, use `DESeqDataSetFromTximport`, which handles the offset properly.
