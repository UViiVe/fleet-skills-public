# LASSO Parameter Reference

Quick reference for the knobs in `run_lasso_panel`. All defaults follow the Ali et al. (*Nat. Med.* 2025) design unless noted.

## Alpha — the mixing parameter

| Value | Behaviour | When to pick |
|-------|-----------|--------------|
| `1.0` | Pure LASSO (L1) | Sparsest panel. Clinical assay development, proteomics |
| `0.5` | Elastic net (mixed L1+L2) | Default for transcriptomics; retains correlated features |
| `0.0` | Ridge (L2) | No feature selection — not appropriate for panel discovery |

Practical rule: start with `1.0`. If stability selection retains fewer than three features, retry with `0.5` and check whether the AUC improves.

## Lambda — the regularization strength

Selected per-iteration by the inner `cv.glmnet` 10-fold loop.

- `lambda.min` (default) — minimises CV error.
- `lambda.1se` — largest lambda within 1 standard error of the minimum; more conservative, fewer features.

Switch to `lambda.1se` only when the panel is too large to be clinically practical and you can afford to trade a small AUC drop for a smaller panel.

## Stability threshold

| Value | Behaviour | Best for |
|-------|-----------|----------|
| `0.8` (default) | Feature must appear in ≥80% of the 50 outer folds | Standard publication |
| `0.9` | Stringent — only features in ≥90% of folds | Confirmatory work |
| `0.6` | Relaxed | Tiny cohorts where 0.8 yields nothing |

If you tune this knob, document it in the methods section.

## `n_repeats`

| Value | Trade-off |
|-------|-----------|
| `20` | Fast screen, noisy frequency estimates |
| `50` (default) | Stable enough for the published design |
| `100` | Smoother AUC CI at 2× wall-time |

## `train_fraction`

70/30 split is the published default. Larger training fractions inflate the CV AUC; smaller (60/40) shrink it but tighten the variance estimate.

## Interpreting outputs

### Selection frequency

Fraction of folds where the feature received a non-zero coefficient.

- ≥ 0.8 — highly stable, include in panel.
- 0.5–0.8 — moderately stable; consider including with a stability annotation.
- < 0.5 — unstable; likely spurious or redundant.

### Panel coefficients

- Positive coefficient — higher feature value pushes the probability of `outcome = 1` upward.
- Negative coefficient — opposite.
- Magnitude — relative importance. Sign matters more for interpretation; magnitude is harder to compare across features on different scales.
- 95% CI crossing zero — the direction of effect isn't reliable; treat the feature as a candidate, not a finding.

### AUC

| AUC range | Clinical interpretation |
|-----------|-------------------------|
| 0.9–1.0 | Excellent discrimination |
| 0.8–0.9 | Good — clinically useful |
| 0.7–0.8 | Moderate — useful alongside other factors |
| 0.6–0.7 | Weak above chance |
| 0.5–0.6 | Near-random |

The CV AUC carries an optimism bias of about 0.02–0.05 AUC units; external validation is the only way to drop it.

## Why nested CV matters

Plain `cv.glmnet` selects lambda **and** evaluates AUC on the same folds — the AUC is biased upward. The nested design:

1. Outer loop: 50 random 70/30 splits.
2. Inner loop on each training set: 10-fold `cv.glmnet` selects lambda.
3. AUC is computed on the outer test set, which never saw the inner lambda choice.

This is the Ali et al. (2025) design. The cost is modest (~50× CV) but the bias correction is essential for any publishable estimate.

## Cascading from WGCNA hub genes

Shen et al. (*Cell* 2024) pre-filter LASSO inputs to WGCNA hub genes:

```r
hub_genes <- readRDS("coexpression_results/hub_genes.rds")
expr_filtered <- expression[rownames(expression) %in% hub_genes, ]
features <- prepare_feature_matrix(expr_filtered, metadata, outcome_col)
```

This shrinks thousands of genes to dozens of biologically coherent candidates, often raising stability without losing AUC.
