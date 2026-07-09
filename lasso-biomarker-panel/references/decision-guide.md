# Decision Guide: which knob, which skill, which fix

## LASSO vs elastic net vs ridge

| Situation | Recommendation | Alpha |
|-----------|----------------|-------|
| Biomarker panel for clinical assay (sparsest) | Pure LASSO | `1.0` |
| Features highly correlated (pathway co-expression) | Elastic net | `0.5` |
| Many weak signals, no dominant features | Elastic net | `0.5` |
| Bulk RNA-seq (some correlation, default starting point) | LASSO; switch to elastic net if <3 stable | `1.0 → 0.5` |
| Proteomics (SOMAscan, Olink) | LASSO — typically lower feature correlation | `1.0` |
| Multi-omics with cross-omics correlation | Elastic net | `0.5` |

## When to use this skill (vs. alternatives)

| Goal | Use this skill? | If not, use |
|------|-----------------|-------------|
| Select a minimal predictive panel | yes | — |
| Predict a binary clinical outcome from omics | yes | — |
| Unsupervised patient clustering | no | `bulk-omics-clustering` |
| Find differentially expressed genes | no | `bulk-rnaseq-counts-to-de-deseq2` |
| Multi-omics factor analysis | no | `multi-omics-integration` |
| Pathway enrichment of selected panel | yes — after panel discovery | `pathway-enrichment` / `functional-enrichment-from-degs` |
| Co-expression context for the panel | either side of LASSO | `coexpression-network` |

## Pre-filtering strategy

### When to pre-filter

- `>10,000` features — pre-filter to top 5,000 by variance (`prepare_features.R` does this by default).
- `>20,000` features — consider DE pre-filtering, but be careful about circularity.
- `<5,000` features — no pre-filter needed.

### Pre-filtering methods, best to worst

1. **DE-based.** Statistically principled, but introduces circularity if the DE used the same outcome you're trying to predict. Only safe when the DE was run on an independent slice.
2. **Variance-based (default).** No circularity. Picks up informative features in practice.
3. **WGCNA hub-gene-based.** Biologically informed — reduces to pathway-representative features.
4. **Random subset.** Don't.

## Troubleshooting decision tree

```
No stable features (all < 80% selection frequency)
├── Try alpha = 0.5 (elastic net)
│   ├── Still none?
│   │   ├── Lower stability threshold to 0.6
│   │   └── If still none: signal probably too weak for this cohort
│   └── Features found → proceed with elastic net panel

Too many stable features (>20)
├── Raise stability threshold to 0.9
├── Set top_n_features = 10 per iteration
└── Switch from elastic net back to LASSO (alpha = 1.0)

Low AUC (<0.6)
├── Check class balance (need ≥30% minority class)
├── Try DE pre-filter (on an independent slice to avoid circularity)
├── Raise n_repeats to 100 for more stable CV estimates
└── Accept: effect size may be too small for this cohort
```

## Class balance and small cohorts

- Aim for ≥30% minority class. Below that, LASSO will lean on majority-class features and AUC inflates.
- The pipeline does class-balanced sampling per outer fold automatically.
- For severely imbalanced data (>90/10), consider SMOTE upstream or `glmnet(..., weights = ...)` with inverse-class-frequency weights.

## When to skip external validation

Never. External validation is what separates an exploratory panel from a deployable test. If you absolutely cannot run external validation, label every reported AUC as "discovery-cohort estimate, expected to be optimistic" and don't use the word "validated".
