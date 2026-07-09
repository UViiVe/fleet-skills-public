# Multi-Cohort Validation Guide

A panel only earns the label "validated" after it has been tested in an independent cohort. This note describes the validation hierarchy, the cross-platform pitfalls, and the bundled cross-cohort examples.

## Validation hierarchy (strongest to weakest)

1. **Independent cohort, different platform** — e.g. RNA-seq discovery validated on microarray. The most demanding test.
2. **Independent cohort, same platform** — different study site, same technology.
3. **Temporal split** — train on earlier samples, validate on later samples.
4. **Nested cross-validation only** — the internal floor, not external validation.

## Cross-platform validation

When discovery and validation use different platforms (RNA-seq vs. Affymetrix vs. Olink, …):

### Gene symbol mapping

1. Map platform-specific IDs to a common gene symbol set (HGNC for human).
2. For multi-probe genes, keep the probe with the highest variance.
3. Intersect discovery and validation gene-symbol sets.
4. Only the matched features are usable for validation.

### Expected feature loss

| Source → Target | Typical overlap |
|-----------------|-----------------|
| RNA-seq → Affymetrix | 70–85% |
| RNA-seq → Olink proteomics | 0–5% (different molecular types) |
| Same platform | 95–100% |

### Handling missing panel features

If some panel features are absent in the validation matrix:

- Report how many panel features survived mapping.
- The locked model still produces predictions — `predict_biomarker_panel` skips missing features.
- AUC drops because information is lost; report the cross-platform AUC alongside the discovery AUC.
- If >50% of features are missing, the validation result should be interpreted with caution and labelled as such.

## Sample-size targets

| Design | Discovery | Validation | Total |
|--------|-----------|------------|-------|
| Exploratory | ≥50 | ≥30 | ≥80 |
| Confirmatory | ≥100 | ≥50 | ≥150 |
| Regulatory-grade | ≥200 | ≥100 | ≥300 |

## Class balance

- Aim for ≥30% minority class in both cohorts.
- The pipeline uses class-balanced sampling per outer fold automatically.
- Severely imbalanced cohorts (>90/10) may need SMOTE upstream or weighted LASSO (`glmnet(..., weights = ...)`).

## Bundled cross-cohort examples

### Breast cancer discovery → GSE32646 validation

Same platform, different cohort, related-but-distinct endpoint.

| Property | Discovery (GSE25055) | Validation (GSE32646) |
|----------|----------------------|------------------------|
| Disease | Breast cancer | Breast cancer |
| Platform | Affymetrix U133A | Affymetrix |
| Samples | ~218 | ~115 |
| Outcome | Basal vs. Luminal A | Pathologic complete response |
| Treatment | Mixed chemo | Paclitaxel + FEC |

Useful because the same platform keeps gene-overlap high, the cohorts are independent, and the related-but-distinct endpoint tests the generalizability of the molecular signal.

### UNIFI → PURSUIT UC validation

A more demanding test: different drug, same disease.

| Property | Discovery (GSE206285) | Validation (GSE92415) |
|----------|------------------------|-----------------------|
| Disease | Ulcerative colitis | Ulcerative colitis |
| Platform | Affymetrix HT HG-U133+ PM | Same family |
| Samples | ~542 | ~87 |
| Outcome | Week-8 mucosal healing | Week-6 mucosal healing |
| Drug | Ustekinumab (anti-IL-12/23) | Golimumab (anti-TNF) |

A cross-drug validation tests whether the response biology is drug-independent. If a panel discovered on ustekinumab-treated patients also predicts response to golimumab, it's a stronger signal than a same-drug validation.

## Multi-cohort framework

Cruchaga, Western, et al. (2024 preprint) propose a six-cohort, 3,107-sample design:

1. **Discovery** — largest, most homogeneous. Train the locked panel here.
2. **Internal validation** — same-study subset. Nested CV provides this in our pipeline.
3. **External validation 1** — same disease, different study.
4. **External validation 2** — related disease (cross-disease).
5. **External validation 3** — different platform (cross-platform).
6. **External validation 4** — different demographic / geographic population.

Each layer adds a different generalizability claim. A panel that passes all six is the closest thing to "ready for prospective trial design".

## Reporting checklist

Before publishing a "validated panel" claim:

- [ ] CV AUC from nested cross-validation (with 95% CI from the fold distribution).
- [ ] Validation AUC on at least one independent cohort.
- [ ] Number of panel features surviving cross-platform mapping.
- [ ] Per-feature direction of effect consistent between cohorts (sign of LASSO coefficient).
- [ ] Calibration plot for both cohorts.
- [ ] Decision-curve analysis if the panel is intended for clinical use.
- [ ] Statement of optimism bias if external validation is missing.
