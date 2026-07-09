---
name: lasso-biomarker-panel
id: "skill_166819affcab4340a91a894872e75315"
display-name: "LASSO Biomarker Panel Selection"
short-description: "Select a minimal (5–15 feature) biomarker panel from high-dimensional omics with nested-CV penalized logistic regression and stability selection; validate on an independent cohort. Use when you nee..."
category: clinical
keywords: [LASSO, biomarker, panel, feature selection, clinical]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Select a minimal (5–15 feature) biomarker panel from high-dimensional omics with nested-CV penalized logistic regression and stability selection; validate on an independent cohort. Use when you need a parsimonious binary classifier with AUC/ROC and feature-stability evidence."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "glmnet (R) / scikit-learn (Python)"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Friedman J et al., J Stat Softw 2010, doi:10.18637/jss.v033.i01; Pedregosa F et al., JMLR 2011"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
  - version: 1.0.0
    date: 2026-05-20
    event: initial-import
    summary: "Initial public import."
status: active
promoted: true
---

# LASSO Biomarker Panel Discovery & Validation

Build a small, interpretable biomarker panel from a high-dimensional omics matrix by combining penalized logistic regression (LASSO / elastic net) with repeated nested cross-validation and stability selection. Produces the artefacts a clinical-exploratory submission needs: ROC curves with discovery and validation, feature-stability bar plot, calibration curve, and a forest plot of LASSO coefficients.

The discovery design follows Ali et al. (*Nat. Med.* 2025) — 50 repeated 70/30 splits, inner 10-fold `cv.glmnet` for lambda, stability threshold of 0.8 by default.

## When to Use This Skill

Use this skill when you need to:

- ✅ Select a minimal biomarker panel (5–15 features) from thousands of candidates
- ✅ Build a predictive model for a binary clinical outcome (responder vs. non-responder, disease vs. control, subtype A vs. B)
- ✅ Pair discovery with an independent validation cohort
- ✅ Generate regulatory-grade outputs: ROC/AUC with 95% CI, calibration, decision curves
- ✅ Provide a locked panel that can be projected onto new patients with `predict_biomarker_panel`

Skip this skill for:

- ❌ Unsupervised clustering → `bulk-omics-clustering`
- ❌ Differential expression of all genes → `bulk-rnaseq-counts-to-de-deseq2`
- ❌ Multi-omics factor discovery → `multi-omics-integration`
- ❌ Continuous outcomes — this skill assumes binary classification (extend to `glmnet(family = "gaussian")` yourself if needed)

## Quick Start

Around five to ten minutes end-to-end on the sepsis demo (479 ICU patients, MARS consortium):

```r
source("scripts/load_example_data.R")
data <- load_sepsis_data(outcome = "endotype")           # or other demos below

source("scripts/prepare_features.R")
features <- prepare_feature_matrix(data$expression, data$metadata,
                                    data$outcome_col, top_n_variable = 2000)

source("scripts/lasso_workflow.R")
model <- run_lasso_panel(features$X, features$y, alpha = 0.5)

source("scripts/biomarker_plots.R")
generate_all_plots(model, X = features$X, y = features$y, output_dir = "results")

source("scripts/export_results.R")
export_all(model, output_dir = "results", data = data, features = features)
```

## Installation

| Package | Version | License | Commercial use | Install command |
|---------|---------|---------|----------------|-----------------|
| glmnet | ≥4.1 | GPL-2 | yes | `install.packages("glmnet")` |
| pROC | ≥1.18 | GPL (≥3) | yes | `install.packages("pROC")` |
| ggplot2 | ≥3.4 | MIT | yes | `install.packages("ggplot2")` |
| ggprism | ≥1.0.3 | GPL (≥3) | yes | `install.packages("ggprism")` |
| ggrepel | ≥0.9 | GPL-3 | yes | `install.packages("ggrepel")` |
| ComplexHeatmap | ≥2.10 | MIT | yes | `BiocManager::install("ComplexHeatmap")` |
| circlize | ≥0.4.15 | MIT | yes | `install.packages("circlize")` |
| GEOquery | ≥2.62 | Artistic-2.0 | yes | `BiocManager::install("GEOquery")` (demos) |
| Biobase | ≥2.54 | Artistic-2.0 | yes | `BiocManager::install("Biobase")` (demos) |
| msigdbr | ≥7.5 | MIT | yes | `install.packages("msigdbr")` (optional, enrichment) |
| clusterProfiler | ≥4.2 | Artistic-2.0 | yes | `BiocManager::install("clusterProfiler")` (optional) |

```r
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages(c("glmnet", "pROC", "ggplot2", "ggprism", "ggrepel", "circlize"))
BiocManager::install(c("ComplexHeatmap", "GEOquery", "Biobase"))
```

## Inputs

**Required.**

- A numeric expression matrix with rows = features (genes/proteins) and columns = samples. Normalized counts, TPM, FPKM, log-protein abundances all fine.
- A sample metadata `data.frame` with a binary outcome column (0/1 or two-level factor). Row names must match the column names of the expression matrix.
- Minimum sample sizes: 20 per group for exploratory work, 40+ recommended.

**Optional.**

- A validation cohort: expression matrix + metadata with the same outcome encoding.
- A pre-filtered feature list (e.g. DE results, WGCNA hub genes) to seed LASSO.

## Outputs

**CSV tables.**

- `biomarker_panel.csv` — final stable features with LASSO coefficients and selection frequencies.
- `all_feature_stability.csv` — every feature ranked by selection frequency.
- `discovery_performance.csv` — per-fold AUC, sensitivity, specificity from the nested CV.
- `validation_performance.csv` — AUC, sensitivity, specificity in the independent cohort (if validation step was run).

**Plots (PNG + SVG at 300 DPI).**

- `roc_curve` — ROC for the locked panel, discovery + validation overlaid when available.
- `stability_barplot` — selection frequency per candidate feature.
- `coefficient_forest` — LASSO coefficient with 95% CI for each panel feature.
- `calibration_curve` — predicted vs. observed event rate.
- `auc_distribution` — distribution of CV-AUC across the 50 outer folds.
- `feature_heatmap` — clustered heatmap of panel-feature expression across samples.

**Analysis object.**

- `lasso_model.rds` — list containing the final model, stable features, feature importance table, fold-wise AUCs, CV predictions, parameters. Load with `model <- readRDS("results/lasso_model.rds")`. Project new samples with `predict_biomarker_panel(model, new_X)`.

**Reports.**

- `summary_report.md` — methods, dataset description, performance table, panel composition.
- `summary_report.pdf` — same content as PDF (requires `tinytex` or a LaTeX installation).

## Clarification Questions

⚠️ **CRITICAL: ASK THIS FIRST.**

1. **Demo or your own data?** The bundled demos: sepsis blood transcriptomics (MARS consortium, ~15–25 gene panel), breast cancer subtypes (GSE25055), IMvigor210 bladder cancer immuno-oncology, or UNIFI UC ustekinumab trial. If your own, what's the expression matrix and the binary outcome column?

2. **Outcome.** What does outcome `= 1` mean (responder/non-responder, alive/dead, subtype A vs. B)? What organism (human / mouse / rat — drives downstream enrichment).

3. **Regularization.** Pure LASSO (`alpha = 1.0`, sparsest) or elastic net (`alpha = 0.5`, retains correlated features)? Default is elastic net for transcriptomic data; switch to pure LASSO for proteomics.

4. **Pre-filter.** Top 2,000 most-variable features (default for >5,000 inputs), top 500 (smaller datasets), or use all features?

5. **Validation cohort.** Will you provide an independent expression matrix + metadata? Strongly recommended — without external validation the reported AUC is discovery-cohort only and expected to be optimistic.

## Standard Workflow

Five scripts, four mandatory steps plus a strongly recommended fifth.

**Step 1 — load data.**

```r
source("scripts/load_example_data.R")
data <- load_sepsis_data(outcome = "endotype")            # 479 patients, MARS endotype
# or: load_sepsis_data(outcome = "mortality")
# or: load_breast_cancer_pcr_data(outcome = "subtype")   # 218 tumours, Basal vs LumA
# or: load_imvigor210_data()                              # 190 bladder cancer, atezolizumab
# or: load_unifi_data()                                   # 542 UC patients, ustekinumab
```

Verification: `okay <dataset> loaded successfully!`.

**Step 2 — prepare features.**

```r
source("scripts/prepare_features.R")
features <- prepare_feature_matrix(data$expression, data$metadata,
                                    data$outcome_col, top_n_variable = 2000)
```

`features$X` is a samples × features matrix on a standardized scale; `features$y` is the 0/1 outcome vector aligned to the rows.

**Step 3 — fit LASSO panel.**

```r
source("scripts/lasso_workflow.R")
model <- run_lasso_panel(features$X, features$y, alpha = 0.5)
```

Verification: `okay LASSO panel selection completed successfully!`.

**Step 4 — plots.**

```r
source("scripts/biomarker_plots.R")
generate_all_plots(model, X = features$X, y = features$y, output_dir = "results")
```

Verification: `okay All biomarker plots generated successfully!`.

**Step 5 — exports.**

```r
source("scripts/export_results.R")
export_all(model, output_dir = "results", data = data, features = features)
```

Verification: `=== Export Complete ===`.

**Step 6 (strongly recommended) — external validation.**

```r
source("scripts/load_example_data.R")
val_data <- load_breast_cancer_validation_data("GSE32646")   # or another demo
source("scripts/prepare_features.R")
val_features <- prepare_validation_features(val_data$expression, features$feature_names)
source("scripts/validate_external.R")
val_result <- validate_panel(model, val_features$X,
                              val_data$metadata[[val_data$outcome_col]],
                              cohort_name = "GSE32646")
export_all(model, validation_result = val_result, output_dir = "results",
           data = data, features = features)
```

Without this step the discovery AUC is an in-cohort estimate and almost certainly optimistic by 0.02–0.05 AUC units.

## When Scripts Fail

1. **Fix and retry** — install the missing package, re-source.
2. **Modify in place** — adjust `n_repeats`, `alpha`, or `stability_threshold` in `lasso_workflow.R`.
3. **Read as a reference** — copy a helper for a one-off analysis.
4. **Write from scratch** — last resort, document why.

## Decision Points

- **`alpha = 1.0` vs `alpha = 0.5`.** Pure LASSO produces sparser panels; elastic net keeps correlated features (useful when multiple co-regulated genes carry the same signal). Try LASSO first; switch to elastic net if the stability step yields fewer than three features.
- **Stability threshold.** 0.8 default. Relax to 0.6 for noisy small-cohort data; raise to 0.9 only for confirmatory publications.
- **`n_repeats`.** 50 is the published default. 100 stabilizes the AUC CI at twice the cost. 20 is fine for screening.
- **Pre-filter.** Variance-based default. DE-based pre-filter introduces circularity if the DE used the outcome being predicted; only use it on an independent slice.
- **Cite the CV AUC** (from `discovery_performance.csv`), not the final-model AUC printed on the ROC plot. The final-model AUC is in-sample.

## Common Issues

| Symptom | Likely cause | Resolution |
|---------|--------------|------------|
| `No shared samples` | Column names of expression don't match metadata row names | `colnames(expression)` must match `rownames(metadata)`; reconcile and re-run |
| `Outcome must be binary` | Multi-level or numeric outcome | Recode to 0/1 or a two-level factor before calling `prepare_feature_matrix` |
| `cv.glmnet convergence warning` | Too few samples for 10-fold CV in fold | Reduce `n_inner_folds` to 5 |
| No features pass stability | Weak signal or small effect | Script auto-relaxes to top 10 by frequency; consider `alpha = 0.5` |
| Low validation AUC | Cross-platform gene loss | Inspect `val_features$n_matched`; if <50% matched, use same-platform validation |
| `SVG export warning` | `svglite` not installed | Ignore — script falls back to base R `svg()` |
| GEOquery download fails | Bioconductor mirror / firewall | Retry; or download manually and pass the file path |
| PDF report missing | tinytex / LaTeX absent | `tinytex::install_tinytex()`; markdown report is generated regardless |
| Final-model AUC > CV-AUC by 0.1+ | Overfitting flag | Verify reported AUC is the **CV** AUC, not the final-model AUC |
| Validation AUC ≪ discovery AUC | Real out-of-sample drop or platform mismatch | Quantify cross-platform overlap; report both with provenance |

## Agent Summary Guidelines

When presenting results, the agent must:

1. Cite the CV AUC from `discovery_performance.csv` — never the final-model AUC on the ROC plot.
2. State explicitly whether external validation was performed. If not, append: "These results are from the discovery cohort only. No external validation was performed. AUC estimates are expected to be optimistic."
3. Not describe panel-gene biology unless `pathway-enrichment` or `functional-enrichment-from-degs` has been run on the panel in this session.
4. Report PDF status — if it failed, say so and link to the markdown report.
5. Reserve the word "validated" for panels that have passed an independent cohort.

## Interpretation Guidelines

- **AUC > 0.8** — strong discrimination, clinically useful as a standalone test.
- **AUC 0.7–0.8** — moderate; useful alongside clinical factors.
- **AUC 0.6–0.7** — weak above-chance, typical of transcriptomic response prediction.
- **Stability ≥ 80%** — robust feature, include in the panel.
- **Positive coefficient** — higher expression raises the probability of `outcome = 1`.
- **Calibration on the diagonal** — well-calibrated; deviation indicates risk over- or under-estimation.

### Which AUC to quote

Two AUCs come out of the workflow:

| Figure | Value | What it means | Cite it? |
|--------|-------|---------------|----------|
| Mean CV AUC (`discovery_performance.csv`) | e.g. 0.986 | Average over 50 outer-fold test sets | yes |
| Final-model AUC (ROC plot annotation) | e.g. 0.996 | Locked panel applied to its own training data | no (in-sample) |

The CV AUC is the right headline. The final-model AUC is in-sample and will always be at least as optimistic.

### Optimism bias on the CV AUC

The reported CV AUC is still a discovery-cohort estimate, not a true out-of-sample number. Two sources of residual bias:

- Feature stability is computed on the same samples used for AUC estimation.
- The locked panel is re-fit on all samples after selection.

Empirically the bias is 0.02–0.05 AUC units for cohorts of this size — a discovery AUC of 0.986 typically lands at ≥0.93 on an independent cohort when the signal is real. The only way to remove the bias is external validation.

## Suggested Next Steps

1. **Pathway enrichment of the panel** → `pathway-enrichment` or `functional-enrichment-from-degs`. (Mandatory before interpreting biology.)
2. **Co-expression context** — which network modules contain panel genes → `coexpression-network`.
3. **Patient stratification** using the panel as input → `multi-omics-integration`.
4. **Literature validation** — are the panel genes known disease or therapy markers?
5. **Independent cohort replication** — apply the locked panel to other datasets.

## Related Skills

| Skill | Relationship |
|-------|--------------|
| `bulk-rnaseq-counts-to-de-deseq2` | Upstream — DE results can pre-filter LASSO candidates |
| `coexpression-network` | Upstream — WGCNA hub genes seed LASSO candidates |
| `pathway-enrichment` | Downstream — enrichment of panel genes |
| `bulk-omics-clustering` | Alternative — unsupervised patient stratification |
| `multi-omics-integration` | Downstream — multi-omics subtype discovery |

## References

- **[glmnet]**: Friedman J, Hastie T, Tibshirani R. "Regularization Paths for Generalized Linear Models via Coordinate Descent." *Journal of Statistical Software* 2010;33(1):1–22. doi:10.18637/jss.v033.i01 — [link](https://doi.org/10.18637/jss.v033.i01)
- **[scikit-learn]**: Pedregosa F, Varoquaux G, Gramfort A, et al. "Scikit-learn: Machine Learning in Python." *Journal of Machine Learning Research* 2011;12:2825–2830. — [link](https://jmlr.org/papers/v12/pedregosa11a.html)
- **[pROC]**: Robin X, Turck N, Hainard A, et al. "pROC: an open-source package for R and S+ to analyze and compare ROC curves." *BMC Bioinformatics* 2011;12:77.
- **[Discovery design / application]**: Ali M, et al. Cross-disease LASSO proteomics panels. *Nature Medicine* 2025.
- **[UNIFI trial data]**: Sands BE, et al. *NEJM* 2019.
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)
- Local references:
  - [`references/lasso-reference.md`](references/lasso-reference.md) — alpha, lambda, stability thresholds.
  - [`references/decision-guide.md`](references/decision-guide.md) — pre-filter and method choices.
  - [`references/validation-guide.md`](references/validation-guide.md) — multi-cohort study design.