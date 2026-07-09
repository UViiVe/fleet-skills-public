# Batch Effects — Detection and Adjustment at Analysis Time

This is the analysis-side companion to the design-side document in [experimental-design-statistics](../../experimental-design-statistics/references/batch_effect_mitigation.md). At analysis time you can no longer redesign — but you can detect, quantify, and (carefully) correct.

---

## What counts as a batch

In tabular data: any grouping variable that reflects a non-biological processing or collection unit. Common columns:

- `batch`, `run_id`, `library_prep_batch`, `sequencing_run`.
- `collection_date`, `processing_date`.
- `site`, `center`, `institution`.
- `operator`, `technician`.
- `instrument_id`, `plate`, `flowcell`.

If the data has none of these but was collected over weeks or months, the missingness of batch metadata is itself a quality problem.

---

## Detection workflow

### Step 1 — Confounding test

For each batch-like column, test independence with the outcome variable.

```python
from scripts.batch_check import test_confounding
result = test_confounding(df, batch_col="batch", outcome_col="condition")
```

A chi-square p < 0.05 between batch and outcome means batch is confounded with outcome — there is *no* purely statistical fix.

### Step 2 — Variance partition

Quantify how much variance in the outcome is explained by batch versus biology.

```python
# Example via statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
model_full = smf.ols("outcome ~ batch + condition", data=df).fit()
model_bio  = smf.ols("outcome ~ condition",         data=df).fit()
batch_r2 = model_full.rsquared - model_bio.rsquared
```

| Variance from batch | Interpretation |
|---|---|
| < 5% | minor; can be ignored or modeled |
| 5–15% | adjust with batch term |
| 15–30% | strong; use ComBat or surrogate variables |
| > 30% | concerning; some samples may need exclusion |

### Step 3 — Visualize

PCA colored by batch versus colored by condition. The "ideal" picture is samples clustering by biology, not by batch. Heatmaps and per-sample boxplots also reveal batch effects.

---

## Adjustment methods

### Linear model with batch covariate

The simplest, most defensible approach. Include batch as a fixed effect in the model.

```python
model = smf.ols("outcome ~ batch + condition", data=df).fit()
```

R / DESeq2: `design = ~ batch + condition`.

Pros: transparent, well-understood. Works when batch effects are linear and additive.
Cons: limited if effects are non-linear or interact with condition.

### ComBat (sva package, R)

Empirical Bayes adjustment that removes batch shifts in mean and variance.

```r
library(sva)
adjusted <- ComBat(dat = log_counts, batch = meta$batch,
                   mod = model.matrix(~ condition, data = meta))
```

Pros: strong batch removal; widely used.
Cons: can over-correct and remove real biology. Don't use on count data directly — use on log-transformed normalized values.

### RUVSeq

Removes Unwanted Variation using either negative-control genes or replicates.

Pros: data-driven; doesn't require batch labels.
Cons: needs control specification; effect of `k` (number of factors) is meaningful and must be tuned.

### SVA

Discovers surrogate variables (hidden batch-like factors). Add them to the model.

```r
library(sva)
svobj <- sva(dat, mod, mod0)
# include svobj$sv in downstream design matrix
```

Pros: discovers unobserved batch structure.
Cons: aggressive; can absorb biological variation.

### Mixed-effects model with batch as random effect

```r
library(lme4)
lmer(outcome ~ condition + (1 | batch), data = df)
```

Use when batches are a random sample of a population of batches (e.g., multi-site studies). Don't use when batches are the entire population (typical case).

---

## When statistical correction will not work

- **Complete confounding** (batch = condition): unrecoverable.
- **Batch-specific biology**: different protocols, different ages of samples — these are not corrected by any of the above methods.
- **Very strong effects (> 50% variance)**: even ComBat struggles.

In these cases:
- Exclude problematic batches if you can afford the n loss.
- Acknowledge as a major limitation in the manuscript.
- Repeat with proper design (the cheapest long-term fix).

---

## Pitfalls

- **Correcting before exploring.** Always inspect raw data first — sometimes the "batch effect" is something biologically interesting.
- **Adjusting on the test statistic.** Use covariates that are independent of the outcome under the null.
- **Forgetting batch in some downstream models.** If you adjust at one stage, you have to keep adjusting consistently.
- **Over-correcting.** ComBat can erase legitimate biology if `mod` doesn't include condition. Always provide the `mod` argument.

---

## Reporting

> Batch effects were detected via PCA visualization and quantified by variance-partition analysis (batch contribution: 12.4%). Differential expression was performed with DESeq2 using the model formula `~ batch + condition`. Sensitivity analyses with ComBat-corrected counts gave concordant results (Supplementary Figure S5).

---

## References

- Leek JT et al. *Tackling the widespread and critical impact of batch effects.* Nat Rev Genet 2010;11:733–739.
- Johnson WE, Li C, Rabinovic A. *Adjusting batch effects in microarray expression data using empirical Bayes methods.* Biostatistics 2007;8:118–127.
- Risso D, Ngai J, Speed TP, Dudoit S. *Normalization of RNA-seq data using factor analysis of control genes or samples.* Nat Biotechnol 2014;32:896–902.
- Goh WWB, Wang W, Wong L. *Why batch effects matter in omics data.* Trends Biotechnol 2017;35:498–507.
