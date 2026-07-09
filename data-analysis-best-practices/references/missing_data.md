# Missing Data — Mechanisms, Diagnostics, and Imputation

Missing values are not all the same. The mechanism determines what you're allowed to do; the diagnostics tell you which mechanism you have.

---

## Three mechanisms

### Missing Completely At Random (MCAR)

The probability of missingness is independent of all variables, observed or unobserved. The missing values are a random subset.

Test: Little's MCAR test. If you can't reject MCAR, listwise deletion is unbiased (just inefficient).

In practice: pure MCAR is rare. Real-world data usually has at least some MAR structure.

### Missing At Random (MAR)

Missingness depends on observed variables but not on the missing values themselves. After conditioning on observed data, the missingness is random.

Example: men decline to report weight more often than women. If you condition on sex, the missingness is unrelated to the actual weight value.

In practice: most realistic case. Multiple imputation, EM-based approaches, and many machine learning imputers work under MAR.

### Missing Not At Random (MNAR)

Missingness depends on the missing value itself, even after conditioning on observed data.

Example: people with high incomes are less likely to report their income. The probability of missingness depends on the income value, not on observable covariates.

In practice: detection is hard (you can't observe what you're missing). Requires explicit modeling of the missingness mechanism, sensitivity analysis, or accepting bias.

---

## Diagnostic workflow

1. **Quantify missingness per column.** % missing, count, pattern.
2. **Visualize missingness.** `missingno` matrix and heatmap — patterns suggest mechanism.
3. **Test for MAR proxies.** Add a `column_x_missing` indicator and regress against observed covariates. If the indicator is predicted by observed data, you have MAR (not MCAR).
4. **Examine missingness co-occurrence.** When column A is missing, is column B also missing? Joint patterns suggest structural missingness (e.g., entire visits skipped).
5. **Tell the user.** Surface findings before deciding on a method.

```python
from scripts.missing_data import diagnose_missingness
report = diagnose_missingness(df)
```

---

## Decision rules

| Missingness | Default action |
|---|---|
| 0% | nothing |
| < 5% | impute with median (numeric) or mode (categorical), document |
| 5–20% | impute with care; consider multiple imputation; add missingness indicator |
| 20–50% | drop the column unless it's load-bearing; if kept, use multiple imputation |
| > 50% | drop the column |

Caveats:
- Lower thresholds for outcomes than for predictors.
- If multiple imputation, generate ≥ 5 imputed datasets and pool with Rubin's rules.
- For survey data with structured non-response (e.g., questions only asked of subgroup), recode rather than impute.

---

## Imputation methods

### Single imputation

- **Mean / median:** simple, biased, underestimates variance. Acceptable for < 5% missingness on non-critical columns.
- **Mode (categorical):** equivalent to median for unordered factors.
- **KNN imputation:** matches on k nearest neighbors. Reasonable for moderate missingness in continuous data.
- **Iterative imputer (MICE-style, scikit-learn):** chains regressors. Works under MAR.
- **Random sampling from the observed distribution:** preserves marginal distribution; ignores correlations.

### Multiple imputation

Generate m imputed datasets (m ≥ 5, often 20–50). Run the analysis on each. Pool using Rubin's rules. Captures imputation uncertainty.

Tools: `mice` (R), `IterativeImputer` (sklearn — single by default but supports multiple by changing the random state).

```python
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer

imp = IterativeImputer(random_state=42, max_iter=10)
df_imputed = imp.fit_transform(df)
```

### Indicator columns

Add a binary column `was_missing_x` for each imputed column. Lets downstream models account for whether a value was originally observed.

### Carry-forward / interpolation

For time series only. Last-observation-carried-forward (LOCF), interpolation, or model-based smoothing.

---

## When NOT to impute

- The mechanism is MNAR — imputation under MAR gives biased estimates.
- The column is a treatment or outcome variable in a causal analysis.
- Missingness itself is informative (e.g., a clinical lab not ordered may mean the patient was too sick or too healthy).

In these cases:
- Use complete-case analysis with explicit discussion of bias.
- Model missingness explicitly (pattern-mixture models, selection models).
- Run sensitivity analyses across plausible MNAR scenarios.

---

## Pitfalls

- **Imputing the outcome.** Almost always wrong unless you have a very strong reason.
- **Imputing then using cross-validation.** Imputation must happen inside each CV fold, not before, to avoid data leakage.
- **Single imputation with downstream confidence intervals.** The CIs ignore imputation uncertainty and will be too narrow. Use multiple imputation or bootstrap.
- **Treating zeros as missing or vice versa.** Always check the raw file for the distinction between an explicit `NA` and a structural zero.
- **Dropping rows in a sequential pipeline.** Each drop changes the effective sample size. Track at every step.

---

## Documentation requirements

For every imputation decision:

- Method used (median, KNN, MICE, etc.) with parameters.
- Number of imputations if multiple.
- Random seed.
- Per-column count of imputed values.
- Sensitivity analysis if MNAR is suspected.

Capture in the cleaning log:

```
column: serum_creatinine
original_missing: 28 (8.4%)
mechanism: assumed MAR after conditioning on age, sex, eGFR
method: IterativeImputer (sklearn 1.3.0), max_iter=10, random_state=42
indicator_added: yes (was_missing_serum_creatinine)
```

---

## References

- Little RJA, Rubin DB. *Statistical Analysis with Missing Data*, 3rd ed., Wiley, 2019.
- van Buuren S. *Flexible Imputation of Missing Data*, 2nd ed., CRC Press, 2018.
- Sterne JAC et al. *Multiple imputation for missing data in epidemiological and clinical research.* BMJ 2009;338:b2393.
