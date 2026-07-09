# Data Validation Checklist

A complete pre-analysis checklist. Walk through every item before running the first model.

---

## File-level

- [ ] File opens with the expected reader (CSV → pandas, parquet → pyarrow, RDS → readRDS).
- [ ] Character encoding matches expectations (`utf-8` unless declared otherwise).
- [ ] Decimal separator and thousand separator behave as expected (US `.`/`,` vs EU `,`/`.`).
- [ ] Date columns parse with the right format; no silent month/day swaps.
- [ ] No surprise type coercion on read (numeric → object due to one bad cell).

---

## Shape

- [ ] Number of rows matches expectation (or the user has confirmed the count).
- [ ] Number of columns matches expectation.
- [ ] Column names are unique (no `x.1`, `x.2` artifacts).
- [ ] No completely empty rows or columns (unless intentional padding).

---

## IDs

- [ ] ID column declared and present.
- [ ] IDs are unique (no duplicates).
- [ ] IDs have no leading/trailing whitespace.
- [ ] IDs follow a consistent format (case, prefix, length).
- [ ] If joining datasets, the ID intersection has been checked; mismatches logged.

---

## Types

- [ ] Numeric columns are numeric (`int64` or `float64`, not `object`).
- [ ] Categorical columns are categorical (or string with known cardinality).
- [ ] Boolean columns are `bool`, not 0/1 or "yes"/"no".
- [ ] Date/time columns are parsed.
- [ ] No silent demotion of integer to float due to NA presence.

---

## Ranges and domain checks

- [ ] Numeric columns within plausible ranges (age ≥ 0; concentrations ≥ 0; percentages 0–100; correlations -1 to 1).
- [ ] Categorical columns have only expected levels (no typos creating spurious extra factors).
- [ ] Outliers identified, examined, and either kept-with-note or removed-with-reason.
- [ ] Units are consistent across rows (no mix of "kg" and "g" in the same column).

---

## Missingness

- [ ] Per-column missingness counted and reported.
- [ ] Pattern of co-missingness examined (do certain columns go missing together?).
- [ ] Missingness mechanism assessed (MCAR / MAR / MNAR).
- [ ] Imputation method (if any) chosen and documented.
- [ ] Indicator columns added for imputed values if material.

---

## Metadata alignment

- [ ] Sample table and metadata table share the same ID column.
- [ ] After join, no IDs were silently dropped.
- [ ] Metadata fields used in modeling are non-missing or imputed.

---

## Batch / cohort structure

- [ ] Batch column present (or absence documented).
- [ ] Batch × outcome chi-square not confounded (p > 0.05).
- [ ] Sex / age / site / collection date balanced as expected.
- [ ] PCA colored by batch reviewed.

---

## Distribution sanity

- [ ] Min, max, mean, median, SD of each numeric column inspected.
- [ ] Histograms or violin plots reviewed for unexpected distributions.
- [ ] Skewness or excess zeros noted (may inform transformation choice).
- [ ] Sample-sample correlations reasonable within groups.

---

## Reproducibility setup

- [ ] Random seeds set at the top of the script.
- [ ] Package versions pinned (requirements.txt or equivalent).
- [ ] Working directory clearly defined.
- [ ] Raw data is read-only; processed data is in a separate location.

---

## Documentation

- [ ] Cleaning log written for every row/column removed.
- [ ] Validation report saved alongside the dataset.
- [ ] Decisions about missingness / outliers / batch documented.
- [ ] Software versions captured.

---

## Sign-off

Before moving to modeling:

- [ ] All checks above passed or explicit exceptions logged.
- [ ] User has reviewed the validation summary.
- [ ] Cleaned data has been hashed; the hash is in the analysis manifest.

If any check failed and was waived, the waiver must be in the cleaning log with a reason. Silent exceptions are not acceptable.

---

## Cross-references

- Missing data — [missing_data.md](missing_data.md).
- Multiple testing — [multiple_testing.md](multiple_testing.md).
- Batch effects — [batch_effects.md](batch_effects.md).
- Reproducibility — [reproducibility.md](reproducibility.md).
