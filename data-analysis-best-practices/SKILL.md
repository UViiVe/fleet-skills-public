---
name: data-analysis-best-practices
id: "skill_5cbca10881eb4a97abc9cda42c1303b1"
display-name: "Data Analysis Best Practices"
short-description: "Pre-analysis validation, missing-data handling, multiple-testing correction, batch-effect checks, reproducibility, and seed/version pinning for any data analysis. Use whenever you're handed a datas..."
category: methodology
keywords: [data analysis, statistics, best practices, QC, methodology]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Pre-analysis validation, missing-data handling, multiple-testing correction, batch-effect checks, reproducibility, and seed/version pinning for any data analysis. Use whenever you're handed a dataset to analyze — before the first model is fit."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "methodology guide (best-practice validation and reproducibility framework)"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Sandve GK, et al. Ten Simple Rules for Reproducible Computational Research. PLOS Computational Biology. 2013;9:e1003285. doi:10.1371/journal.pcbi.1003285 | Wilson G, et al. Good enough practices in scientific computing. PLOS Computational Biology. 2017;13:e1005510. doi:10.1371/journal.pcbi.1005510"
changelog:
  - version: "1.0.1"
    date: "2026-07-07"
    event: citation-provenance
    summary: "Added provenance block and References section."
status: active
promoted: true
---

# Data Analysis Best Practices

A pre-analysis discipline. Before running any model, fitting any regression, or generating any figure on user-supplied data, work through this skill's checks. The goal is to catch the boring, predictable failures that destroy entire analyses — duplicate columns, mis-aligned IDs, silent NAs, an uncorrected p-value table — while they're still cheap to fix.

This is a companion skill: it doesn't do science on its own, it makes other skills' science trustworthy.

---

## When to Use This Skill

✅ **Use when:**
- A user hands you a CSV, TSV, parquet, or RDS and asks for "an analysis".
- You're about to fit a model to data that has not been through QC.
- You need to decide how to handle missing values, duplicates, or ID mismatches.
- You're about to report p-values across many tests.
- You're starting work in a new repo and want to lock in seeds, versions, and environments.
- A teammate's "the numbers look weird" is escalated to you.

❌ **Do not use when:**
- The dataset has already passed a documented QC pipeline (cite that pipeline instead).
- You're prototyping on a synthetic toy dataset that explicitly doesn't need validation.
- The task is purely visual — design a figure on already-clean data.

**Core principle:** trust nothing about a dataset until you have looked at its shape, ID set, missingness, and type assumptions yourself. Most "weird results" are upstream data hygiene, not novel biology or modeling failure.

---

## Quick Start

A 30-second pre-flight on any tabular dataset:

```python
import pandas as pd, numpy as np
from scripts.validate_data import validate_dataframe

df = pd.read_csv("data.csv")
report = validate_dataframe(df, id_col="sample_id")
print(report.summary())
```

`validate_dataframe` returns a structured report with duplicate columns, duplicate IDs, missingness per column, type inference, and outlier flags. Read it before doing anything else.

---

## Installation

| Package | Version | License | Commercial use | Install |
|---|---|---|---|---|
| pandas | ≥ 2.0 | BSD-3 | yes | `pip install pandas` |
| numpy | ≥ 1.24 | BSD-3 | yes | `pip install numpy` |
| scipy | ≥ 1.10 | BSD-3 | yes | `pip install scipy` |
| statsmodels | ≥ 0.14 | BSD-3 | yes | `pip install statsmodels` |
| scikit-learn | ≥ 1.3 | BSD-3 | yes | `pip install scikit-learn` |
| missingno | ≥ 0.5 | MIT | yes | `pip install missingno` |
| pyyaml | ≥ 6.0 | MIT | yes | `pip install pyyaml` |

For a reproducible environment:

```bash
pip install pandas numpy scipy statsmodels scikit-learn missingno pyyaml
pip freeze > requirements.txt
```

For R workflows:

| Package | Purpose |
|---|---|
| `mice` | Multiple imputation |
| `naniar` | Missing-data visualization |
| `renv` | Environment pinning |

---

## Inputs

**Required:**
- A tabular dataset (CSV, TSV, parquet, RDS, h5ad colData, or pandas/data.frame in memory).
- The intended ID column (or "row order is the ID" with confirmation).

**Optional:**
- A metadata table to be joined.
- An a priori list of expected columns / types / value ranges.
- Domain knowledge about valid value ranges (e.g., age 0–120, expression ≥ 0).

---

## Outputs

| File | Description |
|---|---|
| `validation_report.json` | Structured report: shape, duplicates, missingness, type checks, outliers, ID mismatches. |
| `validation_summary.md` | Human-readable summary suitable for sharing. |
| `cleaning_log.txt` | Audit trail of every row/column removed, with reason. |
| `cleaned_data.csv` (optional) | Cleaned dataset with documented exclusions. |
| `environment.yml` / `requirements.txt` | Pinned package versions. |
| `seeds.json` | Random seeds used for every stochastic step. |

The cleaning log is critical: every removal must be auditable. "I dropped 12 samples" is not acceptable; "I dropped 12 samples with > 20% missingness — see cleaning_log.txt rows 1-12" is.

---

## Clarification Questions

### 1. ⚠️ CRITICAL: ASK THIS FIRST — What does this dataset claim to be?

Before validating, get explicit answers:
- What is the unit of observation (a sample? a patient? a cell? a gene?).
- What is the intended ID column? (If "the index", confirm uniqueness.)
- What columns are expected to be present? What types?
- Are there known data-quality issues already?

Without this, validation is just descriptive statistics. With it, validation becomes a contract check.

### 2. How should missing data be handled?

Options:
- Listwise deletion (drop rows with any missing).
- Column-level: drop columns above a missingness threshold.
- Imputation (mean / median / KNN / model-based / multiple imputation).
- Carry through with explicit `NA` and downstream-tool support.

The right answer depends on the mechanism (MCAR / MAR / MNAR) and the downstream method. See [references/missing_data.md](references/missing_data.md).

### 3. Multiple-testing context

If results will involve many comparisons:
- How many tests? (Drives method choice.)
- What's the false-positive tolerance — FDR or FWER?
- Standard threshold (0.05) or stricter (0.01) or looser (0.10 for exploratory)?

See [references/multiple_testing.md](references/multiple_testing.md).

### 4. Reproducibility expectations

- Should this analysis be exactly re-runnable from raw data?
- Where do random seeds get set?
- What's the package-pinning strategy (lockfile, container, both)?

See [references/reproducibility.md](references/reproducibility.md).

### 5. Batch / cohort structure

- Were samples collected over time, at different sites, or by different operators?
- Is batch information recorded in the metadata?

Even if the immediate task isn't a batch-effect study, batch structure changes which checks to run.

---

## Standard Workflow

### Step 1 — Load and shape-check

```python
from scripts.validate_data import load_and_check
df = load_and_check("data.csv", id_col="sample_id")
```

What this checks:
- File loads without ambiguous type coercion (no surprise string-to-numeric or vice versa).
- Shape matches expectations (or asks if expectations weren't provided).
- ID column is present and unique.

Verify: `✓ Data loaded: <n_rows> rows, <n_cols> cols`. If shape is unexpected, stop and ask.

### Step 2 — Duplicate detection

```python
from scripts.validate_data import find_duplicates
dup_report = find_duplicates(df, id_col="sample_id")
```

Three kinds of duplicate to catch:
- Duplicate column names. Pandas appends `.1`, `.2`, etc. — these are usually a sign of an upstream join mistake. Drop or rename.
- Duplicate IDs. Same sample appearing twice. Investigate before deduping — sometimes it's a legitimate technical replicate.
- Fully duplicated rows. Often a CSV-paste accident.

### Step 3 — Missingness audit

```python
from scripts.validate_data import missingness_report
miss = missingness_report(df)
```

Produces per-column missing counts, percentages, and missingness pattern co-occurrence. Decide column-by-column:

| Missingness | Default action |
|---|---|
| 0% | nothing |
| < 5% | impute with median / mode and document |
| 5–20% | impute with care; consider multiple imputation |
| 20–50% | drop the column or carefully model the missingness |
| > 50% | drop the column unless it's load-bearing |

Mechanism matters: missing-at-random (MAR) and missing-completely-at-random (MCAR) can be imputed; missing-not-at-random (MNAR) usually requires modeling the missingness explicitly. See [references/missing_data.md](references/missing_data.md).

### Step 4 — Metadata alignment

If a separate metadata file is provided:

```python
from scripts.validate_data import align_metadata
df, meta, overlap = align_metadata(df, meta, id_col="sample_id")
```

Reports IDs only in data, only in metadata, and the intersection. Always reduce to the intersection (unless explicitly told otherwise) and log what was dropped.

### Step 5 — Type and range checks

```python
from scripts.validate_data import check_types_and_ranges
check_types_and_ranges(df, expected_types={"age": "int", "expression": "float"},
                       ranges={"age": (0, 120)})
```

Catches: ages > 120, negative concentrations, scientific-notation strings parsed as object, columns silently demoted to object due to one bad cell.

### Step 6 — Batch / cohort structure check

```python
from scripts.validate_data import batch_summary
batch_summary(df, batch_col="batch", outcome_col="condition")
```

Runs a chi-square test of independence between batch and condition. If p < 0.05, you have a confounding problem — see [experimental-design-statistics](../experimental-design-statistics/SKILL.md) before proceeding.

### Step 7 — Apply multiple-testing correction (if relevant)

```python
from scripts.multiple_testing import adjust_pvalues
df["q_value"] = adjust_pvalues(df["p_value"], method="bh")
```

Pre-specify the method. Don't shop methods after seeing the count of significant hits. See [references/multiple_testing.md](references/multiple_testing.md).

### Step 8 — Reproducibility lock-in

```python
from scripts.repro import write_environment, set_seeds
set_seeds(seed=42)
write_environment("requirements.txt")
```

Writes pinned versions and seeds to disk. Re-running with the same seeds and pins should reproduce the analysis exactly.

### Step 9 — Export the cleaning audit

```python
from scripts.validate_data import write_cleaning_log
write_cleaning_log("cleaning_log.txt", removed_rows, removed_cols, reasons)
```

---

## When Scripts Fail

1. **Fix and retry (90%).** Install the missing package or adjust the input format. Re-run the script.
2. **Modify the script (5%).** If a threshold needs tuning (e.g., 20% missingness is too strict for survey data), edit the script and record the change in the cleaning log.
3. **Use as reference (4%).** If the dataset shape genuinely doesn't fit the assumed schema (e.g., long-form when wide was assumed), read the script, adapt, and cite.
4. **Write from scratch (1%).** Only if the existing scripts genuinely cannot handle the data type.

Never silently bypass a validation check because it's "in the way" — fix the underlying issue and document.

---

## Decision Points

| Decision | Default | Override when |
|---|---|---|
| Missingness column drop threshold | 20% | Survey data with structured missingness — higher; biomarker data — lower |
| Imputation method | median + indicator | Many missing patterns → multiple imputation (`mice`) |
| ID mismatch policy | use intersection | Sometimes asymmetric joins are intentional |
| Multiple-testing method | BH-FDR | < 100 tests → Bonferroni; GWAS → standard 5e-8 |
| Random seed | 42 (set everywhere) | Replication studies use the seed from the prior run |
| Outlier handling | flag, don't drop | Drop only after biological justification |
| Numeric tolerance | float64 | float32 for memory; document loss |
| Type coercion | strict (`dtype` per column) | Exploratory only — `low_memory=False` |

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| Duplicate columns named `x.1`, `x.2` | pandas auto-renamed duplicates on read | Drop all but the first; verify upstream join logic |
| IDs in data not in metadata | Upstream filter dropped samples | Use the intersection; log the missing IDs |
| Numeric column read as object | One non-numeric cell ("NA" as string, "—", "<LOD") | Specify `na_values=["NA","—","<LOD"]` on read |
| Missing fraction is 0 but values are 0 | Sentinel zeros being treated as data | Check raw file for explicit nulls vs. zeros |
| p-values all < 0.05 after BH | Test calibration issue | Inspect QQ plot; suspect inflation from batch effects |
| Reproducible runs give different numbers | Stochastic step missing a seed | Set seeds at every `random`/`np.random`/`torch.manual_seed` |
| Tests pass locally, fail in CI | Unpinned dependency drift | Use a lockfile (`pip-compile`, `poetry.lock`, `renv.lock`) |
| Different results across float dtypes | Loss-of-precision in float32 | Pin dtype; convert only at known-safe places |
| Date column parsed as string | Format not detected | Specify `parse_dates=["date"]` with `date_format` |
| Outlier z-score > 5 | Could be data-entry error or true biology | Investigate; only drop with biological justification |
| ID column has whitespace differences | "S01 " vs "S01" treated as different | `.str.strip()` before any join |
| Encoding mismatch (UTF-8 vs latin-1) | Mixed-source data | Specify `encoding=` explicitly on read |

Deeper: [references/missing_data.md](references/missing_data.md), [references/reproducibility.md](references/reproducibility.md), [references/multiple_testing.md](references/multiple_testing.md), [references/batch_effects.md](references/batch_effects.md).

---

## Best Practices

- **Validate at the boundary, trust internally.** When data enters your pipeline, validate exhaustively. After that, downstream code can trust the schema.
- **Never silently drop data.** Every removal goes in the cleaning log with a reason.
- **Lock seeds at the top of the script.** `numpy.random.seed`, `random.seed`, `torch.manual_seed` — all of them, before any stochastic call.
- **Pin versions to a lockfile.** Track what dependency drift can break your numbers.
- **Run cleanings as a function, not inline.** Easier to test; easier to re-run.
- **Pre-specify multiple-testing method before seeing results.** Method-shopping is a form of p-hacking.
- **Re-check after every join.** Joins are the highest-incident point for ID mismatches and accidental row inflation.
- **Save intermediates.** A pickle of the post-validation data lets you re-enter the analysis without re-validating.
- **Treat warnings as errors in CI.** Especially `pandas` `SettingWithCopyWarning` and `numpy` runtime warnings.

---

## Suggested Next Steps

After validation, depending on the analysis type:

1. **experimental-design-statistics** — if the validation surfaced batch confounding or sample-size concerns, escalate to the design skill before continuing.
2. **bulk-rnaseq-counts-to-de-deseq2** / similar assay-specific DE skills — once data is clean and properly aligned with metadata.
3. **functional-enrichment-from-degs** — downstream of DE, with the pre-specified multiple-testing method.
4. **survival-analysis-clinical** — once event/censoring columns have been validated.
5. **lasso-biomarker-panel** — when feature selection follows validation.

---

## Related Skills

- **Companion:** [experimental-design-statistics](../experimental-design-statistics/SKILL.md) for design-time decisions on multiple testing and batch.
- **Upstream of:** every assay-specific DE skill (RNA-seq, ATAC-seq, methylation, proteomics).
- **Reference:** [pathway-enrichment](../pathway-enrichment/SKILL.md) for the multiple-testing discipline downstream of enrichment.

---

## References

**Bundled:**
- [references/missing_data.md](references/missing_data.md) — mechanisms (MCAR, MAR, MNAR), imputation methods, when to use each.
- [references/multiple_testing.md](references/multiple_testing.md) — BH-FDR vs Bonferroni vs q-value vs IHW; how to choose.
- [references/batch_effects.md](references/batch_effects.md) — detection and adjustment when validation reveals batch structure.
- [references/reproducibility.md](references/reproducibility.md) — seeds, lockfiles, environment management.
- [references/data_validation.md](references/data_validation.md) — the validation checklist in long form.

**Scripts:**
- `scripts/validate_data.py` — main validation entry point.
- `scripts/missing_data.py` — missingness audit and imputation.
- `scripts/multiple_testing.py` — p-value adjustment with method selection.
- `scripts/repro.py` — seeds + environment pinning utilities.
- `scripts/batch_check.py` — batch-by-condition confounding detection.

**Primary citations:**
- **Ten Simple Rules**: Sandve GK, Nekrutenko A, Taylor J, Hovig E. "Ten Simple Rules for Reproducible Computational Research." *PLOS Computational Biology* 2013;9:e1003285. doi:10.1371/journal.pcbi.1003285 — [link](https://doi.org/10.1371/journal.pcbi.1003285)
- **Good enough practices**: Wilson G, Bryan J, Cranston K, et al. "Good enough practices in scientific computing." *PLOS Computational Biology* 2017;13:e1005510. doi:10.1371/journal.pcbi.1005510 — [link](https://doi.org/10.1371/journal.pcbi.1005510)
- Wickham H. *Tidy Data.* J Stat Softw 2014;59(10).
- Little RJA, Rubin DB. *Statistical Analysis with Missing Data*, 3rd ed., Wiley, 2019.
- Benjamini Y, Hochberg Y. *Controlling the false discovery rate.* J R Stat Soc B 1995;57:289–300.
- Leek JT et al. *Tackling the widespread and critical impact of batch effects.* Nat Rev Genet 2010;11:733–739.
- Sandve GK, Nekrutenko A, Taylor J, Hovig E. *Ten simple rules for reproducible computational research.* PLoS Comp Biol 2013;9:e1003285.
