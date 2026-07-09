"""Check whether a batch-like column confounds with an outcome.

Inputs:
    --input PATH         CSV/TSV with both batch and outcome columns.
    --batch-col COL      Name of the batch column.
    --outcome-col COL    Name of the outcome column.
    --alpha FLOAT        Significance threshold (default 0.05).

Outputs:
    Prints chi-square test of independence; status PASS or FAIL.

Example:
    python batch_check.py --input data.csv \
        --batch-col batch --outcome-col condition
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("Install pandas with: pip install pandas")

try:
    from scipy import stats
except ImportError:
    sys.exit("Install scipy with: pip install scipy")


def test_confounding(
    df: pd.DataFrame,
    batch_col: str,
    outcome_col: str,
    alpha: float = 0.05,
) -> dict[str, object]:
    """Chi-square test of independence between batch and outcome.

    Returns a dict with status PASS/FAIL, p-value, contingency table.
    """
    if batch_col not in df.columns:
        raise KeyError(f"batch column '{batch_col}' not in dataframe")
    if outcome_col not in df.columns:
        raise KeyError(f"outcome column '{outcome_col}' not in dataframe")

    table = pd.crosstab(df[batch_col], df[outcome_col])
    chi2, p, dof, expected = stats.chi2_contingency(table)
    status = "FAIL" if p < alpha else "PASS"

    return {
        "status": status,
        "confounded": p < alpha,
        "chi2": float(chi2),
        "p_value": float(p),
        "dof": int(dof),
        "alpha": alpha,
        "contingency_table": table.to_dict(),
        "message": (
            f"Batch is CONFOUNDED with {outcome_col} (chi2 p = {p:.4g} < {alpha})"
            if p < alpha
            else f"No confounding between batch and {outcome_col} "
                 f"(chi2 p = {p:.4g} > {alpha})"
        ),
    }


def variance_partition(
    df: pd.DataFrame,
    numeric_outcome: str,
    batch_col: str,
    biology_col: str,
) -> dict[str, float]:
    """Variance explained by batch vs biology in a numeric outcome (OLS R^2)."""
    try:
        import statsmodels.formula.api as smf
    except ImportError:
        sys.exit("Install statsmodels with: pip install statsmodels")

    cols = [numeric_outcome, batch_col, biology_col]
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"columns not found: {missing}")

    full = smf.ols(f"{numeric_outcome} ~ C({batch_col}) + C({biology_col})",
                   data=df).fit()
    bio  = smf.ols(f"{numeric_outcome} ~ C({biology_col})", data=df).fit()
    batch = full.rsquared - bio.rsquared
    return {
        "full_r2":     float(full.rsquared),
        "biology_r2":  float(bio.rsquared),
        "batch_share": float(batch),
    }


def _read_table(path: Path) -> pd.DataFrame:
    ext = path.suffix.lower()
    if ext in (".csv", ".txt"):
        return pd.read_csv(path)
    if ext == ".tsv":
        return pd.read_csv(path, sep="\t")
    if ext == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported extension: {ext}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Batch / outcome confounding check.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--batch-col", required=True)
    parser.add_argument("--outcome-col", required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    args = parser.parse_args()

    df = _read_table(Path(args.input))
    result = test_confounding(df, args.batch_col, args.outcome_col, args.alpha)
    print(result["message"])
    print(f"✓ Confounding check completed: status={result['status']} "
          f"(p={result['p_value']:.4g})")


if __name__ == "__main__":
    main()
