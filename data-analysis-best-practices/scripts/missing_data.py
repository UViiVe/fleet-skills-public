"""Missing-data audit and imputation utilities.

Inputs:
    --input PATH              CSV/TSV/parquet file with missing values.
    --output PATH             Optional path to write imputed dataset.
    --method {median,knn,iterative,drop}
    --drop-threshold FLOAT    Drop columns missing more than this fraction (default 0.5).
    --seed INT                Random seed (default 42).

Outputs:
    JSON missingness report (printed to stdout).
    Imputed CSV at --output (if provided).

Example:
    python missing_data.py --input data.csv --method iterative --output imputed.csv
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Literal

try:
    import pandas as pd
    import numpy as np
except ImportError:
    sys.exit("Install pandas and numpy with: pip install pandas numpy")


def diagnose_missingness(df: pd.DataFrame) -> dict[str, object]:
    """Per-column and co-occurrence diagnostics."""
    pct = (df.isna().mean() * 100).round(2)
    counts = df.isna().sum()

    # Co-occurrence: which columns tend to be missing together?
    missing_mask = df.isna().astype(int)
    co_occurrence: dict[str, dict[str, float]] = {}
    cols_with_missing = [c for c in df.columns if counts[c] > 0]
    for ci in cols_with_missing:
        denom = int(counts[ci])
        if denom == 0:
            continue
        row = {}
        for cj in cols_with_missing:
            if ci == cj:
                continue
            joint = int(((missing_mask[ci] == 1) & (missing_mask[cj] == 1)).sum())
            row[cj] = round(joint / denom, 3)
        co_occurrence[ci] = row

    return {
        "n_rows": len(df),
        "n_cols": df.shape[1],
        "missing_pct_per_col": pct.to_dict(),
        "missing_count_per_col": counts.to_dict(),
        "co_occurrence": co_occurrence,
        "columns_above_20pct": [c for c in df.columns if pct[c] > 20],
        "columns_above_50pct": [c for c in df.columns if pct[c] > 50],
    }


def impute_median(df: pd.DataFrame) -> pd.DataFrame:
    """Median imputation for numeric columns, mode for object columns."""
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_numeric_dtype(out[c]):
            out[c] = out[c].fillna(out[c].median())
        else:
            mode = out[c].mode(dropna=True)
            if not mode.empty:
                out[c] = out[c].fillna(mode.iloc[0])
    return out


def impute_knn(df: pd.DataFrame, n_neighbors: int = 5) -> pd.DataFrame:
    """KNN imputation on numeric columns; categorical columns left unchanged."""
    try:
        from sklearn.impute import KNNImputer
    except ImportError:
        sys.exit("Install scikit-learn with: pip install scikit-learn")

    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    imp = KNNImputer(n_neighbors=n_neighbors)
    out = df.copy()
    if num_cols:
        out[num_cols] = imp.fit_transform(df[num_cols])
    return out


def impute_iterative(df: pd.DataFrame, max_iter: int = 10, seed: int = 42) -> pd.DataFrame:
    """Iterative (MICE-like) imputation on numeric columns."""
    try:
        from sklearn.experimental import enable_iterative_imputer  # noqa
        from sklearn.impute import IterativeImputer
    except ImportError:
        sys.exit("Install scikit-learn with: pip install scikit-learn")

    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    imp = IterativeImputer(max_iter=max_iter, random_state=seed)
    out = df.copy()
    if num_cols:
        out[num_cols] = imp.fit_transform(df[num_cols])
    return out


def drop_high_missing(df: pd.DataFrame, threshold: float = 0.5) -> tuple[pd.DataFrame, list[str]]:
    """Drop columns whose missing fraction exceeds threshold. Return (df, dropped)."""
    fracs = df.isna().mean()
    drop = list(fracs[fracs > threshold].index)
    return df.drop(columns=drop), drop


def add_missingness_indicator(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    """Append a was_missing_<c> boolean column for each column in `columns`."""
    out = df.copy()
    for c in columns:
        out[f"was_missing_{c}"] = df[c].isna()
    return out


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
    parser = argparse.ArgumentParser(description="Missing-data audit and imputation.")
    parser.add_argument("--input", required=True, help="CSV/TSV/parquet file")
    parser.add_argument("--output", default=None, help="Imputed output path (CSV)")
    parser.add_argument(
        "--method",
        choices=["median", "knn", "iterative", "drop", "none"],
        default="none",
    )
    parser.add_argument("--drop-threshold", type=float, default=0.5,
                        help="Drop columns missing > this fraction (default 0.5)")
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    np.random.seed(args.seed)

    df = _read_table(Path(args.input))
    diag = diagnose_missingness(df)
    print(json.dumps(diag, indent=2, default=str))

    if args.method == "none":
        print(f"✓ Diagnosis completed: {df.shape[0]} rows / no imputation requested")
        return

    if args.method == "drop":
        df_out, dropped = drop_high_missing(df, args.drop_threshold)
        print(f"✓ Drop-high-missing completed: removed {len(dropped)} columns: {dropped}")
    elif args.method == "median":
        df_out = impute_median(df)
        print(f"✓ Median imputation completed: {df_out.shape[0]} rows")
    elif args.method == "knn":
        df_out = impute_knn(df)
        print(f"✓ KNN imputation completed: {df_out.shape[0]} rows")
    elif args.method == "iterative":
        df_out = impute_iterative(df, seed=args.seed)
        print(f"✓ Iterative imputation completed: {df_out.shape[0]} rows")
    else:
        raise ValueError(args.method)

    if args.output:
        df_out.to_csv(args.output, index=False)
        print(f"✓ Imputed data written: {args.output}")


if __name__ == "__main__":
    main()
