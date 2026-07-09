"""Validate a tabular dataset before any analysis.

Inputs:
    --input PATH            Path to a CSV/TSV/parquet file.
    --id-col COL            Name of the ID column (default: first column).
    --metadata PATH         Optional metadata file to align against.
    --metadata-id COL       ID column in the metadata file (default matches --id-col).
    --output-dir DIR        Where to write the validation report / cleaning log.

Outputs:
    <output_dir>/validation_report.json   structured validation report
    <output_dir>/validation_summary.md    human-readable summary
    <output_dir>/cleaning_log.txt         audit trail

Example:
    python validate_data.py --input data.csv --id-col sample_id \
        --metadata meta.csv --output-dir out/
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Iterable

try:
    import pandas as pd
except ImportError:
    sys.exit("Install pandas with: pip install pandas")

try:
    import numpy as np
except ImportError:
    sys.exit("Install numpy with: pip install numpy")


@dataclass
class ValidationReport:
    """Structured report from validate_dataframe()."""

    n_rows: int
    n_cols: int
    id_col: str | None
    duplicate_columns: list[str] = field(default_factory=list)
    duplicate_id_count: int = 0
    missingness_pct: dict[str, float] = field(default_factory=dict)
    dtype_summary: dict[str, str] = field(default_factory=dict)
    outliers_per_col: dict[str, int] = field(default_factory=dict)
    id_overlap: dict[str, int] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"shape: {self.n_rows} rows x {self.n_cols} cols",
            f"id column: {self.id_col}",
            f"duplicate columns: {len(self.duplicate_columns)}",
            f"duplicate ids: {self.duplicate_id_count}",
            f"columns with > 20% missing: "
            f"{sum(1 for v in self.missingness_pct.values() if v > 20)}",
            f"warnings: {len(self.warnings)}",
        ]
        return "\n".join(lines)


def _read_table(path: Path) -> pd.DataFrame:
    """Read CSV/TSV/parquet by extension."""
    ext = path.suffix.lower()
    if ext in (".csv", ".txt"):
        return pd.read_csv(path)
    if ext in (".tsv",):
        return pd.read_csv(path, sep="\t")
    if ext in (".parquet",):
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            sys.exit("Install pyarrow with: pip install pyarrow")
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported file extension: {ext}")


def find_duplicate_columns(df: pd.DataFrame) -> list[str]:
    """Return column names that are duplicated (including pandas-renamed .1 suffixes)."""
    cols = df.columns
    dup_mask = cols.duplicated(keep=False) | cols.str.contains(r"\.\d+$")
    return list(set(cols[dup_mask]))


def find_duplicate_ids(df: pd.DataFrame, id_col: str) -> int:
    """Count duplicate IDs in id_col."""
    if id_col not in df.columns:
        return 0
    return int(df[id_col].duplicated().sum())


def find_duplicates(df: pd.DataFrame, id_col: str | None = None) -> dict[str, Any]:
    """Return duplicate-column, duplicate-ID, and duplicate-row counts."""
    return {
        "duplicate_columns": find_duplicate_columns(df),
        "duplicate_id_count": find_duplicate_ids(df, id_col) if id_col else 0,
        "duplicate_row_count": int(df.duplicated().sum()),
    }


def missingness_report(df: pd.DataFrame) -> dict[str, float]:
    """Per-column missingness as percentages, sorted descending."""
    miss = (df.isna().mean() * 100).round(2).to_dict()
    return dict(sorted(miss.items(), key=lambda kv: -kv[1]))


def check_types_and_ranges(
    df: pd.DataFrame,
    expected_types: dict[str, str] | None = None,
    ranges: dict[str, tuple[float, float]] | None = None,
) -> list[str]:
    """Validate expected dtypes and numeric ranges; return warning strings."""
    warnings: list[str] = []
    expected_types = expected_types or {}
    ranges = ranges or {}

    for col, expected in expected_types.items():
        if col not in df.columns:
            warnings.append(f"{col}: missing column")
            continue
        series = df[col]
        if expected in {"int", "integer"} and not pd.api.types.is_integer_dtype(series):
            warnings.append(f"{col}: expected integer, found {series.dtype}")
        elif expected in {"float", "numeric"} and not pd.api.types.is_numeric_dtype(series):
            warnings.append(f"{col}: expected numeric, found {series.dtype}")
        elif expected in {"str", "string"} and not pd.api.types.is_string_dtype(series):
            warnings.append(f"{col}: expected string, found {series.dtype}")

    for col, (lo, hi) in ranges.items():
        if col not in df.columns:
            warnings.append(f"{col}: missing column for range check")
            continue
        numeric = pd.to_numeric(df[col], errors="coerce")
        bad = numeric.notna() & ((numeric < lo) | (numeric > hi))
        if bad.any():
            warnings.append(f"{col}: {int(bad.sum())} value(s) outside [{lo}, {hi}]")

    return warnings


def batch_summary(
    df: pd.DataFrame,
    batch_col: str,
    outcome_col: str,
) -> dict[str, Any]:
    """Summarize batch/outcome balance with an optional chi-square test."""
    if batch_col not in df.columns or outcome_col not in df.columns:
        missing = [c for c in (batch_col, outcome_col) if c not in df.columns]
        raise ValueError(f"missing columns: {missing}")

    table = pd.crosstab(df[batch_col], df[outcome_col])
    result: dict[str, Any] = {
        "table": table,
        "batch_col": batch_col,
        "outcome_col": outcome_col,
    }
    try:
        from scipy.stats import chi2_contingency
    except ImportError:
        result["warning"] = "Install scipy for chi-square p-value"
        return result

    chi2, p_value, dof, expected = chi2_contingency(table)
    result.update({
        "chi2": float(chi2),
        "p_value": float(p_value),
        "degrees_of_freedom": int(dof),
        "expected": expected,
    })
    return result


def outlier_count(series: pd.Series, z_threshold: float = 5.0) -> int:
    """Count |z| > threshold for numeric series; 0 otherwise."""
    if not pd.api.types.is_numeric_dtype(series):
        return 0
    s = series.dropna()
    if len(s) < 3 or s.std(ddof=0) == 0:
        return 0
    z = (s - s.mean()) / s.std(ddof=0)
    return int((z.abs() > z_threshold).sum())


def validate_dataframe(
    df: pd.DataFrame,
    id_col: str | None = None,
    z_threshold: float = 5.0,
) -> ValidationReport:
    """Run all validation checks on a DataFrame and return a ValidationReport."""
    warnings: list[str] = []

    n_rows, n_cols = df.shape
    if id_col is not None and id_col not in df.columns:
        warnings.append(f"id_col '{id_col}' not found; first column will be used")
        id_col = df.columns[0]
    if id_col is None:
        id_col = df.columns[0]

    dup_cols = find_duplicate_columns(df)
    if dup_cols:
        warnings.append(
            f"duplicate columns detected: {dup_cols} — investigate upstream joins"
        )

    dup_ids = find_duplicate_ids(df, id_col)
    if dup_ids > 0:
        warnings.append(f"{dup_ids} duplicate ids in '{id_col}'")

    miss = missingness_report(df)
    high_miss = {c: v for c, v in miss.items() if v > 20}
    if high_miss:
        warnings.append(
            f"columns with > 20% missing: {list(high_miss.keys())}"
        )

    dtypes = {c: str(t) for c, t in df.dtypes.items()}
    object_cols_with_numbers = []
    for c in df.columns:
        if df[c].dtype == object:
            try:
                pd.to_numeric(df[c].dropna(), errors="raise")
                object_cols_with_numbers.append(c)
            except (ValueError, TypeError):
                pass
    if object_cols_with_numbers:
        warnings.append(
            f"object-dtype columns that look numeric: "
            f"{object_cols_with_numbers}; check for stray non-numeric values"
        )

    outliers = {c: outlier_count(df[c], z_threshold) for c in df.columns}
    outliers = {c: v for c, v in outliers.items() if v > 0}

    return ValidationReport(
        n_rows=n_rows,
        n_cols=n_cols,
        id_col=id_col,
        duplicate_columns=dup_cols,
        duplicate_id_count=dup_ids,
        missingness_pct=miss,
        dtype_summary=dtypes,
        outliers_per_col=outliers,
        warnings=warnings,
    )


def align_metadata(
    df: pd.DataFrame,
    meta: pd.DataFrame,
    id_col: str,
    meta_id_col: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Reduce both dataframes to their ID intersection. Return (df, meta, overlap_report)."""
    if meta_id_col is None:
        meta_id_col = id_col

    data_ids = set(df[id_col].dropna().astype(str))
    meta_ids = set(meta[meta_id_col].dropna().astype(str))

    common = data_ids & meta_ids
    only_data = data_ids - meta_ids
    only_meta = meta_ids - data_ids

    overlap = {
        "common": len(common),
        "only_in_data": len(only_data),
        "only_in_metadata": len(only_meta),
    }

    df_out = df[df[id_col].astype(str).isin(common)].copy()
    meta_out = meta[meta[meta_id_col].astype(str).isin(common)].copy()
    return df_out, meta_out, overlap


def load_and_check(path: str | Path, id_col: str | None = None) -> pd.DataFrame:
    """Load a tabular file with a shape check."""
    df = _read_table(Path(path))
    if id_col and id_col not in df.columns:
        sys.exit(f"id_col '{id_col}' not in file columns: {list(df.columns)[:8]}...")
    print(f"✓ Data loaded: {len(df)} rows x {df.shape[1]} cols / {path}")
    return df


def write_report(report: ValidationReport, output_dir: Path) -> None:
    """Write JSON and Markdown reports."""
    output_dir.mkdir(parents=True, exist_ok=True)

    json_path = output_dir / "validation_report.json"
    with json_path.open("w") as fh:
        json.dump(asdict(report), fh, indent=2, default=str)

    md_path = output_dir / "validation_summary.md"
    with md_path.open("w") as fh:
        fh.write("# Validation Summary\n\n")
        fh.write(f"- rows: {report.n_rows}\n")
        fh.write(f"- cols: {report.n_cols}\n")
        fh.write(f"- id column: {report.id_col}\n")
        fh.write(f"- duplicate columns: {report.duplicate_columns}\n")
        fh.write(f"- duplicate ids: {report.duplicate_id_count}\n\n")
        fh.write("## Missingness (top 10)\n\n")
        for c, v in list(report.missingness_pct.items())[:10]:
            fh.write(f"- {c}: {v}%\n")
        fh.write("\n## Warnings\n\n")
        for w in report.warnings:
            fh.write(f"- {w}\n")

    print(f"✓ Validation reports written: {json_path}, {md_path}")


def write_cleaning_log(
    output_path: str | Path,
    removed_rows: Iterable[Any],
    removed_cols: Iterable[str],
    reasons: dict[str, str] | None = None,
) -> None:
    """Append a cleaning audit trail to disk."""
    reasons = reasons or {}
    p = Path(output_path)
    with p.open("a") as fh:
        fh.write(f"# Cleaning log entry — {pd.Timestamp.now().isoformat()}\n")
        for r in removed_rows:
            fh.write(f"removed_row: {r}  reason: {reasons.get(str(r), 'n/a')}\n")
        for c in removed_cols:
            fh.write(f"removed_col: {c}  reason: {reasons.get(c, 'n/a')}\n")
    print(f"✓ Cleaning log appended: {p} ({sum(1 for _ in [removed_rows, removed_cols])} entries)")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate a tabular dataset before analysis."
    )
    parser.add_argument("--input", required=True, help="CSV/TSV/parquet file")
    parser.add_argument("--id-col", default=None, help="ID column name")
    parser.add_argument("--metadata", default=None, help="Optional metadata file")
    parser.add_argument("--metadata-id", default=None, help="ID column in metadata")
    parser.add_argument("--output-dir", default="validation_out", help="Output directory")
    args = parser.parse_args()

    df = load_and_check(args.input, id_col=args.id_col)
    report = validate_dataframe(df, id_col=args.id_col)

    if args.metadata:
        meta = _read_table(Path(args.metadata))
        df, meta, overlap = align_metadata(
            df, meta,
            id_col=args.id_col or df.columns[0],
            meta_id_col=args.metadata_id,
        )
        report.id_overlap = overlap

    write_report(report, Path(args.output_dir))
    print("\n" + report.summary())


if __name__ == "__main__":
    main()
