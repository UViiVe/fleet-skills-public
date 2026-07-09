#!/usr/bin/env python3
"""Filter the parsed variant table by consequence, frequency, and impact.

Public functions:
  - ``filter_by_consequence(df, consequence_levels=...)`` — keep rows whose impact
    level is in the provided set.
  - ``filter_by_frequency(df, max_af=0.01, column='gnomad_af')`` — keep rare
    variants below ``max_af`` (NA treated as rare unless ``strict=True``).
  - ``filter_by_quality(df, min_qual=30, max_dp=None)`` — quality cutoffs.
  - ``filter_by_consequence_term(df, terms=...)`` — substring match on the
    SO consequence term (e.g. ``missense_variant``).
  - ``filter_clinically_significant(df)`` — keep variants with ClinVar
    pathogenic / likely-pathogenic annotation.

Example:
    python scripts/filter_variants.py --input all_variants.csv \\
        --output high_impact.csv --levels HIGH MODERATE
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("Install required packages with: pip install pandas")


_DEFAULT_IMPACT = ("HIGH", "MODERATE")


def filter_by_consequence(
    df: pd.DataFrame, *,
    consequence_levels: tuple[str, ...] = _DEFAULT_IMPACT,
    column: str = "impact",
) -> pd.DataFrame:
    if column not in df.columns:
        print(f"  warning: no '{column}' column — returning input unchanged")
        return df
    levels = {l.upper() for l in consequence_levels}
    keep = df[column].astype(str).str.upper().isin(levels)
    print(f"  filter_by_consequence: {keep.sum():,} / {len(df):,} kept "
          f"(levels = {sorted(levels)})")
    return df.loc[keep].copy()


def filter_by_consequence_term(
    df: pd.DataFrame, terms: list[str],
    *, column: str = "consequence",
) -> pd.DataFrame:
    if column not in df.columns:
        print(f"  warning: no '{column}' column — returning input unchanged")
        return df
    pattern = "|".join(terms)
    keep = df[column].astype(str).str.contains(pattern, regex=True, na=False)
    print(f"  filter_by_consequence_term: {keep.sum():,} / {len(df):,} kept "
          f"(terms = {terms})")
    return df.loc[keep].copy()


def filter_by_frequency(
    df: pd.DataFrame, *,
    max_af: float = 0.01,
    column: str = "gnomad_af",
    strict: bool = False,
) -> pd.DataFrame:
    """Keep rare variants. NA frequencies are treated as rare unless ``strict``."""
    if column not in df.columns:
        # Fall back to a common alternative naming.
        alt = "af" if "af" in df.columns else None
        if alt is None:
            print(f"  warning: no '{column}' or 'af' column — returning input")
            return df
        column = alt
    af = pd.to_numeric(df[column], errors="coerce")
    if strict:
        keep = af < max_af
    else:
        keep = (af.isna()) | (af < max_af)
    print(f"  filter_by_frequency: {int(keep.sum()):,} / {len(df):,} kept "
          f"(AF < {max_af} on '{column}')")
    return df.loc[keep].copy()


def filter_by_quality(
    df: pd.DataFrame, *,
    min_qual: float | None = 30.0,
    max_dp: float | None = None,
    min_dp: float | None = None,
) -> pd.DataFrame:
    keep = pd.Series(True, index=df.index)
    if min_qual is not None and "qual" in df.columns:
        keep &= pd.to_numeric(df["qual"], errors="coerce").fillna(0) >= min_qual
    if min_dp is not None and "dp" in df.columns:
        keep &= pd.to_numeric(df["dp"], errors="coerce").fillna(0) >= min_dp
    if max_dp is not None and "dp" in df.columns:
        keep &= pd.to_numeric(df["dp"], errors="coerce").fillna(0) <= max_dp
    print(f"  filter_by_quality: {int(keep.sum()):,} / {len(df):,} kept "
          f"(qual >= {min_qual})")
    return df.loc[keep].copy()


def filter_clinically_significant(
    df: pd.DataFrame, column: str = "clinvar_sig",
) -> pd.DataFrame:
    if column not in df.columns:
        print(f"  warning: no '{column}' column — returning input unchanged")
        return df
    sig = df[column].astype(str).str.lower()
    keep = sig.str.contains("pathogenic", na=False) & ~sig.str.contains(
        "benign", na=False)
    print(f"  filter_clinically_significant: {int(keep.sum()):,} / {len(df):,} kept")
    return df.loc[keep].copy()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--levels", nargs="+", default=list(_DEFAULT_IMPACT))
    parser.add_argument("--max-af", type=float, default=None)
    parser.add_argument("--af-column", default="gnomad_af")
    parser.add_argument("--min-qual", type=float, default=None)
    parser.add_argument("--clinically-significant", action="store_true")
    args = parser.parse_args()

    print("=== Filtering variants ===")
    df = pd.read_csv(args.input)
    print(f"  loaded {len(df):,} rows")

    df = filter_by_consequence(df, consequence_levels=tuple(args.levels))
    if args.max_af is not None:
        df = filter_by_frequency(df, max_af=args.max_af,
                                 column=args.af_column)
    if args.min_qual is not None:
        df = filter_by_quality(df, min_qual=args.min_qual)
    if args.clinically_significant:
        df = filter_clinically_significant(df)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"✓ Filtering completed: {len(df):,} rows / {args.output}")


if __name__ == "__main__":
    main()
