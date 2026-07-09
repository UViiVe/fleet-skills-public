#!/usr/bin/env python3
"""Aggregate variant-level annotations into a gene-level burden table.

For every gene symbol represented in the parsed variant table, this script
counts variants, breaks them down by impact (HIGH / MODERATE / LOW / MODIFIER),
counts ClinVar pathogenic / likely-pathogenic alleles, computes median CADD if
present, and lists the top consequence.

Public function: ``annotate_gene_summary(df, output_path=None)``.

Example:
    python scripts/annotate_genes.py --input all_variants.csv \\
        --output gene_summary.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install required packages with: pip install pandas numpy")


_IMPACT_LEVELS = ("HIGH", "MODERATE", "LOW", "MODIFIER")


def annotate_gene_summary(
    df: pd.DataFrame, output_path: str | Path | None = None,
) -> pd.DataFrame:
    """Return one row per gene_symbol with burden + impact breakdown."""
    if df.empty or "gene_symbol" not in df.columns:
        return pd.DataFrame()

    work = df.copy()
    work["impact"] = work.get("impact", "").astype(str).str.upper()

    grouped = work.groupby("gene_symbol", dropna=True)
    rows = []
    for gene, g in grouped:
        if not gene:
            continue
        impact_counts = g["impact"].value_counts().to_dict()
        record = {
            "gene_symbol": gene,
            "n_variants":  int(len(g)),
            "n_HIGH":      int(impact_counts.get("HIGH", 0)),
            "n_MODERATE":  int(impact_counts.get("MODERATE", 0)),
            "n_LOW":       int(impact_counts.get("LOW", 0)),
            "n_MODIFIER":  int(impact_counts.get("MODIFIER", 0)),
        }
        if "clinvar_sig" in g.columns:
            sig = g["clinvar_sig"].astype(str).str.lower()
            record["n_pathogenic"] = int(
                (sig.str.contains("pathogenic", na=False)
                 & ~sig.str.contains("benign", na=False)).sum()
            )
        if "cadd_phred" in g.columns:
            cadd = pd.to_numeric(g["cadd_phred"], errors="coerce")
            record["median_cadd"] = float(cadd.median()) if cadd.notna().any() else np.nan
            record["max_cadd"]    = float(cadd.max())    if cadd.notna().any() else np.nan
        if "revel" in g.columns:
            revel = pd.to_numeric(g["revel"], errors="coerce")
            record["max_revel"] = float(revel.max()) if revel.notna().any() else np.nan
        if "consequence" in g.columns:
            top_csq = g["consequence"].astype(str).value_counts().idxmax()
            record["top_consequence"] = top_csq
        rows.append(record)

    out = pd.DataFrame(rows).sort_values("n_HIGH", ascending=False, kind="mergesort")
    out = out.sort_values("n_variants",  ascending=False, kind="mergesort")

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(output_path, index=False)
        print(f"  -> {output_path}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    print("=== Gene-level summary ===")
    df = pd.read_csv(args.input)
    print(f"  loaded {len(df):,} rows")
    summary = annotate_gene_summary(df, args.output)
    print(f"✓ Gene summary completed: {len(summary):,} genes / {args.output}")


if __name__ == "__main__":
    main()
