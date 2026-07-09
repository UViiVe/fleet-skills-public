"""Identify pan-essential genes (broadly required across the panel).

Filters a gene-effect matrix down to genes that meet both a median-effect
threshold and a dependent-line-fraction threshold. These are the housekeeping
fitness genes — translation, splicing, proteasome — that should be excluded
before cohort-selectivity ranking.

Example:
    python find_pan_essentials.py --in gene_effect.parquet --out pan_essentials.csv \\
        --median-thr -1.0 --dep-frac-thr 0.8
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas pyarrow")


def pan_essential_genes(
    gene_effect: "pd.DataFrame",
    median_thr: float = -1.0,
    dep_frac_thr: float = 0.80,
    dep_score_thr: float = -0.5,
) -> "pd.DataFrame":
    """Score every gene; return rows satisfying the pan-essential definition.

    Parameters
    ----------
    gene_effect : DataFrame
        Genes (index) × cell lines (columns) matrix of gene-effect scores.
    median_thr : float
        Maximum allowed median effect (default -1.0).
    dep_frac_thr : float
        Minimum required fraction of cell lines below dep_score_thr.
    dep_score_thr : float
        Threshold defining "dependent" per cell line (default -0.5).
    """
    median = gene_effect.median(axis=1)
    dep_frac = (gene_effect < dep_score_thr).mean(axis=1)
    df = pd.DataFrame({
        "median_effect": median,
        "dependent_fraction": dep_frac,
        "is_pan_essential": (median < median_thr) & (dep_frac >= dep_frac_thr),
    })
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True,
                        help="Gene-effect matrix (parquet or csv)")
    parser.add_argument("--out", required=True, help="Output CSV")
    parser.add_argument("--median-thr", type=float, default=-1.0)
    parser.add_argument("--dep-frac-thr", type=float, default=0.80)
    parser.add_argument("--dep-score-thr", type=float, default=-0.5)
    args = parser.parse_args()

    p = Path(args.in_path)
    if p.suffix == ".parquet":
        ge = pd.read_parquet(p)
    else:
        ge = pd.read_csv(p, index_col=0)
    print(f"Loaded gene-effect matrix: {ge.shape[0]} genes × {ge.shape[1]} cell lines")

    scored = pan_essential_genes(
        ge,
        median_thr=args.median_thr,
        dep_frac_thr=args.dep_frac_thr,
        dep_score_thr=args.dep_score_thr,
    )
    n_pan = int(scored["is_pan_essential"].sum())
    print(f"  pan-essentials (median<{args.median_thr}, dep_frac≥{args.dep_frac_thr}): {n_pan}")

    out = scored.sort_values("median_effect")
    out.to_csv(args.out)
    print(f"✓ Wrote {args.out} — {n_pan} pan-essentials, {len(out) - n_pan} other")


if __name__ == "__main__":
    main()
