"""Export a long-form dependency parquet for downstream skills.

Pivots a (genes × cell lines) gene-effect matrix into a tidy
(cell_line, gene, gene_effect, scoring_algorithm, depmap_release) table —
the shape that ``multi-omics-integration`` and ``lasso-biomarker-panel``
expect.

Example:
    python export_dependencies.py --in gene_effect.parquet \\
        --release 24Q2 --algorithm Chronos --out dependencies.parquet
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas pyarrow")


def export_dependencies(
    gene_effect: "pd.DataFrame",
    release: str,
    scoring_algorithm: str,
    output_path: str,
    include_dependent_only: bool = False,
    dependency_threshold: float = -0.5,
) -> "pd.DataFrame":
    """Return a long-form dependency table and persist it to parquet.

    Parameters
    ----------
    gene_effect : DataFrame
        Genes × cell lines.
    release : str
        DepMap release tag (e.g. '24Q2').
    scoring_algorithm : str
        'CERES' or 'Chronos'.
    output_path : str
        Destination (parquet).
    include_dependent_only : bool
        If True, keep only rows with gene_effect < dependency_threshold.
    """
    df = (
        gene_effect.stack()
        .rename("gene_effect")
        .rename_axis(["gene", "cell_line"])
        .reset_index()
    )
    if include_dependent_only:
        df = df[df["gene_effect"] < dependency_threshold].copy()

    df["depmap_release"] = release
    df["scoring_algorithm"] = scoring_algorithm

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out, index=False)
    print(
        f"✓ Wrote {out} — {len(df):,} (cell_line, gene) rows "
        f"({df['gene'].nunique()} genes × {df['cell_line'].nunique()} cell lines)"
    )
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--release", required=True, help="DepMap release tag, e.g. 24Q2")
    parser.add_argument("--algorithm", default="Chronos", choices=["CERES", "Chronos", "unknown"])
    parser.add_argument("--out", required=True)
    parser.add_argument("--dependent-only", action="store_true",
                        help="Only emit rows where gene_effect < dependency-threshold")
    parser.add_argument("--dependency-threshold", type=float, default=-0.5)
    args = parser.parse_args()

    p = Path(args.in_path)
    ge = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, index_col=0)
    print(f"Loaded gene-effect: {ge.shape}")

    export_dependencies(
        gene_effect=ge,
        release=args.release,
        scoring_algorithm=args.algorithm,
        output_path=args.out,
        include_dependent_only=args.dependent_only,
        dependency_threshold=args.dependency_threshold,
    )


if __name__ == "__main__":
    main()
