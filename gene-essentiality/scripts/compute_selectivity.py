"""Cohort-vs-background selectivity analysis.

For each gene, compares its gene-effect distribution in a cohort of cell
lines to a background set, reporting median delta plus a Wilcoxon p-value
and BH-FDR. Pan-essentials are flagged but not automatically removed —
filter them upstream if desired.

Example:
    python compute_selectivity.py --in gene_effect.parquet \\
        --cohort AML_lines.txt --out selectivity.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Iterable, Optional

try:
    import numpy as np
    import pandas as pd
    from scipy.stats import mannwhitneyu
    from statsmodels.stats.multitest import multipletests
except ImportError:
    sys.exit("Install with: pip install pandas numpy scipy statsmodels pyarrow")


def selectivity_score(
    gene_effect: "pd.DataFrame",
    subset_lines: Iterable[str],
    background_lines: Optional[Iterable[str]] = None,
    min_lines_per_group: int = 3,
    direction: str = "less",
) -> "pd.DataFrame":
    """Per-gene cohort-vs-background selectivity table.

    Parameters
    ----------
    gene_effect : DataFrame
        Genes × cell lines (gene effect; negative = essential).
    subset_lines : iterable of str
        Cell-line IDs forming the cohort.
    background_lines : iterable of str or None
        If None, defaults to all columns not in subset.
    min_lines_per_group : int
        Minimum cohort and background size to run a test on a given gene.
    direction : 'less' or 'two-sided'
        Wilcoxon alternative. 'less' means "cohort more negative (more essential)".
    """
    subset = [c for c in subset_lines if c in gene_effect.columns]
    if background_lines is None:
        background = [c for c in gene_effect.columns if c not in set(subset)]
    else:
        background = [c for c in background_lines if c in gene_effect.columns]
    if not subset or not background:
        raise ValueError("Empty subset or background after intersecting with gene_effect columns")

    sub_med = gene_effect[subset].median(axis=1)
    bg_med = gene_effect[background].median(axis=1)
    delta = sub_med - bg_med

    cohort_dep_frac = (gene_effect[subset] < -0.5).mean(axis=1)
    bg_dep_frac = (gene_effect[background] < -0.5).mean(axis=1)

    pvals = pd.Series(np.nan, index=gene_effect.index)
    for gene in gene_effect.index:
        sub_vals = gene_effect.loc[gene, subset].dropna().values
        bg_vals = gene_effect.loc[gene, background].dropna().values
        if len(sub_vals) >= min_lines_per_group and len(bg_vals) >= min_lines_per_group:
            try:
                _, p = mannwhitneyu(sub_vals, bg_vals, alternative=direction)
                pvals.loc[gene] = p
            except ValueError:
                continue

    mask = pvals.notna()
    fdr = pd.Series(np.nan, index=pvals.index)
    if mask.any():
        _, q, _, _ = multipletests(pvals[mask].values, method="fdr_bh")
        fdr.loc[mask] = q

    df = pd.DataFrame({
        "cohort_median": sub_med,
        "background_median": bg_med,
        "delta": delta,
        "cohort_dep_frac": cohort_dep_frac,
        "background_dep_frac": bg_dep_frac,
        "pval": pvals,
        "fdr": fdr,
    }).sort_values("delta")

    return df


def permutation_selectivity(
    gene_effect: "pd.DataFrame",
    gene: str,
    subset_lines: Iterable[str],
    n_iter: int = 10000,
    seed: int = 0,
) -> dict:
    """Permutation-based selectivity p-value for one gene; use for small cohorts."""
    rng = np.random.default_rng(seed)
    subset = [c for c in subset_lines if c in gene_effect.columns]
    n_sub = len(subset)
    all_lines = list(gene_effect.columns)
    sub_vals = gene_effect.loc[gene, subset].dropna()
    bg_vals = gene_effect.loc[gene, [c for c in all_lines if c not in subset]].dropna()
    observed = float(sub_vals.median() - bg_vals.median())

    count = 0
    for _ in range(n_iter):
        sample = rng.choice(all_lines, size=n_sub, replace=False)
        bg = [c for c in all_lines if c not in set(sample)]
        s = gene_effect.loc[gene, list(sample)].dropna()
        b = gene_effect.loc[gene, bg].dropna()
        if len(s) and len(b) and (s.median() - b.median()) <= observed:
            count += 1
    return {"gene": gene, "delta": observed, "pval_permutation": count / n_iter, "n_iter": n_iter}


def _read_lines(path: str) -> list[str]:
    with open(path) as f:
        return [ln.strip() for ln in f if ln.strip() and not ln.startswith("#")]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True,
                        help="Gene-effect matrix (parquet or CSV)")
    parser.add_argument("--cohort", required=True,
                        help="Text file with one cell-line ID per line")
    parser.add_argument("--background", default=None,
                        help="Optional cell-line IDs file; defaults to all not in cohort")
    parser.add_argument("--out", required=True)
    parser.add_argument("--direction", default="less", choices=["less", "two-sided"])
    parser.add_argument("--min-lines", type=int, default=3)
    args = parser.parse_args()

    p = Path(args.in_path)
    if p.suffix == ".parquet":
        ge = pd.read_parquet(p)
    else:
        ge = pd.read_csv(p, index_col=0)
    print(f"Loaded gene-effect: {ge.shape}")

    cohort = _read_lines(args.cohort)
    background = _read_lines(args.background) if args.background else None
    print(f"Cohort lines: {len(cohort)}; background: {'all-other' if background is None else len(background)}")

    df = selectivity_score(
        ge,
        subset_lines=cohort,
        background_lines=background,
        min_lines_per_group=args.min_lines,
        direction=args.direction,
    )
    df.to_csv(args.out)
    n_sig = int((df["fdr"] < 0.05).sum())
    print(f"✓ Wrote {args.out} — {len(df)} genes, {n_sig} with FDR<0.05")


if __name__ == "__main__":
    main()
