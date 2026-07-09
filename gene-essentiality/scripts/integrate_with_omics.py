"""Correlate gene-effect scores with another omics layer.

Computes Spearman correlations between a gene-effect column and a paired
feature matrix (expression, mutation, proteomics, copy number). Crucially,
inverts the gene-effect score before correlating so that positive
correlation == higher feature → higher essentiality.

Example:
    python integrate_with_omics.py --gene-effect gene_effect.parquet \\
        --feature-matrix expression.parquet --gene KMT2A --out KMT2A_expr.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

try:
    import numpy as np
    import pandas as pd
    from scipy import stats
except ImportError:
    sys.exit("Install with: pip install pandas numpy scipy pyarrow")


def _find_row(matrix: "pd.DataFrame", symbol: str) -> Optional[str]:
    """Find the matrix row whose composite ID starts with '<symbol> ('."""
    candidates = [i for i in matrix.index if i.startswith(f"{symbol} (")]
    if len(candidates) == 0:
        return None
    if len(candidates) > 1:
        print(f"⚠ multiple matches for {symbol!r}: {candidates} — using first")
    return candidates[0]


def expression_essentiality_corr(
    gene_effect: "pd.DataFrame",
    feature_matrix: "pd.DataFrame",
    dep_gene: str,
    feature_gene: Optional[str] = None,
) -> dict:
    """Spearman correlation between feature_gene's value and dep_gene's essentiality.

    If feature_gene is None, defaults to dep_gene (same-gene expression vs
    same-gene essentiality — the standard "lineage-marker dependency" test).

    The gene-effect score is inverted before correlation, so the returned
    rho follows the intuitive sign:

      positive rho → higher feature → higher essentiality of dep_gene
      negative rho → higher feature → lower essentiality
    """
    feature_gene = feature_gene or dep_gene
    ge_id = _find_row(gene_effect, dep_gene)
    feat_id = _find_row(feature_matrix, feature_gene)
    if ge_id is None or feat_id is None:
        return {"error": f"missing row: dep={ge_id!r} feat={feat_id!r}"}

    shared = gene_effect.columns.intersection(feature_matrix.columns)
    ge = gene_effect.loc[ge_id, shared].dropna()
    feat = feature_matrix.loc[feat_id, shared].dropna()
    common = ge.index.intersection(feat.index)
    if len(common) < 10:
        return {"error": f"only {len(common)} cell lines in common"}

    essentiality = -ge.loc[common]  # invert
    rho, pval = stats.spearmanr(feat.loc[common], essentiality)
    return {
        "dep_gene": dep_gene,
        "feature_gene": feature_gene,
        "n_lines": int(len(common)),
        "spearman_rho": float(rho),
        "pval": float(pval),
        "sign_convention": "gene_effect inverted before correlation; positive rho = feature up → essentiality up",
    }


def mut_dependency(
    gene_effect: "pd.DataFrame",
    mutation_matrix: "pd.DataFrame",
    dep_gene: str,
    mut_gene: str,
    min_per_group: int = 3,
) -> dict:
    """Wilcoxon comparing dep_gene essentiality in mut_gene-mutant vs WT lines."""
    ge_id = _find_row(gene_effect, dep_gene)
    mut_id = _find_row(mutation_matrix, mut_gene)
    if ge_id is None or mut_id is None:
        return {"error": f"missing row: dep={ge_id!r} mut={mut_id!r}"}

    shared = gene_effect.columns.intersection(mutation_matrix.columns)
    mut_status = mutation_matrix.loc[mut_id, shared].fillna(0).astype(bool)
    ge = gene_effect.loc[ge_id, shared].dropna()
    mut = mut_status.reindex(ge.index).fillna(False)
    sub_mut = ge[mut]
    sub_wt = ge[~mut]
    if len(sub_mut) < min_per_group or len(sub_wt) < min_per_group:
        return {"error": f"too few: mut={len(sub_mut)}, wt={len(sub_wt)}"}
    delta = float(sub_mut.median() - sub_wt.median())
    _, p = stats.mannwhitneyu(sub_mut, sub_wt, alternative="less")
    return {
        "dep_gene": dep_gene,
        "mut_gene": mut_gene,
        "n_mut": int(len(sub_mut)),
        "n_wt": int(len(sub_wt)),
        "median_mut": float(sub_mut.median()),
        "median_wt": float(sub_wt.median()),
        "delta": delta,
        "pval": float(p),
    }


def genome_wide_correlation(
    gene_effect: "pd.DataFrame",
    feature_matrix: "pd.DataFrame",
    feature_gene: Optional[str] = None,
) -> "pd.DataFrame":
    """For each gene, correlate (feature value across lines) vs (its inverted essentiality).

    If feature_gene is None, this performs the per-gene self-correlation
    (each gene's expression vs its own essentiality across lines). Otherwise
    it uses the same fixed feature for every dep gene.
    """
    shared = gene_effect.columns.intersection(feature_matrix.columns)
    print(f"  {len(shared)} cell lines in common")

    ge = -gene_effect[shared]                     # invert once
    if feature_gene is None:
        common_genes = ge.index.intersection(feature_matrix.index)
        rho = pd.Series(np.nan, index=common_genes)
        pval = pd.Series(np.nan, index=common_genes)
        for i, g in enumerate(common_genes, 1):
            if i % 1000 == 0:
                print(f"  [{i}/{len(common_genes)}]")
            x = feature_matrix.loc[g, shared].dropna()
            y = ge.loc[g, x.index].dropna()
            xy = x.loc[y.index]
            if len(xy) >= 10:
                r, p = stats.spearmanr(xy.values, y.values)
                rho.loc[g] = r
                pval.loc[g] = p
        return pd.DataFrame({"spearman_rho": rho, "pval": pval}).dropna()
    else:
        feat_id = _find_row(feature_matrix, feature_gene)
        feat = feature_matrix.loc[feat_id, shared].dropna()
        rho = pd.Series(np.nan, index=ge.index)
        pval = pd.Series(np.nan, index=ge.index)
        for g in ge.index:
            y = ge.loc[g, feat.index].dropna()
            common = feat.index.intersection(y.index)
            if len(common) >= 10:
                r, p = stats.spearmanr(feat.loc[common].values, y.loc[common].values)
                rho.loc[g] = r
                pval.loc[g] = p
        return pd.DataFrame({"spearman_rho": rho, "pval": pval}).dropna()


def _read_matrix(path: Path) -> "pd.DataFrame":
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path, index_col=0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument("--feature-matrix", required=True)
    parser.add_argument("--mode", choices=["pair", "mut", "genome-wide-self", "genome-wide-feature"],
                        default="pair")
    parser.add_argument("--gene", required=True, help="Dependency gene symbol")
    parser.add_argument("--feature-gene", default=None,
                        help="Feature gene symbol (defaults to --gene for pair/mut)")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    ge = _read_matrix(Path(args.gene_effect))
    feat = _read_matrix(Path(args.feature_matrix))
    print(f"Loaded gene-effect ({ge.shape}) and feature ({feat.shape})")

    if args.mode == "pair":
        result = expression_essentiality_corr(ge, feat, args.gene, args.feature_gene)
        pd.Series(result).to_csv(args.out)
    elif args.mode == "mut":
        result = mut_dependency(ge, feat, args.gene, args.feature_gene or args.gene)
        pd.Series(result).to_csv(args.out)
    elif args.mode == "genome-wide-self":
        df = genome_wide_correlation(ge, feat, feature_gene=None)
        df.to_csv(args.out)
    elif args.mode == "genome-wide-feature":
        df = genome_wide_correlation(ge, feat, feature_gene=args.feature_gene or args.gene)
        df.to_csv(args.out)

    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
