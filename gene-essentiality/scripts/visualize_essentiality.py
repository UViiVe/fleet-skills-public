"""Standard visualizations for DepMap essentiality analyses.

- ``waterfall_plot``: gene-effect distribution across cell lines for one gene
- ``heatmap_selective``: z-scored essentiality heatmap for a gene set × lines
- ``scatter_feature_vs_essentiality``: paired scatter with inverted gene effect

Example:
    python visualize_essentiality.py --gene-effect gene_effect.parquet \\
        --gene KMT2A --sample-info Model.parquet --out-dir figures/
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Iterable, Optional

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
except ImportError:
    sys.exit("Install with: pip install matplotlib seaborn pandas numpy pyarrow")

sns.set_style("ticks")


def _find_row(matrix: "pd.DataFrame", symbol: str) -> Optional[str]:
    cands = [i for i in matrix.index if i.startswith(f"{symbol} (")]
    return cands[0] if cands else None


def _save(fig, path: str, dpi: int = 200) -> None:
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    fig.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  → {path}")


def waterfall_plot(
    gene_effect: "pd.DataFrame",
    gene: str,
    sample_info: Optional["pd.DataFrame"] = None,
    color_col: str = "OncotreeLineage",
    output_path: str = "waterfall.png",
) -> None:
    """Sorted gene-effect distribution across cell lines for one gene."""
    row_id = _find_row(gene_effect, gene)
    if row_id is None:
        raise KeyError(f"{gene!r} not in gene_effect matrix")
    ge = gene_effect.loc[row_id].dropna().sort_values()

    fig, ax = plt.subplots(figsize=(14, 5))
    x = np.arange(len(ge))
    if sample_info is not None and color_col in sample_info.columns:
        lineage = sample_info[color_col].reindex(ge.index).fillna("other")
        cats = sorted(lineage.unique())
        palette = sns.color_palette("tab20", n_colors=len(cats))
        cat2col = dict(zip(cats, palette))
        ax.bar(x, ge.values, color=[cat2col[c] for c in lineage], width=1.0, edgecolor="none")
        handles = [plt.Rectangle((0, 0), 1, 1, color=cat2col[c]) for c in cats]
        ax.legend(handles, cats, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
    else:
        ax.bar(x, ge.values, color="steelblue", width=1.0, edgecolor="none")

    ax.axhline(-1.0, c="gray", ls="--", lw=0.8, label="pan-essential floor")
    ax.axhline(-0.5, c="gray", ls=":", lw=0.8, label="dependency threshold")
    ax.axhline(0,    c="black", lw=0.5)
    ax.set_xlabel("Cell line (sorted by gene effect)")
    ax.set_ylabel("Gene effect (more negative = more essential)")
    ax.set_title(f"Waterfall: {gene}")
    ax.set_xlim(-1, len(ge))
    plt.tight_layout()
    _save(fig, output_path)


def heatmap_selective(
    gene_effect: "pd.DataFrame",
    gene_list: Iterable[str],
    cell_lines: Optional[Iterable[str]] = None,
    cluster: bool = True,
    output_path: str = "heatmap.png",
) -> None:
    """Z-scored within-gene heatmap of selective essentials × cell lines."""
    rows = []
    matched = []
    for g in gene_list:
        rid = _find_row(gene_effect, g) or (g if g in gene_effect.index else None)
        if rid is None:
            print(f"  ⚠ skip {g!r}: not in matrix")
            continue
        rows.append(rid)
        matched.append(g)
    if not rows:
        raise ValueError("None of the requested genes are in the matrix")

    cols = list(cell_lines) if cell_lines is not None else list(gene_effect.columns)
    sub = gene_effect.loc[rows, cols].dropna(axis=1, how="all")
    # z-score within each gene (row)
    z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1).replace(0, 1), axis=0)
    z.index = matched

    g = sns.clustermap(
        z,
        cmap="RdBu_r", center=0, vmin=-2, vmax=2,
        figsize=(max(8, 0.04 * z.shape[1]), max(4, 0.3 * z.shape[0])),
        col_cluster=cluster, row_cluster=cluster,
        xticklabels=False,
        cbar_kws={"label": "z(gene effect)"},
    )
    g.fig.suptitle("Selective essentials (z-scored per gene)", y=1.02)
    _save(g.fig, output_path)


def scatter_feature_vs_essentiality(
    gene_effect: "pd.DataFrame",
    feature_matrix: "pd.DataFrame",
    dep_gene: str,
    feature_gene: Optional[str] = None,
    output_path: str = "scatter.png",
) -> None:
    """Paired scatter; the y-axis is the **inverted** gene effect (essentiality)."""
    feature_gene = feature_gene or dep_gene
    ge_id = _find_row(gene_effect, dep_gene)
    feat_id = _find_row(feature_matrix, feature_gene)
    if ge_id is None or feat_id is None:
        raise KeyError(f"row not found: dep={ge_id!r} feat={feat_id!r}")
    shared = gene_effect.columns.intersection(feature_matrix.columns)
    ge = gene_effect.loc[ge_id, shared].dropna()
    feat = feature_matrix.loc[feat_id, shared].dropna()
    common = ge.index.intersection(feat.index)
    essentiality = -ge.loc[common]
    feat = feat.loc[common]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(feat, essentiality, s=14, alpha=0.5, c="steelblue", edgecolor="white", linewidth=0.5)
    ax.axhline(0, c="gray", lw=0.5, ls="--")
    ax.set_xlabel(f"{feature_gene} (feature)")
    ax.set_ylabel(f"{dep_gene} essentiality (-gene_effect; higher = more essential)")
    ax.set_title(f"{feature_gene} vs {dep_gene} essentiality  (n={len(common)})")
    plt.tight_layout()
    _save(fig, output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gene-effect", required=True)
    parser.add_argument("--mode", choices=["waterfall", "heatmap", "scatter"], default="waterfall")
    parser.add_argument("--gene", required=True)
    parser.add_argument("--gene-list", default=None, help="comma-separated for --mode heatmap")
    parser.add_argument("--sample-info", default=None)
    parser.add_argument("--feature-matrix", default=None)
    parser.add_argument("--feature-gene", default=None)
    parser.add_argument("--color-col", default="OncotreeLineage")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    def _read(path: str):
        p = Path(path)
        return pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, index_col=0)

    ge = _read(args.gene_effect)
    si = _read(args.sample_info) if args.sample_info else None

    if args.mode == "waterfall":
        waterfall_plot(ge, args.gene, sample_info=si, color_col=args.color_col, output_path=args.out)
    elif args.mode == "heatmap":
        genes = [g.strip() for g in (args.gene_list or args.gene).split(",")]
        heatmap_selective(ge, genes, output_path=args.out)
    elif args.mode == "scatter":
        feat = _read(args.feature_matrix)
        scatter_feature_vs_essentiality(ge, feat, args.gene, args.feature_gene, output_path=args.out)

    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
