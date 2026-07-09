"""Per-perturbation UMAP/PCA plots and screen-level summary visualizations.

Saves PNG (always) and SVG (best-effort). The per-perturbation views colour by
classification (outlier vs control), perturbation label, sgRNA, and — if the
target gene was measured — the gene's own expression.

Example:
    python visualize_perturbations.py --de-csv results/initial_screening/TREM2_ttest.csv \\
        --perturbation TREM2 --out-dir figures/
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, Iterable, List, Optional

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import scanpy as sc
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install matplotlib seaborn numpy pandas scanpy anndata")

sns.set_style("ticks")
plt.rcParams["font.family"] = "sans-serif"


def _save(fig, base_path: str, dpi: int = 300) -> None:
    """Write PNG always, SVG if the backend allows; close the figure."""
    parent = os.path.dirname(base_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    stem = os.path.splitext(base_path)[0]
    png = f"{stem}.png"
    fig.savefig(png, dpi=dpi, bbox_inches="tight", facecolor="white")
    print(f"  → {png}")
    try:
        svg = f"{stem}.svg"
        fig.savefig(svg, format="svg", bbox_inches="tight", facecolor="white")
        print(f"  → {svg}")
    except Exception:
        print("  (SVG export unavailable; PNG saved)")
    plt.close(fig)


def visualize_single_perturbation(
    adata: "ad.AnnData",
    gene_name: str,
    gene_col: str = "gene",
    sgrna_col: str = "sgRNA",
    color_by: Optional[List[str]] = None,
    plot_types: Iterable[str] = ("umap", "pca"),
    output_dir: str = "figures/",
    show: bool = False,
) -> None:
    """UMAP and/or PCA for a single perturbation."""
    os.makedirs(output_dir, exist_ok=True)
    if color_by is None:
        color_by = ["classification", gene_col, sgrna_col]
    colors = list(color_by)
    if gene_name in adata.var_names and gene_name not in colors:
        colors.insert(2, gene_name)
    colors = [c for c in colors if c in adata.obs.columns or c in adata.var_names]

    if "umap" in plot_types and "X_umap" in adata.obsm:
        sc.pl.umap(adata, color=colors, save=f"_{gene_name}.pdf", show=show)
    if "pca" in plot_types and "X_pca" in adata.obsm:
        sc.pl.pca(adata, color=colors, save=f"_{gene_name}.pdf", show=show)


def visualize_all_perturbations(
    adata_dict: Dict[str, "ad.AnnData"],
    gene_col: str = "gene",
    sgrna_col: str = "sgRNA",
    color_by: Optional[List[str]] = None,
    plot_types: Iterable[str] = ("umap", "pca"),
    output_dir: str = "figures/",
    show: bool = False,
) -> None:
    """Iterate `visualize_single_perturbation` over every perturbation."""
    print(f"Visualizing {len(adata_dict)} perturbations → {output_dir}")
    for i, (gene, adata) in enumerate(adata_dict.items(), 1):
        if i % 10 == 0:
            print(f"  [{i}/{len(adata_dict)}] {gene}")
        visualize_single_perturbation(
            adata, gene, gene_col=gene_col, sgrna_col=sgrna_col,
            color_by=color_by, plot_types=plot_types, output_dir=output_dir,
            show=show,
        )
    print(f"✓ Visualization complete")


def plot_hit_summary(
    summary_df: "pd.DataFrame",
    output_path: str = "figures/hit_summary.png",
    top_n: int = 20,
) -> None:
    """Two-panel screen overview: top outlier-fractions + DE-vs-outlier scatter."""
    s = summary_df.sort_values("outlier_fraction", ascending=False).head(top_n)
    fig, ax = plt.subplots(1, 2, figsize=(14, 6))

    sns.barplot(data=s, x="outlier_fraction", y="gene", hue="is_hit", ax=ax[0], dodge=False)
    ax[0].set_xlabel("Outlier fraction")
    ax[0].set_ylabel("Perturbation")
    ax[0].set_title(f"Top {top_n} perturbations by outlier fraction")

    sns.scatterplot(
        data=summary_df, x="n_de_genes", y="n_outliers",
        hue="is_hit", size="n_cells", alpha=0.6, ax=ax[1],
    )
    ax[1].set_xlabel("DE genes")
    ax[1].set_ylabel("Outlier cells")
    ax[1].set_title("DE genes vs outliers")

    plt.tight_layout()
    _save(fig, output_path)
    print("✓ Hit summary saved")


def create_volcano_plots(
    de_results: "pd.DataFrame",
    perturbation: str,
    output_dir: str = "figures/",
    fdr_threshold: float = 0.05,
    log2fc_threshold: float = 0.5,
) -> None:
    """Volcano plot for a single perturbation."""
    os.makedirs(output_dir, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 6))

    sig_col = "qval" if "qval" in de_results.columns else "adj_pval"
    sig = (de_results[sig_col] < fdr_threshold) & (de_results["log2fc"].abs() > log2fc_threshold)
    neglog10 = -np.log10(de_results["pval"] + 1e-300)

    ax.scatter(
        de_results.loc[~sig, "log2fc"], neglog10[~sig],
        c="gray", alpha=0.5, s=10, label="ns",
    )
    ax.scatter(
        de_results.loc[sig, "log2fc"], neglog10[sig],
        c="crimson", alpha=0.8, s=12, label=f"FDR<{fdr_threshold}",
    )
    ax.axhline(-np.log10(fdr_threshold), c="steelblue", ls="--", lw=1, alpha=0.6)
    ax.axvline(log2fc_threshold, c="steelblue", ls="--", lw=1, alpha=0.6)
    ax.axvline(-log2fc_threshold, c="steelblue", ls="--", lw=1, alpha=0.6)
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10(p)")
    ax.set_title(f"Volcano — {perturbation}")
    ax.legend()
    plt.tight_layout()
    _save(fig, os.path.join(output_dir, f"volcano_{perturbation}.png"))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--de-csv", required=True, help="Per-perturbation DE CSV")
    parser.add_argument("--perturbation", required=True)
    parser.add_argument("--out-dir", default="figures/")
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--log2fc", type=float, default=0.5)
    args = parser.parse_args()

    df = pd.read_csv(args.de_csv)
    create_volcano_plots(
        df, args.perturbation, args.out_dir,
        fdr_threshold=args.fdr, log2fc_threshold=args.log2fc,
    )
    print("✓ Volcano plot generated")


if __name__ == "__main__":
    main()
