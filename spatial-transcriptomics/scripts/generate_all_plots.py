"""Visualizations for the spatial transcriptomics workflow.

Emits PNG (always) and SVG (best-effort) for the following figures:

    1. QC violins
    2. Spatial clusters overlaid on tissue
    3. Spatial marker-gene panels
    4. UMAP coloured by cluster
    5. Neighborhood enrichment clustermap
    6. Co-occurrence curves
    7. Top SVGs bar chart
    8. Spatial map of the top SVG

Inputs : processed AnnData from spatial_workflow.py
Outputs: PNG/SVG files in `output_dir`

Example:
    python generate_all_plots.py --input visium_results/adata_processed.h5ad \\
        --output-dir visium_results --markers MYH7 TNNI3 MYH6
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Union

try:
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns
    from plotnine import (
        ggplot, aes, geom_violin, geom_point, geom_col, geom_line,
        labs, theme, coord_flip, facet_wrap, element_text, element_blank,
        guides, guide_legend,
    )
    from plotnine_prism import theme_prism
except ImportError as exc:  # pragma: no cover - graceful CLI exit
    sys.exit(
        "Install with: pip install scanpy squidpy plotnine plotnine-prism "
        "seaborn matplotlib  (missing: %s)" % exc.name
    )


# ---------------------------------------------------------------------------
# Save helpers
# ---------------------------------------------------------------------------

def _save_figure(
    fig,
    base_path: Union[str, Path],
    *,
    width: float = 8,
    height: float = 6,
    dpi: int = 300,
) -> None:
    """Write `fig` as PNG and SVG with graceful SVG fallback.

    Works for both matplotlib figures (has ``savefig``) and plotnine plots
    (has ``save``).
    """
    base = Path(base_path).with_suffix("")
    png = base.with_suffix(".png")
    svg = base.with_suffix(".svg")

    try:
        if hasattr(fig, "save"):
            fig.save(str(png), width=width, height=height, dpi=dpi, verbose=False)
        else:
            fig.savefig(str(png), dpi=dpi, bbox_inches="tight", format="png")
        print(f"   Saved: {png}")
    except Exception as exc:  # noqa: BLE001
        print(f"   PNG export failed: {exc}")

    try:
        if hasattr(fig, "save"):
            fig.save(str(svg), width=width, height=height, verbose=False)
        else:
            fig.savefig(str(svg), bbox_inches="tight", format="svg")
        print(f"   Saved: {svg}")
    except Exception:  # noqa: BLE001
        print("   (SVG export failed)")


def _save_and_close(fig, base_path, **kwargs) -> None:
    _save_figure(fig, base_path, **kwargs)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Individual plots
# ---------------------------------------------------------------------------

def _plot_qc_violins(adata, output_dir: Path) -> None:
    qc = pd.DataFrame({
        "Total Counts":    adata.obs["total_counts"].values,
        "Genes Detected":  adata.obs["n_genes_by_counts"].values,
        "Mitochondrial %": adata.obs["pct_counts_mt"].values,
    }).melt(var_name="Metric", value_name="Value")

    p = (
        ggplot(qc, aes(x="Metric", y="Value"))
        + geom_violin(fill="#4C72B0", alpha=0.7)
        + facet_wrap("~Metric", scales="free")
        + labs(title="QC Metrics", x="", y="Value")
        + theme_prism(base_size=12)
        + theme(
            plot_title=element_text(hjust=0.5, face="bold", size=14),
            axis_text_x=element_blank(),
            strip_text=element_text(size=11, face="bold"),
        )
    )
    _save_figure(p, output_dir / "qc_violins.png", width=12, height=5)
    print("  ✓ QC violin plots saved")


def _plot_spatial_clusters(adata, output_dir: Path) -> None:
    import scanpy as sc
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.spatial(
        adata, color="leiden", ax=ax, show=False,
        title="Spatial Clusters (Leiden)", frameon=False, spot_size=1.5,
    )
    _save_and_close(fig, output_dir / "spatial_clusters.png", width=8, height=8)
    print("  ✓ spatial cluster plot saved")


def _plot_spatial_markers(adata, marker_genes: List[str], output_dir: Path) -> None:
    import scanpy as sc

    n_genes = len(marker_genes)
    ncols = min(3, n_genes)
    nrows = (n_genes + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5.5 * nrows))
    if n_genes == 1:
        axes = np.array([axes])
    axes = np.asarray(axes).flatten()

    for i, gene in enumerate(marker_genes):
        sc.pl.spatial(
            adata, color=gene, ax=axes[i], show=False, title=gene,
            frameon=False, spot_size=1.5, cmap="Reds",
        )
    for j in range(n_genes, len(axes)):
        axes[j].set_visible(False)

    fig.suptitle("Marker Gene Expression", fontsize=14, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save_and_close(
        fig,
        output_dir / "spatial_markers.png",
        width=6 * ncols,
        height=5.5 * nrows,
    )
    print(f"  ✓ marker plots saved ({n_genes} genes)")


def _plot_umap_clusters(adata, output_dir: Path) -> None:
    umap = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
    umap["Cluster"] = adata.obs["leiden"].values

    p = (
        ggplot(umap, aes(x="UMAP1", y="UMAP2", color="Cluster"))
        + geom_point(size=0.5, alpha=0.7)
        + labs(title="UMAP (Leiden Clusters)", color="Cluster")
        + theme_prism(base_size=12)
        + theme(plot_title=element_text(hjust=0.5, face="bold", size=14))
        + guides(color=guide_legend(override_aes={"size": 3}))
    )
    _save_figure(p, output_dir / "umap_clusters.png", width=9, height=7)
    print("  ✓ UMAP cluster plot saved")


def _plot_nhood_enrichment(adata, output_dir: Path) -> None:
    zscore = adata.uns["leiden_nhood_enrichment"]["zscore"]
    clusters = sorted(adata.obs["leiden"].unique(), key=int)
    labels = [str(c) for c in clusters]
    z_df = pd.DataFrame(zscore, index=labels, columns=labels)

    g = sns.clustermap(
        z_df,
        cmap="RdBu_r",
        center=0,
        figsize=(10, 9),
        cbar_kws={"label": "Enrichment Z-score"},
        dendrogram_ratio=0.1,
        linewidths=0.5,
        linecolor="white",
    )
    g.ax_heatmap.set_xlabel("Cluster", fontsize=12)
    g.ax_heatmap.set_ylabel("Cluster", fontsize=12)
    g.fig.suptitle("Neighborhood Enrichment", fontsize=14, fontweight="bold", y=1.02)

    _save_and_close(g.fig, output_dir / "neighborhood_enrichment.png", width=10, height=9)
    print("  ✓ neighborhood enrichment heatmap saved")


def _plot_co_occurrence(adata, output_dir: Path) -> None:
    info = adata.uns["leiden_co_occurrence"]
    occ = info["occ"]
    interval = info["interval"]
    midpoints = [(interval[i] + interval[i + 1]) / 2 for i in range(len(interval) - 1)]
    clusters = sorted(adata.obs["leiden"].unique(), key=int)

    rows = []
    for i in range(min(len(clusters), 8)):
        for j, d in enumerate(midpoints):
            rows.append({
                "Distance": d,
                "Co-occurrence": occ[i, i, j],
                "Cluster": str(clusters[i]),
            })
    df = pd.DataFrame(rows)

    p = (
        ggplot(df, aes(x="Distance", y="Co-occurrence", color="Cluster"))
        + geom_line(size=1.2, alpha=0.8)
        + labs(
            title="Cluster Co-occurrence by Distance",
            x="Distance",
            y="Co-occurrence Probability",
            color="Cluster",
        )
        + theme_prism(base_size=12)
        + theme(plot_title=element_text(hjust=0.5, face="bold", size=14))
    )
    _save_figure(p, output_dir / "co_occurrence.png", width=10, height=6)
    print("  ✓ co-occurrence plot saved")


def _plot_top_svgs(svg_results: pd.DataFrame, top_n: int, output_dir: Path) -> None:
    top = svg_results.head(top_n).copy()
    top["Gene"] = pd.Categorical(top.index, categories=top.index[::-1], ordered=True)

    p = (
        ggplot(top, aes(x="Gene", y="I"))
        + geom_col(fill="#4C72B0", alpha=0.85)
        + coord_flip()
        + labs(title=f"Top {top_n} Spatially Variable Genes", x="", y="Moran's I")
        + theme_prism(base_size=12)
        + theme(plot_title=element_text(hjust=0.5, face="bold", size=14))
    )
    _save_figure(p, output_dir / "top_svgs.png", width=8, height=6)
    print("  ✓ top SVGs bar plot saved")


def _plot_spatial_svg(adata, gene: str, output_dir: Path) -> None:
    import scanpy as sc
    fig, ax = plt.subplots(figsize=(8, 8))
    sc.pl.spatial(
        adata, color=gene, ax=ax, show=False, title=f"{gene} (Top SVG)",
        frameon=False, spot_size=1.5, cmap="magma",
    )
    _save_and_close(fig, output_dir / f"spatial_svg_{gene}.png", width=8, height=8)
    print(f"  ✓ spatial SVG plot saved ({gene})")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def generate_all_plots(
    adata: "anndata.AnnData",
    output_dir: str = "visium_results",
    marker_genes: Optional[List[str]] = None,
    top_n_svgs: int = 10,
) -> None:
    """Generate all spatial transcriptomics visualizations.

    Parameters
    ----------
    adata
        Processed AnnData from ``spatial_workflow.run_spatial_analysis``.
    output_dir
        Directory for PNG/SVG output.
    marker_genes
        Genes to overlay on tissue (defaults to cardiac panel).
    top_n_svgs
        Number of top SVGs for the bar chart.
    """
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("\n=== Generating visualizations ===\n")

    if marker_genes is None:
        marker_genes = ["MYH7", "TNNI3", "MYH6", "NPPA", "TTN", "ACTN2"]
    available = [g for g in marker_genes if g in adata.var_names]

    stages = [
        ("QC violins",                   lambda: _plot_qc_violins(adata, out)),
        ("spatial clusters",              lambda: _plot_spatial_clusters(adata, out)),
        ("spatial markers",
            (lambda: _plot_spatial_markers(adata, available, out))
            if available else None),
        ("UMAP clusters",                 lambda: _plot_umap_clusters(adata, out)),
        ("neighborhood enrichment",       lambda: _maybe_nhood(adata, out)),
        ("co-occurrence",                 lambda: _maybe_cooccurrence(adata, out)),
        ("top SVGs / spatial SVG",
            lambda: _maybe_svg_plots(adata, out, top_n_svgs)),
    ]

    for label, fn in stages:
        if fn is None:
            print(f"  - skipping {label} (no matching markers)")
            continue
        try:
            fn()
        except Exception as exc:  # noqa: BLE001
            print(f"  ✗ {label} failed: {exc}")

    print("\n" + "=" * 50)
    print("✓ All visualizations generated successfully!")
    print("=" * 50)


def _maybe_nhood(adata, out):
    if "leiden_nhood_enrichment" in adata.uns:
        _plot_nhood_enrichment(adata, out)
    else:
        print("  - neighborhood enrichment skipped (results not present)")


def _maybe_cooccurrence(adata, out):
    if "leiden_co_occurrence" in adata.uns:
        _plot_co_occurrence(adata, out)
    else:
        print("  - co-occurrence skipped (results not present)")


def _maybe_svg_plots(adata, out, top_n):
    svg = adata.uns.get("svg_results")
    if svg is None or len(svg) == 0:
        print("  - SVG plots skipped (no significant SVGs)")
        return
    _plot_top_svgs(svg, top_n, out)
    top_gene = svg.index[0]
    if top_gene in adata.var_names:
        _plot_spatial_svg(adata, top_gene, out)
    else:
        print(f"  - top SVG spatial plot skipped ('{top_gene}' missing from var_names)")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate spatial transcriptomics visualizations."
    )
    parser.add_argument("--input", required=True, help="Processed AnnData (.h5ad)")
    parser.add_argument("--output-dir", default="visium_results")
    parser.add_argument("--markers", nargs="*", default=None,
                        help="Optional list of marker genes to overlay on tissue.")
    parser.add_argument("--top-n-svgs", type=int, default=10)
    args = parser.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    adata = sc.read_h5ad(args.input)
    generate_all_plots(
        adata,
        output_dir=args.output_dir,
        marker_genes=args.markers,
        top_n_svgs=args.top_n_svgs,
    )


if __name__ == "__main__":
    main()
