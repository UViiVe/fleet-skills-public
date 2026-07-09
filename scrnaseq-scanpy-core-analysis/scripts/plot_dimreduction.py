"""UMAP/PCA plots: cluster overlay, feature plots, batch overlay.

Inputs:  AnnData with `obsm['X_umap']` and clustering columns.
Outputs: PNG figures.
Example:
    python plot_dimreduction.py --input adata.h5ad --cluster-key leiden_0.8 \\
        --output-dir results/umap
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def _save_fig(fig, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"✓ Saved {out_path}")


def plot_umap_clusters(
    adata: "AnnData",
    cluster_key: str = "leiden_0.8",
    output_dir: str = "results/umap",
) -> None:
    """UMAP coloured by cluster."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")
    fig = sc.pl.umap(adata, color=cluster_key, legend_loc="on data",
                     legend_fontsize="x-small", legend_fontoutline=2,
                     show=False, return_fig=True)
    _save_fig(fig, Path(output_dir) / f"umap_{cluster_key}.png")
    plt.close(fig)


def plot_clustering_comparison(
    adata: "AnnData",
    resolutions: Optional[List[float]] = None,
    output_dir: str = "results/umap",
) -> None:
    """Side-by-side UMAP at multiple resolutions."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")
    resolutions = resolutions or [0.4, 0.6, 0.8, 1.0]
    columns = [f"leiden_{r}" for r in resolutions if f"leiden_{r}" in adata.obs.columns]
    if not columns:
        print("  no leiden_<res> columns found")
        return
    fig = sc.pl.umap(adata, color=columns, ncols=2, show=False, return_fig=True)
    _save_fig(fig, Path(output_dir) / "umap_resolution_comparison.png")
    plt.close(fig)


def plot_feature_umap(
    adata: "AnnData",
    features: List[str],
    output_dir: str = "results/umap",
    cmap: str = "Reds",
) -> None:
    """UMAP for each feature (gene or obs column) in `features`."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")
    valid = [f for f in features if f in adata.var_names or f in adata.obs.columns]
    if not valid:
        print("  none of the requested features were found")
        return
    fig = sc.pl.umap(adata, color=valid, ncols=3, cmap=cmap, show=False, return_fig=True)
    _save_fig(fig, Path(output_dir) / "umap_features.png")
    plt.close(fig)


def plot_umap_styled(
    adata: "AnnData",
    color_key: str,
    output_dir: str = "results/umap",
) -> None:
    """Custom styled UMAP via matplotlib (legend outside, no frame)."""
    try:
        import matplotlib.pyplot as plt
        import pandas as pd
        import seaborn as sns
    except ImportError:
        sys.exit("Install with: pip install matplotlib pandas seaborn")
    if "X_umap" not in adata.obsm:
        print("  no X_umap in obsm")
        return
    df = pd.DataFrame(adata.obsm["X_umap"], columns=["UMAP1", "UMAP2"])
    df[color_key] = adata.obs[color_key].astype(str).values
    fig, ax = plt.subplots(figsize=(8, 6))
    for v, sub in df.groupby(color_key):
        ax.scatter(sub["UMAP1"], sub["UMAP2"], s=2, alpha=0.7, label=v)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, markerscale=4, fontsize=8)
    sns.despine(ax=ax)
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    _save_fig(fig, Path(output_dir) / f"umap_styled_{color_key}.png")
    plt.close(fig)


def plot_embedding_density(
    adata: "AnnData",
    groupby: str,
    output_dir: str = "results/umap",
) -> None:
    """Density plot of UMAP coordinates split by `groupby`."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")
    if "X_umap" not in adata.obsm:
        return
    sc.tl.embedding_density(adata, basis="umap", groupby=groupby)
    fig = sc.pl.embedding_density(adata, basis="umap", key=f"umap_density_{groupby}",
                                  show=False, return_fig=True)
    _save_fig(fig, Path(output_dir) / f"umap_density_{groupby}.png")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--cluster-key", default="leiden_0.8")
    parser.add_argument("--features", nargs="*", default=[])
    parser.add_argument("--output-dir", default="results/umap")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    plot_umap_clusters(adata, cluster_key=args.cluster_key, output_dir=args.output_dir)
    if args.features:
        plot_feature_umap(adata, args.features, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
