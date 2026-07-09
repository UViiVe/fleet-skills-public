"""Find cluster markers (Wilcoxon by default) and produce heatmaps/dot plots.

These statistics are exploratory cluster-defining tests — not valid for inferring
condition effects (use `pseudobulk_de.py` for that).

Inputs:  AnnData with `.X` = log-normalized counts and a cluster column.
Outputs: A DataFrame of ranked markers per cluster (long form), CSV exports, plots.
Example:
    python find_markers.py --input adata_clust.h5ad --output markers.csv \\
        --cluster-key leiden_0.8 --method wilcoxon
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def find_all_cluster_markers(
    adata: "AnnData",
    cluster_key: str = "leiden_0.8",
    method: str = "wilcoxon",
    min_logfoldchange: float = 0.25,
    min_pct: float = 0.1,
    n_genes: int = 200,
) -> pd.DataFrame:
    """Rank genes across all clusters and return a long-form DataFrame."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    if cluster_key not in adata.obs.columns:
        raise KeyError(f"`{cluster_key}` not in adata.obs")

    print(f"rank_genes_groups: groupby={cluster_key}, method={method}")
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method=method, n_genes=n_genes)

    rec = adata.uns["rank_genes_groups"]
    cats = rec["names"].dtype.names
    rows = []
    for cat in cats:
        for i in range(len(rec["names"][cat])):
            rows.append({
                "cluster": cat,
                "gene": rec["names"][cat][i],
                "score": float(rec["scores"][cat][i]),
                "logfoldchange": float(rec["logfoldchanges"][cat][i]),
                "pval": float(rec["pvals"][cat][i]),
                "pval_adj": float(rec["pvals_adj"][cat][i]),
            })
    df = pd.DataFrame(rows)

    df = df[df["logfoldchange"] >= min_logfoldchange].copy()
    print(f"✓ Markers identified: {df.shape[0]} rows across {df['cluster'].nunique()} clusters")
    return df


def find_markers_for_cluster(
    adata: "AnnData",
    cluster: str,
    cluster_key: str = "leiden_0.8",
    method: str = "wilcoxon",
    n_genes: int = 50,
) -> pd.DataFrame:
    """Find markers of `cluster` vs the rest."""
    df = find_all_cluster_markers(adata, cluster_key=cluster_key, method=method, n_genes=n_genes)
    sub = df[df["cluster"] == str(cluster)].copy().reset_index(drop=True)
    return sub


def export_marker_tables(markers: pd.DataFrame, output_dir: str = "results/markers") -> None:
    """Write one CSV with all markers and per-cluster CSVs."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    all_path = out / "cluster_markers_all.csv"
    markers.to_csv(all_path, index=False)
    for cluster, sub in markers.groupby("cluster"):
        sub.to_csv(out / f"markers_cluster_{cluster}.csv", index=False)
    print(f"✓ Exported markers to {out} ({markers.shape[0]} rows)")


def plot_top_markers_heatmap(
    adata: "AnnData",
    markers: pd.DataFrame,
    n_top: int = 10,
    cluster_key: str = "leiden_0.8",
    output_dir: str = "results/markers",
) -> None:
    """Plot a heatmap of top-N markers per cluster."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    top = (
        markers.sort_values(["cluster", "logfoldchange"], ascending=[True, False])
        .groupby("cluster")
        .head(n_top)
    )
    gene_list = top["gene"].drop_duplicates().tolist()

    sc.pl.heatmap(
        adata,
        var_names=gene_list,
        groupby=cluster_key,
        show=False,
        save=False,
        figsize=(12, 8),
    )
    out = Path(output_dir) / "marker_heatmap.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved {out}")


def plot_markers_dotplot(
    adata: "AnnData",
    markers: pd.DataFrame,
    n_top: int = 5,
    cluster_key: str = "leiden_0.8",
    output_dir: str = "results/markers",
) -> None:
    """Plot a dot plot of top-N markers per cluster."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    top = (
        markers.sort_values(["cluster", "logfoldchange"], ascending=[True, False])
        .groupby("cluster")
        .head(n_top)
    )
    gene_list = top["gene"].drop_duplicates().tolist()
    sc.pl.dotplot(adata, var_names=gene_list, groupby=cluster_key,
                  show=False, save=False, figsize=(10, 6))
    out = Path(output_dir) / "marker_dotplot.png"
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"✓ Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", default=None, help="CSV output for markers")
    parser.add_argument("--cluster-key", default="leiden_0.8")
    parser.add_argument("--method", default="wilcoxon")
    parser.add_argument("--plots", action="store_true", help="Also save heatmap + dot plot")
    parser.add_argument("--plot-dir", default="results/markers")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    markers = find_all_cluster_markers(adata, cluster_key=args.cluster_key, method=args.method)
    if args.output:
        markers.to_csv(args.output, index=False)
        print(f"✓ Wrote {args.output}: {markers.shape[0]} rows")
    if args.plots:
        plot_top_markers_heatmap(adata, markers, cluster_key=args.cluster_key,
                                  output_dir=args.plot_dir)
        plot_markers_dotplot(adata, markers, cluster_key=args.cluster_key,
                              output_dir=args.plot_dir)


if __name__ == "__main__":
    main()
