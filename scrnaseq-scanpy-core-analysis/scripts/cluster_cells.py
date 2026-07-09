"""Build a neighbor graph and cluster with Leiden (or Louvain).

Inputs:  AnnData with PCA or an integrated representation in `obsm`.
Outputs: `adata.obs['leiden_<res>']` columns (one per resolution).
Example:
    python cluster_cells.py --input adata_pca.h5ad --output adata_clust.h5ad \\
        --use-rep X_pca --n-pcs 30 --resolutions 0.4 0.6 0.8 1.0
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def build_neighbor_graph(
    adata: "AnnData",
    n_neighbors: int = 10,
    n_pcs: int = 30,
    use_rep: Optional[str] = "X_pca",
    metric: str = "euclidean",
    random_state: int = 0,
) -> "AnnData":
    """Compute the KNN graph used by clustering and UMAP."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    if use_rep == "X_pca":
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric=metric,
                        random_state=random_state)
    else:
        if use_rep not in adata.obsm:
            raise KeyError(f"`{use_rep}` not found in adata.obsm")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep, metric=metric,
                        random_state=random_state)
    print(f"✓ Neighbor graph: n_neighbors={n_neighbors}, use_rep={use_rep}")
    return adata


def cluster_leiden(
    adata: "AnnData",
    resolution: float = 0.8,
    key_added: Optional[str] = None,
    random_state: int = 42,
) -> "AnnData":
    """Single-resolution Leiden clustering."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")
    key_added = key_added or f"leiden_{resolution}"
    sc.tl.leiden(adata, resolution=resolution, key_added=key_added, random_state=random_state)
    n = adata.obs[key_added].nunique()
    print(f"✓ Leiden {resolution}: {n} clusters → {key_added}")
    return adata


def cluster_leiden_multiple_resolutions(
    adata: "AnnData",
    resolutions: Optional[List[float]] = None,
    random_state: int = 42,
) -> "AnnData":
    """Sweep multiple Leiden resolutions; store each as `leiden_<res>`."""
    resolutions = resolutions or [0.4, 0.6, 0.8, 1.0]
    print(f"Sweeping Leiden resolutions: {resolutions}")
    for r in resolutions:
        adata = cluster_leiden(adata, resolution=r, random_state=random_state)
    return adata


def cluster_louvain(
    adata: "AnnData",
    resolution: float = 0.8,
    key_added: Optional[str] = None,
    random_state: int = 42,
) -> "AnnData":
    """Single-resolution Louvain clustering (legacy)."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")
    key_added = key_added or f"louvain_{resolution}"
    sc.tl.louvain(adata, resolution=resolution, key_added=key_added, random_state=random_state)
    n = adata.obs[key_added].nunique()
    print(f"✓ Louvain {resolution}: {n} clusters → {key_added}")
    return adata


def calculate_cluster_qc_stats(adata: "AnnData", cluster_key: str = "leiden_0.8") -> "AnnData":
    """Summarize cluster sizes and median QC values."""
    try:
        import pandas as pd
    except ImportError:
        sys.exit("Install with: pip install pandas")
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"`{cluster_key}` not in adata.obs")
    metric_cols = [c for c in ("n_genes_by_counts", "total_counts", "pct_counts_mt") if c in adata.obs.columns]
    grouped = adata.obs.groupby(cluster_key)
    stats = grouped[metric_cols].median()
    stats["n_cells"] = grouped.size()
    print(stats.to_string())
    adata.uns[f"{cluster_key}_qc_stats"] = stats.to_dict()
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--use-rep", default="X_pca")
    parser.add_argument("--n-neighbors", type=int, default=10)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--resolutions", type=float, nargs="+", default=[0.4, 0.6, 0.8, 1.0])
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    adata = build_neighbor_graph(adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs,
                                 use_rep=args.use_rep)
    adata = cluster_leiden_multiple_resolutions(adata, resolutions=args.resolutions)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
