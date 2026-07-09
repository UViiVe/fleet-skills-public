"""Compute UMAP (and optionally t-SNE / diffusion map) embeddings.

Inputs:  AnnData with a neighbor graph (see `cluster_cells.py`).
Outputs: `adata.obsm['X_umap']` (and `X_tsne` / `X_diffmap` if requested).
Example:
    python run_umap.py --input adata_clust.h5ad --output adata_umap.h5ad --min-dist 0.5
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def run_umap_reduction(
    adata: "AnnData",
    n_neighbors: Optional[int] = None,
    min_dist: float = 0.5,
    spread: float = 1.0,
    random_state: int = 42,
) -> "AnnData":
    """Compute UMAP. Always set `random_state` — UMAP is stochastic."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    if "neighbors" not in adata.uns:
        if n_neighbors is None:
            raise ValueError("No neighbor graph; pass n_neighbors or build one first.")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors)

    sc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=random_state)
    print(f"✓ UMAP computed (min_dist={min_dist}, spread={spread})")
    return adata


def run_tsne_reduction(
    adata: "AnnData",
    n_pcs: int = 30,
    perplexity: float = 30,
    random_state: int = 42,
) -> "AnnData":
    """Compute t-SNE (slower than UMAP; use for comparison only)."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")
    sc.tl.tsne(adata, n_pcs=n_pcs, perplexity=perplexity, random_state=random_state)
    print(f"✓ t-SNE computed (perplexity={perplexity})")
    return adata


def run_diffmap(adata: "AnnData", n_comps: int = 15) -> "AnnData":
    """Compute diffusion map (input for some trajectory inference methods)."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")
    sc.tl.diffmap(adata, n_comps=n_comps)
    print(f"✓ Diffusion map computed (n_comps={n_comps})")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--min-dist", type=float, default=0.5)
    parser.add_argument("--spread", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--tsne", action="store_true", help="Also compute t-SNE")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    adata = run_umap_reduction(adata, min_dist=args.min_dist, spread=args.spread,
                                random_state=args.seed)
    if args.tsne:
        adata = run_tsne_reduction(adata, random_state=args.seed)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
