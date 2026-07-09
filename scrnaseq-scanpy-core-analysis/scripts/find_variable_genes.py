"""Identify highly variable genes (HVGs) and optionally subset to them.

Inputs:  AnnData with raw counts in `layers['counts']`; log-normalized `.X`.
Outputs: `adata.var['highly_variable']` flag and dispersion stats.
Example:
    python find_variable_genes.py --input adata_norm.h5ad --output adata_hvg.h5ad \\
        --n-top-genes 2000 --flavor seurat_v3
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def find_highly_variable_genes(
    adata: "AnnData",
    n_top_genes: int = 2000,
    flavor: str = "seurat_v3",
    batch_key: Optional[str] = None,
    min_mean: float = 0.0125,
    max_mean: float = 3,
    min_disp: float = 0.5,
) -> "AnnData":
    """Select HVGs. `seurat_v3` is preferred; runs on raw counts in 'counts' layer."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    print(f"HVG selection: flavor={flavor}, n_top_genes={n_top_genes}, batch_key={batch_key}")

    kwargs = dict(n_top_genes=n_top_genes, flavor=flavor)
    if batch_key and batch_key in adata.obs.columns:
        kwargs["batch_key"] = batch_key

    if flavor == "seurat_v3":
        if "counts" not in adata.layers:
            print("  warning: 'counts' layer missing — Seurat v3 expects raw counts; using `.X`.")
        else:
            kwargs["layer"] = "counts"
    else:
        kwargs.update(min_mean=min_mean, max_mean=max_mean, min_disp=min_disp)

    sc.pp.highly_variable_genes(adata, **kwargs)
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"✓ Selected {n_hvg} HVGs")
    return adata


def filter_to_hvgs(adata: "AnnData") -> "AnnData":
    """Subset to highly-variable genes in place."""
    if "highly_variable" not in adata.var.columns:
        raise ValueError("Call find_highly_variable_genes first.")
    n0 = adata.n_vars
    adata = adata[:, adata.var["highly_variable"]].copy()
    print(f"✓ Subset to HVGs: {adata.n_vars}/{n0}")
    return adata


def plot_variable_genes(adata: "AnnData", output_dir: str = "results/hvg") -> None:
    """Visualize mean–dispersion HVG selection."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fig = sc.pl.highly_variable_genes(adata, show=False, save=False, return_fig=True)
    out = Path(output_dir) / "hvg_selection.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--flavor", default="seurat_v3")
    parser.add_argument("--batch-key", default=None)
    parser.add_argument("--subset", action="store_true", help="Subset adata to HVGs")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    adata = find_highly_variable_genes(adata, n_top_genes=args.n_top_genes,
                                        flavor=args.flavor, batch_key=args.batch_key)
    if args.subset:
        adata = filter_to_hvgs(adata)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
