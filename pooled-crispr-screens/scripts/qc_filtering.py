"""QC metrics and per-library filtering for pooled CRISPR screens.

Calculates gene count, total UMI, and mitochondrial fraction per cell, then
applies cell-type-configurable thresholds. Both dense and sparse matrices are
handled.

Example:
    python qc_filtering.py --in mapped.h5ad --out filtered.h5ad \\
        --min-genes 2000 --min-counts 8000 --max-mito 0.15
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

try:
    import numpy as np
    import scanpy as sc
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install numpy scanpy anndata")


def calculate_qc_metrics(adata: "ad.AnnData", mito_prefix: str = "MT-") -> "ad.AnnData":
    """Add n_genes (if absent), n_counts, and percent_mito to obs."""
    mito_mask = adata.var_names.str.startswith(mito_prefix)
    print(f"Mitochondrial genes ({mito_prefix}*): {int(mito_mask.sum())}")

    total = np.asarray(adata.X.sum(axis=1)).flatten()
    adata.obs["n_counts"] = total
    if mito_mask.sum() > 0:
        mito_counts = np.asarray(adata[:, mito_mask].X.sum(axis=1)).flatten()
        adata.obs["percent_mito"] = mito_counts / np.where(total > 0, total, 1)
    else:
        adata.obs["percent_mito"] = 0.0

    if "n_genes" not in adata.obs.columns:
        adata.obs["n_genes"] = np.asarray((adata.X > 0).sum(axis=1)).flatten()

    print(
        f"  n_genes:      mean {adata.obs['n_genes'].mean():.0f}, median {adata.obs['n_genes'].median():.0f}\n"
        f"  n_counts:     mean {adata.obs['n_counts'].mean():.0f}, median {adata.obs['n_counts'].median():.0f}\n"
        f"  percent_mito: mean {adata.obs['percent_mito'].mean():.3f}, median {adata.obs['percent_mito'].median():.3f}"
    )
    return adata


def apply_qc_filters(
    adata: "ad.AnnData",
    min_genes: int = 2000,
    min_counts: int = 8000,
    max_mito_pct: float = 0.15,
    calculate_metrics: bool = True,
    mito_prefix: str = "MT-",
) -> "ad.AnnData":
    """Drop low-quality cells using cell-count, UMI, and mito thresholds."""
    n_before = adata.n_obs
    if calculate_metrics:
        adata = calculate_qc_metrics(adata, mito_prefix=mito_prefix)

    print(f"QC thresholds: min_genes={min_genes}, min_counts={min_counts}, max_mito_pct={max_mito_pct}")

    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print(f"  dropped {n0 - adata.n_obs} cells with <{min_genes} genes")

    n0 = adata.n_obs
    sc.pp.filter_cells(adata, min_counts=min_counts)
    print(f"  dropped {n0 - adata.n_obs} cells with <{min_counts} counts")

    n0 = adata.n_obs
    adata = adata[adata.obs["percent_mito"] < max_mito_pct, :].copy()
    print(f"  dropped {n0 - adata.n_obs} cells with mito% >= {max_mito_pct}")

    retention = adata.n_obs / n_before * 100
    print(
        f"✓ QC filtering complete: {adata.n_obs}/{n_before} cells retained ({retention:.1f}%)"
    )
    return adata


def plot_qc_metrics(adata: "ad.AnnData", save_prefix: Optional[str] = None) -> None:
    """Violin + scatter plots for the three QC metrics."""
    sc.pl.violin(
        adata, ["n_genes", "n_counts", "percent_mito"],
        jitter=0.4, multi_panel=True,
        save=f"_{save_prefix}_violin.pdf" if save_prefix else None,
    )
    sc.pl.scatter(
        adata, x="n_counts", y="percent_mito",
        save=f"_{save_prefix}_counts_mito.pdf" if save_prefix else None,
    )
    sc.pl.scatter(
        adata, x="n_counts", y="n_genes",
        save=f"_{save_prefix}_counts_genes.pdf" if save_prefix else None,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--min-genes", type=int, default=2000)
    parser.add_argument("--min-counts", type=int, default=8000)
    parser.add_argument("--max-mito", type=float, default=0.15)
    parser.add_argument("--mito-prefix", default="MT-",
                        help="Gene-symbol prefix for mitochondrial genes (default MT-, mouse mt-)")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.in_path)
    adata = apply_qc_filters(
        adata,
        min_genes=args.min_genes,
        min_counts=args.min_counts,
        max_mito_pct=args.max_mito,
        mito_prefix=args.mito_prefix,
    )
    adata.write(args.out)
    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
