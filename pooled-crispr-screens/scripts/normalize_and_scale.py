"""Normalize, log-transform, and optionally scale a CRISPR-screen AnnData.

`normalize_only` produces log1p(CPM) suitable for the t-test screening tier.
`normalize_and_scale_data` additionally regresses out covariates and scales —
needed for PCA / UMAP.

Example:
    python normalize_and_scale.py --in filtered.h5ad --out normalized.h5ad \\
        --target-sum 1e6 --no-scale
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional

try:
    import scanpy as sc
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install scanpy anndata")


def normalize_only(
    adata: "ad.AnnData",
    target_sum: float = 1e6,
    exclude_highly_expressed: bool = True,
    save_raw: bool = False,
) -> "ad.AnnData":
    """CPM-style normalize and log1p; no regression, no scaling."""
    if save_raw:
        adata.raw = adata.copy()
    sc.pp.normalize_total(
        adata, target_sum=target_sum, exclude_highly_expressed=exclude_highly_expressed
    )
    sc.pp.log1p(adata)
    print(f"✓ Normalization complete: target_sum={target_sum:g}, log1p applied")
    return adata


def normalize_and_scale_data(
    adata: "ad.AnnData",
    target_sum: float = 1e6,
    exclude_highly_expressed: bool = True,
    regress_out: Optional[List[str]] = None,
    max_value: float = 10.0,
    save_raw: bool = True,
    n_jobs: int = 4,
) -> "ad.AnnData":
    """Full pipeline: normalize → log → save raw → regress → scale.

    After this, ``adata.X`` is scaled (for PCA/UMAP) and ``adata.raw`` holds
    the log-normalized layer (for DE and visualization).
    """
    print(f"[1/5] Normalize to {target_sum:g} per cell (exclude_highly_expressed={exclude_highly_expressed})")
    sc.pp.normalize_total(
        adata, target_sum=target_sum, exclude_highly_expressed=exclude_highly_expressed
    )
    print("[2/5] log1p")
    sc.pp.log1p(adata)
    if save_raw:
        print("[3/5] Snapshot log-normalized layer into adata.raw")
        adata.raw = adata
    else:
        print("[3/5] Skip raw snapshot (save_raw=False)")
    if regress_out:
        print(f"[4/5] Regress out: {regress_out}")
        sc.pp.regress_out(adata, regress_out, n_jobs=n_jobs)
    else:
        print("[4/5] Skip regression (regress_out=None)")
    print(f"[5/5] Scale with max_value={max_value}")
    sc.pp.scale(adata, max_value=max_value)
    print("✓ Normalization + scaling complete")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--target-sum", type=float, default=1e6)
    parser.add_argument("--no-scale", action="store_true",
                        help="Only normalize+log (skip regression and scaling)")
    parser.add_argument("--regress-out", nargs="*", default=None,
                        help="Columns to regress out (e.g. n_counts percent_mito)")
    parser.add_argument("--max-value", type=float, default=10.0)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.in_path)
    if args.no_scale:
        adata = normalize_only(adata, target_sum=args.target_sum)
    else:
        adata = normalize_and_scale_data(
            adata, target_sum=args.target_sum,
            regress_out=args.regress_out, max_value=args.max_value,
        )
    adata.write(args.out)
    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
