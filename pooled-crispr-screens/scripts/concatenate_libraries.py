"""Concatenate multiple per-library AnnData objects into one screen-level matrix.

Outer-joins genes across libraries so that minor manifest differences do not
silently truncate the matrix. Optionally adds batch labels.

Example:
    python concatenate_libraries.py --in lib1.h5ad lib2.h5ad lib3.h5ad \\
        --labels lib1 lib2 lib3 --out combined.h5ad
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional

try:
    import pandas as pd
    import scanpy as sc
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install pandas scanpy anndata")


def concatenate_libraries(
    adata_list: List["ad.AnnData"],
    batch_labels: Optional[List[str]] = None,
    batch_key: str = "batch",
) -> "ad.AnnData":
    """Outer-join AnnData objects with optional batch labels."""
    print(f"Joining {len(adata_list)} libraries (outer join)")
    for i, a in enumerate(adata_list):
        label = batch_labels[i] if batch_labels else f"batch_{i}"
        msg = f"  {label}: {a.n_obs} cells × {a.n_vars} genes"
        if "gene" in a.obs.columns:
            msg += f" | {a.obs['gene'].nunique()} perturbations"
        print(msg)

    if batch_labels:
        for a, label in zip(adata_list, batch_labels):
            a.obs[batch_key] = label

    combined = ad.concat(adata_list, join="outer")

    print(f"  combined: {combined.n_obs} cells × {combined.n_vars} genes")
    if "gene" in combined.obs.columns:
        cells_per = combined.obs.groupby("gene").size()
        print(
            f"  perturbations: {combined.obs['gene'].nunique()} | "
            f"cells/pert: mean {cells_per.mean():.1f}, median {cells_per.median():.1f}, "
            f"min {cells_per.min()}, max {cells_per.max()}"
        )
    print(f"✓ Libraries concatenated: {combined.n_obs} cells")
    return combined


def check_batch_balance(
    adata: "ad.AnnData", batch_key: str = "batch", gene_key: str = "gene"
) -> None:
    """Print a cross-batch coverage diagnostic per perturbation."""
    crosstab = pd.crosstab(adata.obs[gene_key], adata.obs[batch_key])
    n_all = (crosstab > 0).all(axis=1).sum()
    n_one = (crosstab > 0).sum(axis=1).eq(1).sum()
    print(
        f"Batch balance — {len(crosstab)} perturbations, "
        f"present in every batch: {n_all}, in only one: {n_one}"
    )
    cv = crosstab.std(axis=1) / crosstab.mean(axis=1).replace(0, 1)
    print("Top 10 batch-imbalanced perturbations (coefficient of variation):")
    print(cv.nlargest(10))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_paths", nargs="+", required=True,
                        help="Per-library h5ad paths")
    parser.add_argument("--labels", nargs="+", default=None,
                        help="Batch labels matching --in (default: numeric)")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    if args.labels and len(args.labels) != len(args.in_paths):
        sys.exit("Number of --labels must match number of --in libraries")

    adata_list = [sc.read_h5ad(p) for p in args.in_paths]
    combined = concatenate_libraries(adata_list, args.labels, args.batch_key)
    if "gene" in combined.obs.columns:
        check_batch_balance(combined, args.batch_key, "gene")
    combined.write(args.out)
    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
