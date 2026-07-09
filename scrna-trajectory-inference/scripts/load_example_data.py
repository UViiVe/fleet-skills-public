"""Load the pancreatic endocrinogenesis demo dataset for trajectory inference.

Inputs:  none (the dataset is downloaded by scvelo if not cached).
Output:  AnnData with spliced/unspliced layers and a `clusters` column.
Example:
    python load_example_data.py --output pancreas.h5ad
"""

from __future__ import annotations

import argparse
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def load_example_data() -> "AnnData":
    """Return the pancreas endocrinogenesis demo dataset (3,696 cells)."""
    try:
        import scvelo as scv
    except ImportError:
        sys.exit("Install with: pip install scvelo")

    print("Fetching pancreatic endocrinogenesis demo dataset...")
    adata = scv.datasets.pancreas()
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    print(f"✓ Data loaded successfully! {adata.n_obs} cells × {adata.n_vars} genes")
    print(f"  clusters: {sorted(adata.obs['clusters'].unique().tolist())}")
    return adata


def load_user_data(path: str) -> "AnnData":
    """Load a user-provided `.h5ad`."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")
    adata = sc.read_h5ad(path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    print(f"✓ Data loaded successfully! {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", help="Path to write the demo .h5ad")
    args = parser.parse_args()

    adata = load_example_data()
    if args.output:
        adata.write_h5ad(args.output)
        print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
