"""Load a built-in example dataset for testing the scanpy pipeline.

Inputs (CLI): dataset name (`pbmc3k` or `pbmc68k_reduced`).
Output: an in-memory `AnnData`; optional `.h5ad` snapshot.
Example:
    python load_example_data.py --dataset pbmc3k --output pbmc3k.h5ad
"""

from __future__ import annotations

import argparse
import sys
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


_AVAILABLE = ("pbmc3k", "pbmc68k_reduced")


def load_example_data(dataset: str = "pbmc3k") -> "AnnData":
    """Fetch a small example dataset distributed with scanpy."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    if dataset not in _AVAILABLE:
        raise ValueError(f"Unknown dataset {dataset!r}. Available: {_AVAILABLE}")

    print(f"Fetching example dataset: {dataset}")
    if dataset == "pbmc3k":
        adata = sc.datasets.pbmc3k()
        descr = "2,700 PBMCs (10X Chromium v1, healthy donor)"
    else:
        adata = sc.datasets.pbmc68k_reduced()
        descr = "Subsampled 68k PBMC reference"

    print(f"  {descr}")
    print(f"✓ Data loaded successfully! {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", choices=_AVAILABLE, default="pbmc3k")
    parser.add_argument("--output", help="Optional: write .h5ad")
    args = parser.parse_args()

    adata = load_example_data(args.dataset)
    if args.output:
        adata.write_h5ad(args.output)
        print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
