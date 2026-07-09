"""Load 10X Feature-Barcode matrices for pooled CRISPR screens.

Reads one or more `raw_feature_bc_matrix.h5` files produced by CellRanger with
Feature Barcoding enabled (gene expression + sgRNA UMI counts in the same
matrix) and returns AnnData objects with deduplicated gene names.

Example:
    python load_10x_libraries.py --h5 lib1.h5 lib2.h5 --out adata_list.h5ad
"""

from __future__ import annotations

import argparse
import sys
from typing import List

try:
    import scanpy as sc
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install scanpy anndata")


def load_library(h5_path: str) -> "ad.AnnData":
    """Load a single 10X h5 file and dedupe variable names."""
    print(f"Reading {h5_path}")
    adata = sc.read_10x_h5(h5_path)
    adata.var_names_make_unique()
    print(f"  → {adata.n_obs} cells × {adata.n_vars} features")
    return adata


def load_multiple_libraries(h5_paths: List[str]) -> List["ad.AnnData"]:
    """Load a list of 10X h5 files, returning one AnnData per library."""
    out = [load_library(p) for p in h5_paths]
    print(f"✓ Libraries loaded: {len(out)}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--h5", nargs="+", required=True,
                        help="Paths to one or more raw_feature_bc_matrix.h5 files")
    parser.add_argument("--out-prefix", default="library",
                        help="Output filename prefix; one h5ad per library")
    args = parser.parse_args()

    adata_list = load_multiple_libraries(args.h5)
    for i, adata in enumerate(adata_list, 1):
        path = f"{args.out_prefix}_{i}.h5ad"
        adata.write(path)
        print(f"✓ Wrote {path}")


if __name__ == "__main__":
    main()
