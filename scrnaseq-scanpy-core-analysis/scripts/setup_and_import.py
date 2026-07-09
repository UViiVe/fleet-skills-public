"""Set up scanpy and load scRNA-seq data from common formats.

Inputs (CLI): a path to a 10X directory, .h5, .h5ad, or CSV/TSV matrix.
Outputs: an `AnnData` object in memory; CLI mode writes a `.h5ad` snapshot.
Example:
    python setup_and_import.py --input data/filtered_feature_bc_matrix/ \\
        --format 10x --output adata.h5ad
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path
from typing import Optional, Union

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")


def setup_scanpy_environment(verbose: bool = True) -> None:
    """Configure scanpy global settings and suppress noisy warnings."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    sc.settings.verbosity = 3
    sc.settings.set_figure_params(dpi=300, facecolor="white", frameon=False)

    warnings.filterwarnings("ignore", category=FutureWarning)
    warnings.filterwarnings("ignore", category=UserWarning)

    if verbose:
        print(f"scanpy   {sc.__version__}")
        print(f"numpy    {np.__version__}")
        print(f"pandas   {pd.__version__}")
        print("Environment ready.")


def import_10x_data(
    data_dir: Union[str, Path],
    var_names: str = "gene_symbols",
    cache: bool = False,
    min_cells: int = 3,
    min_genes: int = 200,
) -> "anndata.AnnData":
    """Load a CellRanger filtered/raw feature_bc_matrix directory."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    data_dir = Path(data_dir)
    if not data_dir.exists():
        raise FileNotFoundError(f"No such directory: {data_dir}")

    print(f"Reading 10X directory: {data_dir}")
    adata = sc.read_10x_mtx(data_dir, var_names=var_names, cache=cache)
    adata.var_names_make_unique()

    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=min_genes)

    print(f"✓ Data loaded successfully! {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def import_h5_data(
    h5_file: Union[str, Path],
    genome: Optional[str] = None,
    min_cells: int = 3,
    min_genes: int = 200,
) -> "anndata.AnnData":
    """Load a 10X-formatted `.h5` file."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    h5_file = Path(h5_file)
    if not h5_file.exists():
        raise FileNotFoundError(f"No such file: {h5_file}")

    print(f"Reading H5 file: {h5_file}")
    adata = sc.read_10x_h5(h5_file, genome=genome)
    adata.var_names_make_unique()

    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=min_genes)

    print(f"✓ Data loaded successfully! {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def import_count_matrix(
    file_path: Union[str, Path],
    transpose: bool = False,
    sep: Optional[str] = None,
    min_cells: int = 3,
    min_genes: int = 200,
) -> "anndata.AnnData":
    """Load a CSV/TSV count matrix. Default orientation is genes × cells."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"No such file: {file_path}")

    if sep is None:
        sep = "\t" if file_path.suffix.lower() == ".tsv" else ","

    print(f"Reading count matrix: {file_path} (sep={sep!r})")
    counts = pd.read_csv(file_path, sep=sep, index_col=0)
    if transpose:
        counts = counts.T

    # scanpy expects cells as rows, genes as columns.
    adata = sc.AnnData(counts.T)
    adata.var_names_make_unique()

    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=min_genes)

    print(f"✓ Data loaded successfully! {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def import_loom_data(
    loom_file: Union[str, Path],
    sparse: bool = True,
    cleanup: bool = True,
) -> "anndata.AnnData":
    """Load a `.loom` file (e.g. velocyto output)."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    loom_file = Path(loom_file)
    if not loom_file.exists():
        raise FileNotFoundError(f"No such file: {loom_file}")

    print(f"Reading Loom file: {loom_file}")
    adata = sc.read_loom(loom_file, sparse=sparse, cleanup=cleanup)
    print(f"✓ Data loaded successfully! {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def add_metadata(
    adata: "anndata.AnnData",
    metadata_file: Union[str, Path],
    merge_on: str = "index",
) -> "anndata.AnnData":
    """Join per-cell metadata from a CSV onto `adata.obs`."""
    metadata_file = Path(metadata_file)
    if not metadata_file.exists():
        raise FileNotFoundError(f"No such file: {metadata_file}")

    meta = pd.read_csv(metadata_file)
    if merge_on == "index" and meta.index.name is None:
        meta.set_index(meta.columns[0], inplace=True)

    common = set(adata.obs_names) & set(meta.index)
    if not common:
        raise ValueError("No matching barcodes between AnnData and metadata file.")

    for col in meta.columns:
        adata.obs[col] = meta.loc[adata.obs_names, col]
    print(f"✓ Added {len(meta.columns)} metadata columns to {len(common)} cells")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help="Path to input")
    parser.add_argument(
        "--format",
        choices=["10x", "h5", "csv", "tsv", "loom", "h5ad"],
        required=True,
    )
    parser.add_argument("--transpose", action="store_true", help="CSV/TSV: transpose if cells are rows")
    parser.add_argument("--output", help="Optional: write .h5ad snapshot")
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--min-genes", type=int, default=200)
    args = parser.parse_args()

    setup_scanpy_environment()

    if args.format == "10x":
        adata = import_10x_data(args.input, min_cells=args.min_cells, min_genes=args.min_genes)
    elif args.format == "h5":
        adata = import_h5_data(args.input, min_cells=args.min_cells, min_genes=args.min_genes)
    elif args.format in ("csv", "tsv"):
        adata = import_count_matrix(args.input, transpose=args.transpose,
                                    min_cells=args.min_cells, min_genes=args.min_genes)
    elif args.format == "loom":
        adata = import_loom_data(args.input)
    elif args.format == "h5ad":
        import scanpy as sc
        adata = sc.read_h5ad(args.input)

    if args.output:
        adata.write_h5ad(args.output)
        print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
