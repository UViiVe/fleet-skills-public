"""Map sgRNA assignments to cells.

Joins a TSV of (cell_barcode, sgRNA_id) onto an AnnData, parses the target
gene out of the sgRNA ID, and keeps only mapped cells. Doublets should be
filtered upstream — this script does not attempt to call them.

Example:
    python map_sgrna_to_cells.py --in adata.h5ad --mapping mapping.tsv --out mapped.h5ad
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

try:
    import pandas as pd
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install pandas anndata")


def map_sgrna_to_adata(
    adata: "ad.AnnData",
    sgrna_mapping_file: str,
    sgrna_delimiter: str = "_",
    gene_position: int = 0,
) -> "ad.AnnData":
    """Attach sgRNA and target-gene labels to cells; drop unmapped cells.

    Parameters
    ----------
    adata : AnnData
    sgrna_mapping_file : str
        Path to a 2-column TSV (no header): cell_barcode<TAB>sgRNA_id
    sgrna_delimiter : str
        Character used to split sgRNA_id (default "_", giving e.g.
        ``GENE_sgRNA1`` -> ``GENE``).
    gene_position : int
        Index after split where the gene symbol lives (default 0).

    Returns
    -------
    AnnData with ``obs['sgRNA']`` and ``obs['gene']`` set, restricted to
    cells that appeared in the mapping file.
    """
    n_before = adata.n_obs
    mapping = pd.read_table(sgrna_mapping_file, header=None, index_col=0)
    print(f"Loaded {len(mapping)} cell-guide rows from {sgrna_mapping_file}")

    joined = adata.obs.merge(mapping, left_index=True, right_index=True, how="inner")
    adata = adata[adata.obs.index.isin(joined.index), :].copy()
    adata.obs["sgRNA"] = joined[1]
    adata.obs["gene"] = adata.obs["sgRNA"].apply(
        lambda s: s.split(sgrna_delimiter)[gene_position]
    )

    rate = adata.n_obs / n_before * 100
    print(
        f"  Mapping rate: {adata.n_obs}/{n_before} cells ({rate:.1f}%) | "
        f"{adata.obs['sgRNA'].nunique()} sgRNAs, "
        f"{adata.obs['gene'].nunique()} target genes"
    )
    print(f"✓ sgRNA mapping complete: {adata.n_obs} cells")
    return adata


def check_mapping_quality(adata: "ad.AnnData") -> "pd.DataFrame":
    """Per-gene summary: cells assigned and number of unique sgRNAs."""
    stats = pd.DataFrame(
        {
            "cells_per_gene": adata.obs.groupby("gene").size(),
            "sgrnas_per_gene": adata.obs.groupby("gene")["sgRNA"].nunique(),
        }
    )
    print(
        f"Cells per gene — mean {stats['cells_per_gene'].mean():.1f} | "
        f"median {stats['cells_per_gene'].median():.1f} | "
        f"<20 cells: {(stats['cells_per_gene'] < 20).sum()} | "
        f"<10 cells: {(stats['cells_per_gene'] < 10).sum()}"
    )
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True, help="Input h5ad")
    parser.add_argument("--mapping", required=True, help="TSV: barcode<TAB>sgRNA")
    parser.add_argument("--delimiter", default="_", help="sgRNA-ID delimiter (default '_')")
    parser.add_argument("--gene-position", type=int, default=0,
                        help="Index of gene symbol after splitting (default 0)")
    parser.add_argument("--out", required=True, help="Output h5ad")
    args = parser.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    adata = sc.read_h5ad(args.in_path)
    adata = map_sgrna_to_adata(adata, args.mapping, args.delimiter, args.gene_position)
    check_mapping_quality(adata)
    adata.write(args.out)
    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
