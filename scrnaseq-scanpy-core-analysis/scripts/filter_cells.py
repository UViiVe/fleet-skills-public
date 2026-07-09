"""Filter cells by QC criteria, run doublet detection, and remove genes/cells.

Inputs:  AnnData with QC metrics already calculated (see `qc_metrics.py`).
Outputs: A filtered AnnData with `predicted_doublet` and `doublet_score` columns
         added when Scrublet runs.
Example:
    python filter_cells.py --input adata_qc.h5ad --output adata_filt.h5ad \\
        --remove-doublets --doublet-threshold 0.25
"""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional, TYPE_CHECKING

try:
    import numpy as np
except ImportError:
    sys.exit("Install with: pip install numpy")

if TYPE_CHECKING:
    from anndata import AnnData


def _estimate_doublet_rate(n_cells: int) -> float:
    """10X cell-loading expected doublet rate ~ 0.8% per 1000 cells."""
    return min(0.10, max(0.005, n_cells / 1000 * 0.008))


def filter_cells_by_qc(
    adata: "AnnData",
    min_genes: int = 200,
    max_genes: Optional[int] = None,
    max_pct_mt: Optional[float] = None,
) -> "AnnData":
    """Apply fixed gene-count and MT% thresholds."""
    n0 = adata.n_obs
    mask = adata.obs["n_genes_by_counts"] >= min_genes
    if max_genes is not None:
        mask &= adata.obs["n_genes_by_counts"] <= max_genes
    if max_pct_mt is not None:
        mask &= adata.obs["pct_counts_mt"] <= max_pct_mt
    adata = adata[mask, :].copy()
    print(f"✓ filter_cells_by_qc: kept {adata.n_obs}/{n0} ({100*adata.n_obs/n0:.1f}%)")
    return adata


def filter_cells_by_tissue(
    adata: "AnnData",
    tissue: str = "pbmc",
    remove_doublets: bool = True,
) -> "AnnData":
    """Apply tissue-appropriate thresholds and optionally remove doublets."""
    from qc_metrics import get_tissue_qc_thresholds  # local import to avoid cycle

    th = get_tissue_qc_thresholds(tissue)
    adata = filter_cells_by_qc(
        adata,
        min_genes=int(th["min_genes"]),
        max_genes=int(th["max_genes"]),
        max_pct_mt=float(th["max_mt"]),
    )
    if remove_doublets and "predicted_doublet" in adata.obs.columns:
        n0 = adata.n_obs
        adata = adata[~adata.obs["predicted_doublet"].astype(bool), :].copy()
        print(f"✓ removed doublets: kept {adata.n_obs}/{n0} ({100*adata.n_obs/n0:.1f}%)")
    return adata


def filter_doublets(adata: "AnnData", score_threshold: Optional[float] = None) -> "AnnData":
    """Remove cells flagged by Scrublet."""
    if "predicted_doublet" not in adata.obs.columns:
        print("  no `predicted_doublet` column — run Scrublet first.")
        return adata
    n0 = adata.n_obs
    if score_threshold is not None and "doublet_score" in adata.obs.columns:
        keep = adata.obs["doublet_score"] < score_threshold
    else:
        keep = ~adata.obs["predicted_doublet"].astype(bool)
    adata = adata[keep, :].copy()
    print(f"✓ filter_doublets: kept {adata.n_obs}/{n0} ({100*adata.n_obs/n0:.1f}%)")
    return adata


def run_scrublet_detection(
    adata: "AnnData",
    batch_key: Optional[str] = "batch",
    expected_doublet_rate: Optional[float] = None,
    min_counts: int = 2,
    min_cells: int = 3,
    n_prin_comps: int = 30,
    random_state: int = 0,
) -> "AnnData":
    """Run Scrublet, optionally per-batch."""
    try:
        import scrublet as scr
    except ImportError:
        sys.exit("Install with: pip install scrublet")

    if batch_key is None or batch_key not in adata.obs.columns:
        rate = expected_doublet_rate if expected_doublet_rate is not None else _estimate_doublet_rate(adata.n_obs)
        print(f"Scrublet (single run, expected rate={rate:.3f})")
        s = scr.Scrublet(adata.X, expected_doublet_rate=rate, random_state=random_state)
        scores, predicted = s.scrub_doublets(min_counts=min_counts, min_cells=min_cells,
                                              min_gene_variability_pctl=85,
                                              n_prin_comps=n_prin_comps)
        adata.obs["doublet_score"] = scores
        adata.obs["predicted_doublet"] = predicted
    else:
        scores = np.zeros(adata.n_obs, dtype=float)
        predicted = np.zeros(adata.n_obs, dtype=bool)
        for batch in adata.obs[batch_key].unique():
            mask = (adata.obs[batch_key] == batch).values
            sub = adata[mask, :]
            rate = expected_doublet_rate if expected_doublet_rate is not None else _estimate_doublet_rate(sub.n_obs)
            print(f"  Scrublet [{batch}, n={sub.n_obs}, expected rate={rate:.3f}]")
            try:
                s = scr.Scrublet(sub.X, expected_doublet_rate=rate, random_state=random_state)
                sc_scores, sc_pred = s.scrub_doublets(min_counts=min_counts, min_cells=min_cells,
                                                       min_gene_variability_pctl=85,
                                                       n_prin_comps=n_prin_comps)
                scores[mask] = sc_scores
                predicted[mask] = sc_pred
            except Exception as exc:
                print(f"    Scrublet failed on batch {batch}: {exc}")
        adata.obs["doublet_score"] = scores
        adata.obs["predicted_doublet"] = predicted

    n = int(adata.obs["predicted_doublet"].sum())
    pct = 100 * n / adata.n_obs
    print(f"✓ Scrublet detected {n} doublets ({pct:.1f}%)")
    return adata


def filter_by_mad_outliers(
    adata: "AnnData",
    remove_doublets: bool = True,
    doublet_score_threshold: Optional[float] = None,
) -> "AnnData":
    """Drop cells flagged by `batch_mad_outlier_detection` and optionally doublets."""
    n0 = adata.n_obs
    if "outlier" not in adata.obs.columns:
        print("  no `outlier` column — run batch_mad_outlier_detection first.")
        return adata

    keep = ~adata.obs["outlier"].astype(bool)
    adata = adata[keep, :].copy()
    print(f"  removed outliers: {n0 - adata.n_obs}")

    if remove_doublets and "predicted_doublet" in adata.obs.columns:
        if doublet_score_threshold is not None:
            keep = adata.obs["doublet_score"] < doublet_score_threshold
        else:
            keep = ~adata.obs["predicted_doublet"].astype(bool)
        m0 = adata.n_obs
        adata = adata[keep, :].copy()
        print(f"  removed doublets: {m0 - adata.n_obs}")

    print(f"✓ filter_by_mad_outliers: kept {adata.n_obs}/{n0} ({100*adata.n_obs/n0:.1f}%)")
    return adata


def combine_filters_and_apply(adata: "AnnData", masks: List[str]) -> "AnnData":
    """Apply the OR of multiple boolean columns in `adata.obs` (kept = NOT any flag)."""
    n0 = adata.n_obs
    bad = np.zeros(adata.n_obs, dtype=bool)
    for col in masks:
        if col not in adata.obs.columns:
            print(f"  skipping missing column {col!r}")
            continue
        bad |= adata.obs[col].astype(bool).values
    adata = adata[~bad, :].copy()
    print(f"✓ combined filters: kept {adata.n_obs}/{n0}")
    return adata


def filter_genes(adata: "AnnData", min_cells: int = 3, min_counts: int = 0) -> "AnnData":
    """Drop genes expressed in fewer than `min_cells` cells."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    n0 = adata.n_vars
    if min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if min_counts > 0:
        sc.pp.filter_genes(adata, min_counts=min_counts)
    print(f"✓ filter_genes: kept {adata.n_vars}/{n0}")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad input")
    parser.add_argument("--output", required=True, help=".h5ad output")
    parser.add_argument("--scrublet", action="store_true", help="Run Scrublet first")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--remove-doublets", action="store_true")
    parser.add_argument("--doublet-threshold", type=float, default=None)
    parser.add_argument("--mode", choices=["mad", "tissue"], default="mad")
    parser.add_argument("--tissue", default="pbmc")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    if args.scrublet:
        adata = run_scrublet_detection(adata, batch_key=args.batch_key)
    if args.mode == "mad":
        adata = filter_by_mad_outliers(adata, remove_doublets=args.remove_doublets,
                                        doublet_score_threshold=args.doublet_threshold)
    else:
        adata = filter_cells_by_tissue(adata, tissue=args.tissue, remove_doublets=args.remove_doublets)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
