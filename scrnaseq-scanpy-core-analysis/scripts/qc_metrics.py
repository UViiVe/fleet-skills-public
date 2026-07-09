"""Compute QC metrics for scRNA-seq data and flag outliers.

Adds per-cell metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`,
`pct_counts_ribo`, log-transformed versions) to `adata.obs`. Optionally flags
outlier cells using either fixed tissue thresholds or batch-aware MAD.

Inputs:  AnnData with raw counts.
Outputs: AnnData with `.obs` columns added; CLI mode writes the result to `.h5ad`.
Example:
    python qc_metrics.py --input adata.h5ad --output adata_qc.h5ad \\
        --species human --tissue pbmc --strategy mad --nmads 5
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List, Optional, TYPE_CHECKING

try:
    import numpy as np
except ImportError:
    sys.exit("Install with: pip install numpy")

if TYPE_CHECKING:
    from anndata import AnnData


_MITO_PATTERN: Dict[str, str] = {"human": "MT-", "mouse": "mt-"}
_RIBO_PATTERN: Dict[str, str] = {"human": "^RP[SL]", "mouse": "^Rp[sl]"}


def get_species_mito_pattern(species: str) -> str:
    species = species.lower()
    if species not in _MITO_PATTERN:
        raise ValueError(f"species must be one of {list(_MITO_PATTERN)}")
    return _MITO_PATTERN[species]


def get_species_ribo_pattern(species: str) -> str:
    species = species.lower()
    if species not in _RIBO_PATTERN:
        raise ValueError(f"species must be one of {list(_RIBO_PATTERN)}")
    return _RIBO_PATTERN[species]


def calculate_qc_metrics(
    adata: "AnnData",
    species: str = "human",
    calculate_ribo: bool = True,
    inplace: bool = True,
) -> "AnnData":
    """Add basic QC metrics to `adata.obs`."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    if not inplace:
        adata = adata.copy()

    mito = get_species_mito_pattern(species)
    adata.var["mt"] = adata.var_names.str.startswith(mito)
    qc_vars = ["mt"]
    print(f"  Mitochondrial genes: {int(adata.var['mt'].sum())}")

    if calculate_ribo:
        ribo = get_species_ribo_pattern(species)
        adata.var["ribo"] = adata.var_names.str.match(ribo)
        qc_vars.append("ribo")
        print(f"  Ribosomal genes:     {int(adata.var['ribo'].sum())}")

    sc.pp.calculate_qc_metrics(adata, qc_vars=qc_vars, percent_top=None, log1p=False, inplace=True)

    adata.obs["log10_total_counts"] = np.log10(adata.obs["total_counts"] + 1)
    adata.obs["log10_n_genes_by_counts"] = np.log10(adata.obs["n_genes_by_counts"] + 1)
    adata.obs["log1p_total_counts"] = np.log1p(adata.obs["total_counts"])
    adata.obs["log1p_n_genes_by_counts"] = np.log1p(adata.obs["n_genes_by_counts"])

    print(f"✓ QC metrics computed for {adata.n_obs} cells")
    print(f"  median genes/cell: {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"  median UMIs/cell:  {adata.obs['total_counts'].median():.0f}")
    print(f"  median pct_mt:     {adata.obs['pct_counts_mt'].median():.2f}%")
    return adata


_TISSUE_THRESHOLDS: Dict[str, Dict[str, float]] = {
    "pbmc":   {"min_genes": 200, "max_genes": 2500, "max_mt": 5,  "description": "PBMC"},
    "brain":  {"min_genes": 200, "max_genes": 6000, "max_mt": 10, "description": "Brain"},
    "tumor":  {"min_genes": 200, "max_genes": 5000, "max_mt": 20, "description": "Tumor"},
    "kidney": {"min_genes": 200, "max_genes": 4000, "max_mt": 15, "description": "Kidney"},
    "liver":  {"min_genes": 200, "max_genes": 4000, "max_mt": 15, "description": "Liver"},
    "heart":  {"min_genes": 200, "max_genes": 4000, "max_mt": 18, "description": "Heart"},
    "lung":   {"min_genes": 200, "max_genes": 4000, "max_mt": 10, "description": "Lung"},
    "pancreas": {"min_genes": 200, "max_genes": 5000, "max_mt": 10, "description": "Pancreas"},
    "default": {"min_genes": 200, "max_genes": 4000, "max_mt": 10, "description": "Default"},
}


def get_tissue_qc_thresholds(tissue: str) -> Dict[str, float]:
    t = tissue.lower()
    if t not in _TISSUE_THRESHOLDS:
        print(f"Tissue {tissue!r} not recognized; using default thresholds.")
        t = "default"
    th = _TISSUE_THRESHOLDS[t]
    print(f"Tissue: {th['description']} — min_genes={th['min_genes']}, "
          f"max_genes={th['max_genes']}, max_mt={th['max_mt']}%")
    return th


def mark_outliers_fixed(
    adata: "AnnData",
    tissue: str = "pbmc",
    min_genes: Optional[int] = None,
    max_genes: Optional[int] = None,
    max_mt: Optional[float] = None,
    inplace: bool = True,
) -> "AnnData":
    """Flag outliers using tissue-specific fixed thresholds."""
    if not inplace:
        adata = adata.copy()

    th = get_tissue_qc_thresholds(tissue)
    min_g = min_genes if min_genes is not None else th["min_genes"]
    max_g = max_genes if max_genes is not None else th["max_genes"]
    mt_th = max_mt if max_mt is not None else th["max_mt"]

    low = adata.obs["n_genes_by_counts"] < min_g
    high = adata.obs["n_genes_by_counts"] > max_g
    mt = adata.obs["pct_counts_mt"] > mt_th
    adata.obs["outlier_low_genes"] = low
    adata.obs["outlier_high_genes"] = high
    adata.obs["outlier_high_mt"] = mt
    adata.obs["outlier"] = low | high | mt

    n_out = int(adata.obs["outlier"].sum())
    pct = 100 * n_out / adata.n_obs
    print(f"  low_genes:  {int(low.sum())}")
    print(f"  high_genes: {int(high.sum())}")
    print(f"  high_mt:    {int(mt.sum())}")
    print(f"✓ Fixed-threshold outliers: {n_out} ({pct:.1f}%) — keeping {adata.n_obs - n_out} cells")
    return adata


def batch_mad_outlier_detection(
    adata: "AnnData",
    batch_key: str = "batch",
    metrics: Optional[List[str]] = None,
    nmads: float = 5,
    inplace: bool = True,
) -> "AnnData":
    """Batch-aware MAD-based outlier detection on QC metrics."""
    if not inplace:
        adata = adata.copy()

    if metrics is None:
        metrics = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt"]

    if batch_key not in adata.obs.columns:
        print(f"  No '{batch_key}' column — treating as a single batch.")
        adata.obs[batch_key] = "batch_1"

    for m in ("log1p_total_counts", "log1p_n_genes_by_counts"):
        if m not in adata.obs.columns:
            adata.obs[m] = np.log1p(adata.obs[m.replace("log1p_", "")])

    print(f"Running MAD outlier detection (nmads={nmads}) on {len(metrics)} metrics...")
    adata.obs["outlier"] = False
    for metric in metrics:
        if metric not in adata.obs.columns:
            print(f"  skipping unknown metric {metric!r}")
            continue
        col = f"outlier_{metric}"
        adata.obs[col] = False
        for batch in adata.obs[batch_key].unique():
            mask = adata.obs[batch_key] == batch
            vals = adata.obs.loc[mask, metric].values
            med = float(np.median(vals))
            mad = float(np.median(np.abs(vals - med)))
            if mad == 0:
                continue
            lo, hi = med - nmads * mad, med + nmads * mad
            bad = (vals < lo) | (vals > hi)
            adata.obs.loc[mask, col] = bad
            adata.obs.loc[mask, "outlier"] |= bad
            if bad.any():
                print(f"  {metric} [{batch}]: {int(bad.sum())} outliers ({lo:.2f}, {hi:.2f})")

    n_out = int(adata.obs["outlier"].sum())
    pct = 100 * n_out / adata.n_obs
    print(f"✓ MAD outliers: {n_out} ({pct:.1f}%) — keeping {adata.n_obs - n_out} cells")
    return adata


def calculate_doublet_scores(
    adata: "AnnData",
    expected_doublet_rate: float = 0.06,
    random_state: int = 0,
) -> "AnnData":
    """Compute Scrublet doublet scores in place on `adata.X`."""
    try:
        import scrublet as scr
    except ImportError:
        sys.exit("Install with: pip install scrublet")

    print("Running Scrublet...")
    sc = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate, random_state=random_state)
    scores, predicted = sc.scrub_doublets(min_counts=2, min_cells=3,
                                           min_gene_variability_pctl=85, n_prin_comps=30)
    adata.obs["doublet_score"] = scores
    adata.obs["predicted_doublet"] = predicted
    n = int(predicted.sum())
    pct = 100 * n / len(predicted)
    print(f"✓ Predicted doublets: {n} ({pct:.1f}%)")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad input")
    parser.add_argument("--output", required=True, help=".h5ad output")
    parser.add_argument("--species", choices=["human", "mouse"], default="human")
    parser.add_argument("--tissue", default="pbmc")
    parser.add_argument("--strategy", choices=["mad", "fixed"], default="mad")
    parser.add_argument("--nmads", type=float, default=5)
    parser.add_argument("--batch-key", default="batch")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    adata = calculate_qc_metrics(adata, species=args.species)
    if args.strategy == "mad":
        adata = batch_mad_outlier_detection(adata, batch_key=args.batch_key, nmads=args.nmads)
    else:
        adata = mark_outliers_fixed(adata, tissue=args.tissue)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
