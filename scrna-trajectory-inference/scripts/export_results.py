"""Export trajectory results: .h5ad, CSV tables, pickle, metadata.

Inputs:  AnnData + results dict from `run_trajectory_analysis`.
Outputs: `adata_trajectory.h5ad`, several CSVs, `trajectory_results.pkl`,
         `analysis_metadata.json`.
Example:
    python export_results.py --input adata.h5ad --results results.pkl \\
        --output-dir trajectory_results
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def export_all(
    adata: "AnnData",
    results: Dict[str, Any],
    output_dir: str = "trajectory_results",
) -> None:
    """One-shot export of every artifact."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("Exporting AnnData...")
    adata.write_h5ad(out / "adata_trajectory.h5ad", compression="lzf")
    print(f"  ✓ {out / 'adata_trajectory.h5ad'}")

    print("Exporting pseudotime assignments...")
    pt_df = pd.DataFrame({
        "barcode": adata.obs_names,
        "pseudotime": adata.obs.get("dpt_pseudotime", pd.Series([np.nan] * adata.n_obs)),
    })
    pt_col = results["parameters"]["cluster_key"]
    if pt_col in adata.obs.columns:
        pt_df["cell_type"] = adata.obs[pt_col].astype(str).values
    pt_df.to_csv(out / "pseudotime_assignments.csv", index=False)
    print(f"  ✓ {out / 'pseudotime_assignments.csv'} ({pt_df.shape[0]} cells)")

    print("Exporting trajectory genes...")
    if isinstance(results.get("trajectory_genes"), pd.DataFrame):
        results["trajectory_genes"].to_csv(out / "trajectory_genes.csv", index=False)
        print(f"  ✓ {out / 'trajectory_genes.csv'} ({results['trajectory_genes'].shape[0]} rows)")

    if results.get("velocity_results"):
        print("Exporting velocity genes...")
        try:
            vg = results["velocity_results"]["velocity_genes"]
            vg.to_csv(out / "velocity_genes.csv", index=False)
            print(f"  ✓ {out / 'velocity_genes.csv'}")
        except Exception as exc:
            print(f"  warning: velocity export failed: {exc}")

    if results.get("cellrank_results"):
        print("Exporting CellRank results...")
        try:
            fates = results["cellrank_results"]["fate_probabilities"]
            if fates is not None and not fates.empty:
                fates.to_csv(out / "fate_probabilities.csv")
                print(f"  ✓ {out / 'fate_probabilities.csv'}")
        except Exception as exc:
            print(f"  warning: fate probabilities export failed: {exc}")
        drivers = results["cellrank_results"].get("drivers", {})
        for fate, drv in drivers.items():
            try:
                safe = str(fate).replace("/", "_").replace(" ", "_")
                drv.to_csv(out / f"driver_genes_{safe}.csv")
                print(f"  ✓ {out / f'driver_genes_{safe}.csv'}")
            except Exception as exc:
                print(f"  warning: driver gene export failed for {fate}: {exc}")

    print("Writing trajectory_results.pkl...")
    with open(out / "trajectory_results.pkl", "wb") as f:
        pickle.dump(results, f)
    print(f"  ✓ {out / 'trajectory_results.pkl'}")

    print("Writing analysis_metadata.json...")
    metadata = {
        "timestamp": datetime.utcnow().isoformat() + "Z",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "parameters": results["parameters"],
        "trajectory_genes_count": int(len(results.get("trajectory_genes", []))),
        "has_velocity": results.get("velocity_results") is not None,
        "has_cellrank": results.get("cellrank_results") is not None,
    }
    if isinstance(results.get("pseudotime"), pd.Series):
        finite = results["pseudotime"][np.isfinite(results["pseudotime"])]
        metadata["pseudotime"] = {
            "min": float(finite.min()),
            "max": float(finite.max()),
            "mean": float(finite.mean()),
            "n_finite": int(finite.shape[0]),
        }
    with open(out / "analysis_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2, default=str)
    print(f"  ✓ {out / 'analysis_metadata.json'}")

    print("\n=== Export Complete ===")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad input")
    parser.add_argument("--results", required=True, help="results pickle")
    parser.add_argument("--output-dir", default="trajectory_results")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    with open(args.results, "rb") as f:
        results = pickle.load(f)
    export_all(adata, results, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
