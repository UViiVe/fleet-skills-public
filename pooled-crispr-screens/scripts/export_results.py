"""Bundle CRISPR-screen results for downstream skills.

Writes per-perturbation stats, hit lists, the processed h5ad, a pickle for
direct hand-off to the next skill, and triggers PDF report generation.

Example:
    python export_results.py --in adata_processed.h5ad \\
        --summary perturbation_summary.csv --out-dir results/
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
from typing import Dict, List, Optional

try:
    import pandas as pd
    import anndata as ad
    import scanpy as sc
except ImportError:
    sys.exit("Install with: pip install pandas anndata scanpy")


def export_all_results(
    adata: "ad.AnnData",
    perturbation_summary: "pd.DataFrame",
    de_results: Optional[Dict[str, "pd.DataFrame"]] = None,
    outlier_cells: Optional[List[str]] = None,
    output_dir: str = "results/",
    export_formats: Optional[List[str]] = None,
) -> None:
    """Write the standard result bundle."""
    if export_formats is None:
        export_formats = ["csv", "h5ad", "excel"]
    os.makedirs(output_dir, exist_ok=True)
    print(f"Exporting to {output_dir}")

    if perturbation_summary is not None and "csv" in export_formats:
        perturbation_summary.to_csv(os.path.join(output_dir, "hit_summary.csv"), index=False)
        perturbation_summary.to_csv(os.path.join(output_dir, "per_perturbation_stats.csv"), index=False)
        print("  → hit_summary.csv, per_perturbation_stats.csv")

    if de_results and "excel" in export_formats:
        try:
            with pd.ExcelWriter(os.path.join(output_dir, "top_de_genes.xlsx"), engine="openpyxl") as w:
                for gene, df in de_results.items():
                    df.head(50).to_excel(w, sheet_name=gene[:31], index=False)
            print(f"  → top_de_genes.xlsx ({len(de_results)} sheets)")
        except ImportError:
            print("  (skipped xlsx — install openpyxl for Excel export)")

    if outlier_cells and "csv" in export_formats:
        pd.DataFrame({"cell_barcode": outlier_cells}).to_csv(
            os.path.join(output_dir, "outlier_cells.csv"), index=False
        )
        print(f"  → outlier_cells.csv ({len(outlier_cells)} rows)")

    if "h5ad" in export_formats:
        adata.write(os.path.join(output_dir, "adata_processed.h5ad"))
        print("  → adata_processed.h5ad")

    if perturbation_summary is not None and "is_hit" in perturbation_summary.columns:
        hits = perturbation_summary[perturbation_summary["is_hit"]]
        cols = [c for c in ["perturbation", "gene", "n_outliers", "outlier_fraction", "n_de_genes"]
                if c in hits.columns]
        if cols:
            hits[cols].to_csv(os.path.join(output_dir, "hit_list.csv"), index=False)
            print(f"  → hit_list.csv ({len(hits)} hits)")

    bundle = {
        "perturbation_summary": perturbation_summary,
        "de_results": de_results,
        "outlier_cells": outlier_cells,
    }
    pkl_path = os.path.join(output_dir, "analysis_objects.pkl")
    with open(pkl_path, "wb") as f:
        pickle.dump(bundle, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"  → analysis_objects.pkl  (load with `pickle.load(open('{pkl_path}','rb'))`)")

    try:
        from .generate_report import generate_report  # type: ignore
    except Exception:
        try:
            from generate_report import generate_report  # type: ignore
        except Exception:
            generate_report = None  # type: ignore
    if generate_report is not None:
        try:
            generate_report(
                adata=adata,
                perturbation_summary=perturbation_summary,
                de_results=de_results,
                output_dir=output_dir,
            )
        except Exception as exc:
            print(f"  PDF report skipped: {exc}")
    else:
        print("  PDF report skipped: generate_report not importable")

    print("✓ Export complete")


def export_hit_list(
    perturbation_summary: "pd.DataFrame",
    output_path: str = "results/hit_list.txt",
    sort_by: str = "outlier_fraction",
) -> None:
    """Simple plain-text hit list, one gene per line."""
    hits = perturbation_summary[perturbation_summary["is_hit"]].sort_values(sort_by, ascending=False)
    col = "gene" if "gene" in hits.columns else "perturbation"
    with open(output_path, "w") as f:
        for g in hits[col]:
            f.write(f"{g}\n")
    print(f"✓ Wrote {output_path} ({len(hits)} hits)")


def create_summary_report(
    adata: "ad.AnnData",
    perturbation_summary: "pd.DataFrame",
    de_results: Optional[Dict] = None,
    output_path: str = "results/analysis_report.txt",
) -> None:
    """Plain-text human-readable summary."""
    with open(output_path, "w") as f:
        f.write("Pooled CRISPR Screen — Analysis Report\n")
        f.write("=" * 60 + "\n\n")
        f.write("Dataset\n")
        f.write("-" * 60 + "\n")
        f.write(f"  cells:               {adata.n_obs:,}\n")
        f.write(f"  genes (expression):  {adata.n_vars:,}\n")
        if "gene" in adata.obs.columns:
            f.write(f"  perturbations:       {adata.obs['gene'].nunique()}\n")
        if "sgRNA" in adata.obs.columns:
            f.write(f"  sgRNAs:              {adata.obs['sgRNA'].nunique()}\n")
        f.write("\n")

        f.write("Hit Calling\n")
        f.write("-" * 60 + "\n")
        n_hits = int(perturbation_summary["is_hit"].sum())
        n_tested = len(perturbation_summary)
        f.write(f"  perturbations tested: {n_tested}\n")
        f.write(f"  hits:                 {n_hits} ({n_hits/n_tested*100:.1f}%)\n")
        f.write(f"  outlier cells total:  {int(perturbation_summary['n_outliers'].sum())}\n\n")

        f.write("Top 10 Hits (by outlier fraction)\n")
        f.write("-" * 60 + "\n")
        top = perturbation_summary.nlargest(10, "outlier_fraction")
        col = "gene" if "gene" in top.columns else "perturbation"
        for _, r in top.iterrows():
            f.write(
                f"  {r[col]:20s} | outliers: {int(r['n_outliers']):4d} "
                f"({r['outlier_fraction']:.1%}) | DE: {int(r['n_de_genes']):4d}\n"
            )
    print(f"✓ Wrote {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True, help="Processed h5ad")
    parser.add_argument("--summary", required=True, help="Perturbation summary CSV")
    parser.add_argument("--out-dir", default="results/")
    parser.add_argument("--outliers-csv", default=None)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.in_path)
    summary = pd.read_csv(args.summary)
    outliers = None
    if args.outliers_csv and os.path.exists(args.outliers_csv):
        outliers = pd.read_csv(args.outliers_csv)["cell_barcode"].tolist()

    export_all_results(
        adata=adata,
        perturbation_summary=summary,
        outlier_cells=outliers,
        output_dir=args.out_dir,
    )
    create_summary_report(adata, summary, output_path=os.path.join(args.out_dir, "analysis_report.txt"))


if __name__ == "__main__":
    main()
