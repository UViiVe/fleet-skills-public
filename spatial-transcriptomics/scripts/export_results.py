"""Export spatial transcriptomics analysis results.

Writes the canonical bundle:

    1. adata_processed.h5ad   — processed AnnData (downstream consumer)
    2. spatially_variable_genes.csv
    3. cluster_assignments.csv
    4. neighborhood_enrichment.csv
    5. spot_metadata.csv
    6. analysis_summary.txt

Inputs : processed AnnData from spatial_workflow.run_spatial_analysis
Outputs: files listed above in `output_dir`

Example:
    python export_results.py --input visium_results/adata_processed.h5ad \\
        --output-dir visium_results
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path
from typing import Any

try:
    import numpy as np
    import pandas as pd
except ImportError:  # pragma: no cover
    sys.exit("Install with: pip install numpy pandas")


def export_all(
    adata: "anndata.AnnData",
    output_dir: str = "visium_results",
) -> None:
    """Write all spatial analysis exports under ``output_dir``."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("\n=== Exporting results ===\n")

    _export_h5ad(adata, out)
    _export_svgs(adata, out)
    _export_cluster_assignments(adata, out)
    _export_neighborhood_enrichment(adata, out)
    _export_spot_metadata(adata, out)
    _export_summary(adata, out)

    print("\n" + "=" * 50)
    print("=== Export Complete ===")
    print("=" * 50)
    print(f"\nAll results saved to: {out}/")
    print("\nKey outputs:")
    print("  - adata_processed.h5ad   (load with sc.read_h5ad)")
    print("  - spatially_variable_genes.csv")
    print("  - cluster_assignments.csv")
    print("  - neighborhood_enrichment.csv")
    print("  - spot_metadata.csv")
    print("  - analysis_summary.txt")


# ---------------------------------------------------------------------------
# Per-file writers
# ---------------------------------------------------------------------------

def _export_h5ad(adata, out: Path) -> None:
    print("1. Processed AnnData")
    path = out / "adata_processed.h5ad"
    adata_save = adata.copy()
    _coerce_uns_for_h5ad(adata_save)
    adata_save.write_h5ad(path, compression="gzip")
    mb = path.stat().st_size / (1024 ** 2)
    print(f"   ✓ {path} ({mb:.1f} MB)")
    print("   (load with: adata = sc.read_h5ad('adata_processed.h5ad'))")


def _export_svgs(adata, out: Path) -> None:
    print("\n2. Spatially variable genes")
    svgs = adata.uns.get("svg_results")
    if svgs is None or len(svgs) == 0:
        print("   (no SVG results to export)")
        return
    path = out / "spatially_variable_genes.csv"
    svgs.to_csv(path)
    print(f"   ✓ {path} ({len(svgs)} genes)")


def _export_cluster_assignments(adata, out: Path) -> None:
    print("\n3. Cluster assignments")
    df = pd.DataFrame({
        "barcode":   adata.obs_names,
        "cluster":   adata.obs["leiden"].values,
        "spatial_x": adata.obsm["spatial"][:, 0],
        "spatial_y": adata.obsm["spatial"][:, 1],
    })
    if "total_counts" in adata.obs.columns:
        df["total_counts"] = adata.obs["total_counts"].values
    if "n_genes_by_counts" in adata.obs.columns:
        df["n_genes"] = adata.obs["n_genes_by_counts"].values
    path = out / "cluster_assignments.csv"
    df.to_csv(path, index=False)
    print(f"   ✓ {path} ({len(df)} spots)")


def _export_neighborhood_enrichment(adata, out: Path) -> None:
    print("\n4. Neighborhood enrichment")
    if "leiden_nhood_enrichment" not in adata.uns:
        print("   (no neighborhood enrichment results)")
        return
    zscore = adata.uns["leiden_nhood_enrichment"]["zscore"]
    clusters = sorted(adata.obs["leiden"].unique(), key=int)
    labels = [str(c) for c in clusters]
    df = pd.DataFrame(zscore, index=labels, columns=labels)
    path = out / "neighborhood_enrichment.csv"
    df.to_csv(path)
    print(f"   ✓ {path}")


def _export_spot_metadata(adata, out: Path) -> None:
    print("\n5. Spot metadata")
    path = out / "spot_metadata.csv"
    adata.obs.to_csv(path)
    print(f"   ✓ {path} ({len(adata.obs)} spots, {len(adata.obs.columns)} columns)")


def _export_summary(adata, out: Path) -> None:
    print("\n6. Analysis summary")
    params = adata.uns.get("spatial_analysis_params", {})
    path = out / "analysis_summary.txt"

    lines = []
    lines.append("=" * 70)
    lines.append("SPATIAL TRANSCRIPTOMICS ANALYSIS SUMMARY")
    lines.append(f"Generated: {datetime.now():%Y-%m-%d %H:%M:%S}")
    lines.append("=" * 70)
    lines.append("")
    lines.append("DATASET")
    lines.append("-" * 70)
    lines.append(f"Spots (after QC):  {adata.n_obs}")
    lines.append(f"Genes:             {adata.n_vars}")
    lines.append(f"Clusters (Leiden): {params.get('n_clusters', 'N/A')}")
    lines.append(f"SVGs (FDR < 0.05): {params.get('n_svgs', 'N/A')}")
    lines.append("")
    lines.append("QC FILTERING")
    lines.append("-" * 70)
    n_init = int(params.get("n_spots_initial", 0))
    n_after = int(params.get("n_spots_after", 0))
    lines.append(f"Initial spots:     {n_init or 'N/A'}")
    lines.append(f"After filtering:   {n_after or 'N/A'}")
    if n_init:
        lines.append(f"Retention:         {100 * n_after / n_init:.1f}%")
    lines.append(f"min_genes:         {params.get('min_genes', 'N/A')}")
    lines.append(f"max_pct_mito:      {params.get('max_pct_mito', 'N/A')}%")
    lines.append("")
    lines.append("ANALYSIS PARAMETERS")
    lines.append("-" * 70)
    lines.append(f"HVGs:              {params.get('n_top_genes', 'N/A')}")
    lines.append(f"PCA components:    {params.get('n_pcs', 'N/A')}")
    lines.append(f"Leiden resolution: {params.get('resolution', 'N/A')}")
    lines.append(f"SVG permutations:  {params.get('svgs_n_perms', 'N/A')}")
    lines.append("")

    if "leiden" in adata.obs.columns:
        lines.append("CLUSTER DISTRIBUTION")
        lines.append("-" * 70)
        for cluster in sorted(adata.obs["leiden"].unique(), key=int):
            count = int((adata.obs["leiden"] == cluster).sum())
            pct = 100 * count / adata.n_obs
            lines.append(f"  Cluster {cluster:>3}: {count:>5} spots ({pct:>5.1f}%)")
        lines.append("")

    svgs = adata.uns.get("svg_results")
    if svgs is not None and len(svgs) > 0:
        lines.append("TOP SPATIALLY VARIABLE GENES")
        lines.append("-" * 70)
        for gene, row in svgs.head(20).iterrows():
            lines.append(
                f"  {gene:<15} Moran's I = {row['I']:.4f}  "
                f"FDR = {row['pval_norm_fdr_bh']:.2e}"
            )
        lines.append("")

    lines.append("=" * 70)
    lines.append("Pipeline: Squidpy + Scanpy spatial transcriptomics analysis")
    lines.append("=" * 70)

    path.write_text("\n".join(lines) + "\n")
    print(f"   ✓ {path}")


# ---------------------------------------------------------------------------
# Internal: clean .uns for h5ad serialization
# ---------------------------------------------------------------------------

def _coerce_uns_for_h5ad(adata) -> None:
    for key in list(adata.uns.keys()):
        val = adata.uns[key]
        if isinstance(val, dict):
            adata.uns[key] = _coerce_dict(val)
        elif isinstance(val, (np.integer, np.floating)):
            adata.uns[key] = val.item()


def _coerce_dict(d: dict) -> dict:
    out: dict = {}
    for k, v in d.items():
        if isinstance(v, np.integer):
            out[k] = int(v)
        elif isinstance(v, np.floating):
            out[k] = float(v)
        elif isinstance(v, dict):
            out[k] = _coerce_dict(v)
        else:
            out[k] = v
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export spatial transcriptomics analysis results."
    )
    parser.add_argument("--input", required=True, help="Processed AnnData (.h5ad)")
    parser.add_argument("--output-dir", default="visium_results")
    args = parser.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    adata = sc.read_h5ad(args.input)
    export_all(adata, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
