"""Export analysis artifacts: .h5ad + CSVs + summary text.

Inputs:  Annotated AnnData.
Outputs: `adata_processed.h5ad`, expression matrices, metadata, embeddings,
         marker tables, and `analysis_summary.txt`.
Example:
    python export_results.py --input adata.h5ad --output-dir results --cluster-key cell_type
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional, TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def _check_multimodal_data(adata: "AnnData") -> None:
    """Warn if non-RNA modalities are detected (e.g. CITE-seq ADT)."""
    obsm_keys = list(adata.obsm.keys())
    layer_keys = list(adata.layers.keys())
    if any("adt" in k.lower() or "protein" in k.lower() for k in obsm_keys + layer_keys):
        print("  multimodal data detected (ADT/protein) — RNA-only export; not analysing other modalities")


def _estimate_h5ad_size_mb(adata: "AnnData") -> float:
    """Rough size estimate in MB."""
    n_obs, n_vars = adata.n_obs, adata.n_vars
    bytes_per_value = 4
    layer_count = 1 + len(adata.layers)
    obsm_size = sum(getattr(v, "nbytes", 0) for v in adata.obsm.values())
    return (n_obs * n_vars * bytes_per_value * layer_count + obsm_size) / 1e6


def save_h5ad(
    adata: "AnnData",
    path: str,
    compression: Optional[str] = "lzf",
) -> None:
    """Write the AnnData with the chosen compression."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    size_mb = _estimate_h5ad_size_mb(adata)
    print(f"  estimated h5ad size: ~{size_mb:.1f} MB (compression={compression!r})")
    adata.write_h5ad(str(p), compression=compression)
    print(f"✓ Wrote {p}")


def export_expression_matrix(
    adata: "AnnData",
    output_dir: str,
    layer: Optional[str] = None,
    hvg_only: bool = False,
) -> None:
    """Write a CSV of the expression matrix."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    X = adata.layers[layer] if layer else adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    if hvg_only and "highly_variable" in adata.var.columns:
        idx = np.where(adata.var["highly_variable"].values)[0]
        X = X[:, idx]
        cols = adata.var_names.values[idx]
    else:
        cols = adata.var_names.values
    df = pd.DataFrame(X, index=adata.obs_names.values, columns=cols)
    name = f"expression_{layer or 'X'}{'_hvg' if hvg_only else ''}.csv"
    df.to_csv(out / name)
    print(f"✓ Saved {out / name} ({df.shape})")


def export_metadata(adata: "AnnData", output_dir: str) -> None:
    """Write `.obs` to CSV."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    adata.obs.to_csv(out / "cell_metadata.csv")
    print(f"✓ Saved {out / 'cell_metadata.csv'} ({adata.obs.shape})")


def export_embeddings(adata: "AnnData", output_dir: str) -> None:
    """Write PCA/UMAP/integrated coordinates to CSVs."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for key in ("X_pca", "X_umap", "X_scVI", "X_scANVI", "X_pca_harmony"):
        if key in adata.obsm:
            df = pd.DataFrame(adata.obsm[key], index=adata.obs_names.values,
                              columns=[f"{key[2:]}{i+1}" for i in range(adata.obsm[key].shape[1])])
            df.to_csv(out / f"{key.lower()}_coordinates.csv")
            print(f"✓ Saved {out / f'{key.lower()}_coordinates.csv'}")


def export_anndata_results(
    adata: "AnnData",
    output_dir: str = "results",
    cluster_key: str = "cell_type",
    export_raw: bool = True,
    export_metadata: bool = True,
    export_embeddings_flag: bool = True,
    export_markers: bool = True,
    compression: str = "lzf",
) -> None:
    """One-shot export of every artifact."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    _check_multimodal_data(adata)

    save_h5ad(adata, str(out / "adata_processed.h5ad"), compression=compression)

    if export_metadata:
        globals()["export_metadata"](adata, str(out))   # local function name collision

    if export_embeddings_flag:
        export_embeddings(adata, str(out))

    if export_raw and "counts" in adata.layers:
        try:
            export_expression_matrix(adata, str(out), layer="counts", hvg_only=True)
        except Exception as exc:
            print(f"  warning: couldn't export counts CSV: {exc}")

    try:
        export_expression_matrix(adata, str(out), layer=None, hvg_only=True)
    except Exception as exc:
        print(f"  warning: couldn't export normalized CSV: {exc}")

    if export_markers and "rank_genes_groups" in adata.uns:
        try:
            _export_markers_uns(adata, str(out))
        except Exception as exc:
            print(f"  warning: couldn't export markers: {exc}")

    create_summary_report(adata, output_dir=str(out), cluster_key=cluster_key)
    print("=== Export Complete ===")


def _export_markers_uns(adata: "AnnData", output_dir: str) -> None:
    """Pull markers out of `uns['rank_genes_groups']`."""
    rec = adata.uns["rank_genes_groups"]
    cats = rec["names"].dtype.names
    rows = []
    for cat in cats:
        for i in range(len(rec["names"][cat])):
            rows.append({
                "cluster": cat,
                "gene": rec["names"][cat][i],
                "score": float(rec["scores"][cat][i]),
                "logfoldchange": float(rec["logfoldchanges"][cat][i]),
                "pval": float(rec["pvals"][cat][i]),
                "pval_adj": float(rec["pvals_adj"][cat][i]),
            })
    df = pd.DataFrame(rows)
    df.to_csv(Path(output_dir) / "cluster_markers_all.csv", index=False)
    print(f"✓ Saved {Path(output_dir) / 'cluster_markers_all.csv'}")


def create_summary_report(
    adata: "AnnData",
    output_dir: str = "results",
    cluster_key: str = "cell_type",
) -> None:
    """Write `analysis_summary.txt`."""
    lines = [
        "Single-cell RNA-seq Analysis Summary",
        "=" * 40,
        "",
        f"Cells: {adata.n_obs}",
        f"Genes: {adata.n_vars}",
    ]
    if "highly_variable" in adata.var.columns:
        lines.append(f"HVGs: {int(adata.var['highly_variable'].sum())}")
    for m in ("n_genes_by_counts", "total_counts", "pct_counts_mt"):
        if m in adata.obs.columns:
            lines.append(f"  median {m}: {adata.obs[m].median():.2f}")
    if "predicted_doublet" in adata.obs.columns:
        lines.append(f"  predicted doublets: {int(adata.obs['predicted_doublet'].sum())}")
    if "integration" in adata.uns:
        lines.append("")
        lines.append("Integration:")
        for k, v in adata.uns["integration"].items():
            lines.append(f"  {k}: {v}")
    if cluster_key in adata.obs.columns:
        lines.append("")
        lines.append(f"Cluster counts ({cluster_key}):")
        counts = adata.obs[cluster_key].value_counts()
        for k, v in counts.items():
            lines.append(f"  {k}: {int(v)}")

    out = Path(output_dir) / "analysis_summary.txt"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"✓ Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad")
    parser.add_argument("--output-dir", default="results")
    parser.add_argument("--cluster-key", default="cell_type")
    parser.add_argument("--compression", default="lzf",
                        choices=["lzf", "gzip", "none"])
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    compression = None if args.compression == "none" else args.compression
    export_anndata_results(
        adata,
        output_dir=args.output_dir,
        cluster_key=args.cluster_key,
        compression=compression,
    )


if __name__ == "__main__":
    main()
