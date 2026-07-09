"""Assign cell type labels (manual mapping or CellTypist).

Inputs:  AnnData with a cluster column; optionally a CellTypist model name.
Outputs: `adata.obs['cell_type']` (manual) or `predicted_labels` (CellTypist).
Example:
    python annotate_celltypes.py --input adata.h5ad --output annotated.h5ad \\
        --mode celltypist --model Immune_All_Low.pkl
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData

# Suspect labels that CellTypist over-calls (used by `_validate_celltypist_annotations`).
_SUSPECT_LABELS = {"ILC", "HSC", "stem-like", "low-quality"}


def annotate_clusters_manual(
    adata: "AnnData",
    annotations: Dict[str, str],
    cluster_key: str = "leiden_0.8",
    celltype_col: str = "cell_type",
) -> "AnnData":
    """Apply a `{cluster_id: cell_type}` dictionary to `adata.obs[celltype_col]`."""
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"`{cluster_key}` not in adata.obs")
    clusters = adata.obs[cluster_key].astype(str)
    mapped = clusters.map(annotations)
    n_missing = int(mapped.isna().sum())
    if n_missing:
        print(f"  warning: {n_missing} cells not covered by mapping — labelled 'Unassigned'")
        mapped = mapped.fillna("Unassigned")
    adata.obs[celltype_col] = mapped.astype("category")
    print(f"✓ Annotated {len(annotations)} clusters → {adata.obs[celltype_col].nunique()} unique labels")
    return adata


def annotate_with_celltypist(
    adata: "AnnData",
    model: str = "Immune_All_Low.pkl",
    majority_voting: bool = True,
    over_clustering: Optional[str] = None,
) -> "AnnData":
    """Run CellTypist with `model` and store predicted labels."""
    try:
        import celltypist
    except ImportError:
        sys.exit("Install with: pip install celltypist")

    print(f"CellTypist: model={model}, majority_voting={majority_voting}")
    predictions = celltypist.annotate(
        adata,
        model=model,
        majority_voting=majority_voting,
        over_clustering=over_clustering,
    )
    adata = predictions.to_adata()
    print(f"✓ CellTypist labels: {adata.obs['predicted_labels'].nunique()} unique")
    _validate_celltypist_annotations(adata)
    return adata


def _validate_celltypist_annotations(adata: "AnnData") -> None:
    """Warn about known CellTypist failure modes."""
    if "predicted_labels" not in adata.obs.columns:
        return
    labels = adata.obs["predicted_labels"].astype(str)
    flagged = labels.apply(lambda x: any(s.lower() in x.lower() for s in _SUSPECT_LABELS))
    n = int(flagged.sum())
    if n:
        print(f"  CellTypist labels suspect for {n} cells (ILC/HSC/stem-like). "
              f"Cross-check against canonical markers.")


def plot_annotated_umap(
    adata: "AnnData",
    celltype_col: str = "cell_type",
    output_dir: str = "results/annotation",
) -> None:
    """UMAP coloured by the annotation column."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fig = sc.pl.umap(adata, color=celltype_col, legend_loc="on data",
                     legend_fontsize="x-small", legend_fontoutline=2,
                     show=False, return_fig=True)
    out = Path(output_dir) / "umap_annotated.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved {out}")


def create_annotation_summary(
    adata: "AnnData",
    celltype_col: str = "cell_type",
    output_dir: str = "results/annotation",
) -> None:
    """Save a JSON of cell counts per cell type."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    counts = adata.obs[celltype_col].value_counts().to_dict()
    out = Path(output_dir) / "cell_type_counts.json"
    with open(out, "w") as f:
        json.dump({k: int(v) for k, v in counts.items()}, f, indent=2, sort_keys=True)
    print(f"✓ Saved {out}")


def plot_annotation_sankey(adata: "AnnData",
                            from_col: str = "leiden_0.8",
                            to_col: str = "cell_type",
                            output_dir: str = "results/annotation") -> None:
    """Save a basic flow table of cluster → cell type counts (no plotly dependency)."""
    try:
        import pandas as pd
    except ImportError:
        sys.exit("Install with: pip install pandas")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    cross = pd.crosstab(adata.obs[from_col], adata.obs[to_col])
    out = Path(output_dir) / "cluster_to_celltype.csv"
    cross.to_csv(out)
    print(f"✓ Saved {out}")


def compare_annotations(
    adata: "AnnData",
    col_a: str,
    col_b: str,
    output_dir: str = "results/annotation",
) -> None:
    """Save a crosstab comparing two annotation columns."""
    try:
        import pandas as pd
    except ImportError:
        sys.exit("Install with: pip install pandas")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    cross = pd.crosstab(adata.obs[col_a], adata.obs[col_b])
    out = Path(output_dir) / f"compare_{col_a}_vs_{col_b}.csv"
    cross.to_csv(out)
    print(f"✓ Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--mode", choices=["manual", "celltypist"], default="manual")
    parser.add_argument("--mapping-json", help="JSON dict mapping cluster_id → cell_type")
    parser.add_argument("--cluster-key", default="leiden_0.8")
    parser.add_argument("--model", default="Immune_All_Low.pkl")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)

    if args.mode == "manual":
        if not args.mapping_json:
            sys.exit("--mapping-json required for manual mode")
        with open(args.mapping_json) as f:
            mapping = json.load(f)
        adata = annotate_clusters_manual(adata, mapping, cluster_key=args.cluster_key)
    else:
        adata = annotate_with_celltypist(adata, model=args.model)

    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
