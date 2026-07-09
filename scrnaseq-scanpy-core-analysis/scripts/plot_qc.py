"""QC plots: violins, scatters, histograms, top-expressed genes.

Inputs:  AnnData with QC metrics in `.obs`.
Outputs: PNG figures under the chosen output dir.
Example:
    python plot_qc.py --input adata_qc.h5ad --output-dir results/qc --batch-key batch
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def _save_fig(fig, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    try:
        out_path.with_suffix(".svg")
        fig.savefig(out_path.with_suffix(".svg"), bbox_inches="tight")
    except Exception:
        pass
    print(f"✓ Saved {out_path}")


def plot_qc_violin(
    adata: "AnnData",
    metrics: Optional[List[str]] = None,
    batch_key: Optional[str] = None,
    output_dir: str = "results/qc",
) -> None:
    """Violins of `n_genes_by_counts`, `total_counts`, `pct_counts_mt`."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")

    metrics = metrics or ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    metrics = [m for m in metrics if m in adata.obs.columns]
    if not metrics:
        print("  no QC metrics found")
        return
    fig = sc.pl.violin(
        adata, metrics,
        jitter=0.4, multi_panel=True,
        groupby=batch_key if batch_key in adata.obs.columns else None,
        show=False, return_fig=True,
    )
    _save_fig(fig, Path(output_dir) / "qc_violin.png")
    plt.close(fig)


def plot_qc_scatter(
    adata: "AnnData",
    output_dir: str = "results/qc",
) -> None:
    """Scatter of total_counts vs n_genes_by_counts, colored by pct_counts_mt."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    if not all(c in adata.obs.columns for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt")):
        print("  missing QC metrics — skip")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                    c=adata.obs["pct_counts_mt"], s=3, alpha=0.7, cmap="viridis")
    plt.colorbar(sc, ax=ax, label="pct_counts_mt")
    ax.set_xlabel("total counts")
    ax.set_ylabel("n_genes_by_counts")
    ax.set_xscale("log")
    ax.set_yscale("log")
    _save_fig(fig, Path(output_dir) / "qc_scatter.png")
    plt.close(fig)


def plot_qc_histograms(
    adata: "AnnData",
    output_dir: str = "results/qc",
) -> None:
    """Histograms of each QC metric."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    metrics = [m for m in ("n_genes_by_counts", "total_counts", "pct_counts_mt") if m in adata.obs.columns]
    if not metrics:
        return
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 3))
    if len(metrics) == 1:
        axes = [axes]
    for ax, m in zip(axes, metrics):
        ax.hist(adata.obs[m], bins=80, color="steelblue", alpha=0.7)
        ax.set_xlabel(m)
        ax.set_ylabel("cells")
    plt.tight_layout()
    _save_fig(fig, Path(output_dir) / "qc_histograms.png")
    plt.close(fig)


def plot_highest_expr_genes(
    adata: "AnnData",
    n_top: int = 20,
    output_dir: str = "results/qc",
) -> None:
    """Top-`n_top` genes by total expression (good for spotting ambient soup)."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")
    fig = sc.pl.highest_expr_genes(adata, n_top=n_top, show=False, return_fig=True)
    _save_fig(fig, Path(output_dir) / "highest_expr_genes.png")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", default="results/qc")
    parser.add_argument("--batch-key", default=None)
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    plot_qc_violin(adata, batch_key=args.batch_key, output_dir=args.output_dir)
    plot_qc_scatter(adata, output_dir=args.output_dir)
    plot_qc_histograms(adata, output_dir=args.output_dir)
    plot_highest_expr_genes(adata, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
