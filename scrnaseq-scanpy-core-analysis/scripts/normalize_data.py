"""Normalize raw counts. Defaults to target-sum + log1p; supports Pearson residuals.

Inputs:  AnnData with raw counts in `.X`.
Outputs: AnnData with `.layers["counts"]` = raw, `.X` = log-normalized (or
         `.layers["pearson_residuals"]` set when using Pearson mode).
Example:
    python normalize_data.py --input adata.h5ad --output adata_norm.h5ad \\
        --method standard --target-sum 1e4
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def run_standard_normalization(
    adata: "AnnData",
    target_sum: float = 1e4,
    exclude_highly_expressed: bool = False,
    log_transform: bool = True,
) -> "AnnData":
    """Library-size + log1p normalization. Stores raw in `layers['counts']`."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    print(f"Standard normalization (target_sum={target_sum}, log={log_transform})")

    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=target_sum,
                          exclude_highly_expressed=exclude_highly_expressed)
    if log_transform:
        sc.pp.log1p(adata)

    adata.layers["log1p_norm"] = adata.X.copy()
    print(f"✓ Normalization complete: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def run_pearson_residuals(
    adata: "AnnData",
    theta: float = 100,
    clip: Optional[float] = None,
) -> "AnnData":
    """Analytic Pearson residuals normalization (experimental in scanpy)."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    print(f"Pearson residuals (theta={theta}, clip={clip})")
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    sc.experimental.pp.normalize_pearson_residuals(adata, theta=theta, clip=clip)
    adata.layers["pearson_residuals"] = adata.X.copy()
    print(f"✓ Pearson residuals computed")
    return adata


def plot_normalization_comparison(adata: "AnnData", output_dir: str = "results/normalization") -> None:
    """Save a side-by-side QC plot before/after normalization."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    from pathlib import Path
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    if "counts" in adata.layers:
        try:
            import numpy as np
            raw_means = np.asarray(adata.layers["counts"].sum(axis=1)).ravel()
            axes[0].hist(raw_means, bins=80, color="steelblue", alpha=0.7)
            axes[0].set_title("Raw total counts (per cell)")
            axes[0].set_xlabel("UMIs")
        except Exception:
            pass
    try:
        import numpy as np
        norm_means = np.asarray(adata.X.sum(axis=1)).ravel()
        axes[1].hist(norm_means, bins=80, color="darkorange", alpha=0.7)
        axes[1].set_title("Normalized total (per cell)")
        axes[1].set_xlabel("normalized counts")
    except Exception:
        pass
    plt.tight_layout()
    out = Path(output_dir) / "normalization_comparison.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"✓ Saved {out}")


def plot_library_size_effect(adata: "AnnData", output_dir: str = "results/normalization") -> None:
    """Plot total_counts vs n_genes_by_counts colored by pct_counts_mt."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    from pathlib import Path
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if not all(c in adata.obs.columns for c in ("total_counts", "n_genes_by_counts", "pct_counts_mt")):
        print("  skipping library-size plot — missing QC metrics")
        return
    fig, ax = plt.subplots(figsize=(6, 5))
    sc = ax.scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                    c=adata.obs["pct_counts_mt"], s=3, alpha=0.7, cmap="viridis")
    plt.colorbar(sc, ax=ax, label="pct_counts_mt")
    ax.set_xlabel("total UMIs")
    ax.set_ylabel("genes detected")
    ax.set_xscale("log")
    ax.set_yscale("log")
    plt.tight_layout()
    out = Path(output_dir) / "library_size_effect.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"✓ Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--method", choices=["standard", "pearson"], default="standard")
    parser.add_argument("--target-sum", type=float, default=1e4)
    parser.add_argument("--theta", type=float, default=100)
    parser.add_argument("--clip", type=float, default=None)
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    if args.method == "standard":
        adata = run_standard_normalization(adata, target_sum=args.target_sum)
    else:
        adata = run_pearson_residuals(adata, theta=args.theta, clip=args.clip)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
