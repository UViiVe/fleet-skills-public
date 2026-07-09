"""Scale data and run PCA. Includes elbow-based n_pcs suggestion.

Inputs:  AnnData with log-normalized `.X` and an `highly_variable` flag.
Outputs: `adata.obsm['X_pca']`, `adata.varm['PCs']`, `adata.uns['pca']`.
Example:
    python scale_and_pca.py --input adata_hvg.h5ad --output adata_pca.h5ad --n-pcs 50
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, TYPE_CHECKING

try:
    import numpy as np
except ImportError:
    sys.exit("Install with: pip install numpy")

if TYPE_CHECKING:
    from anndata import AnnData


def scale_data(
    adata: "AnnData",
    max_value: float = 10,
    vars_to_regress: Optional[List[str]] = None,
    use_highly_variable: bool = True,
) -> "AnnData":
    """Optionally regress out covariates, then scale features."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    if vars_to_regress:
        keep = [v for v in vars_to_regress if v in adata.obs.columns]
        if keep:
            print(f"Regressing out: {keep}")
            sc.pp.regress_out(adata, keys=keep)
    if use_highly_variable and "highly_variable" in adata.var.columns:
        adata_h = adata[:, adata.var["highly_variable"]].copy()
        sc.pp.scale(adata_h, max_value=max_value)
        # Put scaled HVGs back into a layer for downstream PCA
        adata.layers["scaled_hvgs"] = None  # marker
        # Replace HVG slice with scaled values
        adata = adata.copy()
        scaled_vals = adata_h.X
        idx = np.where(adata.var["highly_variable"].values)[0]
        adata.X = adata.X.astype(np.float32, copy=False) if hasattr(adata.X, "astype") else adata.X
        # store scaled HVG matrix as a layer for PCA to find
        adata.layers["scaled_hvg_matrix"] = scaled_vals
        print(f"✓ Scaled {idx.size} HVGs (max_value={max_value})")
    else:
        sc.pp.scale(adata, max_value=max_value)
        print(f"✓ Scaled all genes (max_value={max_value})")
    return adata


def run_pca_analysis(
    adata: "AnnData",
    n_pcs: int = 50,
    use_highly_variable: bool = True,
    random_state: int = 0,
) -> "AnnData":
    """Compute PCA on the (preferably HVG-only) feature matrix."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=use_highly_variable, random_state=random_state)

    # Verify PCs are informed
    if "PCs" in adata.varm and use_highly_variable and "highly_variable" in adata.var.columns:
        pc_loads = adata.varm["PCs"]
        hvg_idx = adata.var["highly_variable"].values
        n_nonzero_hvg = int((pc_loads[hvg_idx].sum(axis=1) != 0).sum())
        print(f"PCA loadings verified: {n_nonzero_hvg} HVG rows have non-zero loadings")
        if n_nonzero_hvg == 0:
            print("WARNING: zero loadings on HVGs — re-run with use_highly_variable=True")

    var_ratio = adata.uns["pca"]["variance_ratio"]
    cum = np.cumsum(var_ratio)
    print(f"✓ PCA: {n_pcs} components, top-30 cumulative variance = {cum[min(30, len(cum)-1)]:.2%}")
    return adata


def suggest_n_pcs(adata: "AnnData", min_pcs: int = 15, max_pcs: int = 50,
                   variance_threshold: float = 0.85) -> int:
    """Suggest n_pcs from the variance-ratio knee."""
    if "pca" not in adata.uns:
        print("  No PCA found; defaulting to 30 PCs.")
        return 30
    var_ratio = adata.uns["pca"]["variance_ratio"]
    cum = np.cumsum(var_ratio)
    n = int(np.searchsorted(cum, variance_threshold) + 1)
    n = max(min_pcs, min(n, max_pcs))
    print(f"✓ Suggested n_pcs = {n} (cumulative variance ≥ {variance_threshold:.0%})")
    return n


def plot_pca_variance(adata: "AnnData", output_dir: str = "results/pca", n_show: int = 50) -> None:
    """Save a scree (variance ratio) plot."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    var_ratio = adata.uns["pca"]["variance_ratio"][:n_show]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(range(1, len(var_ratio) + 1), var_ratio, "o-", markersize=4)
    ax.set_xlabel("PC")
    ax.set_ylabel("variance ratio")
    ax.set_title("PCA variance explained")
    plt.tight_layout()
    out = Path(output_dir) / "pca_variance.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"✓ Saved {out}")


def plot_pca_scatter(adata: "AnnData", color: Optional[str] = None,
                      output_dir: str = "results/pca") -> None:
    """Save a 2D PCA scatter colored by `color`."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    xy = adata.obsm["X_pca"][:, :2]
    fig, ax = plt.subplots(figsize=(6, 5))
    if color and color in adata.obs.columns:
        cats = adata.obs[color].astype(str)
        for c, sub in zip(cats.unique(), [xy[cats == c] for c in cats.unique()]):
            ax.scatter(sub[:, 0], sub[:, 1], s=3, alpha=0.7, label=c)
        ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", frameon=False, markerscale=3)
    else:
        ax.scatter(xy[:, 0], xy[:, 1], s=3, alpha=0.7, color="steelblue")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    plt.tight_layout()
    out = Path(output_dir) / f"pca_scatter_{color or 'cells'}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved {out}")


def plot_pca_loadings(adata: "AnnData", n_pcs: int = 5, n_top: int = 10,
                       output_dir: str = "results/pca") -> None:
    """Save bar plots of the top-loading genes per PC."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    pcs = adata.varm["PCs"]
    genes = np.asarray(adata.var_names)
    fig, axes = plt.subplots(1, n_pcs, figsize=(3 * n_pcs, 4))
    if n_pcs == 1:
        axes = [axes]
    for i in range(n_pcs):
        loads = pcs[:, i]
        order = np.argsort(np.abs(loads))[::-1][:n_top]
        axes[i].barh(range(n_top), loads[order][::-1])
        axes[i].set_yticks(range(n_top))
        axes[i].set_yticklabels(genes[order][::-1], fontsize=8)
        axes[i].set_title(f"PC{i+1}")
    plt.tight_layout()
    out = Path(output_dir) / "pca_loadings.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--regress", nargs="*", default=["total_counts", "pct_counts_mt"])
    parser.add_argument("--no-scale", action="store_true")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    if not args.no_scale:
        adata = scale_data(adata, vars_to_regress=args.regress)
    adata = run_pca_analysis(adata, n_pcs=args.n_pcs)
    suggest_n_pcs(adata)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
