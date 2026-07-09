"""Diagnostic metrics for batch integration quality (LISI, ASW) and plots.

Inputs:  AnnData with one or more low-dim representations in `.obsm`.
Outputs: dictionaries of LISI / ASW scores; comparison plots.
Example:
    python integration_diagnostics.py --input adata.h5ad --batch-key batch \\
        --label-key cell_type --use-rep X_scVI --output-dir results/integration
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def _lisi_one(X: np.ndarray, labels: np.ndarray, perplexity: float = 30) -> np.ndarray:
    """Compute per-cell LISI given a representation and categorical labels.

    Approximation via Gaussian-kernel weights on KNN distances. Mirrors the
    Harmony LISI formulation closely enough for relative comparisons.
    """
    try:
        from sklearn.neighbors import NearestNeighbors
    except ImportError:
        sys.exit("Install with: pip install scikit-learn")
    k = int(min(len(labels) - 1, max(int(perplexity * 3), 30)))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    dists, idx = nn.kneighbors(X)
    dists = dists[:, 1:]   # drop self
    idx = idx[:, 1:]
    # Gaussian kernel with bandwidth ~ perplexity
    sigma = np.maximum(dists.mean(axis=1, keepdims=True), 1e-8)
    w = np.exp(-(dists ** 2) / (2 * sigma ** 2))
    w /= np.maximum(w.sum(axis=1, keepdims=True), 1e-8)

    unique = np.unique(labels)
    lisi = np.zeros(X.shape[0])
    for i in range(X.shape[0]):
        neigh_lbl = labels[idx[i]]
        p_squared = 0.0
        for lbl in unique:
            p = w[i][neigh_lbl == lbl].sum()
            p_squared += p ** 2
        lisi[i] = 1.0 / max(p_squared, 1e-8)
    return lisi


def compute_lisi_scores(
    adata: "AnnData",
    batch_key: str = "batch",
    label_key: Optional[str] = None,
    use_rep: str = "X_scVI",
    perplexity: float = 30,
) -> Dict[str, float]:
    """Mean batch (iLISI) and cell type (cLISI) LISI scores."""
    if use_rep not in adata.obsm:
        raise KeyError(f"`{use_rep}` not in adata.obsm")
    X = np.asarray(adata.obsm[use_rep])
    out: Dict[str, float] = {}

    if batch_key in adata.obs.columns:
        lab = adata.obs[batch_key].astype(str).values
        out["batch_LISI"] = float(_lisi_one(X, lab, perplexity=perplexity).mean())
        print(f"  iLISI ({batch_key}, {use_rep}): {out['batch_LISI']:.2f}")

    if label_key and label_key in adata.obs.columns:
        lab = adata.obs[label_key].astype(str).values
        out["celltype_LISI"] = float(_lisi_one(X, lab, perplexity=perplexity).mean())
        print(f"  cLISI ({label_key}, {use_rep}): {out['celltype_LISI']:.2f}")

    return out


def compute_asw_scores(
    adata: "AnnData",
    batch_key: str = "batch",
    label_key: Optional[str] = None,
    use_rep: str = "X_scVI",
) -> Dict[str, float]:
    """Mean batch ASW and cell-type ASW (silhouette)."""
    try:
        from sklearn.metrics import silhouette_score
    except ImportError:
        sys.exit("Install with: pip install scikit-learn")

    if use_rep not in adata.obsm:
        raise KeyError(f"`{use_rep}` not in adata.obsm")
    X = np.asarray(adata.obsm[use_rep])
    out: Dict[str, float] = {}

    if batch_key in adata.obs.columns and adata.obs[batch_key].nunique() > 1:
        lab = adata.obs[batch_key].astype(str).values
        try:
            out["batch_ASW"] = float(silhouette_score(X, lab, sample_size=min(10000, X.shape[0])))
            print(f"  batch ASW: {out['batch_ASW']:.3f} (target ≈ 0)")
        except Exception as exc:
            print(f"  batch ASW failed: {exc}")

    if label_key and label_key in adata.obs.columns and adata.obs[label_key].nunique() > 1:
        lab = adata.obs[label_key].astype(str).values
        try:
            out["celltype_ASW"] = float(silhouette_score(X, lab, sample_size=min(10000, X.shape[0])))
            print(f"  cell type ASW: {out['celltype_ASW']:.3f} (target > 0.5)")
        except Exception as exc:
            print(f"  cell type ASW failed: {exc}")

    return out


def plot_integration_metrics(
    adata: "AnnData",
    lisi: Dict[str, float],
    asw: Optional[Dict[str, float]] = None,
    output_dir: str = "results/integration",
) -> None:
    """Bar plot of LISI and ASW summary scores."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    df_rows = []
    for k, v in lisi.items():
        df_rows.append({"metric": k, "score": v})
    if asw:
        for k, v in asw.items():
            df_rows.append({"metric": k, "score": v})
    if not df_rows:
        return
    df = pd.DataFrame(df_rows)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(df["metric"], df["score"], color="steelblue")
    ax.set_ylabel("score")
    ax.set_title("Integration metrics")
    plt.xticks(rotation=20, ha="right")
    plt.tight_layout()
    out = Path(output_dir) / "integration_metrics.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved {out}")


def compare_integration_methods(
    adata: "AnnData",
    batch_key: str = "batch",
    label_key: Optional[str] = None,
    representations: Optional[List[str]] = None,
    output_dir: str = "results/integration",
) -> pd.DataFrame:
    """Compute LISI and ASW for each representation and write a comparison CSV."""
    representations = representations or [r for r in ("X_pca", "X_scVI", "X_scANVI", "X_pca_harmony") if r in adata.obsm]
    rows = []
    for rep in representations:
        print(f"-- {rep} --")
        lisi = compute_lisi_scores(adata, batch_key=batch_key, label_key=label_key, use_rep=rep)
        asw = compute_asw_scores(adata, batch_key=batch_key, label_key=label_key, use_rep=rep)
        rows.append({"representation": rep, **lisi, **asw})
    df = pd.DataFrame(rows)
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    out = Path(output_dir) / "integration_comparison.csv"
    df.to_csv(out, index=False)
    print(f"✓ Saved {out}")
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--label-key", default=None)
    parser.add_argument("--use-rep", default="X_scVI")
    parser.add_argument("--compare", action="store_true", help="Compare all available representations")
    parser.add_argument("--output-dir", default="results/integration")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    if args.compare:
        compare_integration_methods(adata, batch_key=args.batch_key,
                                     label_key=args.label_key,
                                     output_dir=args.output_dir)
    else:
        lisi = compute_lisi_scores(adata, batch_key=args.batch_key,
                                    label_key=args.label_key,
                                    use_rep=args.use_rep)
        asw = compute_asw_scores(adata, batch_key=args.batch_key,
                                  label_key=args.label_key,
                                  use_rep=args.use_rep)
        plot_integration_metrics(adata, lisi, asw, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
