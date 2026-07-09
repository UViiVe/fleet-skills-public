"""Per-perturbation outlier detection in PCA-on-DE-genes space.

For every perturbation, run a fast t-test against control to find its DE genes,
PCA those genes, fit an outlier detector (LocalOutlierFactor by default) on
control cells, and classify the perturbed cells. The number of outliers per
perturbation becomes the hit-call.

Example:
    python detect_perturbed_cells.py --in normalized.h5ad --out-dir results/
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Dict, Optional

try:
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata as ad
    import diffxpy.api as de
    from sklearn.neighbors import LocalOutlierFactor
    from sklearn.ensemble import IsolationForest
    from sklearn.svm import OneClassSVM
except ImportError:
    sys.exit(
        "Install with: pip install scanpy anndata diffxpy scikit-learn"
    )


def _make_clf(method: str):
    if method == "LocalOutlierFactor":
        return LocalOutlierFactor(novelty=True, contamination="auto")
    if method == "IsolationForest":
        return IsolationForest(contamination="auto", random_state=42)
    if method == "OneClassSVM":
        return OneClassSVM(nu=0.1)
    raise ValueError(f"Unknown method {method!r}")


def detect_perturbed_cells(
    adata: "ad.AnnData",
    control_group: str = "non-targeting",
    gene_col: str = "gene",
    sgrna_col: str = "sgRNA",
    n_pcs: int = 4,
    method: str = "LocalOutlierFactor",
    min_de_genes: int = 10,
    min_outliers: int = 10,
    pval_threshold: float = 0.05,
    save_plots: bool = True,
    plot_dir: str = "figures/",
) -> Dict:
    """Run the outlier-detection pipeline across every perturbation."""
    os.makedirs(plot_dir, exist_ok=True)

    outliers_global: list[str] = []
    summary_rows: list[dict] = []
    de_results: Dict[str, "pd.DataFrame"] = {}
    adata_per_pert: Dict[str, "ad.AnnData"] = {}

    perturbations = [g for g in adata.obs[gene_col].unique() if g != control_group]
    print(f"Detecting perturbed cells for {len(perturbations)} perturbations")
    print(f"  method={method} n_pcs={n_pcs} min_de_genes={min_de_genes} min_outliers={min_outliers}")

    for i, gene in enumerate(perturbations, 1):
        if i % 10 == 0:
            print(f"  [{i}/{len(perturbations)}] {gene}")

        adata_pair = adata[adata.obs[gene_col].isin([gene, control_group])].copy()
        test = de.test.t_test(data=adata_pair, grouping=gene_col)
        de_df = test.summary()
        de_results[gene] = de_df

        diff_genes = de_df[de_df["pval"] < pval_threshold]["gene"].tolist()
        n_de = len(diff_genes)
        n_pert = int((adata_pair.obs[gene_col] == gene).sum())

        if n_de < min_de_genes:
            summary_rows.append({
                "gene": gene, "n_cells": n_pert, "n_de_genes": n_de,
                "n_outliers": 0, "outlier_fraction": 0.0, "is_hit": False,
                "reason": f"insufficient_DE_genes (<{min_de_genes})",
            })
            continue

        adata_pair = adata_pair[:, diff_genes].copy()
        sc.tl.pca(adata_pair, svd_solver="arpack", n_comps=n_pcs)

        ctrl_mask = (adata_pair.obs[gene_col] == control_group).values
        pert_mask = (adata_pair.obs[gene_col] == gene).values
        control_pca = adata_pair.obsm["X_pca"][ctrl_mask]
        perturbed_pca = adata_pair.obsm["X_pca"][pert_mask]

        clf = _make_clf(method).fit(control_pca)
        preds = clf.predict(perturbed_pca)
        n_out = int((preds == -1).sum())
        frac = n_out / len(perturbed_pca) if len(perturbed_pca) else 0.0

        adata_view = adata[adata_pair.obs_names].copy()
        adata_view.obsm = adata_pair.obsm.copy()
        adata_view.obs["classification"] = 2
        adata_view.obs.loc[pert_mask, "classification"] = preds
        adata_view.obs["classification"] = adata_view.obs["classification"].astype("category")
        sc.pp.neighbors(adata_view, n_pcs=n_pcs)
        sc.tl.umap(adata_view)

        if save_plots:
            colors = ["classification", gene_col, sgrna_col]
            if gene in adata_view.var_names:
                colors.insert(2, gene)
            try:
                sc.pl.umap(adata_view, color=colors, save=f"_{gene}.pdf", show=False)
                sc.pl.pca(adata_view, color=colors, save=f"_{gene}.pdf", show=False)
            except Exception as exc:
                print(f"  plot failure for {gene}: {exc}")

        adata_per_pert[gene] = adata_view
        is_hit = n_out >= min_outliers
        if is_hit:
            outliers_global.extend(adata_view[adata_view.obs["classification"] == -1].obs_names.tolist())

        summary_rows.append({
            "gene": gene, "n_cells": n_pert, "n_de_genes": n_de,
            "n_outliers": n_out, "outlier_fraction": frac, "is_hit": is_hit,
            "reason": "hit" if is_hit else f"insufficient_outliers (<{min_outliers})",
        })

    summary_df = pd.DataFrame(summary_rows)
    control_cells = adata[adata.obs[gene_col] == control_group].obs_names.tolist()
    kept = outliers_global + control_cells

    n_hits = int(summary_df["is_hit"].sum())
    print(
        f"✓ Outlier detection complete: {n_hits}/{len(perturbations)} hits, "
        f"{len(outliers_global)} outlier cells"
    )
    return {
        "outlier_cells": outliers_global,
        "kept_cells": kept,
        "perturbation_summary": summary_df,
        "de_results": de_results,
        "adata_dict": adata_per_pert,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out-dir", default="results/")
    parser.add_argument("--control", default="non-targeting")
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--method", default="LocalOutlierFactor",
                        choices=["LocalOutlierFactor", "IsolationForest", "OneClassSVM"])
    parser.add_argument("--n-pcs", type=int, default=4)
    parser.add_argument("--min-de-genes", type=int, default=10)
    parser.add_argument("--min-outliers", type=int, default=10)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.in_path)
    os.makedirs(args.out_dir, exist_ok=True)
    results = detect_perturbed_cells(
        adata, control_group=args.control, gene_col=args.gene_col,
        method=args.method, n_pcs=args.n_pcs,
        min_de_genes=args.min_de_genes, min_outliers=args.min_outliers,
        plot_dir=os.path.join(args.out_dir, "figures"),
    )
    results["perturbation_summary"].to_csv(
        os.path.join(args.out_dir, "perturbation_summary.csv"), index=False
    )
    pd.DataFrame({"cell_barcode": results["outlier_cells"]}).to_csv(
        os.path.join(args.out_dir, "outlier_cells.csv"), index=False
    )
    print(f"✓ Wrote summary + outlier list to {args.out_dir}")


if __name__ == "__main__":
    main()
