"""Batch integration via scVI / scANVI / Harmony.

Inputs:  AnnData with `.layers['counts']` (raw counts) and a batch column.
Outputs: Latent representation in `adata.obsm['X_scVI']`, `X_scANVI`, or
         `X_pca_harmony`; trained scVI/scANVI model is saved to disk if requested.
Example:
    python integrate_scvi.py --input adata.h5ad --output adata_int.h5ad \\
        --method scvi --batch-key batch --condition-key condition --n-latent 30
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def setup_for_integration(adata: "AnnData") -> "AnnData":
    """Ensure raw counts are in `layers['counts']`."""
    if "counts" not in adata.layers:
        print("  storing current .X into layers['counts'] — make sure it's raw counts")
        adata.layers["counts"] = adata.X.copy()
    return adata


def _check_batch_condition_confounding(
    adata: "AnnData",
    batch_key: str,
    condition_key: Optional[str],
) -> None:
    """Warn when any condition has only one batch (N=1) — pseudobulk DE invalid."""
    if condition_key is None or condition_key not in adata.obs.columns:
        return
    try:
        import pandas as pd
    except ImportError:
        return
    cross = pd.crosstab(adata.obs[condition_key], adata.obs[batch_key])
    n_batches_per_cond = (cross > 0).sum(axis=1)
    singletons = n_batches_per_cond[n_batches_per_cond == 1].index.tolist()
    if singletons:
        print(f"  WARNING: condition(s) {singletons} represented by 1 sample (N=1).")
        print("  Clustering and annotation remain valid; pseudobulk DE will refuse to run.")
        print("  Quantitative composition comparisons across conditions require caveats.")


def run_scvi_integration(
    adata: "AnnData",
    batch_key: str = "batch",
    condition_key: Optional[str] = None,
    n_latent: int = 30,
    n_layers: int = 2,
    n_hidden: int = 128,
    dropout_rate: float = 0.1,
    gene_likelihood: str = "nb",
    max_epochs: int = 400,
    early_stopping: bool = True,
    use_gpu: bool = True,
    save_path: Optional[str] = None,
) -> "AnnData":
    """Fit scVI on raw counts and store the latent representation in `X_scVI`."""
    try:
        import scvi
    except ImportError:
        sys.exit("Install with: pip install scvi-tools")

    adata = setup_for_integration(adata)
    _check_batch_condition_confounding(adata, batch_key, condition_key)

    print(f"scVI: n_latent={n_latent}, n_layers={n_layers}, max_epochs={max_epochs}")
    scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key, layer="counts")
    model = scvi.model.SCVI(
        adata,
        n_latent=n_latent,
        n_layers=n_layers,
        n_hidden=n_hidden,
        dropout_rate=dropout_rate,
        gene_likelihood=gene_likelihood,
    )

    accelerator = "gpu" if use_gpu else "cpu"
    model.train(max_epochs=max_epochs, early_stopping=early_stopping, accelerator=accelerator)

    adata.obsm["X_scVI"] = model.get_latent_representation()
    adata.uns["integration"] = {
        "method": "scVI",
        "n_latent": n_latent,
        "max_epochs": max_epochs,
        "batch_key": batch_key,
    }
    if save_path:
        model.save(save_path, overwrite=True)
        print(f"  saved scVI model to {save_path}")
    print(f"✓ scVI integration complete — X_scVI shape: {adata.obsm['X_scVI'].shape}")
    return adata


def run_scanvi_integration(
    adata: "AnnData",
    batch_key: str = "batch",
    labels_key: str = "cell_type_coarse",
    unlabeled_category: str = "Unknown",
    n_latent: int = 30,
    max_epochs: int = 400,
    use_gpu: bool = True,
) -> "AnnData":
    """Fit scANVI (semi-supervised) with partial cell-type labels."""
    try:
        import scvi
    except ImportError:
        sys.exit("Install with: pip install scvi-tools")

    adata = setup_for_integration(adata)
    if labels_key not in adata.obs.columns:
        raise KeyError(f"`{labels_key}` not in adata.obs")

    print(f"scANVI: labels_key={labels_key}, unlabeled={unlabeled_category!r}")
    scvi.model.SCANVI.setup_anndata(
        adata,
        batch_key=batch_key,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
        layer="counts",
    )
    model = scvi.model.SCANVI(adata, n_latent=n_latent, unlabeled_category=unlabeled_category)
    accelerator = "gpu" if use_gpu else "cpu"
    model.train(max_epochs=max_epochs, accelerator=accelerator)

    adata.obsm["X_scANVI"] = model.get_latent_representation()
    adata.obs["scanvi_pred"] = model.predict()
    adata.uns["integration"] = {"method": "scANVI", "n_latent": n_latent,
                                 "max_epochs": max_epochs, "batch_key": batch_key}
    print(f"✓ scANVI complete — X_scANVI shape: {adata.obsm['X_scANVI'].shape}")
    return adata


def run_harmony_integration(
    adata: "AnnData",
    batch_key: str = "batch",
    basis: str = "X_pca",
    n_pcs: int = 30,
    theta: float = 2.0,
    max_iter: int = 10,
) -> "AnnData":
    """Run Harmony on a PCA representation; store in `X_pca_harmony`."""
    try:
        import harmonypy as hm
    except ImportError:
        sys.exit("Install with: pip install harmonypy")

    if basis not in adata.obsm:
        raise KeyError(f"`{basis}` not in adata.obsm — run PCA first")
    pcs = adata.obsm[basis][:, :n_pcs]
    print(f"Harmony: theta={theta}, n_pcs={n_pcs}, max_iter={max_iter}")
    ho = hm.run_harmony(pcs, adata.obs, batch_key, theta=theta,
                        max_iter_harmony=max_iter)
    adata.obsm["X_pca_harmony"] = ho.Z_corr.T
    adata.uns["integration"] = {"method": "Harmony", "theta": theta,
                                 "n_pcs": n_pcs, "batch_key": batch_key}
    print(f"✓ Harmony complete — X_pca_harmony shape: {adata.obsm['X_pca_harmony'].shape}")
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--method", choices=["scvi", "scanvi", "harmony"], default="scvi")
    parser.add_argument("--batch-key", default="batch")
    parser.add_argument("--condition-key", default=None)
    parser.add_argument("--labels-key", default="cell_type_coarse")
    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--max-epochs", type=int, default=400)
    parser.add_argument("--theta", type=float, default=2.0)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--cpu", action="store_true", help="Force CPU for scVI/scANVI")
    parser.add_argument("--save-model", default=None)
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    if args.method == "scvi":
        adata = run_scvi_integration(adata, batch_key=args.batch_key,
                                      condition_key=args.condition_key,
                                      n_latent=args.n_latent, max_epochs=args.max_epochs,
                                      use_gpu=not args.cpu, save_path=args.save_model)
    elif args.method == "scanvi":
        adata = run_scanvi_integration(adata, batch_key=args.batch_key,
                                        labels_key=args.labels_key,
                                        n_latent=args.n_latent, max_epochs=args.max_epochs,
                                        use_gpu=not args.cpu)
    else:
        adata = run_harmony_integration(adata, batch_key=args.batch_key, n_pcs=args.n_pcs,
                                         theta=args.theta)
    adata.write_h5ad(args.output)
    print(f"✓ Wrote {args.output}")


if __name__ == "__main__":
    main()
