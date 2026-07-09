"""Run a three-level trajectory inference pipeline (PAGA, scVelo, CellRank).

Level 1 (always runs):  PAGA + diffusion pseudotime on scanpy.
Level 2 (optional):     scVelo dynamical model (or stochastic fallback).
Level 3 (optional):     CellRank fate probabilities and driver genes.

Inputs:  preprocessed AnnData with PCA/UMAP and a cluster column.
Output:  a dict with pseudotime, PAGA, gene tables, and optional velocity/fate results.
Example:
    python run_trajectory_analysis.py --input adata.h5ad --output results.pkl \\
        --root-cell-type Ductal --cluster-key clusters
"""

from __future__ import annotations

import argparse
import pickle
import sys
import warnings
from typing import Any, Dict, Optional, TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def run_trajectory(
    adata: "AnnData",
    root_cell_type: str = "Ductal",
    cluster_key: str = "clusters",
    n_neighbors: int = 30,
    n_pcs: int = 30,
    n_dcs: int = 10,
    run_velocity: bool = True,
    run_cellrank: bool = True,
) -> Dict[str, Any]:
    """Run all three levels of trajectory inference."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    results: Dict[str, Any] = {
        "velocity_results": None,
        "cellrank_results": None,
        "parameters": {
            "root_cell_type": root_cell_type,
            "cluster_key": cluster_key,
            "n_neighbors": n_neighbors,
            "n_pcs": n_pcs,
            "n_dcs": n_dcs,
        },
    }

    print("=" * 60)
    print("LEVEL 1: Core trajectory (scanpy)")
    print("=" * 60)

    print("  Computing neighbor graph...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)

    print("  Running PAGA...")
    sc.tl.paga(adata, groups=cluster_key)
    results["paga_connectivities"] = adata.uns["paga"]["connectivities"].toarray()

    print("  Re-initialising UMAP with PAGA layout...")
    try:
        sc.pl.paga(adata, plot=False)
        sc.tl.umap(adata, init_pos="paga")
    except Exception as exc:
        print(f"  warning: PAGA-initialised UMAP failed: {exc}; keeping existing UMAP")

    print("  Computing diffusion map...")
    sc.tl.diffmap(adata, n_comps=n_dcs)
    results["diffmap"] = adata.obsm["X_diffmap"]

    root_cell, root_idx = _select_root_cell(adata, root_cell_type, cluster_key)
    adata.uns["iroot"] = root_idx
    results["root_cell"] = root_cell
    print(f"  Root cell: {root_cell} (type={root_cell_type})")

    print("  Computing diffusion pseudotime...")
    sc.tl.dpt(adata)
    results["pseudotime"] = adata.obs["dpt_pseudotime"].copy()
    finite = results["pseudotime"][np.isfinite(results["pseudotime"])]
    print(f"  Pseudotime range: [{finite.min():.3f}, {finite.max():.3f}]")
    if not finite.shape[0] == adata.n_obs:
        print(f"  warning: {adata.n_obs - finite.shape[0]} cells have non-finite pseudotime")

    print("  Identifying trajectory-associated genes...")
    results["trajectory_genes"] = _find_trajectory_genes(adata)
    print(f"  Found {len(results['trajectory_genes'])} trajectory genes")
    print("\n✓ Core trajectory analysis complete")

    has_velocity_layers = ("spliced" in adata.layers and "unspliced" in adata.layers)
    if run_velocity and has_velocity_layers:
        print("\n" + "=" * 60)
        print("LEVEL 2: RNA velocity (scVelo)")
        print("=" * 60)
        try:
            results["velocity_results"] = _run_rna_velocity(adata)
            print("\n✓ RNA velocity analysis complete")
        except ImportError:
            print("  scVelo not installed — skipping (`pip install scvelo`)")
        except Exception as exc:
            print(f"  RNA velocity failed: {exc}")
    elif run_velocity:
        print("\n  No spliced/unspliced layers — skipping RNA velocity")
    else:
        print("\n  RNA velocity skipped (run_velocity=False)")

    if run_cellrank:
        print("\n" + "=" * 60)
        print("LEVEL 3: CellRank fate mapping")
        print("=" * 60)
        try:
            results["cellrank_results"] = _run_cellrank(
                adata,
                cluster_key=cluster_key,
                has_velocity=results["velocity_results"] is not None,
            )
            print("\n✓ CellRank fate mapping complete")
        except ImportError:
            print("  CellRank not installed — skipping (`pip install cellrank`)")
        except Exception as exc:
            print(f"  CellRank failed: {exc}")
    else:
        print("\n  CellRank skipped (run_cellrank=False)")

    print("\n" + "=" * 60)
    print("Trajectory analysis summary")
    print("=" * 60)
    print(f"  Cells:            {adata.n_obs:,}")
    print(f"  Root cell type:   {root_cell_type}")
    print(f"  Pseudotime range: [{finite.min():.3f}, {finite.max():.3f}]")
    print(f"  PAGA clusters:    {adata.obs[cluster_key].nunique()}")
    print(f"  Trajectory genes: {len(results['trajectory_genes'])}")
    print(f"  RNA velocity:     {'yes' if results['velocity_results'] else 'no'}")
    print(f"  CellRank fates:   {'yes' if results['cellrank_results'] else 'no'}")
    print("\n✓ Trajectory analysis completed successfully!")
    return results


def _select_root_cell(adata, root_cell_type, cluster_key):
    if cluster_key not in adata.obs.columns:
        raise KeyError(f"`{cluster_key}` not in adata.obs")
    mask = adata.obs[cluster_key].astype(str) == root_cell_type
    if not mask.any():
        raise KeyError(
            f"Root cell type {root_cell_type!r} not found in {cluster_key}. "
            f"Available: {sorted(adata.obs[cluster_key].unique().tolist())}"
        )
    idx = int(np.where(mask)[0][0])
    return adata.obs_names[idx], idx


def _find_trajectory_genes(adata,
                            min_pct_cells: float = 0.10,
                            fdr_threshold: float = 0.05,
                            min_abs_corr: float = 0.30) -> pd.DataFrame:
    try:
        from scipy.stats import spearmanr
        from statsmodels.stats.multitest import multipletests
    except ImportError:
        sys.exit("Install with: pip install scipy statsmodels")

    pt = adata.obs["dpt_pseudotime"].values
    finite_mask = np.isfinite(pt)
    pt = pt[finite_mask]
    X = adata.X[finite_mask, :]
    if hasattr(X, "toarray"):
        X = X.toarray()

    min_cells = int(min_pct_cells * X.shape[0])
    expressed = (X > 0).sum(axis=0) >= min_cells
    genes = adata.var_names[expressed]
    X = X[:, expressed]

    corrs = np.zeros(X.shape[1])
    pvals = np.ones(X.shape[1])
    for i in range(X.shape[1]):
        try:
            r, p = spearmanr(X[:, i], pt)
            corrs[i] = r
            pvals[i] = p
        except Exception:
            pass

    _, fdr, _, _ = multipletests(pvals, method="fdr_bh")
    df = pd.DataFrame({
        "gene": genes,
        "correlation": corrs,
        "pval": pvals,
        "fdr": fdr,
    })
    df["direction"] = np.where(df["correlation"] > 0, "up", "down")
    df = df[(df["fdr"] < fdr_threshold) & (df["correlation"].abs() >= min_abs_corr)]
    df = df.sort_values("correlation", key=abs, ascending=False).reset_index(drop=True)
    return df


def _run_rna_velocity(adata) -> Dict[str, Any]:
    try:
        import scvelo as scv
    except ImportError:
        raise ImportError("scvelo")

    scv.settings.verbosity = 1
    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    mode_used = "dynamical"
    try:
        scv.tl.recover_dynamics(adata, n_jobs=4)
        scv.tl.velocity(adata, mode="dynamical")
    except Exception as exc:
        print(f"  dynamical model failed ({exc}); falling back to stochastic")
        scv.tl.velocity(adata, mode="stochastic")
        mode_used = "stochastic"
    scv.tl.velocity_graph(adata)
    if mode_used == "dynamical":
        scv.tl.latent_time(adata)
        scv.tl.velocity_pseudotime(adata)
    scv.tl.rank_velocity_genes(adata, groupby="clusters" if "clusters" in adata.obs.columns else None,
                                min_corr=0.3)
    velocity_genes_df = pd.DataFrame(adata.uns["rank_velocity_genes"]["names"])

    return {
        "mode": mode_used,
        "velocity_genes": velocity_genes_df,
        "has_latent_time": "latent_time" in adata.obs.columns,
    }


def _run_cellrank(adata, cluster_key, has_velocity: bool) -> Dict[str, Any]:
    try:
        import cellrank as cr
    except ImportError:
        raise ImportError("cellrank")

    if has_velocity:
        try:
            from cellrank.kernels import VelocityKernel, ConnectivityKernel
            vk = VelocityKernel(adata).compute_transition_matrix()
            ck = ConnectivityKernel(adata).compute_transition_matrix()
            kernel = 0.8 * vk + 0.2 * ck
        except Exception:
            from cellrank.kernels import PseudotimeKernel
            kernel = PseudotimeKernel(adata, time_key="dpt_pseudotime").compute_transition_matrix()
    else:
        from cellrank.kernels import PseudotimeKernel
        kernel = PseudotimeKernel(adata, time_key="dpt_pseudotime").compute_transition_matrix()

    from cellrank.estimators import GPCCA
    g = GPCCA(kernel)
    g.fit(n_states=10, cluster_key=cluster_key)
    g.predict_terminal_states()
    g.compute_fate_probabilities()

    fate_probs = pd.DataFrame(
        adata.obsm["lineages_fwd"],
        index=adata.obs_names,
        columns=g.terminal_states.cat.categories.tolist()
            if hasattr(g.terminal_states, "cat") else None,
    )

    drivers = {}
    for ts in g.terminal_states.cat.categories:
        try:
            drv = g.compute_lineage_drivers(lineages=ts)
            drv["fate"] = ts
            drivers[ts] = drv
        except Exception as exc:
            print(f"  driver gene computation failed for {ts}: {exc}")

    return {
        "fate_probabilities": fate_probs,
        "drivers": drivers,
        "terminal_states": list(g.terminal_states.cat.categories) if hasattr(g.terminal_states, "cat") else [],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad input")
    parser.add_argument("--output", required=True, help="Pickle output for results")
    parser.add_argument("--root-cell-type", default="Ductal")
    parser.add_argument("--cluster-key", default="clusters")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-dcs", type=int, default=10)
    parser.add_argument("--no-velocity", action="store_true")
    parser.add_argument("--no-cellrank", action="store_true")
    parser.add_argument("--adata-out", help="Optional: write annotated .h5ad")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    results = run_trajectory(
        adata,
        root_cell_type=args.root_cell_type,
        cluster_key=args.cluster_key,
        n_neighbors=args.n_neighbors,
        n_pcs=args.n_pcs,
        n_dcs=args.n_dcs,
        run_velocity=not args.no_velocity,
        run_cellrank=not args.no_cellrank,
    )
    with open(args.output, "wb") as f:
        pickle.dump(results, f)
    print(f"✓ Wrote {args.output}")
    if args.adata_out:
        adata.write_h5ad(args.adata_out)
        print(f"✓ Wrote {args.adata_out}")


if __name__ == "__main__":
    main()
