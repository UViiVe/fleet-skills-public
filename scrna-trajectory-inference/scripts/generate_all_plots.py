"""Render the standard set of trajectory plots.

Plots PAGA, pseudotime UMAP, diffusion components, gene heatmap, gene trends,
plus velocity / fate / driver plots when results are available.

Inputs:  AnnData + results dict from `run_trajectory_analysis.py`.
Outputs: PNG (and SVG when possible) figures under `output_dir`.
Example:
    python generate_all_plots.py --input adata.h5ad --results results.pkl \\
        --output-dir trajectory_results --cluster-key clusters
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def _save(fig, base_path: Path) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base_path, dpi=300, bbox_inches="tight")
    try:
        fig.savefig(base_path.with_suffix(".svg"), bbox_inches="tight")
    except Exception:
        pass
    print(f"  ✓ {base_path}")


def generate_all_plots(
    adata: "AnnData",
    results: Dict[str, Any],
    output_dir: str = "trajectory_results",
    cluster_key: str = "clusters",
) -> None:
    """Generate every plot supported by the results dict."""
    try:
        import matplotlib.pyplot as plt
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install matplotlib scanpy")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("Generating PAGA graph plot...")
    try:
        fig = sc.pl.paga_compare(adata, basis="umap", show=False, return_fig=True)
        _save(fig, out / "paga_graph.png")
        plt.close(fig)
    except Exception as exc:
        print(f"  warning: PAGA plot failed: {exc}")

    print("Generating pseudotime UMAP...")
    try:
        fig = sc.pl.umap(adata, color="dpt_pseudotime", cmap="viridis",
                         show=False, return_fig=True)
        _save(fig, out / "pseudotime_umap.png")
        plt.close(fig)
    except Exception as exc:
        print(f"  warning: pseudotime UMAP failed: {exc}")

    print("Generating pseudotime violin plot...")
    try:
        fig = sc.pl.violin(adata, keys="dpt_pseudotime", groupby=cluster_key,
                           rotation=45, show=False, return_fig=True)
        _save(fig, out / "pseudotime_violin.png")
        plt.close(fig)
    except Exception as exc:
        print(f"  warning: pseudotime violin failed: {exc}")

    print("Generating diffusion components plot...")
    try:
        sc.pl.diffmap(adata, color=cluster_key, show=False)
        fig = plt.gcf()
        _save(fig, out / "diffusion_components.png")
        plt.close(fig)
    except Exception as exc:
        print(f"  warning: diffusion components plot failed: {exc}")

    print("Generating gene trend / heatmap plots...")
    try:
        traj_genes = results.get("trajectory_genes")
        if isinstance(traj_genes, pd.DataFrame) and not traj_genes.empty:
            top = traj_genes.head(30)["gene"].tolist()
            top = [g for g in top if g in adata.var_names]
            if top:
                fig = sc.pl.heatmap(adata, var_names=top, groupby=cluster_key,
                                     show=False, return_fig=True if hasattr(sc.pl.heatmap, "return_fig") else False)
                if fig is None:
                    fig = plt.gcf()
                _save(fig, out / "gene_heatmap.png")
                plt.close(fig)

                _plot_gene_trends(adata, top[:6], out / "gene_trends.png")
    except Exception as exc:
        print(f"  warning: trajectory-gene plot failed: {exc}")

    print("Generating PAGA connectivity heatmap...")
    try:
        conn = results["paga_connectivities"]
        labels = list(adata.obs[cluster_key].astype(str).unique())
        fig, ax = plt.subplots(figsize=(8, 6))
        im = ax.imshow(conn, cmap="viridis")
        ax.set_xticks(range(len(labels)))
        ax.set_yticks(range(len(labels)))
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_title("PAGA connectivity")
        plt.colorbar(im, ax=ax)
        _save(fig, out / "paga_connectivity.png")
        plt.close(fig)
    except Exception as exc:
        print(f"  warning: PAGA connectivity heatmap failed: {exc}")

    if results.get("velocity_results"):
        _plot_velocity(adata, results, out)

    if results.get("cellrank_results"):
        _plot_cellrank(adata, results, out, cluster_key=cluster_key)

    print("\n✓ All plots generated successfully!")


def _plot_gene_trends(adata, genes, output_path: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    pt = adata.obs["dpt_pseudotime"].values
    finite = np.isfinite(pt)
    pt = pt[finite]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, gene in zip(axes.ravel(), genes):
        if gene not in adata.var_names:
            ax.set_visible(False)
            continue
        vals = adata[:, gene].X[finite, :]
        if hasattr(vals, "toarray"):
            vals = vals.toarray()
        vals = np.asarray(vals).ravel()
        order = np.argsort(pt)
        ax.scatter(pt[order], vals[order], s=1, alpha=0.3, color="grey")
        try:
            bins = np.linspace(pt.min(), pt.max(), 30)
            digit = np.digitize(pt, bins)
            means = [vals[digit == i].mean() if (digit == i).any() else np.nan for i in range(1, len(bins))]
            ax.plot(bins[:-1], means, color="crimson", lw=2)
        except Exception:
            pass
        ax.set_xlabel("pseudotime")
        ax.set_ylabel("expression")
        ax.set_title(gene)
    plt.tight_layout()
    _save(fig, output_path)
    plt.close(fig)


def _plot_velocity(adata, results, output_dir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import scvelo as scv
    except ImportError:
        return
    try:
        scv.pl.velocity_embedding_stream(adata, basis="umap", show=False)
        fig = plt.gcf()
        _save(fig, output_dir / "velocity_stream.png")
        plt.close(fig)
    except Exception as exc:
        print(f"  warning: velocity_stream failed: {exc}")
    try:
        scv.pl.velocity_confidence(adata, show=False)
        fig = plt.gcf()
        _save(fig, output_dir / "velocity_confidence.png")
        plt.close(fig)
    except Exception as exc:
        print(f"  warning: velocity_confidence failed: {exc}")
    if results["velocity_results"].get("has_latent_time"):
        try:
            scv.pl.scatter(adata, color="latent_time", color_map="gnuplot", show=False)
            fig = plt.gcf()
            _save(fig, output_dir / "latent_time.png")
            plt.close(fig)
        except Exception as exc:
            print(f"  warning: latent_time plot failed: {exc}")
    try:
        velocity_genes = results["velocity_results"]["velocity_genes"]
        first_col = velocity_genes.columns[0]
        top = velocity_genes[first_col].head(6).tolist()
        top = [g for g in top if g in adata.var_names]
        if top:
            scv.pl.scatter(adata, var_names=top, ncols=3, show=False)
            fig = plt.gcf()
            _save(fig, output_dir / "velocity_top_genes.png")
            plt.close(fig)
    except Exception as exc:
        print(f"  warning: top velocity genes plot failed: {exc}")


def _plot_cellrank(adata, results, output_dir: Path, cluster_key: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return
    fates = results["cellrank_results"]["fate_probabilities"]
    if fates is not None and not fates.empty:
        try:
            import scanpy as sc
            for fate in fates.columns:
                col_name = f"fate_{fate}"
                adata.obs[col_name] = fates[fate].values
                fig = sc.pl.umap(adata, color=col_name, cmap="viridis",
                                  show=False, return_fig=True)
                _save(fig, output_dir / f"fate_{fate}.png")
                plt.close(fig)
        except Exception as exc:
            print(f"  warning: fate UMAP failed: {exc}")
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            im = ax.imshow(fates.values.T, aspect="auto", cmap="viridis")
            ax.set_xlabel("cells (sorted)")
            ax.set_yticks(range(fates.shape[1]))
            ax.set_yticklabels(fates.columns)
            ax.set_title("Fate probabilities")
            plt.colorbar(im, ax=ax)
            _save(fig, output_dir / "fate_heatmap.png")
            plt.close(fig)
        except Exception as exc:
            print(f"  warning: fate heatmap failed: {exc}")

    drivers = results["cellrank_results"].get("drivers", {})
    for fate, drv in drivers.items():
        try:
            top = drv.head(15)
            fig, ax = plt.subplots(figsize=(6, 6))
            corr_col = [c for c in top.columns if "corr" in c.lower()]
            if not corr_col: continue
            ax.barh(range(len(top)), top[corr_col[0]], color="steelblue")
            ax.set_yticks(range(len(top)))
            ax.set_yticklabels(top.index, fontsize=8)
            ax.invert_yaxis()
            ax.set_xlabel("correlation with fate prob.")
            ax.set_title(f"Drivers — {fate}")
            _save(fig, output_dir / f"driver_genes_{fate}.png")
            plt.close(fig)
        except Exception as exc:
            print(f"  warning: drivers plot for {fate} failed: {exc}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad input")
    parser.add_argument("--results", required=True, help="Pickle from run_trajectory_analysis")
    parser.add_argument("--output-dir", default="trajectory_results")
    parser.add_argument("--cluster-key", default="clusters")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    with open(args.results, "rb") as f:
        results = pickle.load(f)
    generate_all_plots(adata, results,
                        output_dir=args.output_dir,
                        cluster_key=args.cluster_key)


if __name__ == "__main__":
    main()
