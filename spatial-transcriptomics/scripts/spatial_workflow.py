"""Spatial transcriptomics analysis workflow (10x Visium).

Runs the full Visium pipeline:

    1. QC metrics (mitochondrial, total counts, genes per spot)
    2. Spot/gene filtering
    3. Total-count normalization + log1p
    4. Highly variable gene selection (seurat_v3 on raw counts)
    5. PCA, neighbor graph, UMAP, Leiden clustering
    6. Spatial neighbors graph (hex grid)
    7. Spatially variable genes (Moran's I)
    8. Neighborhood enrichment
    9. Co-occurrence

Inputs : raw Visium AnnData (from load_example_data.py)
Outputs: processed AnnData with spatial statistics in .uns

Example:
    python spatial_workflow.py --input adata.h5ad --output-dir visium_results
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Optional


def run_spatial_analysis(
    adata: "anndata.AnnData",
    min_genes: int = 200,
    min_cells: int = 10,
    max_pct_mito: float = 50.0,
    n_top_genes: int = 2000,
    n_pcs: int = 30,
    resolution: float = 0.8,
    n_neighbors: int = 15,
    svgs_n_perms: int = 100,
    svgs_n_jobs: int = 1,
    output_dir: str = "visium_results",
) -> "anndata.AnnData":
    """Run the complete spatial transcriptomics pipeline.

    Parameters
    ----------
    adata
        Raw Visium AnnData with ``.obsm['spatial']`` populated.
    min_genes
        Minimum genes per retained spot.
    min_cells
        Minimum spots per retained gene.
    max_pct_mito
        Maximum mitochondrial percent per retained spot.
    n_top_genes
        Number of highly variable genes to use for PCA/clustering.
    n_pcs
        PCA dimensionality.
    resolution
        Leiden resolution.
    n_neighbors
        Expression-graph neighbor count.
    svgs_n_perms
        Permutations for the Moran's I null distribution.
    svgs_n_jobs
        Parallel jobs for Moran's I.
    output_dir
        Where to write artefacts; created if missing.

    Returns
    -------
    anndata.AnnData
        Processed AnnData with clustering, SVG results, enrichment, and
        co-occurrence attached.
    """
    try:
        import scanpy as sc
        import squidpy as sq
    except ImportError:
        sys.exit("Install with: pip install scanpy squidpy")

    os.makedirs(output_dir, exist_ok=True)

    print("\n=== Spatial transcriptomics analysis ===\n")
    n_initial = adata.n_obs

    adata = _qc_metrics(adata, sc)
    adata = _filter_spots(adata, sc, min_genes, min_cells, max_pct_mito, n_initial)
    adata = _normalize(adata, sc)
    adata = _hvg(adata, sc, n_top_genes)
    adata = _dimred_cluster(adata, sc, n_pcs, n_neighbors, resolution)
    adata = _spatial_graph(adata, sq)
    adata, n_svgs = _svgs(adata, sq, svgs_n_perms, svgs_n_jobs)
    n_clusters = adata.obs["leiden"].nunique()
    adata = _enrichment_and_cooccurrence(adata, sq, n_clusters)

    adata.uns["spatial_analysis_params"] = {
        "min_genes": min_genes,
        "min_cells": min_cells,
        "max_pct_mito": max_pct_mito,
        "n_top_genes": n_top_genes,
        "n_pcs": n_pcs,
        "resolution": resolution,
        "n_neighbors": n_neighbors,
        "svgs_n_perms": svgs_n_perms,
        "n_spots_initial": int(n_initial),
        "n_spots_after": int(adata.n_obs),
        "n_clusters": int(n_clusters),
        "n_svgs": int(n_svgs),
    }

    print("\n" + "=" * 50)
    print("✓ Spatial analysis completed successfully!")
    print("=" * 50)
    print(f"  Spots analyzed : {adata.n_obs}")
    print(f"  Clusters found : {n_clusters}")
    print(f"  SVGs identified: {n_svgs}")
    return adata


# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------

def _qc_metrics(adata, sc):
    print("1. QC metrics")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    print(f"   median genes/spot : {adata.obs['n_genes_by_counts'].median():.0f}")
    print(f"   median counts/spot: {adata.obs['total_counts'].median():.0f}")
    print(f"   median %MT        : {adata.obs['pct_counts_mt'].median():.1f}%")
    return adata


def _filter_spots(adata, sc, min_genes, min_cells, max_pct_mito, n_initial):
    print(f"\n2. Filter spots (min_genes={min_genes}, max_mito={max_pct_mito}%)")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs["pct_counts_mt"] < max_pct_mito, :].copy()
    retention = 100 * adata.n_obs / n_initial if n_initial else 0
    print(f"   spots: {n_initial} -> {adata.n_obs} ({retention:.1f}% retained)")
    return adata


def _normalize(adata, sc):
    print("\n3. Normalize (target_sum=1e4 + log1p)")
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["normalized"] = adata.X.copy()
    return adata


def _hvg(adata, sc, n_top_genes):
    print(f"\n4. HVG selection (top {n_top_genes}, seurat_v3 on raw counts)")
    sc.pp.highly_variable_genes(
        adata, n_top_genes=n_top_genes, flavor="seurat_v3", layer="counts"
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    print(f"   highly variable genes: {n_hvg}")
    return adata


def _dimred_cluster(adata, sc, n_pcs, n_neighbors, resolution):
    print("\n5. PCA + neighbors + UMAP + Leiden")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution, key_added="leiden")
    n_clusters = adata.obs["leiden"].nunique()
    print(f"   PCA components : {n_pcs}")
    print(f"   Leiden clusters: {n_clusters} (resolution={resolution})")
    return adata


def _spatial_graph(adata, sq):
    print("\n6. Spatial neighbors graph (grid, n_neighs=6)")
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)
    n_edges = adata.obsp["spatial_connectivities"].nnz
    print(f"   spatial graph edges: {n_edges}")
    return adata


def _svgs(adata, sq, n_perms, n_jobs):
    print(f"\n7. Spatially variable genes (Moran's I, n_perms={n_perms})")
    # Restore normalized expression (scale() left .X scaled).
    adata.X = adata.layers["normalized"].copy()
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=n_perms, n_jobs=n_jobs)

    moran = adata.uns["moranI"].copy()
    svgs = moran[moran["pval_norm_fdr_bh"] < 0.05].sort_values("I", ascending=False)
    n_svgs = len(svgs)
    print(f"   significant SVGs (FDR < 0.05): {n_svgs}")
    if n_svgs:
        top5 = svgs.index[:5].tolist()
        print(f"   top SVGs: {', '.join(top5)}")
    adata.uns["svg_results"] = svgs
    return adata, n_svgs


def _enrichment_and_cooccurrence(adata, sq, n_clusters):
    print("\n8. Neighborhood enrichment")
    sq.gr.nhood_enrichment(adata, cluster_key="leiden")
    print(f"   enrichment matrix: {n_clusters} x {n_clusters}")

    print("\n9. Co-occurrence")
    sq.gr.co_occurrence(adata, cluster_key="leiden")
    print(f"   co-occurrence computed for {n_clusters} clusters")
    return adata


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the full 10x Visium spatial transcriptomics workflow."
    )
    parser.add_argument("--input", required=True, help="Path to input .h5ad file.")
    parser.add_argument("--output-dir", default="visium_results", help="Output directory.")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=10)
    parser.add_argument("--max-pct-mito", type=float, default=50.0)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--resolution", type=float, default=0.8)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--svgs-n-perms", type=int, default=100)
    parser.add_argument("--svgs-n-jobs", type=int, default=1)
    args = parser.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy squidpy")

    adata = sc.read_h5ad(args.input)
    adata = run_spatial_analysis(
        adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_pct_mito=args.max_pct_mito,
        n_top_genes=args.n_top_genes,
        n_pcs=args.n_pcs,
        resolution=args.resolution,
        n_neighbors=args.n_neighbors,
        svgs_n_perms=args.svgs_n_perms,
        svgs_n_jobs=args.svgs_n_jobs,
        output_dir=args.output_dir,
    )

    out_h5ad = os.path.join(args.output_dir, "adata_processed.h5ad")
    adata.write_h5ad(out_h5ad, compression="gzip")
    print(f"\n✓ workflow completed: wrote {out_h5ad}")


if __name__ == "__main__":
    main()
