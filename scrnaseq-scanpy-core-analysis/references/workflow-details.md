# Workflow Details

Expanded code for every step of the standard scanpy pipeline. The `scripts/` directory is the canonical source — this reference shows how the pieces fit together.

## Phase 1: QC and filtering

### Step 1.1 — Optional: ambient RNA correction (raw input only)

Skip this step if you're starting from `filtered_feature_bc_matrix/`. Otherwise, see `ambient_rna_correction.md` for full guidance.

CellBender:

```python
from scripts.remove_ambient_rna import run_cellbender, estimate_contamination

adata = run_cellbender(
    raw_h5="raw_feature_bc_matrix.h5",
    expected_cells=10000,
    total_droplets=30000,
    epochs=200,
    output_dir="results/ambient",
)
print(f"Contamination: {estimate_contamination(adata):.1%}")
```

SoupX (Python wrapper):

```python
from scripts.remove_ambient_rna import run_soupx_python

adata = run_soupx_python(
    raw_matrix_dir="raw_feature_bc_matrix/",
    filtered_matrix_dir="filtered_feature_bc_matrix/",
)
```

### Step 1.2 — Load data

```python
from scripts.setup_and_import import (
    setup_scanpy_environment, import_10x_data, import_h5_data, import_count_matrix
)

setup_scanpy_environment()
adata = import_10x_data("filtered_feature_bc_matrix/")
# or
adata = import_h5_data("filtered_feature_bc_matrix.h5")
# or
adata = import_count_matrix("counts.csv", transpose=True)
# or
import scanpy as sc; adata = sc.read_h5ad("data.h5ad")
```

### Step 1.3 — QC metrics

```python
from scripts.qc_metrics import calculate_qc_metrics, batch_mad_outlier_detection, mark_outliers_fixed

adata = calculate_qc_metrics(adata, species="human")  # adds n_genes_by_counts, total_counts, pct_counts_mt, pct_counts_ribo, log10_*

if "batch" in adata.obs.columns:
    adata = batch_mad_outlier_detection(
        adata,
        batch_key="batch",
        metrics=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt"],
        nmads=5,
    )
else:
    adata = mark_outliers_fixed(adata, tissue="pbmc")  # or "brain", "tumor", etc.
```

### Step 1.4 — QC visualization

```python
from scripts.plot_qc import plot_qc_violin, plot_qc_scatter, plot_qc_histograms

plot_qc_violin(adata, batch_key="batch", output_dir="results/qc")
plot_qc_scatter(adata, output_dir="results/qc")
plot_qc_histograms(adata, output_dir="results/qc")
```

### Step 1.5 — Doublet detection and filtering

```python
from scripts.filter_cells import run_scrublet_detection, filter_by_mad_outliers, filter_cells_by_tissue

adata = run_scrublet_detection(adata, batch_key="batch")
adata_filtered = filter_by_mad_outliers(adata, remove_doublets=True, doublet_score_threshold=0.25)

print(f"Retained: {adata_filtered.n_obs}/{adata.n_obs} ({100 * adata_filtered.n_obs / adata.n_obs:.1f}%)")
```

Retention should be > 70%. If lower, raise `nmads` or revisit tissue-specific thresholds.

## Phase 2: Normalization and feature selection

### Step 2.1 — Normalize

```python
from scripts.normalize_data import run_standard_normalization

adata = run_standard_normalization(adata, target_sum=1e4, log_transform=True)
# stores raw counts in adata.layers["counts"], normalized in adata.X
```

Alternative — Pearson residuals:

```python
from scripts.normalize_data import run_pearson_residuals
adata = run_pearson_residuals(adata, theta=100, clip=None)
```

### Step 2.2 — Highly variable genes

```python
from scripts.find_variable_genes import find_highly_variable_genes

adata = find_highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3", batch_key="batch")
```

### Step 2.3 — Scale and PCA

```python
from scripts.scale_and_pca import scale_data, run_pca_analysis, suggest_n_pcs, plot_pca_variance

adata = scale_data(adata, max_value=10, vars_to_regress=["total_counts", "pct_counts_mt"])
adata = run_pca_analysis(adata, n_pcs=50)
plot_pca_variance(adata, output_dir="results/pca")
n_pcs = suggest_n_pcs(adata)
print(f"Recommended n_pcs: {n_pcs}")
```

## Phase 3: Integration (multi-batch only)

### Step 3.1 — scVI

```python
from scripts.integrate_scvi import run_scvi_integration

adata = run_scvi_integration(
    adata,
    batch_key="batch",
    condition_key="condition",   # optional, enables N=1 confounding check
    n_latent=30,
    n_layers=2,
    max_epochs=400,
    early_stopping=True,
    use_gpu=True,
)
```

### Step 3.2 — scANVI (when partial labels exist)

```python
from scripts.integrate_scvi import run_scanvi_integration

adata = run_scanvi_integration(
    adata,
    batch_key="batch",
    labels_key="cell_type_coarse",
    unlabeled_category="Unknown",
    n_latent=30,
    max_epochs=400,
)
```

### Step 3.3 — Harmony (fast, CPU)

```python
from scripts.integrate_scvi import run_harmony_integration

adata = run_harmony_integration(
    adata,
    batch_key="batch",
    basis="X_pca",
    n_pcs=30,
    theta=2.0,
    max_iter=10,
)
```

### Step 3.4 — Integration diagnostics

```python
from scripts.integration_diagnostics import (
    compute_lisi_scores, compute_asw_scores,
    plot_integration_metrics, compare_integration_methods,
)

lisi = compute_lisi_scores(adata, batch_key="batch", label_key="cell_type", use_rep="X_scVI")
asw = compute_asw_scores(adata, batch_key="batch", label_key="cell_type", use_rep="X_scVI")
plot_integration_metrics(adata, lisi, asw, output_dir="results/integration")
compare_integration_methods(adata, batch_key="batch", representations=["X_pca", "X_scVI"], output_dir="results/integration")
```

Targets: batch iLISI ≈ `n_batches`, cell-type cLISI ≈ 1, batch ASW ≈ 0, cell-type ASW > 0.5.

## Phase 4: Clustering and visualization

### Step 4.1 — Neighbor graph and clustering

```python
from scripts.cluster_cells import build_neighbor_graph, cluster_leiden_multiple_resolutions

use_rep = "X_scVI" if "X_scVI" in adata.obsm else "X_pca"
adata = build_neighbor_graph(adata, n_neighbors=10, n_pcs=n_pcs, use_rep=use_rep, metric="euclidean")
adata = cluster_leiden_multiple_resolutions(adata, resolutions=[0.4, 0.6, 0.8, 1.0], random_state=42)
```

### Step 4.2 — UMAP

```python
from scripts.run_umap import run_umap_reduction

adata = run_umap_reduction(adata, n_neighbors=10, min_dist=0.5, spread=1.0, random_state=42)
```

### Step 4.3 — UMAP plots

```python
from scripts.plot_dimreduction import (
    plot_umap_clusters, plot_clustering_comparison,
    plot_feature_umap, plot_umap_styled,
)

plot_umap_clusters(adata, cluster_key="leiden_0.8", output_dir="results/umap")
plot_clustering_comparison(adata, resolutions=[0.4, 0.6, 0.8, 1.0], output_dir="results/umap")
plot_feature_umap(adata, features=["n_genes_by_counts", "total_counts", "pct_counts_mt"], output_dir="results/umap")
plot_feature_umap(adata, features=["CD3D", "CD8A", "CD4", "CD14", "MS4A1", "NKG7"], output_dir="results/umap")
```

## Phase 5: Marker genes and annotation

### Step 5.1 — Find markers (exploratory)

```python
from scripts.find_markers import (
    find_all_cluster_markers, export_marker_tables,
    plot_top_markers_heatmap, plot_markers_dotplot,
)

markers = find_all_cluster_markers(adata, cluster_key="leiden_0.8", method="wilcoxon",
                                    min_logfoldchange=0.25, min_pct=0.1)
export_marker_tables(markers, output_dir="results/markers")
plot_top_markers_heatmap(adata, markers, n_top=10, cluster_key="leiden_0.8", output_dir="results/markers")
plot_markers_dotplot(adata, markers, n_top=5, output_dir="results/markers")
```

These statistics are for cluster characterization. They are not valid for condition contrasts — use pseudobulk (Phase 6) for that.

### Step 5.2 — Annotate cell types

Manual:

```python
from scripts.annotate_celltypes import annotate_clusters_manual, plot_annotated_umap

annotations = {
    "0": "CD4 T cells",
    "1": "CD14+ Monocytes",
    "2": "CD8 T cells",
    "3": "B cells",
    "4": "NK cells",
    "5": "CD16+ Monocytes",
    "6": "Dendritic cells",
    "7": "Megakaryocytes",
}
adata = annotate_clusters_manual(adata, annotations, cluster_key="leiden_0.8")
plot_annotated_umap(adata, output_dir="results/annotation")
```

CellTypist (automated):

```python
from scripts.annotate_celltypes import annotate_with_celltypist

adata = annotate_with_celltypist(adata, model="Immune_All_Low.pkl", majority_voting=True)
```

Hybrid (recommended):

```python
adata = annotate_with_celltypist(adata, model="Immune_All_Low.pkl")
# Inspect predicted_labels, then refine
refined = {"CD4 T cells": "CD4 T", "CD8 T cells": "CD8 T", ...}
adata = annotate_clusters_manual(adata, refined, cluster_key="predicted_labels")
```

## Phase 6: Inferential analysis (multi-sample, multi-condition)

### Step 6.1 — Pseudobulk DE

```python
from scripts.pseudobulk_de import (
    aggregate_to_pseudobulk, run_deseq2_analysis,
    export_de_results, plot_volcano, plot_ma,
)

pseudobulk = aggregate_to_pseudobulk(
    adata, sample_key="sample_id", celltype_key="cell_type",
    min_cells=10, min_counts=1,
)

de_results = run_deseq2_analysis(
    pseudobulk,
    formula="~ batch + condition",
    contrast=["condition", "treated", "control"],
    celltype_key="cell_type",
    alpha=0.05,
    output_dir="results/pseudobulk_de",
)

export_de_results(de_results, output_dir="results/pseudobulk_de")

for ct, res in de_results.items():
    plot_volcano(res, celltype=ct, padj_threshold=0.05, log2fc_threshold=0.5, output_dir="results/pseudobulk_de")
    plot_ma(res, celltype=ct, output_dir="results/pseudobulk_de")
```

The DE script will refuse to run with N=1 in any group. That's correct behavior — use cell-level Wilcoxon as exploratory analysis only and label it as such.

See `pseudobulk_de_guide.md` for details on design formulas, minimum sample size, and how to interpret results.

## Phase 7: Export

### Step 7.1 — Export everything

```python
from scripts.export_results import export_anndata_results

export_anndata_results(
    adata,
    output_dir="results",
    cluster_key="cell_type",
    export_raw=True,
    export_metadata=True,
    export_embeddings=True,
    export_markers=True,
)
```

Outputs:

- `adata_processed.h5ad` — full AnnData object.
- `expression_matrix_counts.csv`, `expression_matrix_normalized.csv` — full or HVG-only.
- `cell_metadata.csv` — `obs` table.
- `pca_coordinates.csv`, `umap_coordinates.csv` — embeddings.
- `cluster_markers_all.csv` — Wilcoxon ranked genes.
- `analysis_summary.txt` — dataset shape, retention, integration scores.

## Worked example

See `common-patterns.md` for compact end-to-end recipes covering PBMC, multi-batch integration, ambient correction, condition contrasts, automated annotation, and Pearson residuals.
