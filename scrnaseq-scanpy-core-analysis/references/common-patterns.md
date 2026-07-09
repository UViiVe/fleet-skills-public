# Common Patterns

Compact, copy-paste-ready recipes for common analysis scenarios.

## Pattern 1: Standard PBMC, single batch

```python
from scripts.setup_and_import import import_10x_data
from scripts.qc_metrics import calculate_qc_metrics
from scripts.filter_cells import run_scrublet_detection, filter_by_mad_outliers
from scripts.normalize_data import run_standard_normalization
from scripts.find_variable_genes import find_highly_variable_genes
from scripts.scale_and_pca import scale_data, run_pca_analysis
from scripts.cluster_cells import build_neighbor_graph, cluster_leiden_multiple_resolutions
from scripts.run_umap import run_umap_reduction
from scripts.find_markers import find_all_cluster_markers
from scripts.annotate_celltypes import annotate_clusters_manual
from scripts.plot_dimreduction import plot_umap_clusters
from scripts.export_results import export_anndata_results

adata = import_10x_data("filtered_feature_bc_matrix/")
adata = calculate_qc_metrics(adata, species="human")
adata = run_scrublet_detection(adata)
adata = filter_by_mad_outliers(adata, remove_doublets=True)
adata = run_standard_normalization(adata, target_sum=1e4)
adata = find_highly_variable_genes(adata, n_top_genes=2000)
adata = scale_data(adata, vars_to_regress=["total_counts", "pct_counts_mt"])
adata = run_pca_analysis(adata, n_pcs=50)
adata = build_neighbor_graph(adata, use_rep="X_pca", n_neighbors=10, n_pcs=30)
adata = cluster_leiden_multiple_resolutions(adata, resolutions=[0.4, 0.6, 0.8, 1.0])
adata = run_umap_reduction(adata, n_neighbors=10)

markers = find_all_cluster_markers(adata, cluster_key="leiden_0.8")

annotations = {
    "0": "CD4 T cells", "1": "CD14+ Monocytes", "2": "CD8 T cells",
    "3": "B cells", "4": "NK cells", "5": "CD16+ Monocytes",
    "6": "Dendritic cells", "7": "Megakaryocytes",
}
adata = annotate_clusters_manual(adata, annotations, cluster_key="leiden_0.8")
plot_umap_clusters(adata, cluster_key="cell_type", output_dir="results/umap")
export_anndata_results(adata, output_dir="results", cluster_key="cell_type")
```

## Pattern 2: Multi-batch with scVI

```python
import scanpy as sc
from scripts.setup_and_import import import_10x_data
from scripts.qc_metrics import calculate_qc_metrics, batch_mad_outlier_detection
from scripts.filter_cells import run_scrublet_detection, filter_by_mad_outliers
from scripts.normalize_data import run_standard_normalization
from scripts.find_variable_genes import find_highly_variable_genes
from scripts.scale_and_pca import scale_data, run_pca_analysis
from scripts.integrate_scvi import run_scvi_integration
from scripts.integration_diagnostics import compute_lisi_scores, plot_integration_metrics
from scripts.cluster_cells import build_neighbor_graph, cluster_leiden_multiple_resolutions
from scripts.run_umap import run_umap_reduction
from scripts.find_markers import find_all_cluster_markers
from scripts.export_results import export_anndata_results

batches = []
for name in ["batch1", "batch2", "batch3"]:
    ad = import_10x_data(f"data/{name}/filtered_feature_bc_matrix/")
    ad.obs["batch"] = name
    batches.append(ad)
adata = sc.concat(batches, join="outer")

adata = calculate_qc_metrics(adata, species="human")
adata = batch_mad_outlier_detection(adata, batch_key="batch", nmads=5)
adata = run_scrublet_detection(adata, batch_key="batch")
adata = filter_by_mad_outliers(adata, remove_doublets=True)

adata = run_standard_normalization(adata, target_sum=1e4)
adata = find_highly_variable_genes(adata, n_top_genes=2000, batch_key="batch")
adata = scale_data(adata, vars_to_regress=["total_counts", "pct_counts_mt"])
adata = run_pca_analysis(adata, n_pcs=50)

adata = run_scvi_integration(adata, batch_key="batch", n_latent=30, max_epochs=400)
lisi = compute_lisi_scores(adata, batch_key="batch", use_rep="X_scVI")
plot_integration_metrics(adata, lisi, output_dir="results/integration")

adata = build_neighbor_graph(adata, use_rep="X_scVI", n_neighbors=10)
adata = cluster_leiden_multiple_resolutions(adata, resolutions=[0.4, 0.6, 0.8, 1.0])
adata = run_umap_reduction(adata)

markers = find_all_cluster_markers(adata, cluster_key="leiden_0.8")
export_anndata_results(adata, output_dir="results", cluster_key="leiden_0.8")
```

## Pattern 3: Multi-batch with Harmony (CPU)

Same as Pattern 2, but replace the scVI step with:

```python
from scripts.integrate_scvi import run_harmony_integration

adata = run_harmony_integration(adata, batch_key="batch", basis="X_pca", n_pcs=30, theta=2.0)
adata = build_neighbor_graph(adata, use_rep="X_pca_harmony", n_neighbors=10)
```

## Pattern 4: Raw data with ambient correction

```python
from scripts.remove_ambient_rna import run_cellbender, estimate_contamination
from scripts.qc_metrics import calculate_qc_metrics
from scripts.filter_cells import run_scrublet_detection, filter_by_mad_outliers

adata = run_cellbender(
    raw_h5="raw_feature_bc_matrix.h5",
    expected_cells=10000,
    total_droplets=30000,
    output_dir="results/ambient",
)
print(f"Contamination: {estimate_contamination(adata):.1%}")

adata = calculate_qc_metrics(adata, species="human")
adata = run_scrublet_detection(adata)
adata = filter_by_mad_outliers(adata, remove_doublets=True)
# ... continue with normalization, clustering, etc.
```

## Pattern 5: Condition comparison with pseudobulk DE

```python
import scanpy as sc
from scripts.pseudobulk_de import aggregate_to_pseudobulk, run_deseq2_analysis, plot_volcano

# assume adata has 'sample_id', 'condition', 'batch', and 'cell_type' set.

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

for ct, res in de_results.items():
    plot_volcano(res, celltype=ct, padj_threshold=0.05, log2fc_threshold=0.5,
                 output_dir="results/pseudobulk_de")
    sig = res[res["padj"] < 0.05]
    sig.to_csv(f"results/pseudobulk_de/{ct}_sig_genes.csv")
```

Requires N ≥ 2 per condition. The script will refuse to run otherwise.

## Pattern 6: Automated annotation with CellTypist

```python
from scripts.annotate_celltypes import annotate_with_celltypist
from scripts.plot_dimreduction import plot_umap_clusters

adata = annotate_with_celltypist(adata, model="Immune_All_Low.pkl", majority_voting=True)
plot_umap_clusters(adata, cluster_key="predicted_labels", output_dir="results/annotation")
```

## Pattern 7: Hybrid annotation

```python
from scripts.annotate_celltypes import annotate_with_celltypist, annotate_clusters_manual
from scripts.find_markers import find_all_cluster_markers
from scripts.plot_dimreduction import plot_feature_umap

adata = annotate_with_celltypist(adata, model="Immune_All_Low.pkl")
markers = find_all_cluster_markers(adata, cluster_key="leiden_0.8")

for cluster in adata.obs["leiden_0.8"].unique():
    rows = markers[markers["cluster"] == cluster].head(10)
    print(f"\nCluster {cluster}:")
    print(rows[["gene", "logfoldchange", "pval_adj"]])

plot_feature_umap(adata, features=["CD3D", "CD4", "CD8A", "CD14", "MS4A1", "NKG7", "FCGR3A"],
                  output_dir="results/markers")

refined = {
    "0": "CD4 T cells", "1": "CD14+ Monocytes", "2": "CD8 T cells",
    "3": "B cells", "4": "NK cells", "5": "CD16+ Monocytes",
}
adata = annotate_clusters_manual(adata, refined, cluster_key="leiden_0.8")
```

## Pattern 8: Pearson residuals normalization

```python
from scripts.normalize_data import run_pearson_residuals
from scripts.find_variable_genes import find_highly_variable_genes
from scripts.scale_and_pca import run_pca_analysis

adata = run_pearson_residuals(adata, theta=100, clip=None)
adata = find_highly_variable_genes(adata, n_top_genes=2000)
adata = run_pca_analysis(adata, n_pcs=50)  # no scaling needed; Pearson residuals are already scaled
```

## Pattern 9: Tissue-specific QC

```python
from scripts.qc_metrics import calculate_qc_metrics, mark_outliers_fixed
from scripts.filter_cells import run_scrublet_detection, filter_cells_by_tissue
from scripts.plot_qc import plot_qc_violin

adata = calculate_qc_metrics(adata, species="human")
plot_qc_violin(adata, output_dir="results/qc")

adata = mark_outliers_fixed(adata, tissue="brain", min_genes=500, max_genes=6000, max_mt=10.0)
adata = run_scrublet_detection(adata)
adata = filter_cells_by_tissue(adata, tissue="brain", remove_doublets=True)
```

## Pattern 10: scANVI semi-supervised integration

```python
from scripts.integrate_scvi import run_scanvi_integration
from scripts.cluster_cells import build_neighbor_graph, cluster_leiden_multiple_resolutions

# Suppose batch1 was already annotated; mark others as "Unknown".
adata.obs["cell_type_coarse"] = "Unknown"
mask = adata.obs["batch"] == "batch1"
adata.obs.loc[mask, "cell_type_coarse"] = adata.obs.loc[mask, "prior_labels"]

adata = run_scanvi_integration(
    adata, batch_key="batch", labels_key="cell_type_coarse",
    unlabeled_category="Unknown", n_latent=30, max_epochs=400,
)

adata = build_neighbor_graph(adata, use_rep="X_scANVI", n_neighbors=10)
adata = cluster_leiden_multiple_resolutions(adata, resolutions=[0.4, 0.6, 0.8])
```
