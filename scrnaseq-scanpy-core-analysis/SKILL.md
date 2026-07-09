---
name: scrnaseq-scanpy-core-analysis
id: "skill_11a42d9fbcea498ca24a60c4323622e4"
display-name: "scRNA-seq Scanpy Core Analysis"
short-description: "End-to-end scRNA-seq analysis in Python using Scanpy and the scverse stack — QC, normalization, integration, clustering, annotation, pseudobulk DE. Use when analyzing 10X / Drop-seq / Smart-seq dat..."
category: transcriptomics
keywords: [scanpy, scRNA-seq, clustering, umap, python]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "End-to-end scRNA-seq analysis in Python using Scanpy and the scverse stack — QC, normalization, integration, clustering, annotation, pseudobulk DE. Use when analyzing 10X / Drop-seq / Smart-seq data and you want a Python (not R) workflow."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "Scanpy (scverse stack)"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Wolf FA et al., Genome Biology 2018, doi:10.1186/s13059-017-1382-0"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Single-Cell RNA-seq Core Analysis (Scanpy)

A complete Python workflow for single-cell RNA-seq, from a raw count matrix to an annotated `AnnData` object with publication-ready figures. The skill leans on the scverse stack: `scanpy` for the core pipeline, `scvi-tools` or `harmonypy` for integration, `scrublet` for doublets, `celltypist` for reference-based annotation, and `pydeseq2`/DESeq2 for proper inferential statistics on multi-sample designs.

## When to Use This Skill

Use this skill when:

- You have 10X Chromium output (CellRanger `filtered_feature_bc_matrix/`, `.h5`), Drop-seq, inDrop, or Smart-seq2 count data and want to take it through QC, clustering, and cell type annotation.
- You're working in Python and want to stay in the scverse ecosystem (`anndata`, `scanpy`, `scvi-tools`).
- You have multiple samples or batches and need integration (scVI, scANVI, or Harmony) plus an honest assessment of batch mixing.
- You have biological replicates and need pseudobulk + DESeq2 to compare conditions (NOT cell-level Wilcoxon for inference).
- You want curated, opinionated defaults rather than re-deriving them from the docs.

Do NOT use this skill for:

- **Bulk RNA-seq.** Use `bulk-rnaseq-counts-to-de-deseq2` instead.
- **R-based scRNA-seq.** Use `scrnaseq-seurat-core-analysis` if you're committed to Seurat/Bioconductor.
- **Spatial transcriptomics.** Use `spatial-transcriptomics` — the QC and neighbor-graph assumptions differ.
- **Trajectory / pseudotime as the primary question.** This skill produces an annotated object; pass it to `scrna-trajectory-inference` for PAGA, Palantir, scVelo, etc.
- **Cell-cell communication.** Use `cell-cell-communication` after this skill finishes annotating cell types.

## Quick Start

```python
from scripts.load_example_data import load_example_data
from scripts.qc_metrics import calculate_qc_metrics, batch_mad_outlier_detection
from scripts.filter_cells import run_scrublet_detection, filter_by_mad_outliers
from scripts.normalize_data import run_standard_normalization
from scripts.find_variable_genes import find_highly_variable_genes
from scripts.scale_and_pca import scale_data, run_pca_analysis
from scripts.cluster_cells import build_neighbor_graph, cluster_leiden_multiple_resolutions
from scripts.run_umap import run_umap_reduction
from scripts.find_markers import find_all_cluster_markers
from scripts.export_results import export_anndata_results

adata = load_example_data("pbmc3k")
adata = calculate_qc_metrics(adata, species="human")
adata = run_scrublet_detection(adata)
adata = filter_by_mad_outliers(adata, remove_doublets=True)
adata = run_standard_normalization(adata, target_sum=1e4)
adata = find_highly_variable_genes(adata, n_top_genes=2000)
adata = scale_data(adata, vars_to_regress=["total_counts", "pct_counts_mt"])
adata = run_pca_analysis(adata, n_pcs=50)
adata = build_neighbor_graph(adata, use_rep="X_pca", n_neighbors=10, n_pcs=30)
adata = cluster_leiden_multiple_resolutions(adata, resolutions=[0.4, 0.6, 0.8, 1.0])
adata = run_umap_reduction(adata)
markers = find_all_cluster_markers(adata, cluster_key="leiden_0.8")
export_anndata_results(adata, output_dir="results", cluster_key="leiden_0.8")
```

Runtime for PBMC 3k on a laptop CPU: ~2 minutes (scVI adds 5–20 min on a GPU, 30+ on CPU).

## Installation

| Package      | Version | License        | Commercial use | Install command                |
| ------------ | ------- | -------------- | -------------- | ------------------------------ |
| scanpy       | ≥ 1.10  | BSD-3-Clause   | Permitted      | `pip install scanpy`           |
| anndata      | ≥ 0.10  | BSD-3-Clause   | Permitted      | `pip install anndata`          |
| numpy        | ≥ 1.23  | BSD-3-Clause   | Permitted      | `pip install numpy`            |
| pandas       | ≥ 2.0   | BSD-3-Clause   | Permitted      | `pip install pandas`           |
| matplotlib   | ≥ 3.7   | PSF            | Permitted      | `pip install matplotlib`       |
| seaborn      | ≥ 0.13  | BSD-3-Clause   | Permitted      | `pip install seaborn`          |
| scrublet     | ≥ 0.2.3 | MIT            | Permitted      | `pip install scrublet`         |
| scvi-tools   | ≥ 1.1   | BSD-3-Clause   | Permitted      | `pip install scvi-tools`       |
| harmonypy    | ≥ 0.0.9 | GPL-3          | Permitted      | `pip install harmonypy`        |
| celltypist   | ≥ 1.6   | MIT            | Permitted      | `pip install celltypist`       |
| pydeseq2     | ≥ 0.4   | MIT            | Permitted      | `pip install pydeseq2`         |
| leidenalg    | ≥ 0.10  | GPL-3          | Permitted      | `pip install leidenalg`        |
| igraph       | ≥ 0.11  | GPL-2          | Permitted      | `pip install python-igraph`    |

Minimum: `pip install "scanpy[leiden]" anndata scrublet`. Add `scvi-tools harmonypy celltypist pydeseq2` for the full pipeline.

## Inputs

**Required (one of):**

- 10X CellRanger output directory containing `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz` (filtered or raw).
- 10X `.h5` file (filtered or raw).
- Pre-built `.h5ad` AnnData file.
- Count matrix as CSV/TSV (cells × genes or genes × cells; declare which).

**Optional:**

- Sample metadata table (CSV/TSV) keyed on cell barcode or sample ID — columns like `sample_id`, `condition`, `batch`, `donor`, `tissue`.
- Custom marker gene list for annotation (CSV with `gene` and `cell_type` columns).
- Reference atlas for transfer learning (CellTypist `.pkl` model or scVI-pretrained model).

**Data requirements:** Human or mouse data work out of the box (`species="human"|"mouse"` selects `MT-`/`mt-` and `RP[SL]`/`Rp[sl]` patterns). Other species require explicit gene-symbol patterns. At least ~500 cells per sample is recommended; below that, clustering granularity collapses.

## Outputs

**AnnData artifact (primary):**

- `adata_processed.h5ad` — full annotated object with:
  - `adata.X` = log1p-normalized counts (used for visualization and DE).
  - `adata.layers["counts"]` = raw integer UMI counts (used for pseudobulk DE and re-normalization).
  - `adata.obs` with QC metrics (`n_genes_by_counts`, `total_counts`, `pct_counts_mt`, `pct_counts_ribo`, `doublet_score`), clusters at multiple resolutions, and cell type labels.
  - `adata.obsm["X_pca"]`, `adata.obsm["X_umap"]`, optionally `X_scVI` / `X_scANVI` / `X_pca_harmony`.
  - `adata.uns["rank_genes_groups"]` cluster marker tables.

Load with: `import scanpy as sc; adata = sc.read_h5ad("adata_processed.h5ad")`.

**CSV exports (for downstream consumers without AnnData):**

- `cell_metadata.csv` — per-cell `obs` table.
- `expression_matrix_counts.csv` / `expression_matrix_normalized.csv` — full or HVG-only.
- `pca_coordinates.csv`, `umap_coordinates.csv` — embeddings.
- `cluster_markers_all.csv` — top markers per cluster (Wilcoxon).
- `{celltype}_de_results.csv` — pseudobulk DE per cell type when `condition_key` is provided.

**Figures (PNG + SVG at 300 DPI):**

- QC: violin, scatter (counts vs genes colored by MT%), per-batch distributions, doublet score histograms.
- Embedding: PCA scree, UMAP by cluster, UMAP by batch, UMAP by feature (marker genes, QC metrics).
- Marker: clustermap heatmap, dot plot, ranked-gene plots.
- DE: per-cell-type volcano and MA plots.

**PDF report style rules** (when emitting a report alongside outputs):

- US Letter (8.5 × 11 in); set page size explicitly, don't rely on library defaults.
- No Unicode superscripts in numeric output — write `3.4e-06` or `3.4 × 10^(-6)`, never `10⁻⁶` (renders as black squares in most PDF fonts).
- Multi-panel figures must be ≥ 80% of page width — readers can't parse 50%-width UMAPs.
- Page-break only before major sections; keep headings glued to the first paragraph below them.

## Clarification Questions

**⚠️ CRITICAL: ASK THIS FIRST**

**Q1. Where is the data, and is it filtered or raw?**

- Path to CellRanger directory, `.h5`, or `.h5ad`?
- Filtered (CellRanger already called cells) or raw (you need to run cell-calling and ambient correction)?
- Or run on the bundled PBMC 3k demo to validate the pipeline first?

Without this answer, the first step (load + ambient correction decision) is undefined.

**Q2. Species and tissue.**

- Human or mouse? (Sets MT/ribo gene patterns.)
- Tissue? (Drives QC thresholds — PBMC tolerates 5% MT, brain 10%, tumor up to 20%, cardiomyocytes ~ 20% by biology.)

**Q3. Single batch or multi-batch?**

- Single sample → skip integration, fixed QC thresholds OK.
- Multi-batch → MAD-based QC per batch and choose integration: scVI (best, needs GPU), Harmony (fast CPU), scANVI (when you already have partial labels).

**Q4. Are you comparing conditions?**

- If yes: how many biological replicates per condition? Pseudobulk DESeq2 needs N ≥ 2 per group (and ≥ 3 to be honest about variance).
- If no: stop after annotation; pseudobulk step is skipped.

**Q5. Clustering granularity expectation.**

- Major lineages only (T/B/Mono/NK): resolution 0.4–0.5.
- Standard (subtypes like CD4/CD8 and classical/non-classical mono): 0.6–0.8.
- Fine (T-cell substates, transitional populations): 1.0–1.5.
- Default: test 0.4 / 0.6 / 0.8 / 1.0 and pick by marker coherence + stability.

## Standard Workflow

The scripts under `scripts/` are the canonical entry points. **Prefer the script functions over inline scanpy code** — they encode opinionated defaults, verification prints, and the failure modes we've already debugged. Re-write inline only as a last resort (see "When Scripts Fail" below).

### Step 1 — Load and run QC

Scripts: `setup_and_import.py`, `load_example_data.py`, `qc_metrics.py`, `plot_qc.py`, `filter_cells.py`.

```python
from scripts.setup_and_import import setup_scanpy_environment, import_10x_data
from scripts.qc_metrics import calculate_qc_metrics, batch_mad_outlier_detection
from scripts.filter_cells import run_scrublet_detection, filter_by_mad_outliers

setup_scanpy_environment()
adata = import_10x_data("filtered_feature_bc_matrix/")
adata = calculate_qc_metrics(adata, species="human")
adata = batch_mad_outlier_detection(adata, batch_key="batch", nmads=5)
adata = run_scrublet_detection(adata, batch_key="batch")
adata = filter_by_mad_outliers(adata, remove_doublets=True)
```

For raw matrices in high-soup tissues (brain, tumor, lung), prepend CellBender — see `references/ambient_rna_correction.md`.

Target cell retention ≥ 70%. If lower, relax `nmads` from 5 to 7 or fall back to fixed tissue thresholds with `mark_outliers_fixed(tissue="brain")`.

### Step 2 — Normalize, find HVGs, scale, PCA

Scripts: `normalize_data.py`, `find_variable_genes.py`, `scale_and_pca.py`.

```python
from scripts.normalize_data import run_standard_normalization
from scripts.find_variable_genes import find_highly_variable_genes
from scripts.scale_and_pca import scale_data, run_pca_analysis, suggest_n_pcs

adata = run_standard_normalization(adata, target_sum=1e4)
adata = find_highly_variable_genes(adata, n_top_genes=2000, batch_key="batch")
adata = scale_data(adata, vars_to_regress=["total_counts", "pct_counts_mt"])
adata = run_pca_analysis(adata, n_pcs=50)
n_pcs = suggest_n_pcs(adata)
```

The standard normalization (target-sum + log1p) is fine for ~95% of datasets. Use `run_pearson_residuals()` if the variance–mean relationship is wildly heteroscedastic (see `references/scanpy_best_practices.md`).

**Never use < 15 PCs.** Heuristics like "last PC where variance > 5% of PC1" collapse subpopulations. Default to 30; trust `suggest_n_pcs()` over your eye on the scree plot.

### Step 3 — Integrate (multi-batch only)

Scripts: `integrate_scvi.py`, `integration_diagnostics.py`.

```python
from scripts.integrate_scvi import run_scvi_integration
from scripts.integration_diagnostics import compute_lisi_scores

adata = run_scvi_integration(
    adata,
    batch_key="batch",
    condition_key="condition",
    n_latent=30,
    max_epochs=400,
)
lisi = compute_lisi_scores(adata, batch_key="batch", use_rep="X_scVI")
```

`run_scvi_integration` checks for batch–condition confounding (N=1 in any condition) and emits a warning if found — clustering still valid, composition comparisons not. For CPU-only or huge datasets, swap in `run_harmony_integration`.

### Step 4 — Neighbor graph, clustering, UMAP, markers

Scripts: `cluster_cells.py`, `run_umap.py`, `find_markers.py`, `plot_dimreduction.py`.

```python
from scripts.cluster_cells import build_neighbor_graph, cluster_leiden_multiple_resolutions
from scripts.run_umap import run_umap_reduction
from scripts.find_markers import find_all_cluster_markers
from scripts.plot_dimreduction import plot_umap_clusters

use_rep = "X_scVI" if "X_scVI" in adata.obsm else "X_pca"
adata = build_neighbor_graph(adata, use_rep=use_rep, n_neighbors=10, n_pcs=n_pcs)
adata = cluster_leiden_multiple_resolutions(adata, resolutions=[0.4, 0.6, 0.8, 1.0])
adata = run_umap_reduction(adata, random_state=42)
markers = find_all_cluster_markers(adata, cluster_key="leiden_0.8")
plot_umap_clusters(adata, cluster_key="leiden_0.8", output_dir="results/umap")
```

`find_all_cluster_markers` uses Wilcoxon by default. This is **exploratory** — for condition contrasts go to Step 6.

### Step 5 — Annotate cell types

Script: `annotate_celltypes.py`. References: `references/marker_gene_database.md`.

```python
from scripts.annotate_celltypes import annotate_clusters_manual, annotate_with_celltypist

annotations = {
    "0": "CD4 T cells",
    "1": "CD14+ Monocytes",
    "2": "CD8 T cells",
    "3": "B cells",
    "4": "NK cells",
}
adata = annotate_clusters_manual(adata, annotations, cluster_key="leiden_0.8")

# Optional: validate with CellTypist
adata = annotate_with_celltypist(adata, model="Immune_All_Low.pkl", majority_voting=True)
```

The hybrid pattern (CellTypist first pass → manual refinement against your markers) is the most reliable; pure automated annotation has known failure modes (ILC misassigned as NK, HSC over-called, low-complexity clusters dumped into "stem-like").

### Step 6 — Pseudobulk DE (multi-sample, multi-condition)

Script: `pseudobulk_de.py`. Reference: `references/pseudobulk_de_guide.md`.

```python
from scripts.pseudobulk_de import aggregate_to_pseudobulk, run_deseq2_analysis, plot_volcano

pseudobulk = aggregate_to_pseudobulk(
    adata,
    sample_key="sample_id",
    celltype_key="cell_type",
    min_cells=10,
)
de_results = run_deseq2_analysis(
    pseudobulk,
    formula="~ batch + condition",
    contrast=["condition", "treated", "control"],
)
for celltype, res in de_results.items():
    plot_volcano(res, celltype=celltype, output_dir="results/de")
```

The function refuses to run if N=1 in any condition — that's correct, not a bug. Run cell-level Wilcoxon as a candidate-gene exploration in that case, and explicitly label the result as exploratory.

### Step 7 — Export

Script: `export_results.py`.

```python
from scripts.export_results import export_anndata_results
export_anndata_results(adata, output_dir="results", cluster_key="cell_type")
```

Produces `adata_processed.h5ad`, CSV tables, and a `analysis_summary.txt` with dataset shape, retention, integration scores, and cluster counts. The export is idempotent — safe to re-run if it's interrupted.

## When Scripts Fail

Use this hierarchy, in order. Document the level you used.

1. **Fix and retry (90% of failures):** install the missing package, re-run.
2. **Modify the script in place (5%):** edit the script file, comment the change, re-run.
3. **Use the script as reference (4%):** read it, write a one-off inline cell that follows the same algorithm, cite the source.
4. **Write from scratch (1%):** only when the algorithm itself is wrong for your data. Explain why in the report.

Never skip past (1) and write inline scanpy without trying the script.

## Decision Points

| Decision                         | Default                                   | Override when                                                                |
| -------------------------------- | ----------------------------------------- | ---------------------------------------------------------------------------- |
| Ambient RNA correction           | Skip (filtered PBMC)                      | Raw matrix, or high-soup tissue (brain, lung, tumor) → run CellBender        |
| QC thresholding                  | MAD with `nmads=5`                        | Single batch with strong tissue priors → `mark_outliers_fixed(tissue=...)`   |
| Doublet rate (Scrublet)          | Auto-scaled: ~0.8% per 1000 loaded cells  | Plate-based (Smart-seq2) → skip; very high cell loading → 10–15%             |
| Normalization                    | Target-sum 1e4 + log1p                    | Strongly heteroscedastic data → Pearson residuals                            |
| `n_top_genes` (HVG)              | 2000                                      | Small atlases or rare populations → 3000; tiny datasets → 1000               |
| `n_pcs`                          | 30                                        | Use `suggest_n_pcs()` output; never < 15                                     |
| Integration method               | scVI (`n_latent=30`)                      | CPU-only or > 500k cells → Harmony; partial labels available → scANVI        |
| Leiden resolution                | Test 0.4, 0.6, 0.8, 1.0                   | Pick by marker coherence and stability across resolutions                    |
| Annotation                       | Manual based on Step 4 markers            | Quick-look on standard tissues → CellTypist; treat as validation, not truth  |
| Cell-level DE for inference      | **Never**                                 | Always use pseudobulk + DESeq2 for condition tests                           |

## Common Issues

| Issue                                                              | Likely cause                                                       | Fix                                                                                          |
| ------------------------------------------------------------------ | ------------------------------------------------------------------ | -------------------------------------------------------------------------------------------- |
| `ImportError: No module named 'scanpy'`                            | Missing install                                                    | `pip install "scanpy[leiden]" anndata scrublet`                                              |
| `FileNotFoundError: barcodes.tsv.gz`                               | Pointed at the wrong 10X directory                                 | `ls` the dir; it must hold `barcodes`, `features`, `matrix` files (or use `.h5`)             |
| Cell retention < 70% after QC                                      | MAD too strict, or genuinely degraded sample                       | Raise `nmads` to 7; visualize `pct_counts_mt` distribution first                             |
| Mitochondrial genes detected = 0                                   | Wrong species or non-standard gene symbols                         | Verify `MT-` (human) vs `mt-` (mouse); pass `species` arg correctly                          |
| Clusters driven by `pct_counts_mt` / `total_counts`                | Insufficient QC or no regression                                   | Tighten MT filter; pass `vars_to_regress` to `scale_data`                                    |
| Batches still separate in UMAP                                     | Integration under-corrected or skipped                             | Switch from PCA to `X_scVI`; increase `n_latent` to 40 or `theta` (Harmony) to 3             |
| Cell types over-merge after integration                            | Over-correction                                                    | Lower `n_latent` to 20 or `theta` to 1; check `cLISI`                                        |
| `MemoryError` during `scale_data`                                  | Dense matrix on > 50k cells                                        | Subset to HVGs first: `adata = adata[:, adata.var.highly_variable].copy()`                   |
| UMAP changes every run                                             | Stochastic seed                                                    | Always pass `random_state=42` to `run_umap_reduction` and to Leiden                          |
| Top markers are all ribosomal / MT genes                           | Forgot to exclude before ranking                                   | Drop with `~adata.var_names.str.startswith(("MT-", "RPS", "RPL"))` before `find_markers`     |
| Pseudobulk DE "refuses to run"                                     | N=1 in a condition (script guard fires)                            | Don't fight it — report N=1 as a design limitation; use cell-level Wilcoxon as exploratory   |
| `H5AD` write hangs at the end                                      | gzip compression on large object                                   | `save_h5ad(adata, path, compression="lzf")` (default) or `compression=None`                  |
| `CUDA out of memory` during scVI                                   | GPU too small                                                      | Reduce `n_hidden` (128 → 64), drop `n_layers` (2 → 1), or `use_gpu=False`                    |
| CellTypist mislabels ILC as NK or HSC as T                         | Known model failure modes                                          | Always cross-check predicted labels against manual markers; relabel ILC→NK, HSC→Unassigned   |

For deeper triage, see `references/troubleshooting_guide.md`.

## Best Practices

- **Always store raw counts** in `adata.layers["counts"]` before normalizing. You need them later for pseudobulk DE and re-runs.
- **Set the random seed everywhere.** UMAP, Leiden, Scrublet, scVI — all stochastic. `random_state=42` is the convention.
- **Visualize QC distributions before filtering.** A unimodal MT% peak at 12% is biology in cardiomyocytes, not a filtering signal.
- **Document the integration method, parameters, and quality scores** in `adata.uns["integration"]` and in the report. Future-you will thank you.
- **Don't use integrated counts for DE.** Pseudobulk on raw counts (in `adata.layers["counts"]`), include batch as a covariate in the formula.
- **Re-run resolution sweeps after every QC change.** Cluster identity is downstream of QC; you can't reuse old labels after re-filtering.

## Suggested Next Steps

1. **`scrna-trajectory-inference`** — pass the annotated `.h5ad` to PAGA / Palantir / scVelo when your tissue has developmental structure (hematopoiesis, neurogenesis, fibroblast activation, etc.).
2. **`cell-cell-communication`** — once cell types are annotated, run CellChat / LIANA / NicheNet on the same object to score ligand–receptor activity between populations.
3. **`functional-enrichment-from-degs`** — take the pseudobulk DE gene lists and run GSEA / Enrichr / GO over-representation to interpret programs.
4. **`spatial-transcriptomics`** — if you also have Visium / Xenium / Slide-seq from the same tissue, use this object as the dissociated reference for cell-type deconvolution.
5. **`multi-omics-integration`** — combine with paired CITE-seq, scATAC-seq, or bulk profiles via MOFA / WNN / Harmony.

## Related Skills

- **R-based alternative:** `scrnaseq-seurat-core-analysis`.
- **Downstream interpretation:** `functional-enrichment-from-degs`, `coexpression-network`, `upstream-regulator-analysis`.
- **Adjacent omics:** `bulk-rnaseq-differential-expression`, `bulk-omics-clustering`.

## References

Files under `references/`:

- `qc_guidelines.md` — tissue-specific QC thresholds and adaptive filtering recipes.
- `scanpy_best_practices.md` — normalization, HVG, PCA, clustering, UMAP defaults.
- `integration_methods.md` — scVI vs scANVI vs Harmony vs Seurat CCA/RPCA, with LISI/ASW interpretation.
- `marker_gene_database.md` — curated cell-type markers across PBMC, brain, lung, gut, stromal, endothelial.
- `pseudobulk_de_guide.md` — pseudoreplication problem, DESeq2 design formulas, minimum sample size.
- `ambient_rna_correction.md` — CellBender / SoupX, tissue-specific contamination expectations.
- `troubleshooting_guide.md` — error messages, root causes, fixes.
- `workflow-details.md` — fully expanded code per phase.
- `common-patterns.md` — copy-paste recipes for PBMC, multi-batch, ambient correction, pseudobulk, automated annotation.

- **[Scanpy]**: Wolf FA, Angerer P, Theis FJ. "SCANPY: large-scale single-cell gene expression data analysis." *Genome Biology* 2018;19:15. doi:10.1186/s13059-017-1382-0 — [link](https://doi.org/10.1186/s13059-017-1382-0)
- **[scvi-tools]**: Gayoso A, Lopez R, Xing G, et al. "A Python library for probabilistic analysis of single-cell omics data." *Nature Methods* 2022;19:1457–1461. doi:10.1038/s41592-021-01343-9 — [link](https://doi.org/10.1038/s41592-021-01343-9)
- **[Harmony]**: Korsunsky I, Millard N, Fan J, et al. "Fast, sensitive and accurate integration of single-cell data with Harmony." *Nature Methods* 2019;16:1289–1296. doi:10.1038/s41592-019-0619-0 — [link](https://doi.org/10.1038/s41592-019-0619-0)
- **[Best practices]**: Heumos L, Schaar AC, Lance C, et al. "Best practices for single-cell analysis across modalities." *Nature Reviews Genetics* 2023;24:550–572. doi:10.1038/s41576-023-00586-w — [link](https://doi.org/10.1038/s41576-023-00586-w)
- **[Pseudobulk DE]**: Squair JW, Gautier M, Kathe C, et al. "Confronting false discoveries in single-cell differential expression." *Nature Communications* 2021;12:5692. doi:10.1038/s41467-021-22786-0 — [link](https://doi.org/10.1038/s41467-021-22786-0)
- **[CellTypist]**: Domínguez Conde C, Xu C, Jarvis LB, et al. "Cross-tissue immune cell analysis reveals tissue-specific features in humans." *Science* 2022;376:eabl5197. doi:10.1126/science.abl5197 — [link](https://doi.org/10.1126/science.abl5197)
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)