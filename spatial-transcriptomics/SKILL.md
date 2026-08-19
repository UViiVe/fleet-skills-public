---
id: skill_67259cbfdaace62619de8378528a0f4c
name: spatial-transcriptomics
display-name: Spatial Transcriptomics
short-description: Run a complete 10x Visium pipeline — QC, clustering, spatially
  variable genes, neighborhood enrichment, co-occurrence — on Space Ranger output,
  H5AD, or a built-in tissue dataset. Use when you have spot-level gene expression
  with spatial coordinates and need tissue domains, SVGs, and ready-to-use figures.
category: transcriptomics
keywords: spatial transcriptomics, Visium, Squidpy, Scanpy, SVG, tissue domains
version: 1.0.1
last-updated: '2026-05-20'
license: MIT
author: 001TMF (imported by Zevo 2026-05-22)
created: 2026-05-20
source: https://github.com/001TMF/agentic-science
provenance:
  wraps: Squidpy / Scanpy (10x Visium)
  upstream_source: https://github.com/001TMF/agentic-science
  citation: Palla G et al., Nature Methods 2022, doi:10.1038/s41592-021-01358-2
changelog:
- version: 1.0.1
  date: 2026-07-07
  event: citation-provenance
  summary: Added provenance block and References section with upstream tool citations.
    Vienna public subset.
description: Run a complete 10x Visium pipeline — QC, clustering, spatially variable
  genes, neighborhood enrichment, co-occurrence — on Space Ranger output, H5AD, or
  a built-in tissue dataset. Use when you have spot-level gene expression with spatial
  coordinates and need tissue domains, SVGs, and ready-to-use figures.
---

# Spatial Transcriptomics (10x Visium)

This skill takes a 10x Visium dataset from raw spot-level counts to tissue-domain plots, spatially variable genes, and a CSV/figure bundle suitable for a manuscript supplement. It standardizes on Squidpy + Scanpy because that pair handles the Visium hex grid natively and exposes the spatial statistics (Moran's I, neighborhood enrichment, co-occurrence) without extra glue.

## When to Use This Skill

Use when:
- You have **10x Visium** spatial gene expression data — Space Ranger output, an `.h5` filtered matrix, or an `.h5ad` AnnData.
- You want to identify **spatially variable genes** (SVGs) across a tissue section.
- You want to discover **spatial tissue domains** through clustering and overlay them on the H&E image.
- You need a **neighborhood enrichment** matrix that quantifies which clusters spatially co-localize or segregate.
- You want **co-occurrence** curves to see how spatial relationships change with distance.

Do not use when:
- You have single-molecule imaging data (MERFISH/Xenium, MERSCOPE, Stereo-seq) — those need imaging-specific pipelines.
- You have Slide-seq data with bead-level resolution — squidpy can read it but the QC and SVG thresholds differ; use a Slide-seq-specific workflow.
- You only have scRNA-seq without spatial coordinates — use `scrnaseq-scanpy-core-analysis` instead.
- You want cell type **deconvolution** of each spot — that is a separate workflow (cell2location, RCTD, STdeconvolve).

## Quick Start

```python
import sys
sys.path.insert(0, "scripts")

from load_example_data import load_visium_heart
from spatial_workflow import run_spatial_analysis
from generate_all_plots import generate_all_plots
from export_results import export_all

adata = load_visium_heart()                              # cardiac demo
adata = run_spatial_analysis(adata, output_dir="visium_results")
generate_all_plots(adata, output_dir="visium_results")   # auto cardiac markers
export_all(adata, output_dir="visium_results")
```

Wall time on the bundled cardiac dataset: ~5–10 minutes total (first run downloads ~50 MB).

## Installation

| Package | Version | License | Commercial use | Install command |
|---|---|---|---|---|
| squidpy | ≥ 1.4 | BSD-3-Clause | Permitted | `pip install squidpy` |
| scanpy | ≥ 1.9 | BSD-3-Clause | Permitted | `pip install scanpy` |
| anndata | ≥ 0.8 | BSD-3-Clause | Permitted | `pip install anndata` |
| plotnine | ≥ 0.12 | MIT | Permitted | `pip install plotnine` |
| plotnine-prism | ≥ 0.2 | MIT | Permitted | `pip install plotnine-prism` |
| seaborn | ≥ 0.11 | BSD-3-Clause | Permitted | `pip install seaborn` |
| matplotlib | ≥ 3.5 | PSF | Permitted | `pip install matplotlib` |
| scikit-learn | ≥ 1.0 | BSD-3-Clause | Permitted | `pip install scikit-learn` |
| scikit-misc | ≥ 0.1 | BSD-3-Clause | Permitted | `pip install scikit-misc` |
| numpy | ≥ 1.21 | BSD-3-Clause | Permitted | `pip install numpy` |
| pandas | ≥ 1.3 | BSD-3-Clause | Permitted | `pip install pandas` |

One-shot install:

```bash
pip install squidpy scanpy anndata scikit-misc plotnine plotnine-prism \
            seaborn matplotlib numpy pandas scikit-learn
```

All packages use permissive licenses (BSD, MIT, PSF) and are clear for commercial deployment.

## Inputs

**Required**
- Visium expression data in one of:
  - A Space Ranger output directory containing `filtered_feature_bc_matrix/` and a `spatial/` folder.
  - A `filtered_feature_bc_matrix.h5` file (10x HDF5 format).
  - An `.h5ad` AnnData written by Scanpy/Squidpy.

**Optional**
- Tissue H&E image, embedded in the Space Ranger output or carried in `.uns["spatial"]`. Plots fall back to spot-only views if no image is present.
- Marker gene list — defaults to cardiac markers (MYH7, TNNI3, MYH6, NPPA, TTN, ACTN2) when the demo dataset is used; supply your own for other tissues.

**Built-in demo**
- `V1_Human_Heart` from 10x Genomics: ~4,247 spots, ~33,538 genes, includes H&E. Downloaded on first use (~50 MB) and cached by Scanpy.

## Outputs

Written to `output_dir/` (default `visium_results/`):

| File | Description |
|---|---|
| `adata_processed.h5ad` | Processed AnnData (clusters, embeddings, SVG results, spatial graph). Load with `sc.read_h5ad(...)`. |
| `spatially_variable_genes.csv` | SVGs ranked by Moran's I with FDR-corrected p-values. |
| `cluster_assignments.csv` | Spot barcode, Leiden cluster, spatial X/Y, QC summaries. |
| `neighborhood_enrichment.csv` | Cluster × cluster enrichment z-score matrix. |
| `spot_metadata.csv` | Full `adata.obs` metadata for every retained spot. |
| `analysis_summary.txt` | Human-readable run summary (parameters, retention, cluster sizes, top SVGs). |
| `qc_violins.{png,svg}` | QC metric distributions (total counts, genes, % mito). |
| `spatial_clusters.{png,svg}` | Leiden clusters overlaid on tissue. |
| `spatial_markers.{png,svg}` | Selected marker genes overlaid on tissue. |
| `umap_clusters.{png,svg}` | UMAP coloured by Leiden cluster. |
| `neighborhood_enrichment.{png,svg}` | Cluster-cluster enrichment heatmap (clustermap). |
| `co_occurrence.{png,svg}` | Co-occurrence probability vs distance per cluster. |
| `top_svgs.{png,svg}` | Bar chart of top spatially variable genes. |
| `spatial_svg_<GENE>.{png,svg}` | Spatial expression map of the top SVG. |

`adata_processed.h5ad` is the downstream-friendly artefact: any further analysis (deconvolution, multi-sample integration, pathway enrichment of SVGs) can pick up from this single file.

## Clarification Questions

### 1. Input source ⚠️ CRITICAL: ASK THIS FIRST
Do you have Visium data on disk, or should the workflow use the bundled `V1_Human_Heart` demonstration dataset?
- Your data: provide a path to an `.h5ad`, `.h5`, or Space Ranger output directory.
- Demo: parameters are pre-tuned for cardiac tissue; skip the remaining questions.

### 2. QC thresholds (your data only)
- Mitochondrial cutoff: cardiac/muscle tissue is naturally mito-rich (~30–40%); use 50%. For brain, liver, immune, or epithelium, drop to 20%.
- Minimum genes per spot: 200 is the standard floor; lower for sparse tissues.

### 3. Clustering resolution (your data only)
- 0.5 — fewer, broader domains (good first pass).
- 0.8 — default; matches H&E-visible structures in most tissues.
- 1.2 — finer subdivisions; risk of fragmenting biologically meaningful domains.

### 4. Marker genes (your data only)
Provide a small panel (5–10 genes) that you want overlaid on the tissue. If you skip this, the plotting helper falls back to its cardiac default — useful only for the demo dataset.

## Standard Workflow

The four scripts are designed to be called in sequence. Add `scripts/` to `sys.path` (or run from this skill's directory) and call the public function in each module.

1. **Load**
   ```python
   from load_example_data import load_visium_heart, load_visium_data
   adata = load_visium_heart()                        # demo
   # adata = load_visium_data("path/to/visium/")      # your data
   ```
   Verification: `✓ Data loaded successfully!` followed by spot/gene counts.

2. **Analyze**
   ```python
   from spatial_workflow import run_spatial_analysis
   adata = run_spatial_analysis(adata, output_dir="visium_results")
   ```
   Verification: `✓ Spatial analysis completed successfully!` with cluster/SVG counts.

3. **Plot**
   ```python
   from generate_all_plots import generate_all_plots
   generate_all_plots(adata, output_dir="visium_results")
   # generate_all_plots(adata, output_dir="...", marker_genes=["MYC", "TP53"])
   ```
   Verification: `✓ All visualizations generated successfully!`.

4. **Export**
   ```python
   from export_results import export_all
   export_all(adata, output_dir="visium_results")
   ```
   Verification: `=== Export Complete ===` followed by a file listing.

If you do not see those verification lines, you wrote inline scanpy/squidpy code that is bypassing the script — re-import the module and call its public function.

## When Scripts Fail

Apply the hierarchy in order — do not jump straight to writing inline analysis code.

1. **Fix and retry (~90%).** Almost always a missing package; install it and re-run.
2. **Modify the script (~5%).** Edit `scripts/*.py` directly when behavior needs to change. Add a one-line note at the top of the file describing the change.
3. **Use the script as reference (~4%).** Read the script, lift the relevant function, adapt, cite the source.
4. **Write from scratch (~1%).** Only when the existing helper genuinely cannot be coaxed into producing the needed output. Document why.

## Decision Points

| Choice | Default | Override when |
|---|---|---|
| `max_pct_mito` | 50 | Drop to 20 for non-cardiac tissues. |
| `min_genes` | 200 | Drop to 100 for sparse FFPE tissues. |
| `n_top_genes` | 2000 | Raise to 3000 for highly heterogeneous tissue. |
| `resolution` (Leiden) | 0.8 | Lower for broader anatomical regions; raise to 1.2 to split confluent domains. |
| `svgs_n_perms` | 100 | Raise to 1000 for publication-quality FDR. |
| `coord_type` (spatial graph) | `grid` (n_neighs=6) | Switch to `generic` with `radius=` or `n_neighs=` for non-Visium platforms. |
| FDR threshold for SVGs | 0.05 | Tighten to 0.01 when you want only strongly spatial genes. |

## Interpretation

**Spatially variable genes.** Moran's I close to +1 marks tissue-clustered expression (the usual signal of interest). Values near 0 indicate no spatial pattern; near −1 marks a dispersed/checkerboard pattern (rare in real tissue).

**Neighborhood enrichment.** Z-scores beyond ±2 are conventionally interpreted as significant co-localization (+) or segregation (−). Diagonal entries quantify how cohesive a cluster is on its own.

**Co-occurrence curves.** A curve above its expected baseline at a given distance means the two clusters appear together more than chance. Curves that flip sign with distance (above at short range, below at long range) are characteristic of border-zone biology.

**Validate against the H&E.** Always cross-check the cluster map with the histology. Clusters that line up with visible tissue structures are biologically meaningful; clusters that look like noise on the H&E may be technical.

## Common Issues

| Error | Cause | Resolution |
|---|---|---|
| `ModuleNotFoundError: squidpy` | Package missing | `pip install squidpy` |
| `ModuleNotFoundError: plotnine_prism` | Theme package missing | `pip install plotnine-prism` |
| `KeyError: 'spatial'` | AnnData has no `.obsm['spatial']` | Reload with `sc.read_visium()` or attach coordinates manually before running. |
| `0 SVGs found (FDR < 0.05)` | Low signal / too few permutations | Raise `svgs_n_perms=1000`, relax FDR to 0.10 for exploration, or check that normalization ran. |
| `ValueError: coord_type='grid'` | Non-Visium spatial layout | Pass `coord_type='generic'` and set `n_neighs` or `radius`. |
| Memory error during HVG selection | > 30k spots, default in-memory matrix | Subsample with `sc.pp.subsample(adata, n_obs=10000)` or filter to ROI. |
| NaN values in co-occurrence | Known squidpy issue with custom `n_splits` | Stick to the default `n_splits` — the workflow already does. |
| SVG export warns | No `svglite`/cairo backend | Expected — PNG is always written; SVG is best-effort. |
| `library_id` required by `sc.read_visium` | Multiple libraries in one folder | Pass `library_id=` explicitly to `load_visium_data()`. |
| `var_names_make_unique` warning | Duplicate gene symbols | Handled in the loader; safe to ignore. |
| Cardiac markers absent from non-cardiac tissue | Tissue mismatch with default panel | Pass `marker_genes=[…]` appropriate to your tissue. |

## Best Practices

- Run with `svgs_n_perms=100` during iteration and bump to 1000 only for the final report — the speed difference is roughly 10x.
- Always inspect QC violins **before** filtering. If `pct_counts_mt` is bimodal, choose the cutoff in the valley between modes rather than a hardcoded value.
- Store the seed when reporting results: Leiden depends on initialization. The default workflow leaves Scanpy's seed, which is reproducible across runs in the same environment.
- Save the raw counts in `.layers["counts"]` before normalization. The workflow already does this — do not overwrite the layer.
- For multi-sample studies, run this skill per sample and integrate later (Harmony, Scanorama, scVI). Do not concatenate before clustering.

## Suggested Next Steps

- **Functional enrichment on SVGs.** Feed `spatially_variable_genes.csv` into `functional-enrichment-from-degs` for GO / Reactome / KEGG enrichment.
- **Cell-type deconvolution.** Pair the SVG list with a matched scRNA-seq reference using cell2location, RCTD, or STdeconvolve to produce per-spot cell-type probabilities.
- **Spatial cell-cell communication.** Run `cell-cell-communication` (CellChat v2 supports spatial mode) on the deconvolved spots to identify spatially constrained ligand–receptor signaling.
- **Multi-sample comparison.** Repeat this skill per condition, then run differential SVG and differential domain analysis across samples.
- **Regulatory networks.** Hand the cluster-specific gene programs to `grn-pyscenic` for transcription factor inference.

## Related Skills

| Skill | Relationship |
|---|---|
| `scrnaseq-scanpy-core-analysis` | Companion non-spatial workflow; share the same Scanpy idioms. |
| `cell-cell-communication` | Downstream — ligand–receptor inference on cluster-defined regions. |
| `functional-enrichment-from-degs` | Downstream — enrichment on SVG gene lists. |
| `grn-pyscenic` | Downstream — regulatory networks on spatial cluster gene programs. |
| `coexpression-network` | Downstream — co-expression on each spatial domain. |

## References

- **Squidpy**: Palla G, et al. "Squidpy: a scalable framework for spatial omics analysis." *Nature Methods* 2022;19:171–178. doi:10.1038/s41592-021-01358-2 — [link](https://doi.org/10.1038/s41592-021-01358-2)
- **10x Visium platform**: Ståhl PL, et al. "Visualization and analysis of gene expression in tissue sections by spatial transcriptomics." *Science* 2016;353:78–82. doi:10.1126/science.aaf2403 — [link](https://doi.org/10.1126/science.aaf2403)
- **Scanpy**: Wolf FA, Angerer P, Theis FJ. "SCANPY: large-scale single-cell gene expression data analysis." *Genome Biology* 2018;19:15. doi:10.1186/s13059-017-1382-0 — [link](https://doi.org/10.1186/s13059-017-1382-0)
- **Moran's I**: Moran PAP. "Notes on continuous stochastic phenomena." *Biometrika* 1950;37:17–23. doi:10.2307/2332142
- **Upstream skill source**: https://github.com/001TMF/agentic-science (MIT license)
- `references/spatial-analysis-guide.md` — parameter tuning, coordinate types, advanced options.
