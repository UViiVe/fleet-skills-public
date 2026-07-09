---
name: scrna-trajectory-inference
id: "skill_12025f79ff54476aa6c7005a7bf33848"
display-name: "scRNA-seq Trajectory Inference"
short-description: "Reconstruct differentiation trajectories from scRNA-seq using PAGA, diffusion pseudotime, scVelo RNA velocity, and CellRank fate probabilities. Use when ordering cells along a developmental axis, i..."
category: transcriptomics
keywords: [trajectory, pseudotime, scRNA, monocle, velocity]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Reconstruct differentiation trajectories from scRNA-seq using PAGA, diffusion pseudotime, scVelo RNA velocity, and CellRank fate probabilities. Use when ordering cells along a developmental axis, identifying branch points, or computing cell-fate probabilities."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "PAGA, scVelo, CellRank"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Wolf FA et al., Genome Biology 2019, doi:10.1186/s13059-019-1663-x"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Single-Cell Trajectory Inference

Reconstruct developmental and disease trajectories from preprocessed scRNA-seq data. The skill is layered: PAGA + diffusion pseudotime always run on `scanpy` alone; RNA velocity (`scVelo`) and fate probabilities (`CellRank`) layer on when their dependencies and data layers are available.

Designed to chain directly off `scrnaseq-scanpy-core-analysis`: hand in the `.h5ad`, get back trajectory annotations, plots, and gene lists ready for enrichment.

## When to Use This Skill

Use this skill when:

- You have preprocessed scRNA-seq data (PCA, UMAP, clusters in `obs`) and want to order cells along a differentiation or disease axis.
- You want to identify branching points and terminal fates from the topology of the cluster graph.
- You have spliced/unspliced counts (STARsolo, velocyto, kallisto|bustools) and want to add RNA velocity directionality.
- You want CellRank fate probabilities and driver genes per terminal state.
- You're investigating a system with known developmental structure: hematopoiesis, neurogenesis, β-cell differentiation, fibroblast activation, EMT, organoid differentiation, T-cell exhaustion progression, etc.

Do NOT use this skill for:

- **Un-preprocessed data.** Run `scrnaseq-scanpy-core-analysis` first to get a clean `.h5ad` with clusters and a UMAP.
- **Terminally differentiated systems with no real trajectory.** Resting PBMCs from one timepoint don't have one; you'll get a "trajectory" of noise.
- **Bulk RNA-seq over time.** Use `disease-progression-longitudinal` for TimeAx-style pseudotemporal ordering of bulk samples.
- **Very small datasets** (< 200 cells). The neighbor graph is too noisy and pseudotime collapses.
- **R-only environments.** This skill is Python-native. For Seurat, use Monocle3 or Slingshot via the Bioconductor stack (no Seurat equivalent in this skillpack yet).

## Quick Start

```python
from scripts.load_example_data import load_example_data
from scripts.run_trajectory_analysis import run_trajectory
from scripts.generate_all_plots import generate_all_plots
from scripts.export_results import export_all

adata = load_example_data()  # pancreatic endocrinogenesis demo
results = run_trajectory(adata, root_cell_type="Ductal", cluster_key="clusters")
generate_all_plots(adata, results, output_dir="trajectory_results", cluster_key="clusters")
export_all(adata, results, output_dir="trajectory_results")
```

Runtime: ~3 minutes (core trajectory), ~5 minutes (+ scVelo), ~8 minutes (+ CellRank) on the demo dataset.

## Installation

| Package      | Version | License      | Commercial use | Install command            |
| ------------ | ------- | ------------ | -------------- | -------------------------- |
| scanpy       | ≥ 1.10  | BSD-3-Clause | Permitted      | `pip install scanpy`       |
| anndata      | ≥ 0.10  | BSD-3-Clause | Permitted      | `pip install anndata`      |
| numpy        | ≥ 1.23  | BSD-3-Clause | Permitted      | `pip install numpy`        |
| pandas       | ≥ 2.0   | BSD-3-Clause | Permitted      | `pip install pandas`       |
| scipy        | ≥ 1.10  | BSD-3-Clause | Permitted      | `pip install scipy`        |
| matplotlib   | ≥ 3.7   | PSF          | Permitted      | `pip install matplotlib`   |
| seaborn      | ≥ 0.13  | BSD-3-Clause | Permitted      | `pip install seaborn`      |
| scvelo       | ≥ 0.3   | BSD-3-Clause | Permitted      | `pip install scvelo`       |
| cellrank     | ≥ 2.0   | BSD-3-Clause | Permitted      | `pip install cellrank`     |
| statsmodels  | ≥ 0.14  | BSD-3-Clause | Permitted      | `pip install statsmodels`  |
| reportlab    | ≥ 4.0   | BSD          | Permitted      | `pip install reportlab`    |

`pip install scanpy anndata scipy matplotlib seaborn statsmodels` is the minimum (PAGA + DPT). Add `scvelo cellrank reportlab` for the full pipeline. Scripts detect missing optional packages and downgrade gracefully — they don't crash.

## Inputs

**Required:**

- Preprocessed `.h5ad` with:
  - `adata.X` log-normalized counts.
  - `adata.obsm["X_pca"]` and `adata.obsm["X_umap"]`.
  - A cluster or cell-type column in `adata.obs` (e.g. `"clusters"`, `"cell_type"`, `"leiden"`).
  - ≥ 200 cells, ≥ 100 genes.

The output of `scrnaseq-scanpy-core-analysis` (`adata_processed.h5ad`) works directly. Seurat objects need conversion via `SeuratDisk::SaveH5Seurat` + `Convert(dest = "h5ad")`.

**For RNA velocity (optional):**

- `adata.layers["spliced"]` and `adata.layers["unspliced"]` count matrices.
- Produced by `velocyto.py`, STARsolo with `--soloFeatures Gene Velocyto`, or kallisto|bustools `lamanno` workflow.

## Outputs

**AnnData artifact (primary):**

- `adata_trajectory.h5ad` — annotated with:
  - `adata.obs["dpt_pseudotime"]` — diffusion pseudotime per cell.
  - `adata.obsm["X_diffmap"]` — diffusion map coordinates.
  - `adata.uns["paga"]` — PAGA connectivities and positions.
  - `adata.uns["iroot"]` — root cell index.
  - `adata.layers["velocity"]` / `adata.obsm["velocity_umap"]` (if scVelo ran).
  - `adata.obsm["lineages_fwd"]` (if CellRank ran).

Load with: `import scanpy as sc; adata = sc.read_h5ad("adata_trajectory.h5ad")`.

**Tables (CSV):**

- `pseudotime_assignments.csv` — barcode, pseudotime, cell type.
- `trajectory_genes.csv` — genes correlated with pseudotime (correlation, FDR, direction).
- `velocity_genes.csv` — top RNA velocity genes (if scVelo ran).
- `fate_probabilities.csv` — per-cell fate probabilities (if CellRank ran).
- `driver_genes_<fate>.csv` — driver genes per terminal fate (if CellRank ran).

**Other artifacts:**

- `trajectory_results.pkl` — full results dict (pickle).
- `analysis_metadata.json` — parameters and quality metrics.

**Figures (PNG + SVG at 300 DPI):**

- `paga_graph.png` — PAGA cluster graph on UMAP.
- `pseudotime_umap.png` — UMAP coloured by pseudotime.
- `pseudotime_violin.png` — pseudotime distribution per cell type.
- `diffusion_components.png` — first diffusion components.
- `gene_heatmap.png` — top trajectory genes heatmap.
- `gene_trends.png` — expression trends along pseudotime.
- `paga_connectivity.png` — connectivity heatmap.
- `velocity_stream.png`, `velocity_confidence.png`, `latent_time.png`, `velocity_top_genes.png` — if scVelo ran.
- `fate_probabilities.png`, `fate_heatmap.png`, `driver_genes.png` — if CellRank ran.

**Reports:**

- `trajectory_analysis_report.pdf` — full report (requires `reportlab`).
- `trajectory_analysis_report.md` — markdown fallback.

PDF style rules: US Letter, no Unicode superscripts (use `3.4e-06`), figures ≥ 80% page width, no page breaks splitting a heading from its content.

## Clarification Questions

**⚠️ CRITICAL: ASK THIS FIRST**

**Q1. Where is the preprocessed `.h5ad`, or run the demo?**

- Path to `adata.h5ad` (output of `scrnaseq-scanpy-core-analysis` or equivalent)?
- Or use the bundled demo: pancreatic endocrinogenesis (3,696 cells, ductal → α/β/δ/ε branching)?

If demo, parameters are pre-set — only Q2 applies, then start the workflow.

**Q2. Analysis depth.**

- **(a) Core only:** PAGA + diffusion pseudotime. Always works. ~ 3 min.
- **(b) Core + scVelo:** add RNA velocity if `spliced`/`unspliced` layers exist. ~ 5 min.
- **(c) Full:** add CellRank fate probabilities and driver genes. ~ 8 min.

**Q3. Root cell type.** (User-data only.)

Which cell type is the start of the trajectory? In the demo it's "Ductal" (progenitors). Choose the cluster you have biological reason to believe is upstream.

**Q4. Cluster key.** (User-data only.)

Which column in `adata.obs` holds the cell-type labels? Common: `clusters`, `cell_type`, `leiden`, `louvain`.

**Q5. Expected topology.** (User-data only.)

- **(a) Linear** — one lineage, one fate.
- **(b) Branching** — multiple fates from one progenitor (most common).
- **(c) Not sure** — let the data decide; PAGA will reveal the structure.

## Standard Workflow

Scripts under `scripts/` are the canonical entry points. They encode opinionated defaults, sensible verbose output, and the graceful-degradation logic. Don't open-code inline equivalents.

### Step 1 — Load data

```python
from scripts.load_example_data import load_example_data
adata = load_example_data()    # demo
# or for your own data:
import scanpy as sc
adata = sc.read_h5ad("path/to/adata_processed.h5ad")
```

Expect: `"✓ Data loaded successfully!"`

### Step 2 — Run trajectory

```python
from scripts.run_trajectory_analysis import run_trajectory
results = run_trajectory(
    adata,
    root_cell_type="Ductal",     # or your progenitor cluster
    cluster_key="clusters",
    run_velocity=True,           # auto-skips if scVelo or layers missing
    run_cellrank=True,           # auto-skips if cellrank missing
)
```

Pipeline:

1. **Neighbors + PAGA** — graph abstraction at cluster level.
2. **UMAP re-init from PAGA** — trajectory-aware layout.
3. **Diffusion map + DPT** — diffusion pseudotime from the root cell.
4. **Trajectory genes** — Spearman correlation with pseudotime, BH-corrected.
5. **scVelo** (if available): preprocess velocity, dynamical model with stochastic fallback, latent time.
6. **CellRank** (if available): kernel from velocity (or PAGA pseudotime), terminal states, fate probabilities, driver genes per fate.

Expect: `"✓ Trajectory analysis completed successfully!"`

### Step 3 — Plots

```python
from scripts.generate_all_plots import generate_all_plots
generate_all_plots(
    adata, results,
    output_dir="trajectory_results",
    cluster_key="clusters",
)
```

Generates PAGA, pseudotime UMAP, diffusion components, gene heatmap, gene trends, plus velocity / fate / driver plots when the relevant results are available.

Expect: `"✓ All plots generated successfully!"`

### Step 4 — Export

```python
from scripts.export_results import export_all
export_all(adata, results, output_dir="trajectory_results")
```

Writes `adata_trajectory.h5ad`, all CSVs, the pickle of the results dict, and `analysis_metadata.json`. Export is idempotent — safe to re-run.

Expect: `"=== Export Complete ==="`

### Step 5 — Report (optional)

```python
from scripts.generate_report import generate_report
generate_report(adata, results, output_dir="trajectory_results")
```

PDF if `reportlab` is installed, markdown fallback otherwise.

## When Scripts Fail

In order, document the level used:

1. **Fix and retry (90%):** install the missing package (`pip install scvelo` or `cellrank`).
2. **Modify the script in place (5%):** edit the file, comment the change, re-run.
3. **Use the script as reference (4%):** read it, adapt approach inline, cite source.
4. **Write from scratch (1%):** only if the algorithm itself doesn't fit. Justify in the report.

Never skip past (1) and write inline PAGA / DPT / scVelo code first.

## Decision Points

| Decision                          | Default                                  | Override when                                                              |
| --------------------------------- | ---------------------------------------- | -------------------------------------------------------------------------- |
| `n_neighbors` for trajectory      | 30                                       | Small dataset (< 1000 cells) → 15; very large (> 50k) → 50                  |
| `n_pcs`                           | 30                                       | Use the same value chosen during preprocessing (`suggest_n_pcs`)            |
| `n_dcs` (diffusion components)    | 10                                       | Complex branching → 15; very simple linear → 5                              |
| Root cell selection               | First cell of `root_cell_type` cluster   | Override with a specific barcode if you have biological knowledge          |
| Trajectory gene correlation       | Spearman                                 | Pearson when the relationship is known to be linear                        |
| scVelo model                      | Dynamical (with stochastic fallback)     | Force `stochastic` for very low spliced/unspliced ratios                    |
| CellRank kernel                   | Velocity kernel when scVelo ran, else PAGA pseudotime | Use connectivity kernel for very high-dimensional embeddings    |
| Number of terminal states         | Auto-detected via Schur decomposition    | Manual override when biology dictates (e.g., known 3 mature cell types)    |
| Pseudotime smoothing for gene trends | None (raw GAM fit)                    | Increase `n_bins` to 50–100 for smoother trends in noisy datasets           |

## Common Issues

| Issue                                       | Cause                                                           | Fix                                                                                |
| ------------------------------------------- | --------------------------------------------------------------- | ---------------------------------------------------------------------------------- |
| `"Root cell type not found"`                | Cluster name doesn't match `cluster_key`                        | `adata.obs["clusters"].unique()` to find the exact name                            |
| `ModuleNotFoundError: No module named 'scvelo'` | Optional dependency missing                                  | `pip install scvelo`, or accept that velocity will be skipped                      |
| `ModuleNotFoundError: No module named 'cellrank'` | Optional dependency missing                                | `pip install cellrank`, or accept that fate mapping will be skipped                |
| scVelo dynamical model fails to converge    | Insufficient spliced/unspliced counts                           | Script auto-falls back to stochastic model; check fallback was used                |
| Pseudotime has infinity values              | Disconnected components in the neighbor graph                   | Raise `n_neighbors` to 50; verify no batch effects splitting the manifold          |
| PAGA graph shows no connectivity            | Resolution too high — clusters too small and isolated           | Re-cluster at a coarser resolution before running trajectory                       |
| Wrong terminal states from CellRank         | Auto-detection misled by transient states                       | Manually set terminal states via `g.set_terminal_states({...})`                    |
| `"Too few cells"`                           | < 200 cells in the object                                       | Trajectory inference is unreliable below this; expand the dataset                  |
| Diffusion map is degenerate (all-zero)      | All cells too similar — no biological signal in low-dim space   | Recompute PCA with more variable genes; raise `n_dcs`                              |
| Trajectory genes table empty                | Pseudotime range too narrow                                     | Confirm the root cell is upstream; check `dpt_pseudotime` distribution             |
| `IndexError: invalid root cell index`        | Root cell type column has unexpected categories                 | Cast `adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")`           |
| SVG export fails                            | Missing system library                                          | Normal — PNG always works. The skill writes both when possible                     |
| Gene trends look spiky / non-monotonic      | Few cells per pseudotime bin                                    | Increase `n_bins` or restrict to one branch at a time                              |

For deeper triage, see `references/troubleshooting.md`.

## Best Practices

- **Always justify the root cell biologically.** A wrong root inverts pseudotime and your gene lists go the opposite direction.
- **Inspect the PAGA graph before interpreting pseudotime.** It reveals the topology (linear vs branching) and identifies the cell-type connections that should and shouldn't exist.
- **Don't over-interpret pseudotime as real time.** It's an ordering, not a clock. The mapping to actual hours/days requires orthogonal labels (timepoints, lineage tracing).
- **Use scVelo only when the chemistry supports it.** 10X v3 with intronic reads works; 5'-tag scRNA-seq generally doesn't.
- **Validate CellRank terminal states against known biology.** Auto-detection picks Schur eigenvectors that may not match the expected mature cell types.
- **For multi-branch trajectories, analyse each branch separately downstream.** Trajectory genes computed across the full pseudotime can mix branch-specific signals.

## Suggested Next Steps

1. **`functional-enrichment-from-degs`** — take `trajectory_genes.csv` and run pathway enrichment to interpret the programs driving differentiation.
2. **`grn-pyscenic`** — feed `adata_trajectory.h5ad` into pySCENIC to find transcription factors that regulate fate decisions.
3. **`upstream-regulator-analysis`** — predict upstream regulators of the trajectory-gene set.
4. **`scrnaseq-scanpy-core-analysis`** marker finding — run `find_all_cluster_markers` on cells at branch points to identify fate-determining genes.

## Related Skills

- **Upstream (preprocessing):** `scrnaseq-scanpy-core-analysis`, `scrnaseq-seurat-core-analysis` (convert via SeuratDisk).
- **Downstream (interpretation):** `functional-enrichment-from-degs`, `grn-pyscenic`, `upstream-regulator-analysis`, `cell-cell-communication`.
- **Alternative trajectory methods:** `disease-progression-longitudinal` for bulk longitudinal data (TimeAx).

## References

Files under `references/`:

- `trajectory-methods-guide.md` — PAGA, DPT, Monocle3, Slingshot, comparisons and parameter choices.
- `rna-velocity-guide.md` — scVelo dynamical vs stochastic, when velocity is valid, kinetic parameter interpretation.
- `cellrank-guide.md` — fate probability theory, terminal-state detection, driver gene identification.
- `troubleshooting.md` — common pitfalls and fixes.

- **[PAGA]**: Wolf FA, Hamey FK, Plass M, et al. "PAGA: graph abstraction reconciles clustering with trajectory inference." *Genome Biology* 2019;20:59. doi:10.1186/s13059-019-1663-x — [link](https://doi.org/10.1186/s13059-019-1663-x)
- **[scVelo]**: Bergen V, Lange M, Peidli S, Wolf FA, Theis FJ. "Generalizing RNA velocity to transient cell states through dynamical modeling." *Nature Biotechnology* 2020;38:1408–1414. doi:10.1038/s41587-020-0591-3 — [link](https://doi.org/10.1038/s41587-020-0591-3)
- **[CellRank]**: Lange M, Bergen V, Klein M, et al. "CellRank for directed single-cell fate mapping." *Nature Methods* 2022;19:159–170. doi:10.1038/s41592-021-01346-6 — [link](https://doi.org/10.1038/s41592-021-01346-6)
- **[Diffusion pseudotime]**: Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ. "Diffusion pseudotime robustly reconstructs lineage branching." *Nature Methods* 2016;13:845–848. doi:10.1038/nmeth.3971 — [link](https://doi.org/10.1038/nmeth.3971)
- **[RNA velocity]**: La Manno G, Soldatov R, Zeisel A, et al. "RNA velocity of single cells." *Nature* 2018;560:494–498. doi:10.1038/s41586-018-0414-6 — [link](https://doi.org/10.1038/s41586-018-0414-6)
- **[Trajectory benchmarking]**: Saelens W, Cannoodt R, Todorov H, Saeys Y. "A comparison of single-cell trajectory inference methods." *Nature Biotechnology* 2019;37:547–554. doi:10.1038/s41587-019-0071-9 — [link](https://doi.org/10.1038/s41587-019-0071-9)
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)