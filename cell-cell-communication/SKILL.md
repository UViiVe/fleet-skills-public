---
name: cell-cell-communication
id: "skill_8864ffedb57042a095a59200517b15fd"
display-name: "Cell-Cell Communication Analysis"
short-description: "Infer ligand-receptor signaling networks from annotated single-cell RNA-seq data and produce publication-ready communication plots. Use when you have a Seurat object with cell type labels and want ..."
category: transcriptomics
keywords: [cell communication, ligand, receptor, scRNA, CellChat]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Infer ligand-receptor signaling networks from annotated single-cell RNA-seq data and produce publication-ready communication plots. Use when you have a Seurat object with cell type labels and want to map sender/receiver/mediator roles, pathway activity, or cell-cell interaction matrices."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "CellChat v2 / LIANA+"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Jin S et al., Nature Communications 2021, doi:10.1038/s41467-021-21246-9"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Cell-Cell Communication Analysis (CellChat v2)

CellChat v2 turns an annotated single-cell RNA-seq dataset into a map of intercellular signaling: which cell types talk to which, through what ligand–receptor pairs, and along which pathways. This skill wraps that workflow end-to-end — loading, inference, plotting, export — so you can drop in a Seurat object and get a complete results bundle without writing the analysis by hand.

## When to Use This Skill

Use when:
- You have an annotated Seurat object (cell type labels in `meta.data`) and want to characterize cross-population signaling.
- You need to identify which cell types are dominant senders, receivers, mediators, or influencers in specific pathways.
- You want chord diagrams, bubble plots, or signaling-role heatmaps for a manuscript figure.
- You are chaining from an upstream Seurat clustering / annotation workflow and need a `cellchat_object.rds` for downstream comparison work.

Do not use when:
- The dataset is not annotated. Cluster and label cells first (e.g. with `scrnaseq-seurat-core-analysis`).
- You need transcription-factor regulatory networks rather than ligand–receptor interactions — use `grn-pyscenic`.
- You only have bulk RNA-seq counts (CellChat requires single-cell resolution).
- You want spatially constrained communication only — CellChat v2 supports that mode but it requires spatial coordinates and is not the workflow here.

## Quick Start

```r
source("scripts/load_data.R")
source("scripts/run_cellchat.R")
source("scripts/cellchat_plots.R")
source("scripts/export_results.R")

seurat_obj <- load_cellchat_data()                 # PBMC 3k example
cellchat   <- run_cellchat_analysis(seurat_obj)    # ~2–5 min
generate_all_plots(cellchat, output_dir = "results")
export_all(cellchat, seurat_obj = seurat_obj, output_dir = "results")
```

Typical wall time on the bundled PBMC 3k dataset: 3–6 minutes including plotting and export. Larger datasets (20k+ cells) commonly take 10–20 minutes.

## Installation

| Package | Version | License | Commercial use | Install command |
|---|---|---|---|---|
| CellChat | ≥ 2.0.0 | GPL-3 | Permitted | `devtools::install_github("jinworks/CellChat")` |
| Seurat | ≥ 5.0.0 | MIT | Permitted | `install.packages("Seurat")` |
| SeuratData | ≥ 0.2.1 | GPL-3 | Permitted | `devtools::install_github("satijalab/seurat-data")` |
| NMF | ≥ 0.23 | GPL-2+ | Permitted | `install.packages("NMF")` |
| circlize | ≥ 0.4.12 | MIT | Permitted | `install.packages("circlize")` |
| ComplexHeatmap | ≥ 2.12 | MIT | Permitted | `BiocManager::install("ComplexHeatmap")` |
| ggprism | ≥ 1.0.3 | GPL-3 | Permitted | `install.packages("ggprism")` |
| presto | ≥ 1.0 | GPL-3 | Permitted | `remotes::install_github("immunogenomics/presto")` |
| ggalluvial | ≥ 0.12 | GPL-2 | Permitted | `install.packages("ggalluvial")` |
| rmarkdown | ≥ 2.20 | GPL-3 | Permitted | `install.packages("rmarkdown")` |

> CellChat is **not on CRAN**. Install from the `jinworks/CellChat` fork — the original `sqjin/CellChat` repository is archived and no longer receives updates.

## Inputs

**Required**
- An annotated Seurat (v3/v4/v5) object containing:
  - Normalized expression in `data` (slot/layer). The loader runs `NormalizeData()` if the layer looks empty.
  - A metadata column with cell type labels (the workflow defaults to `celltype`).
  - At least 3 cell types, with ≥ 10 cells per type recommended (CellChat will filter smaller groups by default).

**Accepted formats**
- An `.rds` file written by `saveRDS()` (typical of upstream Seurat workflows).
- A live Seurat object already in the R session.
- The built-in PBMC 3k example (`SeuratData::pbmc3k.final`) is loaded automatically when no path is supplied.

**Optional**
- Species selection (`"human"` or `"mouse"`) to pick the right `CellChatDB` table.
- A signaling category filter: `"all"`, `"Secreted Signaling"`, `"ECM-Receptor"`, or `"Cell-Cell Contact"`.

## Outputs

Written to `output_dir/` (default `results/`):

| File | Description |
|---|---|
| `cellchat_object.rds` | Trained CellChat object. Required for multi-condition comparison and pathway deep dives. |
| `significant_interactions.csv` | All significant ligand–receptor pairs with source, target, pathway, probability, p-value. |
| `pathway_summary.csv` | Pathway-level aggregated communication probabilities. |
| `interaction_count_matrix.csv` | Cell type × cell type interaction counts. |
| `interaction_strength_matrix.csv` | Cell type × cell type communication weights. |
| `signaling_roles.csv` | Per-pathway centrality scores (sender, receiver, mediator, influencer). |
| `top_interactions.csv` | Top 20 interactions ranked by communication probability. |
| `interaction_count_network.{png,svg}` | Circle plot of interaction counts. |
| `interaction_strength_network.{png,svg}` | Circle plot of communication strength. |
| `chord_aggregated.{png,svg}` | Chord diagram of the full network. |
| `bubble_ligand_receptor.{png,svg}` | Bubble plot of L-R pairs across cell type pairs. |
| `signaling_outgoing_heatmap.{png,svg}` | Outgoing signaling pattern heatmap. |
| `signaling_incoming_heatmap.{png,svg}` | Incoming signaling pattern heatmap. |
| `signaling_role_scatter.{png,svg}` | Dominant senders vs receivers scatter. |
| `analysis_report.md` | Markdown summary (always generated). |
| `analysis_report.pdf` | PDF version (only when `rmarkdown` + LaTeX/tinytex are installed). |

**Downstream consumers.** `cellchat_object.rds` is the entry point for multi-condition `mergeCellChat()` analysis, for pathway-specific hierarchy plots, and for any custom `subsetCommunication()` extraction. The CSV tables feed cleanly into downstream functional enrichment, network construction, or biomarker selection skills.

## Clarification Questions

### 1. Input source ⚠️ CRITICAL: ASK THIS FIRST
Do you have an annotated Seurat `.rds` file, or should the workflow run on the bundled PBMC 3k demonstration dataset?
- If your file: provide the path. Confirm the dataset is normalized and has cell type labels.
- If demo: the rest of the parameters are pre-set. Skip directly to running the workflow.

### 2. Species (only for user-supplied data)
- Human (use `CellChatDB.human`) — default.
- Mouse (use `CellChatDB.mouse`).

### 3. Cell type column (only for user-supplied data)
Which metadata column holds the cell type labels? Common choices: `celltype`, `cell_type`, `singler_labels`, `predicted.celltype.l2`. Inspect with `colnames(seurat_obj@meta.data)`.

### 4. Signaling scope
Which interaction categories to include?
- All three (Secreted + ECM-Receptor + Cell-Cell Contact) — recommended.
- Secreted signaling only.
- Cell-Cell Contact only.
- ECM-Receptor only.

## Standard Workflow

The four-step pipeline calls the four scripts in order. Each `source()` brings in one function — call it; do not unroll the inline CellChat API.

1. **Load**
   ```r
   source("scripts/load_data.R")
   seurat_obj <- load_cellchat_data()                       # demo
   # seurat_obj <- load_cellchat_data("path/to/seurat.rds") # own data
   ```
   Verification: `"✓ Data loaded successfully! N cells, M cell types"`.

2. **Run CellChat**
   ```r
   source("scripts/run_cellchat.R")
   cellchat <- run_cellchat_analysis(
       seurat_obj,
       species  = "human",
       group.by = "celltype"
   )
   ```
   Verification: `"✓ CellChat analysis completed: N significant interactions across M pathways"`.

3. **Plot**
   ```r
   source("scripts/cellchat_plots.R")
   generate_all_plots(cellchat, output_dir = "results")
   ```
   Verification: `"✓ All plots generated successfully: 6 visualizations saved"`.

4. **Export**
   ```r
   source("scripts/export_results.R")
   export_all(cellchat, seurat_obj = seurat_obj, output_dir = "results")
   ```
   Verification: `"=== Export Complete ==="` followed by a file listing.

If you do not see those verification lines, you wrote inline code instead of calling the helper. Stop, re-`source()` the script, and call the function.

## When Scripts Fail

Apply this hierarchy in order — do not jump past it.

1. **Fix and retry (~90% of failures).** A package is missing or an argument is wrong. Install the package or correct the call and re-run.
2. **Modify the script (~5%).** Edit `scripts/*.R` directly when the underlying logic needs to change (e.g. you need an extra column in `signaling_roles.csv`). Note your change at the top of the file.
3. **Use the script as reference (~4%).** Read it, lift the relevant snippet, adapt to your needs, and cite the source.
4. **Write from scratch (~1%).** Only if the existing function genuinely cannot be coaxed into producing what's needed. Explain why.

Inline ad-hoc CellChat code without trying the script first is the most common mistake — it usually re-implements something the script already does, with subtle differences in slot names that break the next step.

## Decision Points

| Choice | Default | Override when |
|---|---|---|
| `type` (`computeCommunProb`) | `triMean` | Use `truncatedMean` for exploratory analysis or large datasets; `median` when specificity matters more than sensitivity. |
| `min.cells` (`filterCommunication`) | 10 | Drop to 5 to retain rare cell types; raise to 20+ to suppress small-group noise. |
| `population.size` | `FALSE` | Set `TRUE` when cell type abundance varies sharply and you want population-size-weighted probabilities. |
| `db_category` | `"all"` | Restrict to one category to isolate secreted-only or contact-only signaling. |
| Species | `"human"` | Mouse data → `"mouse"` so the correct `CellChatDB.mouse` is used. |
| `do.fast` (`identifyOverExpressedGenes`) | `TRUE` if presto installed | Fall back to standard Wilcoxon if presto is unavailable; results equivalent, just slower. |

## Common Issues

| Issue | Cause | Resolution |
|---|---|---|
| `CellChat` not found | Not installed from GitHub, or wrong fork | `devtools::install_github("jinworks/CellChat")` — the `sqjin/CellChat` repo is archived. |
| `presto` required for Wilcoxon test | `presto` package missing | `remotes::install_github("immunogenomics/presto")`. The script falls back to the standard Wilcoxon test automatically. |
| Zero significant interactions returned | Too few cells per type or overly strict filter | Reduce `min.cells`, merge rare cell types, or check that `group.by` points to the right column. |
| Out-of-memory error on > 50k cells | Working set exceeds available RAM | Downsample with `subset(seurat_obj, downsample = 5000)` or raise R's memory limit. |
| Chord diagram error | `circlize` not installed | `install.packages("circlize")`. |
| SVG export warns "svglite required" | `svglite` not available | Expected — the helper falls back to base R `svg()` and still produces both PNG and SVG. Do not install svglite manually unless you need higher-fidelity output. |
| `group.by` not found | Wrong metadata column name | Inspect `colnames(seurat_obj@meta.data)` and pass the correct column. |
| "no slot of name 'images'" | Seurat v3/v4 object passed to v5 code | Already handled — `UpdateSeuratObject()` is called during loading. If the error persists, run `seurat_obj <- UpdateSeuratObject(seurat_obj)` manually. |
| `NMF` not available | Package not installed | `install.packages("NMF")`. Required for centrality decomposition. |
| PDF report skipped | No LaTeX installation | `install.packages("tinytex"); tinytex::install_tinytex()`. The markdown report is always generated. |
| `subsetCommunication()` returns 0 rows | Pathway aggregation produced no significant probabilities | Re-check that normalized data is present; consider relaxing the default 5% probability filter inside CellChat. |
| Plot titles overlap with chord | Default `par(mar=…)` too tight | Inside the script, set `par(mar = c(2, 2, 3, 2))` before the call. |

## Best Practices

- Validate that `cellchat@idents` matches the levels you expect before plotting; rare-type merges are best done **before** running `createCellChat()`.
- For comparative analyses (disease vs healthy, treated vs untreated), run the four-step pipeline once per condition, save each `cellchat_object.rds`, then `mergeCellChat()` to compare — see `references/cellchat-guide.md`.
- Re-use the same `CellChatDB` subset across conditions to ensure pathway labels align.
- Inspect `cellchat@netP$pathways` after Step 2 to know which pathways are worth drilling into; pathways with very low contribution to total flow are often noise.
- When publishing, report `n_perms` and `min.cells` — both affect the count of "significant" interactions.

## Suggested Next Steps

- **Multi-condition comparison.** Merge two trained CellChat objects with `mergeCellChat()` and use `compareInteractions()`, `netVisual_diffInteraction()`, `rankNet()` to identify rewired signaling. The recipe lives in `references/cellchat-guide.md`.
- **Pathway deep dive.** For each notable pathway in `cellchat@netP$pathways`, generate hierarchy/circle/chord layouts with `netVisual_aggregate()` and the contribution plot with `netAnalysis_contribution()`.
- **Gene regulatory networks.** Run `grn-pyscenic` on the same scRNA-seq data to find transcription factors driving the observed communication.
- **Functional enrichment.** Feed top sender/receiver gene sets into `functional-enrichment-from-degs` for GO / Reactome / KEGG enrichment.
- **Spatial mapping.** If a matched Visium dataset exists, project the inferred signaling onto tissue with `spatial-transcriptomics` plus CellChat's spatial mode.

## Related Skills

| Skill | Relationship |
|---|---|
| `scrnaseq-seurat-core-analysis` | Upstream — produces the annotated Seurat object. |
| `scrnaseq-scanpy-core-analysis` | Alternative upstream (Python) — convert AnnData → Seurat with `MuDataSeurat` or `sceasy` before this skill. |
| `grn-pyscenic` | Complementary — TF-driven gene networks from the same single-cell data. |
| `spatial-transcriptomics` | Companion — spatial context for the communication CellChat infers. |
| `functional-enrichment-from-degs` | Downstream — enrichment on factor- or pathway-specific gene lists. |

## References

- **[CellChat v1]**: Jin S, Guerrero-Juarez CF, Zhang L, Chang I, Ramos R, Kuan C-H, Myung P, Plikus MV, Nie Q. "Inference and analysis of cell-cell communication using CellChat." *Nature Communications* 2021;12:1088. doi:10.1038/s41467-021-21246-9 — [link](https://doi.org/10.1038/s41467-021-21246-9)
- **[CellChat v2]**: Jin S, Plikus MV, Nie Q. "CellChat for systematic analysis of cell-cell communication from single-cell and spatially resolved transcriptomics." *Nature Protocols* 2024. doi:10.1038/s41596-024-01045-4 — [link](https://doi.org/10.1038/s41596-024-01045-4)
- **[LIANA+]**: Dimitrov D, Schäfer PSL, Farr E, et al. "LIANA+: an all-in-one cell-cell communication framework." *Nature Cell Biology* 2024;26:1613–1622. doi:10.1038/s41556-024-01469-w — [link](https://doi.org/10.1038/s41556-024-01469-w)
- **[CellChat v2 source]**: https://github.com/jinworks/CellChat
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)