---
name: pooled-crispr-screens
id: "skill_d585c170870946c0a5ea2a23b52afd98"
display-name: "Pooled CRISPR Screen Analysis"
short-description: "Tiered analysis of pooled CRISPR screens with single-cell readout: fast t-test screening, on-target validation, and batch-corrected DE. Use when you have Perturb-seq, CROP-seq, or CRISPRi/a/ko scre..."
category: functional_genomics
keywords: [CRISPR, screen, MAGeCK, pooled, essentiality]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Tiered analysis of pooled CRISPR screens with single-cell readout: fast t-test screening, on-target validation, and batch-corrected DE. Use when you have Perturb-seq, CROP-seq, or CRISPRi/a/ko screens with 10X Feature Barcoding output and want to call hit perturbations with transcriptional phenotypes."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "MAGeCK / BAGEL2 (count-based dropout pipeline, referenced as alternative to scRNA-seq readout)"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Li W et al., Genome Biology 2014, doi:10.1186/s13059-014-0554-4; Hart T & Bhatt DL, G3 2021, doi:10.1093/g3journal/jkab201"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Pooled CRISPR Screen Analysis

This skill walks a pooled CRISPR screen with scRNA-seq readout through three stages: rapid t-test screening across every perturbation, target-gene validation to filter out non-functional guides, and donor-corrected differential expression on the surviving hits via glmGamPoi. The output is a ranked hit list backed by per-perturbation DE tables, volcano plots, and a publication-ready PDF.

## When to Use This Skill

Use when:
- You have a **pooled CRISPRi/CRISPRa/CRISPRko screen** read out by single-cell RNA-seq (Perturb-seq, CROP-seq, ECCITE-seq, direct-capture variants)
- Input is **10X Feature Barcoding `.h5`** plus a sgRNA-to-cell mapping table, or an equivalent assignment file from a CROP-seq pipeline
- You want to **call hit perturbations** by transcriptional phenotype and hand them off to pathway / network / clustering tools
- You have **multi-library experiments** with biological replicates and want batch correction at the DE step

Do not use when:
- The screen was **arrayed** (one perturbation per well) — use bulk DE skills instead
- The readout is **not transcriptomic** (FACS, growth/dropout, protein) — those need MAGeCK, BAGEL, or DrugZ-style count-based pipelines
- The sgRNA has not yet been assigned to cells — run CellRanger Feature Barcoding or a CROP-seq aligner first
- You want to model perturbation-perturbation interactions in combinatorial screens — that is a separate genetic-interaction-manifold problem

## Quick Start

```python
from scripts.load_example_data import load_example_data
data = load_example_data()                  # Papalexi 2021, ~140 MB, cached
adata_list   = data["adata_list"]           # one AnnData per batch
mapping_files = data["mapping_files"]
meta         = data["metadata"]
```

Dataset: Papalexi & Satija 2021 ECCITE-seq CRISPRi in THP-1 monocytic leukemia cells (~20,700 cells x 18,649 genes across 4 retained batches, 25 target genes targeting immune checkpoints with ~4 guides each plus non-targeting controls). End-to-end runtime on a laptop: ~10 minutes.

For offline runs use `dataset="demo"` — a synthetic 1k-cell, 5k-gene set that completes in ~2 minutes.

## Installation

```bash
conda create -n crispr-screen python=3.10
conda activate crispr-screen
pip install "scanpy>=1.9" "anndata>=0.8" pandas numpy scipy scikit-learn diffxpy matplotlib seaborn
pip install reportlab        # PDF report (optional but recommended)
pip install rpy2 mygene      # glmGamPoi bridge + gene-name synonyms (optional)
```

| Package | Version | License | Commercial use | Install |
|---|---|---|---|---|
| scanpy | >=1.9 | BSD-3 | yes | `pip install scanpy` |
| anndata | >=0.8 | BSD-3 | yes | `pip install anndata` |
| diffxpy | >=0.7 | BSD-3 | yes | `pip install diffxpy` |
| scikit-learn | >=1.0 | BSD-3 | yes | `pip install scikit-learn` |
| matplotlib | >=3.5 | PSF | yes | `pip install matplotlib` |
| seaborn | >=0.12 | BSD-3 | yes | `pip install seaborn` |
| reportlab | >=3.6 | BSD | yes | `pip install reportlab` |
| rpy2 | >=3.5 | GPL-2 (linking) | yes (compliance with GPL) | `pip install rpy2` |
| glmGamPoi (R) | >=1.10 | GPL-3 | yes (compliance with GPL) | `BiocManager::install("glmGamPoi")` |
| mygene | >=3.2 | Apache-2.0 | yes | `pip install mygene` |

glmGamPoi additionally needs R >=4.0 and the `zellkonverter` R package for h5ad I/O.

## Inputs

Required:
- **10X Feature-Barcode matrices** (`.h5`, one per library) — output of CellRanger with Feature Barcoding enabled. Contains gene expression and sgRNA UMI counts in a single matrix.
- **sgRNA mapping files** (TSV, one per library) — two columns, no header: `cell_barcode<TAB>sgRNA_id`. Doublets should already be filtered out upstream.

Optional:
- **Gene-name correction table** (`dict[str, str]`) — for legacy sgRNA-library names that have been re-symbolized in Ensembl (e.g. `TMEM55A` → `PIP4P2`).
- **Donor / batch labels** — needed for the glmGamPoi tier to model nested random effects.

sgRNA naming convention assumed by default: `GENE_sgRNAn` (underscore-delimited, gene first). The delimiter and position are CLI-configurable.

See `references/best_practices.md` for library-design and capture-prep details.

## Outputs

Per-perturbation results:
- `initial_screening/<gene>_ttest.csv` — fast t-test results sorted by p-value
- `glmgampoi/<gene>_glmgampoi.csv` — batch-corrected NB GLM results for surviving hits
- `figures/volcano_<gene>.png|.svg` — volcano plot per hit

Summary artefacts:
- `hit_list.csv` — validated hits with DE-gene counts and outlier fractions
- `hit_summary.csv` / `per_perturbation_stats.csv` — full perturbation table including non-hits
- `validation_results.csv` — target-gene effect direction, fold-change, status per perturbation
- `screening_summary.txt` — hit-rate distribution as plain text
- `crispr_screen_report.pdf` — Introduction / Methods / Results / Conclusions report with embedded figures (falls back to text-only if reportlab is missing)

Processed objects (downstream consumers):
- `adata_processed.h5ad` — normalized AnnData with `obs["gene"]`, `obs["sgRNA"]`, `obs["batch"]`. Loadable by **functional-enrichment-from-degs**, **coexpression-network**, **bulk-omics-clustering**.
- `analysis_objects.pkl` — pickled `{perturbation_summary, de_results, outlier_cells}` dict for direct hand-off.

## Clarification Questions

⚠️ **CRITICAL: ASK THIS FIRST**

1. **Do you want to test with the example Papalexi 2021 dataset, or analyze your own data?** If example data: parameters are fixed (CRISPRi, THP-1 immune cells, expected direction = down) and questions 2–6 should be skipped. If user data: proceed.
2. **Screen modality**: CRISPRi (knockdown), CRISPRa (activation), or CRISPRko (knockout)? This determines the expected direction of the target-gene effect at the validation step.
3. **Library scale and design**: number of biological replicates, number of unique perturbations (~10, 100, 1000+?), sgRNAs per gene, and whether non-targeting controls are present.
4. **Cell type**: primary neurons, iPSC-derived, immune cells, cancer cell lines, organoids — this determines the QC thresholds you should adopt (see Decision Points).
5. **sgRNA naming**: does the sgRNA ID embed the target gene symbol (`GENE_sgRNAn`)? If not, a gene-lookup table is required.
6. **Downstream goals**: hit identification only, or pathway enrichment / network inference / clustering as well?

## Standard Workflow

1. **Load libraries.** `load_multiple_libraries()` reads each `.h5` into an AnnData and de-duplicates gene names.
2. **Map sgRNAs to cells.** `map_sgrna_to_adata()` does an inner join of cell barcodes with the mapping TSV, parses the gene name out of the sgRNA ID, and drops unmapped cells.
3. **QC filter per library.** `apply_qc_filters()` calculates `n_genes`, `n_counts`, `percent_mito` and applies cell-type-appropriate cutoffs.
4. **Concatenate.** `concatenate_libraries()` joins libraries with batch labels using outer-joined gene sets.
5. **Correct gene-symbol mismatches.** `correct_gene_names()` (and `suggest_corrections()` if mygene is available) renames legacy aliases that no longer match the expression matrix.
6. **Filter low-expression genes.** `filter_by_expression()` drops genes whose group-wise mean falls below threshold — pure noise reduction.
7. **Normalize.** `normalize_only()` runs counts-per-million normalization with high-expressor exclusion and log1p. Keep raw counts available for the glmGamPoi tier.
8. **Fast t-test screen.** `screen_all_perturbations()` runs Welch's t-test per perturbation vs control, FDR-adjusts, and writes per-gene CSVs.
9. **Call preliminary hits.** `call_hits(min_de_genes=10)` filters by DE-gene count and |log2FC|.
10. **Validate target effect.** `validate_target_effect()` checks that the target gene's own expression moved in the expected direction at the expected magnitude.
11. **Run glmGamPoi on validated hits.** `run_glmgampoi_batch()` fits a negative-binomial GLM with batch as a fixed effect for the top ~50 hits.
12. **Visualize.** `create_volcano_plots()` per hit, `plot_hit_summary()` for the whole screen.
13. **Export.** `export_all_results()` writes CSVs, the h5ad, the pickle bundle, and the PDF report.

The scripts are designed to compose — pass the output of one into the next. Avoid writing inline replacements for the helpers; the verification messages they print (`✓ ...`) are how the workflow signals progress.

## When Scripts Fail

Resolve in this order:
1. **Fix the obvious cause** — missing package, wrong delimiter, wrong column name. Re-run the same script.
2. **Edit the script in place** — if a parameter needs to differ from the default in a way the CLI does not expose, change the script and note the change.
3. **Reference, then adapt** — copy the logic into your notebook only if subprocess/IO is the actual blocker (e.g. an HPC environment without R for glmGamPoi). The t-test tier alone is publishable.
4. **Write from scratch** — last resort; document why.

## Decision Points

QC thresholds by cell type:

| Cell type | min_genes | min_counts | max_mito_pct |
|---|---|---|---|
| Primary neurons | 2500–3500 | 8000–12000 | 0.10–0.15 |
| iPSC-derived | 2000–2500 | 5000–8000 | 0.15–0.20 |
| Immune cells | 1000–1500 | 3000–5000 | 0.15–0.25 |
| Cancer cell lines | 1500–2500 | 5000–8000 | 0.15–0.25 |
| Macrophages (activated) | 1000 | 3000–5000 | 0.18–0.20 |

Hit-calling defaults:
- `min_de_genes = 10` and `|log2FC| > 0.5` for the screening tier — generous, designed to retain hits for validation
- Validation thresholds: log2FC < -0.5 for CRISPRi, > 0.5 for CRISPRa, target p < 0.05
- Outlier-fraction call: use LocalOutlierFactor with `contamination='auto'`; switch to IsolationForest if the control outlier rate exceeds 5%

When to override:
- **Genome-wide screens (>1000 perturbations)**: raise `min_de_genes` to 20 and require ≥2 sgRNAs per gene to agree
- **Targeted screens (<100 perturbations)**: use Bonferroni rather than FDR for the screening tier
- **Lethal / essential genes**: skip validation (target is depleted, not just modulated); flag in the hit list manually

See `references/statistical_methods.md` for derivations and `references/qc_guidelines.md` for the full QC reasoning.

## Common Issues

| Symptom | Cause | Resolution |
|---|---|---|
| sgRNA mapping rate < 30% | Mapping file mismatch, low capture, MOI < 0.2 | Verify barcode format alignment; check Feature Barcode library on Bioanalyzer; re-titrate virus |
| Doublet rate > 15% | MOI > 0.8, ambient sgRNA RNA | Tighten singlet threshold to >70% UMI from top sgRNA; in future runs lower MOI to 0.3–0.5 |
| Target gene not differentially expressed | Incomplete knockdown, inactive sgRNA, paralog buffering | Confirm dCas9 expression; validate sgRNA by qPCR; check for paralogs that compensate |
| Validation rate < 50% | Weak perturbation system or wrong `expected_direction` | Re-check CRISPRi vs CRISPRa setting; loosen log2FC threshold to -0.3 |
| Control outlier rate > 10% | Batch effects driving false positives | Run Harmony on `X_pca` before outlier detection; or raise `n_neighbors` |
| Too many hits (>50% of perturbations) | Permissive contamination or strong batch signal | Drop contamination to 0.03; check if PC1 is batch-driven |
| glmGamPoi subprocess fails | R, glmGamPoi, or zellkonverter not installed | Install R 4.0+ then `BiocManager::install(c("glmGamPoi","zellkonverter"))`; or stop at the t-test tier — it is still defensible |
| MemoryError during DE | >100k cells in a single AnnData | Switch to backed-mode `sc.read_h5ad(..., backed="r")`; or process perturbations in batches of 100 |
| KeyError on target gene during validation | Gene-symbol mismatch (legacy alias) | Run `suggest_corrections(use_synonyms=True)` and apply via `correct_gene_names()` |
| Inconsistent hits across replicates | Population drift, batch-specific dropout | Compute hit-overlap CV; analyze replicates separately; check growth bottlenecks |
| Excel export refuses ≥256 sheets | Hitting Excel's sheet limit on genome-wide screens | Split into multiple xlsx files of ~100 perturbations each, or emit CSV |
| SVG export silently absent | Matplotlib SVG backend missing | Harmless — PNG is always saved as fallback |

See `references/troubleshooting.md` for full diagnostics with code snippets.

## Best Practices

- **Inspect the controls before celebrating hits.** A clean screen has a non-targeting outlier rate below 5% and a tight control distribution in PCA.
- **Validate before enriching.** Pathway enrichment on unvalidated hits gives confidently wrong answers.
- **Keep the raw counts.** glmGamPoi needs raw counts, not the normalized layer. Save `adata.raw` before scaling.
- **Document the cutoffs.** Hit calling is a chain of thresholds; record each one (min_de_genes, log2FC, FDR, contamination) in the run config.
- **Run two-stage QC if you expect lethal phenotypes.** Strict QC on controls, relaxed on perturbed populations — see `references/qc_guidelines.md`.

## Suggested Next Steps

- **functional-enrichment-from-degs** — feed the per-hit DE tables (top 200 genes by FDR) into GSEA / ORA to assign each hit a pathway label. Direct consumer of `glmgampoi/*.csv`.
- **coexpression-network** — build WGCNA modules from the perturbed-cell pseudobulk and project the hit set onto modules to recover regulatory neighborhoods.
- **bulk-omics-clustering** — cluster perturbations by their transcriptional signature to identify functional groups (a perturbation manifold), useful for combinatorial follow-up.
- **upstream-regulator-analysis** — given a DE list, infer the transcription-factor regulons that explain it; complements per-hit pathway calls.

## Related Skills

- **scrnaseq-scanpy-core-analysis** — for upstream QC and clustering if your data is not yet AnnData
- **functional-enrichment-from-degs** — primary downstream consumer
- **coexpression-network** — module-based interpretation
- **bulk-omics-clustering** — perturbation-similarity grouping

## References

Files under `references/`:
- `best_practices.md` — experimental design, MOI titration, cell-number planning
- `qc_guidelines.md` — cell-type QC thresholds, sgRNA assignment criteria, batch-balance checks
- `statistical_methods.md` — t-test, Wilcoxon, negative-binomial GLM, MAST; outlier-detection comparison; FDR
- `troubleshooting.md` — diagnostics for every common failure mode
- `umi_optimization.md` — UMI-threshold scan for sgRNA assignment in low-capture screens

Primary citations:
- Dixit A, et al. (2016) Perturb-Seq: Dissecting molecular circuits with scalable single-cell RNA profiling of pooled genetic screens. *Cell* 167(7):1853–1866.
- Datlinger P, et al. (2017) Pooled CRISPR screening with single-cell transcriptome readout. *Nature Methods* 14(3):297–301.
- Replogle JM, et al. (2020) Combinatorial single-cell CRISPR screens by direct guide RNA capture and targeted sequencing. *Nature Biotechnology* 38(8):954–961.
- Papalexi E, et al. (2021) Characterizing the molecular regulation of inhibitory immune checkpoints with multimodal single-cell screens. *Nature Genetics* 53(3):322–331.
- Ahlmann-Eltze C, Huber W (2021) glmGamPoi: fitting Gamma-Poisson generalized linear models on single-cell count data. *Bioinformatics* 36(24):5701–5702.
- Peidli S, et al. (2024) scPerturb: harmonized single-cell perturbation data. *Nature Methods* 21(3):531–540.

External resources:
- 10X Feature Barcoding: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/algorithms/crispr
- scPerturb data portal: https://scperturb.org

## References

- **[MAGeCK]**: Li W, et al. "MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens." *Genome Biology* 2014;15:554. doi:10.1186/s13059-014-0554-4 — [link](https://doi.org/10.1186/s13059-014-0554-4)
- **[BAGEL2]**: Hart T, Bhatt DL. "BAGEL2: a computational framework for identifying essential genes from pooled library screens." *G3* 2021;11:jkab201. doi:10.1093/g3journal/jkab201 — [link](https://doi.org/10.1093/g3journal/jkab201)
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)