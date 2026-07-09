---
name: gene-essentiality
id: "skill_761e46870fea4c83b1ebf6b04eb964d2"
display-name: "Gene Essentiality Analysis"
short-description: "Query, interpret, and integrate DepMap (Cancer Dependency Map) CRISPR knockout essentiality scores — including CERES, Chronos, gene-effect inversion, pan-essential vs selective-essential separation..."
category: functional_genomics
keywords: [essentiality, CRISPR, screen, depmap, fitness]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Query, interpret, and integrate DepMap (Cancer Dependency Map) CRISPR knockout essentiality scores — including CERES, Chronos, gene-effect inversion, pan-essential vs selective-essential separation, and integration with expression / mutation / proteomics. Use when you need to know whether a gene is a dependency in a tumour type, find selective vulnerabilities, or correlate essentiality with another omics layer."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "DepMap portal (Cancer Dependency Map, Broad Institute)"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Meyers RM et al., Nature Genetics 2017, doi:10.1038/ng.3984; Dempster JM et al., Genome Biology 2021, doi:10.1186/s13059-021-02540-7"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
  - version: 1.0.0
    date: 2026-05-20
    event: initial-import
    summary: "Initial public import."
status: active
promoted: true
---

# Gene Essentiality Analysis

The Cancer Dependency Map (DepMap) at the Broad Institute publishes genome-wide CRISPR-Cas9 knockout screens across ~1,100+ cancer cell lines (Achilles project) plus orthogonal RNAi data (DEMETER2). Each gene gets a per-cell-line *gene effect score* — a real-valued number indicating how strongly knocking that gene out affects cell fitness.

This skill covers the three things that go wrong most often when people use these data:
1. **Sign convention** — gene-effect scores are *negative for essential* genes. Forgetting to invert breaks every downstream correlation.
2. **Scoring algorithm** — DepMap shipped CERES through release 20Q4 and switched to **Chronos** from 21Q1 onward. They are not interchangeable.
3. **Pan-essential vs selective** — a gene with a strong knockout phenotype in every cell line is not a useful drug target. Selective essentiality is what you want to identify.

## When to Use This Skill

Use when:
- You want to know **whether a candidate gene is a dependency** in a tumour type, cell line, or cohort
- You need to **rank genes by selective essentiality** for therapeutic prioritization
- You are **correlating gene-effect scores with another layer** (expression, mutation, methylation, proteomics) and need the sign convention and confounders explained
- You are building a **classifier of "context-essential" genes** conditioned on cell state
- You want to **distinguish CERES from Chronos** data when reproducing an older analysis or pooling across releases

Do not use when:
- You are running a **new CRISPR screen** in your own lab — that is `pooled-crispr-screens` (single-cell readout) or a bulk readcount-based pipeline (MAGeCK, BAGEL, drugZ)
- You only need **patient-level mutation / expression**, not perturbation data — use TCGA / GTEx tooling instead
- You want **RNAi knockdown data without CRISPR** — DEMETER2 has its own quirks; this skill mentions but does not specialize in it

## Quick Start

```python
from scripts.load_depmap import load_gene_effect, load_sample_info
gene_effect = load_gene_effect(release="24Q2")          # genes × cell lines, Chronos
samples     = load_sample_info(release="24Q2")          # cell-line metadata

# Single-gene query across all cell lines
mll_dep = gene_effect.loc["MLL1 (4297)"]                # one row, ~1100 lines
print(mll_dep.sort_values().head(10))                    # most dependent lines first

# Rank selective essentials in a subset (e.g. AML lines)
aml_lines = samples.query("OncotreeLineage == 'Myeloid' and OncotreeSubtype == 'AML'").index
from scripts.compute_selectivity import selectivity_score
sel = selectivity_score(gene_effect, subset_lines=aml_lines)
print(sel.head(20))                                      # top AML-selective dependencies
```

Runtime: 1–3 minutes to download a DepMap release (~500 MB cached); milliseconds per query after that.

## Installation

```bash
pip install "pandas>=2.0" numpy scipy scikit-learn matplotlib seaborn
pip install pyarrow                                      # parquet I/O for the release files
pip install depmap-utils                                  # optional convenience client (see Notes)
```

| Package | Version | License | Commercial use | Install |
|---|---|---|---|---|
| pandas | >=2.0 | BSD-3 | yes | `pip install pandas` |
| numpy | >=1.24 | BSD-3 | yes | `pip install numpy` |
| scipy | >=1.10 | BSD-3 | yes | `pip install scipy` |
| scikit-learn | >=1.3 | BSD-3 | yes | `pip install scikit-learn` |
| pyarrow | >=14 | Apache-2.0 | yes | `pip install pyarrow` |
| matplotlib | >=3.5 | PSF | yes | `pip install matplotlib` |
| seaborn | >=0.12 | BSD-3 | yes | `pip install seaborn` |

DepMap data itself is under the [Terra Data Use Agreement](https://depmap.org/portal/terms/) — academic + commercial use permitted with attribution.

## Inputs

Required (the loader scripts fetch these by release tag from `https://depmap.org/portal/api/download/`):

- **`CRISPRGeneEffect.csv`** — primary gene-effect matrix (genes × cell lines), Chronos-scored from 21Q1 onward, CERES-scored before
- **`Model.csv`** (or `sample_info.csv` in older releases) — cell-line metadata: lineage, primary disease, mutations summary
- **`OmicsExpressionProteinCodingGenesTPMLogp1.csv`** (optional) — for correlating with expression
- **`OmicsSomaticMutationsMatrixHotspot.csv`** (optional) — for stratifying essentiality by mutation status
- **`ProteinMatrixAveragedNormalizedAbundance.csv`** (optional) — CCLE proteomics for protein-level integration

Optional user inputs:
- Cell-line cohort definition (CSV with `cell_line_id` and `group` columns)
- Gene list (one symbol per line, or Entrez IDs) to restrict the analysis

## Outputs

Per analysis:
- **`gene_effect_ranked.csv`** — long-form table: gene, lineage, median gene effect, % dependent lines, selectivity score
- **`selective_essentials.csv`** — genes ranked by between-group selectivity (subset-vs-rest), with FDR-corrected p-values
- **`pan_essentials.csv`** — the "useless for therapy" reference set (high dependency in nearly every line)
- **`correlation_matrix.csv`** — gene-vs-feature correlation (Spearman) with **scores inverted before correlation**

Plots:
- **`waterfall_<gene>.png`** — gene-effect distribution across lines with disease colouring
- **`heatmap_selective.png`** — top selective essentials × cohort lines
- **`scatter_<gene>_vs_<feature>.png`** — paired-feature scatter with inversion correctly applied

For downstream skills:
- **`dependencies.parquet`** — long-form `{cell_line, gene, gene_effect, scoring_algorithm, depmap_release}` table. Consumable by **multi-omics-integration**, **lasso-biomarker-panel**.

## Clarification Questions

⚠️ **CRITICAL: ASK THIS FIRST**

1. **Are you using gene-effect scores from a DepMap release published 21Q1 or later (Chronos) or earlier (CERES)?** If the user is reproducing an old paper, they must use the matching release; if they are starting fresh, default to the latest Chronos-scored release. CERES and Chronos scores are not directly comparable on the same scale.
2. **What is the question?** (a) "Is gene X essential here?" → single-gene look-up. (b) "What is selective in cohort Y?" → group-vs-rest selectivity. (c) "What predicts essentiality of gene X?" → multivariate / biomarker. (d) "Are dependencies and expression correlated for gene X?" → bivariate.
3. **Define the cohort.** Lineage (e.g. lung), subtype (e.g. NSCLC adenocarcinoma), mutation status (e.g. KRAS G12C+), or an explicit cell-line list? The selectivity calculation depends entirely on this.
4. **Compare against what?** All other cell lines, a lineage-matched cohort, or a specific control set? (This is the *background* against which selectivity is measured.)
5. **Need to integrate with another omics layer?** Expression, mutation, methylation, proteomics, drug response — affects which auxiliary files to load.

## Standard Workflow

1. **Pick a release.** Use `list_releases()` to see available tags; default to the latest. Note in any output which release was used — DepMap revises annotations on every release.
2. **Load the gene-effect matrix** and the sample metadata. The default loader caches on disk.
3. **Verify the sign convention** by checking known pan-essentials: RPS14, RPL5, POLR2A, EIF3A should have median scores < -1.
4. **Compute pan-essential set** (`pan_essential_genes()` — genes with median effect < -1 in >=80% of lines). These are typically removed before downstream selectivity analysis.
5. **Define the cohort** and the background. Selectivity is `(median_effect_in_cohort) - (median_effect_in_background)` with a Wilcoxon p-value.
6. **Score selectivity** (`selectivity_score()`) — returns a per-gene table sorted by `(magnitude, significance, dependency_fraction)`.
7. **Integrate**, optionally, with another layer:
   - Expression: invert gene-effect first, then Spearman
   - Mutation: Wilcoxon between mutant and wildtype lines
   - Proteomics: same as expression; expect noisier signal
8. **Visualize**: waterfall plot per top hit, heatmap across the cohort, scatter for the chosen integration.
9. **Export** the ranked table, the dependency parquet for downstream skills, and the PDF report.

## When Scripts Fail

1. **Fix the obvious.** Wrong release tag, missing pyarrow, missing internet.
2. **Edit the script.** If your cohort definition needs a column not in the sample-metadata index, add it to the loader.
3. **Reference and adapt.** The selectivity calculation is a few lines of pandas — read the script, write your own variant if you need a different test (e.g. limma's trended t-test).
4. **From scratch.** Only if `depmap.org` is unreachable and you have a local snapshot in a non-standard format.

## Decision Points

Pan-essential definition:
- **Default**: median gene effect < -1.0 AND fraction of dependent lines (effect < -0.5) >= 0.8
- More inclusive: median < -0.5 — risks bleed into "selective in many but not all"
- More exclusive: median < -1.5 — limits to ribosome / proteasome core

Dependency call (per cell line):
- **Default**: gene effect < -0.5
- DepMap publishes a probability-of-dependency table (`CRISPRGeneDependency.csv`) using `(prior * 0.5)` — use it directly if you want probabilities rather than thresholds

Selectivity test:
- **Default**: Wilcoxon rank-sum, BH-FDR-corrected
- For small cohorts (<5 lines): permutation test, 10k iterations
- For large cohorts with continuous covariates: limma-trend on gene effects

Sign convention checklist (run on first load):
- RPS14 median effect < -1.0 ✓
- POLR2A median effect < -1.0 ✓
- GFP / luciferase: not present, but if present should sit ~0

CERES vs Chronos:
- Use Chronos for any release tagged 21Q1 or later. Chronos has tighter dynamic range, less noise from copy-number, and reports posterior probabilities for dependency.
- For multi-release comparisons, the script normalizes both to a common scale by re-aligning to a panel of pan-essentials.

## Common Issues

| Symptom | Cause | Fix |
|---|---|---|
| Correlations come out backwards (e.g. "more expression = less essential" when biology says opposite) | Forgot to invert gene-effect scores before correlating | Apply `essentiality_inverted = -gene_effect` (always) before Pearson/Spearman; see `references/sign_conventions.md` |
| Gene appears essential everywhere; "selectivity" makes no sense | Pan-essential gene | Filter against the pan-essential set with `pan_essential_genes()` before computing selectivity |
| Two analyses on "the same gene" disagree | DepMap symbol vs Entrez resolution differs across releases | Use the `gene_id (entrez_id)` composite key as provided; cross-check with `references/sign_conventions.md#gene-identifiers` |
| Older paper says "high CERES = essential" — but data shows opposite | The paper used a *renamed* score (`gene_effect_normalized`) or showed `-gene_effect` | Read the methods carefully; both directions are out there in publications |
| Score difference across releases for the same gene is large | Genuine algorithm change (CERES → Chronos), or new cell lines added to the panel | Document the release tag in every figure caption; do not pool across CERES/Chronos boundaries |
| Spearman p-values absurdly small (<1e-50) for many genes | Sample-size inflation — n ≈ 1100, so any tiny effect is "significant" | Report effect sizes (rho) alongside p-values; impose magnitude thresholds |
| Cohort selectivity hit list is dominated by genes with very few cells expressing them | Confounding by expression: gene is "essential" only because the few expressing lines belong to the cohort | Require minimum mean expression in both groups before testing; or condition the test on expression |
| Heatmap rows are all yellow or all blue | Forgot to z-score within each gene before plotting | `df.sub(df.mean(axis=1), axis=0).div(df.std(axis=1), axis=0)` then plot |
| Cell-line IDs do not match across files | Older releases use ACH-XXXXXX, newer use ModelID | Use the latest `Model.csv` mapping table; `load_depmap.py` does this for you |
| Probability-of-dependency table missing | Older CERES release | Compute manually with `gene_effect < -0.5` as a coarse boolean; or upgrade release |
| `ImportError: pyarrow` | Optional dependency missing | `pip install pyarrow` — parquet I/O is much faster than CSV for the 17k-gene matrix |

## Best Practices

- **State the release.** Every figure or table should carry the DepMap release tag — `21Q1`, `24Q2`, etc. Algorithms and the panel change.
- **Pan-essentials first.** Always pre-compute and remove the pan-essential set before claiming selectivity. Otherwise your "top hit" is a ribosomal subunit.
- **Invert once, document where.** Either invert in the loader and never again, or never invert and remember to flip every figure caption. Pick one rule and write it down.
- **Cohort size matters.** Selectivity with <5 lines in the cohort is fragile — use permutation, not the asymptotic Wilcoxon p-value.
- **Use the dependency-probability table when available.** Chronos releases publish probabilities directly; thresholding is more robust than thresholding the raw score.
- **Integrate with the right denominator.** A gene-effect Spearman vs expression across all 1,100 lines is not the same question as vs expression within one lineage. State which you are doing.

## Suggested Next Steps

- **multi-omics-integration** — combine essentiality with TCGA expression or proteomics, looking for context-dependent vulnerabilities
- **lasso-biomarker-panel** — given a target gene's essentiality, train a sparse regression on baseline expression / mutation features to predict which patient tumours are likely to respond
- **functional-enrichment-from-degs** — feed top selective essentials into pathway enrichment for biological interpretation
- **pooled-crispr-screens** — if you have your own scCRISPR dataset, this skill is the public-data baseline to compare against
- **survival-analysis-clinical** — once you have a candidate dependency, check whether its expression correlates with patient survival in TCGA

## Related Skills

- **multi-omics-integration** — for joining DepMap to TCGA / GTEx / CCLE in one matrix
- **lasso-biomarker-panel** — for predicting essentiality from baseline features
- **pooled-crispr-screens** — for analysing your own scCRISPR-readout screen
- **pathway-enrichment** — for biological interpretation of dependency rankings

## References

- **[CERES algorithm]**: Meyers RM, Bryan JG, McFarland JM, et al. "Computational correction of copy number effect improves specificity of single guide RNA essentiality screens in cancer cells." *Nature Genetics* 2017;49(12):1779–1784. doi:10.1038/ng.3984 — [link](https://doi.org/10.1038/ng.3984)
- **[Chronos algorithm]**: Dempster JM, Boyle I, Vazquez F, et al. "Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects." *Genome Biology* 2021;22:343. doi:10.1186/s13059-021-02540-7 — [link](https://doi.org/10.1186/s13059-021-02540-7)
- **[DepMap / Achilles]**: Tsherniak A, Vazquez F, Montgomery PG, et al. "Defining a Cancer Dependency Map." *Cell* 2017;170(3):564–576.
- **[DEMETER2 / RNAi]**: McFarland JM, Ho ZV, Kugener G, et al. "Improved estimation of cancer dependencies from large-scale RNAi screens using model-based normalization and data integration." *Nature Communications* 2018;9:4610.
- **[Project Score integration]**: Pacini C, Dempster JM, Boyle I, et al. "A comprehensive clinically informed map of dependencies in cancer cells and framework for target prioritization." *Cancer Cell* 2024;42(2):301–316.
- **[DepMap portal]**: https://depmap.org/portal/
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)

Files under `references/`:
- `depmap_overview.md` — what DepMap is, what is in each release, how to fetch it
- `sign_conventions.md` — gene-effect sign rules, gene-identifier formats, common foot-guns
- `ceres_vs_chronos.md` — algorithm differences, when each applies, how to reconcile across releases
- `selectivity_analysis.md` — pan-essential filter, statistical tests for cohort selectivity, effect-size thresholds
- `integration_with_other_omics.md` — joining DepMap to expression, mutation, methylation, proteomics, drug-response