---
name: functional-enrichment-from-degs
id: "skill_b96b1468316043dfb6c86ab5a6616a0f"
display-name: "Functional Enrichment from DEGs"
short-description: "Pathway / gene-set enrichment on differential-expression results using clusterProfiler. GSEA (rank-based, all genes) is the default; ORA (cutoff-based, significant genes only) complements it. Use a..."
category: functional_genomics
keywords: [enrichment, GO, KEGG, pathways, DEGs]
last-updated: "2026-07-07"
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Pathway / gene-set enrichment on differential-expression results using clusterProfiler. GSEA (rank-based, all genes) is the default; ORA (cutoff-based, significant genes only) complements it. Use after running DE to identify which pathways or biological processes drive the signal."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
version: 1.0.1
provenance:
  wraps: "clusterProfiler"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Wu T et al., The Innovation 2021, doi:10.1016/j.xinn.2021.100141"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Functional Enrichment from DE Results

GSEA (rank-based) and ORA (cutoff-based) enrichment using `clusterProfiler` and MSigDB. Reads a DE results table from any of DESeq2 / edgeR / limma, builds the ranked list or filtered hit list, runs the analysis, and produces publication-ready outputs.

## When to Use This Skill

✅ **Use this skill when you have:**
- A differential-expression results table with log2 fold changes and (ideally) test statistics or p-values.
- A specific question: "What pathways / biological processes are activated or suppressed in my comparison?"
- A need to translate a gene-level result into a biological-process interpretation suitable for a figure or a paragraph.

❌ **Don't use this skill for:**
- A gene list without context. ORA can still work, but choose a sensible background; otherwise you'll get arbitrary "enrichment" of common pathways.
- Transcription-factor target enrichment specifically — that has dedicated tools (DoRothEA, ChIP-Atlas, CollecTRI). See `upstream-regulator-analysis` and `chip-atlas-target-genes`.
- Direct network analysis — for protein-protein interaction or co-expression network modules, use the coexpression-network skill, not pathway enrichment.
- Drug-perturbation signatures — use the L1000 / Connectivity Map skills if you have one available.

## Quick Start

```r
source("scripts/load_example_data.R"); de <- load_example_de_results()

source("scripts/prepare_gene_lists.R")
source("scripts/get_msigdb_genesets.R")
source("scripts/run_gsea.R")
source("scripts/generate_plots.R")
source("scripts/export_results.R")

ranked    <- create_ranked_list(de)
term2gene <- get_msigdb_genesets(species = "human", categories = c("H"))   # Hallmark
gsea      <- run_gsea(ranked, term2gene)
generate_all_plots(gsea, output_dir = "enrich_out/figures")
export_all(gsea, ranked, output_dir = "enrich_out")
```

For the bundled airway-dexamethasone example, expect 5-15 significant Hallmark pathways (inflammation, TNF-α signaling, interferon response) in under two minutes.

## Installation

| Package | Version | License | Commercial | Install |
|---------|---------|---------|------------|---------|
| clusterProfiler | ≥ 4.6 | Artistic-2.0 | ✅ | `BiocManager::install("clusterProfiler")` |
| enrichplot | ≥ 1.18 | Artistic-2.0 | ✅ | `BiocManager::install("enrichplot")` |
| msigdbr | ≥ 7.5 | MIT | ✅ | `install.packages("msigdbr")` |
| org.Hs.eg.db | latest | Artistic-2.0 | ✅ | `BiocManager::install("org.Hs.eg.db")` |
| org.Mm.eg.db | latest | Artistic-2.0 | ✅ | `BiocManager::install("org.Mm.eg.db")` |
| ggplot2 | ≥ 3.4 | MIT | ✅ | `install.packages("ggplot2")` |
| ggrepel | ≥ 0.9 | GPL-3 | ✅ | `install.packages("ggrepel")` |
| dplyr | ≥ 1.1 | MIT | ✅ | `install.packages("dplyr")` |

**Database licensing.** KEGG requires a commercial license for commercial reuse. MSigDB Hallmark and GO are Creative Commons Attribution 4.0 (commercial use with attribution). Reactome is Creative Commons Attribution 4.0. For commercial work, prefer Hallmark + Reactome over KEGG.

## Inputs

**Required:**

- **DE results table** with columns:
  - `gene` (HGNC symbols preferred; Ensembl IDs need mapping).
  - `log2FC` (or `log2FoldChange`, `logFC`).
  - `padj` (or `FDR`, `adj.P.Val`).
  - Optional but recommended: `stat` (or `t`, `F`) — the test statistic. Best ranking metric for GSEA.
  - Optional: `pvalue` — used as a fallback when no `stat` available.

**Supported sources** (all auto-detected by `load_de_results.R`):

- DESeq2 `results()` output (CSV).
- edgeR `topTags()` output (CSV).
- limma `topTable()` output (CSV).
- A custom CSV with the column names above.

**Optional:**

- **Species** — human (default), mouse, or "custom" with user-supplied gene sets.
- **Database choice** — Hallmark (default), KEGG, Reactome, GO:BP, C6 (oncogenic), C7 (immunologic).
- **ORA thresholds** — `padj ≤ 0.05`, `|log2FC| ≥ 1` (default).

## Outputs

| File | Contents |
|------|----------|
| `gsea_results.csv` | All GSEA results: pathway, setSize, NES, p-value, p.adjust, leading-edge genes |
| `ora_up_results.csv` | (if ORA run) ORA on upregulated genes |
| `ora_down_results.csv` | (if ORA run) ORA on downregulated genes |
| `gsea_result.rds` | Full clusterProfiler `gseaResult` (for further plotting) |
| `ora_up_result.rds` / `ora_down_result.rds` | `enrichResult` objects |
| `ranked_genes.rds` | The named numeric vector used for GSEA |
| `figures/gsea_dotplot.{png,svg}` | Activated/suppressed pathway summary |
| `figures/gsea_running_score.{png,svg}` | Running-enrichment plot for the top pathway |
| `figures/ora_up_barplot.{png,svg}` | Upregulated ORA barplot |
| `figures/ora_down_barplot.{png,svg}` | Downregulated ORA barplot |
| `enrichment_summary.md` | Human-readable summary with top pathways and interpretation |

## Clarification Questions

**⚠️ CRITICAL: ASK THIS FIRST.**

1. **Is the DE results table available, and does it have a test statistic or just log2FC / padj?** A test statistic (DESeq2 `stat`, limma `t`, edgeR `F`) is the best GSEA ranking metric. If you only have log2FC + padj, we can fall back to `sign(log2FC) × -log10(padj)` — workable but not ideal.

2. **GSEA, ORA, or both?** GSEA is the default — it uses the whole ranked list and detects coordinated changes. ORA is good when you have only a gene list (no fold changes) or want a sanity check with hard cutoffs. Running both adds confidence.

3. **Which databases?** Hallmark (50 well-curated sets, good default), KEGG (metabolic + signaling), Reactome (detailed signaling), GO:BP (comprehensive, large, redundant), C6 (oncogenic, for cancer), C7 (immunology). Hallmark + Reactome is a defensible default; Hallmark alone is fine for a quick first look.

4. **Species?** Human (default) and mouse have direct MSigDB support. Other species need either ortholog mapping or custom gene sets.

5. **For ORA: thresholds and direction?** `padj ≤ 0.05`, `|log2FC| ≥ 1` are standard. **Crucially**: ORA should be run separately on upregulated and downregulated gene lists — mixing them washes out direction.

## Standard Workflow

```r
# 1. Load DE results (auto-detects source tool, standardizes columns)
source("scripts/load_de_results.R")
de <- load_de_results("your_de_results.csv")

# 2. Prepare inputs
source("scripts/prepare_gene_lists.R")
ranked <- create_ranked_list(de)             # for GSEA
sig    <- filter_significant_genes(de,
                                   padj_thresh   = 0.05,
                                   log2fc_thresh = 1.0)  # for ORA

# 3. Get gene sets
source("scripts/get_msigdb_genesets.R")
term2gene <- get_msigdb_genesets(species = "human",
                                 categories = c("H", "C2:CP:REACTOME"))

# 4. Run analyses
source("scripts/run_gsea.R")
gsea <- run_gsea(ranked, term2gene)

source("scripts/run_ora.R")
ora_up   <- run_ora(sig$up,   term2gene, background = sig$background,
                    direction = "upregulated")
ora_down <- run_ora(sig$down, term2gene, background = sig$background,
                    direction = "downregulated")

# 5. Plot + export
source("scripts/generate_plots.R")
generate_all_plots(gsea, ora_up = ora_up, ora_down = ora_down,
                   output_dir = "enrich_out/figures")

source("scripts/export_results.R")
export_all(gsea, ranked,
           ora_up = ora_up, ora_down = ora_down,
           output_dir = "enrich_out")
```

## When Scripts Fail

1. **Fix and retry** — install the missing Bioconductor packages, then re-source.
2. **Modify the script** — for non-standard gene-set categories or custom thresholds, edit the call site.
3. **Use as reference** — for species without MSigDB support, adapt `get_msigdb_genesets.R` to read a custom gene-set file.
4. **Write from scratch** — only for exotic enrichment methods (e.g., `topGO`-style hierarchical pruning, fgsea on transcript-level data).

## Decision Points

**Method.** GSEA when you have full DE results with rankable values. ORA when you only have a gene list. Both for a robust analysis with cross-validation.

**Ranking metric (GSEA).** Test statistic > signed -log10(p) > log2FC. The script will pick the best available and tell you which.

**Database.** Hallmark for first-pass exploration. Hallmark + Reactome for the standard publication-grade analysis. Add GO:BP only if you specifically need granular processes — the size makes results harder to interpret. Avoid KEGG for commercial work without a license.

**Gene-set size filter.** `min_size = 15`, `max_size = 500` is the standard. Smaller sets give unstable statistics; larger sets lose biological specificity.

**ORA background.** Always the set of all *tested* genes, not "all known genes". The wrong background inflates p-values massively.

**ORA direction.** Always separate up vs down. Combining directions hides biology — a pathway enriched for upregulated AND downregulated members usually appears non-significant in a combined test even when both directions tell stories worth telling.

## Common Issues

| Symptom | Cause | Fix |
|---------|-------|-----|
| Zero significant GSEA results | Weak signal, mismatched gene IDs, too few permutations | Verify HGNC symbols match the database; increase `n_perm` to 10,000; try a larger database like GO:BP |
| Zero significant ORA results | Thresholds too stringent; too few significant genes; wrong background | Relax to `padj ≤ 0.1`; if hit list has < 50 genes, GSEA is more appropriate |
| ORA hits look generic / non-specific | Wrong background (used "all known genes" instead of all tested genes) | Pass the actual tested-gene set as `background = `; the function defaults to all genes in your DE table |
| Different tools say different things | Different ranking metrics, different cutoffs | Standardize on one ranking metric; report both methods explicitly |
| Gene symbols not recognized | Wrong species, Ensembl IDs without mapping, outdated symbols | Confirm species; map Ensembl → symbol via `org.Hs.eg.db`; check for legacy symbol aliases |
| Memory issues with GO:BP | 7000+ gene sets | Apply size filter `min_size = 30, max_size = 300`; or switch to Hallmark + Reactome |
| GSEA running slow | High permutation count, many gene sets | `n_perm = 1000` for testing; `n_perm = 10000` for final results |
| KEGG license warning | KEGG requires commercial license | Use Hallmark + Reactome for commercial work |
| GSEA dotplot crowded | Too many significant pathways | Filter top 15-20 by NES magnitude |
| Leading-edge genes don't appear in DE table | Outdated gene-set mappings | `msigdbr::msigdbr()` uses latest; check `msigdb_version()` |
| SVG export fails | `svglite` missing or version conflict | Plotting script falls back to base `svg()` device automatically; PNG is always written |
| Duplicate gene names in input | Multiple Ensembl IDs map to same symbol | `prepare_gene_lists.R` keeps the smallest padj per symbol; documented |

## Best Practices

- **Always rank by test statistic** when running GSEA. It accounts for both effect size and statistical confidence.
- **Always specify a background** when running ORA. Default = all tested genes from your DE table.
- **Always run ORA separately on up- and downregulated genes.** Direction matters.
- **Always cite the database version.** `msigdbr_version()` and the MSigDB version (e.g., 7.5.1) need to appear in methods.
- **Always inspect leading-edge genes** for top GSEA hits. A pathway driven by 3-4 outlier genes is fragile; one driven by 20+ coordinated genes is robust.
- **Always check whether top pathways are biologically plausible.** Inflammatory hits in a dexamethasone treatment? Expected. Inflammatory hits in a metabolism knockout? Investigate.
- **Avoid running every database simultaneously.** Hallmark + one other gives focused, interpretable results. Adding GO:BP + KEGG + Reactome + C2 + C5 + C6 + C7 produces 15,000 pathways and 200 "significant" hits, most of which are slight variations on the same theme.

## Suggested Next Steps

1. **Manuscript figure prep** — `generate_plots.R` produces publication-ready PNG + SVG. The GSEA dotplot is the standard headline figure.
2. **Pathway visualization** — `enrichplot::cnetplot` shows gene-pathway connections; `enrichplot::emapplot` shows pathway-pathway similarity. Both consume the `gseaResult` object directly.
3. **Targeted follow-up** — once you have top pathways, examine specific gene members in your DE results. Verify the direction of leading-edge genes makes biological sense.
4. **Compare to public data** — query Enrichr or g:Profiler with the same gene list as a sanity check. Disagreement is informative.
5. **Combine with upstream regulator analysis** — see `upstream-regulator-analysis` for transcription-factor-level interpretation.

## Related Skills

- [bulk-rnaseq-counts-to-de-deseq2](../bulk-rnaseq-counts-to-de-deseq2/SKILL.md) — produces the DE results this skill consumes.
- [bulk-rnaseq-differential-expression](../bulk-rnaseq-differential-expression/SKILL.md) — methodology on padj, LFC thresholds; read before applying ORA cutoffs.
- [bulk-omics-clustering](../bulk-omics-clustering/SKILL.md) — cluster samples / genes; this skill can then characterize each cluster's enrichment.

## References

- `references/gsea_ora_comparison.md` — when to use which method, with worked decision tree.
- `references/database_guide.md` — MSigDB collection breakdown and selection guidance.
- `references/decision-guide.md` — quick decision tree across all parameters.
- `references/interpretation_guidelines.md` — how to read a GSEA dotplot, leading-edge analysis, common interpretation traps.
- `references/method-reference.md` — extended code patterns and adaptation examples.
- `references/gsea_ora_validation_framework.md` — what counts as a valid claim, how to write the methods paragraph.

**Primary citations:**

- Subramanian A et al. (2005). *Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.* PNAS 102: 15545-15550. DOI: 10.1073/pnas.0506580102
- Yu G et al. (2012). *clusterProfiler: an R package for comparing biological themes among gene clusters.* OMICS 16: 284-287. DOI: 10.1089/omi.2011.0118
- Liberzon A et al. (2015). *The Molecular Signatures Database Hallmark Gene Set Collection.* Cell Systems 1: 417-425. DOI: 10.1016/j.cels.2015.12.004
- Reimand J et al. (2019). *Pathway enrichment analysis and visualization of omics data using g:Profiler, GSEA, Cytoscape and EnrichmentMap.* Nat Protoc 14: 482-517. DOI: 10.1038/s41596-018-0103-9
- Geistlinger L et al. (2021). *Toward a gold standard for benchmarking gene set enrichment analysis.* Brief Bioinform 22: bbz158. DOI: 10.1093/bib/bbz158
