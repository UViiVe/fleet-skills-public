# Gene Set Database Guide

Which gene-set collection to choose depends on the biology you're after, redundancy tolerance, and licensing constraints.

## Quick map

| Goal | Database |
|------|----------|
| Default exploratory | Hallmark |
| Default publication-quality | Hallmark + Reactome |
| Metabolism focus | KEGG (or Reactome) |
| Detailed signaling | Reactome |
| Comprehensive GO terms | GO:BP (large) |
| Localization | GO:CC |
| Molecular activities | GO:MF |
| Cancer biology | Hallmark + C6 oncogenic signatures |
| Immunology | Hallmark + C7 immunologic signatures |
| Custom hypothesis | User-supplied gene set |

## MSigDB collections (via `msigdbr`)

### Hallmark (H) — DEFAULT

- **50 gene sets.** Curated for coherence and non-redundancy.
- Highly interpretable; ideal first-pass.
- Examples: `HALLMARK_INFLAMMATORY_RESPONSE`, `HALLMARK_E2F_TARGETS`, `HALLMARK_OXIDATIVE_PHOSPHORYLATION`.

```r
term2gene <- get_msigdb_genesets(species = "human", categories = c("H"))
```

When to use: every analysis, as the first lens. If you only run one database, run Hallmark.

### KEGG (C2:CP:KEGG) — license-restricted

- ~186 human gene sets. Strong metabolic coverage.
- Diagrams available on the KEGG website (no other database provides this).
- **Commercial restriction.** KEGG requires a paid license for commercial use. Free for academic.

```r
term2gene <- get_msigdb_genesets(species = "human", categories = c("KEGG"))
```

When to use: academic metabolic-pathway questions; when you specifically want a KEGG pathway diagram.

### Reactome (C2:CP:REACTOME)

- ~1,615 human gene sets. Frequently updated; detailed signaling.
- Hierarchical (parent/child pathways → some redundancy).
- Creative Commons Attribution 4.0 license; commercial use permitted with attribution.

```r
term2gene <- get_msigdb_genesets(species = "human", categories = c("REACTOME"))
```

When to use: detailed signaling pathway analysis; commercial work where KEGG is unavailable.

### GO Biological Process (C5:GO:BP)

- ~7,658 human gene sets. Comprehensive but redundant due to the GO hierarchy.
- Useful for granular process annotation, painful for headline figures.

```r
term2gene <- get_msigdb_genesets(species = "human", categories = c("GO:BP"))
```

When to use: when Hallmark and Reactome are too high-level for your question. Be prepared to apply heavy size filtering (`min_size = 30, max_size = 300`) and to manually curate the top hits.

### GO Molecular Function (C5:GO:MF) and Cellular Component (C5:GO:CC)

- GO:MF (~1,738): enzymatic, binding, catalytic activities.
- GO:CC (~1,006): subcellular localization.

```r
term2gene <- get_msigdb_genesets(species = "human", categories = c("GO:MF"))
term2gene <- get_msigdb_genesets(species = "human", categories = c("GO:CC"))
```

When to use: only if your question is specifically about enzymatic activities (MF) or localization (CC). Otherwise GO:BP usually covers it.

### Oncogenic Signatures (C6)

- ~189 gene sets derived from experimental perturbations of oncogenes and tumor suppressors.
- Useful for cancer biology — e.g., asking whether your tumor RNA-seq resembles a known KRAS-on signature.

```r
term2gene <- get_msigdb_genesets(species = "human", categories = c("H", "C6"))
```

### Immunologic Signatures (C7)

- ~4,872 gene sets from immunology experiments (T-cell activation, etc.).
- Useful for vaccine, infection, or autoimmune work.

```r
term2gene <- get_msigdb_genesets(species = "human", categories = c("H", "C7"))
```

## Combining databases

```r
# Default exploratory: tight, focused
term2gene <- get_msigdb_genesets("human", c("H"))                       # ~50 sets

# Standard publication-quality
term2gene <- get_msigdb_genesets("human", c("H", "REACTOME"))           # ~1,665 sets

# Comprehensive (be ready for many hits)
term2gene <- get_msigdb_genesets("human", c("H", "KEGG", "REACTOME"))   # ~1,851 sets

# Cancer biology
term2gene <- get_msigdb_genesets("human", c("H", "C6"))                 # ~239 sets

# Immunology
term2gene <- get_msigdb_genesets("human", c("H", "C7"))                 # ~4,922 sets
```

More is not better. The bigger the gene-set universe, the more multiple-testing burden, and the more redundant hits you have to dedupe manually.

## Size filtering

```r
gsea_result <- run_gsea(ranked, term2gene, min_size = 15, max_size = 500)
ora_result  <- run_ora(genes, term2gene, background, min_size = 10, max_size = 500)
```

- **< 10-15 genes**: unstable statistics, high false-positive rate.
- **> 500 genes**: too broad to be specific; tends to dominate the rankings without telling you anything.
- **15-500**: standard window.

For GO:BP specifically, tighten to `min_size = 30, max_size = 300` to manage redundancy.

## Species support

MSigDB supports:

- **Human** (Homo sapiens) — most comprehensive.
- **Mouse** (Mus musculus) — human sets mapped to mouse orthologs.

For other species, your options are:

1. Use human/mouse sets with ortholog mapping (works for most mammals; degrades for distant species).
2. Use species-specific tools (e.g., `enrichGO` with the appropriate `OrgDb`).
3. Build a custom gene-set library from PathwayCommons, AnimalTFDB, or BioCyc.

```r
# Mouse
term2gene <- get_msigdb_genesets(species = "mouse", categories = c("H"))

# Custom (read your own two-column term/gene table)
term2gene <- read.csv("my_genesets.csv")  # columns: term, gene
```

## Licensing summary

| Database | License | Commercial use |
|----------|---------|----------------|
| Hallmark (H) | Creative Commons Attribution 4.0 | ✅ With attribution |
| KEGG (C2:CP:KEGG) | KEGG license | ❌ Requires paid license |
| Reactome (C2:CP:REACTOME) | Creative Commons Attribution 4.0 | ✅ With attribution |
| GO (C5:*) | Creative Commons Attribution 4.0 | ✅ With attribution |
| Oncogenic (C6) | Creative Commons Attribution 4.0 | ✅ With attribution |
| Immunologic (C7) | Creative Commons Attribution 4.0 | ✅ With attribution |

For commercial / industry work, default to **Hallmark + Reactome** to stay license-clean.

## Reporting requirements

Methods sections should state:

- MSigDB version (`packageVersion("msigdbr")` and the embedded MSigDB version it carries).
- Which collections were tested.
- Size filter (`min_size`, `max_size`).
- Species.

Example: *"Gene set enrichment analysis was performed using clusterProfiler v4.6.0 against MSigDB v7.5.1 Hallmark and Reactome collections (msigdbr v7.5.1) for Homo sapiens, with gene-set size restricted to 15-500 members."*
