# GSEA vs ORA: When to Use Which

Two methods, two different questions, two different inputs. Both are useful; they're not interchangeable.

## The fundamental difference

**GSEA** asks: *"Among all my genes ranked by effect, do the genes in pathway X cluster near the top or bottom?"* — uses every gene in the experiment.

**ORA** asks: *"Among my significant genes, is pathway X over-represented relative to the background?"* — uses a cutoff-filtered subset.

## GSEA

### What it does

Walks down the ranked gene list. Increments a running score when a gene belongs to the pathway, decrements when it doesn't. The maximum deviation from zero is the enrichment score (ES). Permutation testing produces a p-value; normalization across pathway sizes gives the normalized enrichment score (NES).

### Input

A ranked vector of all genes. The ranking metric in preference order:

1. **Test statistic** (DESeq2 `stat`, limma `t`, edgeR `F`) — balances effect size and statistical confidence.
2. **Signed -log10(p-value)**: `sign(log2FC) × -log10(pvalue)`. Good fallback.
3. **Log2 fold change alone** — last resort; ignores statistical confidence and is dominated by low-count noisy genes.

### Output

| Column | Meaning |
|--------|---------|
| `pathway` | Gene-set name |
| `setSize` | Number of genes in the pathway present in your ranked list |
| `NES` | Normalized enrichment score; positive = activated, negative = suppressed |
| `pvalue` | Raw permutation p-value |
| `p.adjust` | BH-adjusted p-value (use this) |
| `leading_edge` | The core genes driving the enrichment signal |

### When to use

- You have full DE results with log2FC and ideally a test statistic.
- You want to detect coordinated changes across many genes, even if few individually pass significance.
- You want direction (activated vs suppressed) in the output.

### Limitations

- Requires a meaningful ranking metric.
- Computationally more expensive (10,000 permutations is standard).
- Less interpretable to non-bioinformaticians than ORA.
- Sensitive to ties in the ranking; collapse duplicate gene symbols first.

## ORA (Over-Representation Analysis)

### What it does

Asks: of my N significant genes, k are in pathway X. The background of M total genes contains K members of pathway X. Is the observed overlap larger than expected by chance? Fisher's exact / hypergeometric test gives the p-value.

### Input

- A list of significant genes (e.g., `padj ≤ 0.05`, `|log2FC| ≥ 1`).
- A background gene list (the universe — usually all *tested* genes from the DE analysis).

### Output

| Column | Meaning |
|--------|---------|
| `pathway` | Gene-set name |
| `GeneRatio` | k / N (significant genes in pathway / total significant) |
| `BgRatio` | K / M (pathway size / background size) |
| `pvalue` | Hypergeometric raw p |
| `p.adjust` | BH-adjusted |
| `geneID` | The significant genes contributing to the enrichment |

### When to use

- You only have a gene list, no fold changes (e.g., from a ChIP-seq peak-to-gene mapping, or a manually curated list).
- You want a simple, intuitive result you can explain to a wet-lab collaborator.
- You want a sanity check on GSEA findings.

### Limitations

- Requires arbitrary cutoffs to define "significant".
- Loses information about effect size and direction.
- **Must be run separately on up- and downregulated genes** — combining loses directionality and can mask real biology.
- Sensitive to the choice of background gene set. Using "all known human genes" instead of "all tested genes" inflates p-values.
- Underpowered when the significant gene list is small (< 50).

## Side-by-side

| Question | GSEA | ORA |
|----------|------|-----|
| What does the input look like? | All genes, ranked | A gene list + background |
| Need fold changes? | Yes (or test statistic) | No |
| Need a cutoff? | No | Yes (padj + LFC) |
| Captures coordinated subtle effects? | Yes | No |
| Directional? | Yes (positive vs negative NES) | Yes if you split up/down |
| Computationally heavy? | Yes (permutations) | No |
| Easy to explain? | Medium | Easy |
| Sensitive to ties / duplicates? | Yes | Less so |

## Concordance / discordance interpretation

When you run both:

**Concordant** (pathway significant in both methods): the strongest possible signal.

**GSEA yes, ORA no**: the pathway has a coordinated but modest effect — many genes shifted slightly in the same direction, none individually significant. This is the *typical* GSEA-discovers-it-first scenario. Real biology.

**ORA yes, GSEA no**: a few genes with strong individual effects are in the pathway, but the rest of the pathway didn't respond. May indicate a sub-module of the pathway is affected, or the pathway annotation includes too many irrelevant genes.

**Neither**: the pathway is not affected. Or: wrong species, wrong gene-ID type, missing annotation.

## Decision tree

```
Do you have full DE results (gene + log2FC + padj)?
│
├── No, just a gene list ────────────────────────► ORA
│
└── Yes
    │
    ├── Default / no specific request ─────────► GSEA
    │
    ├── Want validation across methods ────────► Both GSEA and ORA
    │
    ├── Few significant genes (< 50) ──────────► GSEA (ORA underpowered)
    │
    └── User explicitly requested ORA ─────────► ORA
        (e.g., for comparison to Enrichr)
```

## Common pitfalls

**GSEA:**

- Ranking by `log2FC` alone — gives noisy genes (low count, high LFC variance) undue weight.
- Filtering before ranking — defeats the point of GSEA. Pass the full table.
- Too few permutations (< 1000) — unstable p-values.
- Duplicate gene names — collapse to a single rank per gene first.

**ORA:**

- Default background = "all known human genes" → inflated p-values. Use the tested-gene set.
- Combining up and down genes in one ORA → masks direction.
- Cutoff too lenient (padj 0.5) → gene list dominated by noise → enrichment of common housekeeping pathways.
- Cutoff too stringent (padj 0.001, |LFC| > 3) → list too small → no power.

## Recommended methods-section phrasing

For GSEA: *"Gene set enrichment analysis was performed using clusterProfiler v[VERSION] with [PERM] permutations, ranking genes by the Wald test statistic. MSigDB Hallmark and Reactome gene sets (v[VERSION]) were tested with size limits 15-500 genes. Significance was defined as BH-adjusted p < 0.05."*

For ORA: *"Over-representation analysis was performed using clusterProfiler v[VERSION] with the hypergeometric test, separately on up- and downregulated genes (DESeq2 padj < 0.05, |log2FC| ≥ 1). All [N] tested genes served as the background universe. Gene sets from MSigDB Hallmark were tested with size limits 10-500 genes. Significance was defined as BH-adjusted p < 0.05."*
