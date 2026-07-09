# Enrichment Decision Guide

The four decisions in order of importance.

## Decision 1: Method

```
Do you have full DE results with log2FC + (test statistic OR pvalue)?
│
├── Yes ──► GSEA is the default. Run ORA as well if you want validation.
│
└── No (just a gene list) ─► ORA only.
```

Refinements:

- **< 50 significant genes in the input**: ORA is underpowered. Use GSEA on the full ranked list.
- **User explicitly asks for ORA** (e.g., to match a colleague's Enrichr output): honor the request.
- **Both methods together**: gold-standard for publication. Run GSEA on the ranked list; run ORA on the up- and down- subsets.

## Decision 2: Ranking metric (GSEA only)

```
Is there a test statistic column? (DESeq2 stat, limma t, edgeR F)
│
├── Yes ──► Use it. Best metric.
│
└── No ──► Is there a raw pvalue?
    │
    ├── Yes ──► sign(log2FC) × -log10(pvalue)
    │
    └── No ──► sign(log2FC) × -log10(padj)   [last resort]
```

The `prepare_gene_lists.R` script picks the best available metric automatically and reports which.

**Never rank by log2FC alone.** Low-count noisy genes have huge LFCs and would dominate the ranks.

## Decision 3: Database

```
What's your research focus?
│
├── Generic exploratory ─► Hallmark (always start here)
│
├── Want richer pathway detail
│   ├── Commercial work ──────► Hallmark + Reactome
│   └── Academic work ────────► Hallmark + KEGG + Reactome
│
├── Cancer ──────────────────► Hallmark + C6
├── Immunology ──────────────► Hallmark + C7
├── Metabolism (academic) ───► KEGG
└── Need ultra-granular GO ──► GO:BP (apply size filter aggressively)
```

Default to the smallest collection that answers your question. Bigger collections introduce redundancy and dilute the headline list.

## Decision 4: Significance thresholds

### GSEA

`p.adjust < 0.05` is standard. If you have many gene sets (GO:BP), consider `p.adjust < 0.01` to control the headline list size.

Always inspect the *direction* of NES alongside significance:
- NES > 0 = activated in the non-reference condition.
- NES < 0 = suppressed.

### ORA

Standard thresholds for the input gene list:

| padj | \|log2FC\| | When |
|------|-----------|------|
| ≤ 0.05 | ≥ 1 | Default (cell line, drug treatment) |
| ≤ 0.05 | ≥ 0.5 | Tissue / in-vivo (modest effect sizes) |
| ≤ 0.1 | none | Screen / hypothesis-generating |
| ≤ 0.01 | ≥ 2 | High-stringency validation |

ORA enrichment significance: `p.adjust < 0.05` on the resulting hypergeometric tests.

## Quick reference card

| Decision | Default | Override when |
|----------|---------|---------------|
| Method | GSEA | Only a gene list → ORA; want validation → both |
| Ranking | Test statistic | Not present → signed -log10(p) |
| Database | Hallmark | Want signaling depth → +Reactome; cancer → +C6 |
| GSEA p threshold | 0.05 | Many sets (GO:BP) → 0.01 |
| ORA gene threshold | padj ≤ 0.05, \|LFC\| ≥ 1 | Tissue → \|LFC\| ≥ 0.5; screening → padj ≤ 0.1 |
| Gene-set size | 15-500 | GO:BP → 30-300 |
| Permutations | 10,000 | Quick test → 1,000 |
| Background (ORA) | All tested genes | Never "all known genes" |

## Decision points by scenario

**"What pathways are perturbed in my drug treatment?"**
GSEA, Hallmark + Reactome, default thresholds. The headline figure is the dotplot of top NES.

**"Are oncogenic signatures activated in tumor vs adjacent normal?"**
GSEA, Hallmark + C6. Inspect running-enrichment plots for the top hits to confirm coordination.

**"I have 200 ChIP-seq peak-associated genes, no fold changes."**
ORA only. Use the genes within X kb of any tested peak as the background, not the whole genome.

**"Few significant DEGs but consistent direction."**
GSEA. ORA will say "nothing's enriched" because the input list is too small; GSEA can pick up coordinated subtle changes.

**"Reviewer asked for Enrichr comparison."**
Run ORA with the same gene list, Hallmark + KEGG, and report the top pathways alongside the original GSEA. Disagreement is informative; concordance is reassuring.

**"My DE was done in mouse, but I want to relate to human pathways."**
Either: (a) use `species = "mouse"` in `msigdbr` (uses human sets mapped to mouse orthologs), or (b) map your mouse symbols to human orthologs via `biomaRt` and run with `species = "human"`. Option (a) is simpler and equally valid for most purposes.
