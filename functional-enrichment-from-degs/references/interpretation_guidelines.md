# Interpreting Enrichment Results

The most common failure mode of pathway analysis isn't producing the wrong p-value — it's reading the output without enough skepticism. This guide walks through how to read a GSEA / ORA result table and the figures that go with it.

## The GSEA result table

| Column | What to look at |
|--------|----------------|
| `pathway` | Is it biologically plausible for this experiment? |
| `setSize` | Should be 15-500. Outside that range → unstable statistics |
| `NES` | Direction + magnitude. > 2 or < -2 is large; 1.5-2 is moderate |
| `p.adjust` | < 0.05 standard. The headline filter |
| `leading_edge` | The core genes driving the score |
| `rank` (in GSEA-format table) | Where in the ranked list the leading edge sits |

A typical sanity check workflow:

1. Sort by `p.adjust` ascending, look at the top 20.
2. Group conceptually: are they all variations on a theme (e.g., interferon response × 6 pathways) or genuinely distinct biology?
3. For each top hit, inspect the leading-edge gene list. Are they genes you've heard of? Are they coherent (e.g., all proteasome subunits)?
4. For each top hit, look at the running-enrichment plot. Sharp peak near the top of the list → robust. Flat / noisy → fragile.

## The ORA result table

| Column | What to look at |
|--------|----------------|
| `pathway` | Plausible? Specific or generic? |
| `GeneRatio` | k/N, the fraction of your significant genes in this pathway |
| `BgRatio` | K/M, the fraction of the background in this pathway |
| `p.adjust` | < 0.05 standard |
| `geneID` | The specific significant genes contributing |

Sanity check:

1. The top ORA pathway shouldn't be a generic catch-all like "GO:BP regulation of biological process" — that's a sign your gene list is too noisy or your background is wrong.
2. For each up-regulated and down-regulated pathway, name the specific genes in `geneID`. Do they match known biology?
3. ORA hits should *concentrate* in your gene list (high GeneRatio relative to BgRatio). A pathway where 40 of your 200 hits live (20%) but accounts for 5% of the background → strong.

## Reading the GSEA dotplot

The standard headline figure: x-axis is NES (or rank-based score), y-axis is the pathway, dot size = setSize, dot color = `p.adjust`.

- **Cluster of pathways with positive NES on one side, negative on the other**: a balanced activation/suppression. Often a clean treatment effect.
- **All NES on one side**: a global expression shift. May reflect real coordinated biology (e.g., proliferative tumor vs quiescent normal) or a normalization artifact. Investigate.
- **Tiny dots dominating**: gene sets at the size floor (15) producing fragile p-values. Apply `min_size = 30` and rerun.

## Reading the running-enrichment plot

For a single pathway, shows the running enrichment score (top panel), gene rank position of pathway members (middle, tick marks), and the ranking metric (bottom).

What to look for:

- **Sharp peak near the top** → pathway robustly activated; the leading-edge genes carry the signal.
- **Sharp valley near the bottom** → pathway robustly suppressed.
- **Symmetric / flat curve** → no real enrichment, despite a "significant" p-value. The pathway may have a few outlier ranks pulling the score up by chance.
- **Many ticks evenly distributed** → pathway members are scattered through the ranking. Not enriched.

## Leading-edge analysis

The leading-edge genes are the core members that drive a positive NES — the genes ranked highest in your list that also belong to the pathway.

For a top GSEA hit, the leading-edge gene list is what you would put in a downstream wet-lab experiment ("we will validate the top 5 leading-edge genes of pathway X by qPCR").

Pull and inspect:

```r
leading <- as.data.frame(gsea_result)["HALLMARK_INFLAMMATORY_RESPONSE", "core_enrichment"]
# Returns a "/"-separated string of gene symbols.
strsplit(leading, "/")[[1]]
```

If the leading edge is 3 genes, the pathway hit is fragile. If it's 30+ genes, it's robust.

## Common interpretation traps

### Trap 1: "p.adjust = 0.04 means it's real"

p.adjust is a threshold, not a confidence rating. A p.adjust of 0.04 with a leading edge of 3 genes is weaker than a p.adjust of 0.08 with a leading edge of 25 genes. Always cross-check with leading edge + running enrichment.

### Trap 2: Top hits are all redundant variations of one theme

GO:BP especially produces this. "Inflammatory response", "innate immune response", "cytokine signaling", "leukocyte activation" — all related, all listing similar genes.

What to do:
- Use `enrichplot::emapplot(gsea_result)` to visualize pathway-pathway similarity; clustering reveals the redundant groups.
- Pick one representative per cluster for the headline figure.
- Don't claim "12 pathways were enriched" if 11 of them are variations on the same thing.

### Trap 3: Pathway names that sound biologically meaningful but aren't specific

`GOBP_REGULATION_OF_RESPONSE_TO_STIMULUS` is true of almost every gene set. If it appears in your top hits, you have a too-noisy gene list or a poorly chosen background.

### Trap 4: Strong NES but unfamiliar pathway

Look up the original publication that defined the pathway (often a stem-cell or cancer cell-line perturbation experiment). The "pathway" may be specific to a context that doesn't apply to your experiment. Use with caution.

### Trap 5: GSEA significant, ORA not, and you stop there

This is the typical case where GSEA is doing its job: detecting a coordinated subtle effect. Report both, explain the discordance, and trust GSEA when the leading-edge gene set is large and coherent.

### Trap 6: All pathways significant

With a very large gene set collection (GO:BP) and a robust signal, you can easily get 100+ significant pathways. The default reporting style is misleading. Instead:

- Apply a stricter p.adjust threshold (0.01 or 0.001).
- Or apply a NES magnitude threshold (`|NES| > 2`).
- Or cluster by similarity (emapplot) and report one per cluster.

### Trap 7: Treating ORA up + down combined as the result

This loses direction. A pathway with 30 upregulated and 25 downregulated genes appears as enriched but isn't "activated" — it's *perturbed* in both directions. Always split.

### Trap 8: No background passed to ORA

`enricher()` defaults to "all genes annotated in the gene set". This means a pathway like "ribosome assembly" (almost every cell has it) appears enriched in nearly any gene list. Pass the actual tested gene set as `universe = `.

### Trap 9: Citing the pathway but not its definition

"HALLMARK_INFLAMMATORY_RESPONSE was activated" is fine. "Inflammation was activated" is over-stating — Hallmark Inflammation is a specific 200-gene curated set, not the entire concept of inflammation in biology. Be precise about what was tested.

### Trap 10: Comparing pathways across studies via the headline list

Different DE methods, different cutoffs, different gene-set versions all change which pathways appear at the top. For cross-study comparison, use rank-based or NES-based comparisons, not "we both found pathway X".

## What to report in methods

For each enrichment claim:

1. Method (GSEA or ORA), tool version (clusterProfiler), permutations (GSEA) or test (ORA).
2. Gene-set source and version (MSigDB version + `msigdbr` version).
3. Collections tested.
4. Size filter (`min_size`, `max_size`).
5. Significance threshold (`p.adjust < 0.05`).
6. Ranking metric (GSEA) or input criteria (ORA: padj and LFC cutoffs, background definition).
7. Number of significant pathways and the headline hits.

## What to report in results

For each headline pathway:

- Pathway name, NES (or GeneRatio for ORA), p.adjust.
- The number of leading-edge genes (or geneID count for ORA).
- A 1-sentence biological interpretation in your domain.

Avoid claiming biology that goes beyond what the test supports. "Inflammation was activated, suggesting an immune response" is fine. "TNF-α drove the response" is a *hypothesis* that needs experimental follow-up.
