# Validation Framework: What Counts as a Defensible Claim

A pathway enrichment result is a *hypothesis* about biology, not a proof. This document defines the bar for a claim that survives review.

## A defensible claim has all of the following

1. **A coherent leading-edge gene set** (for GSEA) or a meaningful `geneID` set (for ORA). At least 10 genes that are *biologically related* in the same way the pathway implies.
2. **A robust statistical signal**: `p.adjust < 0.05` *and* either a moderate-or-larger NES (|NES| > 1.5 for GSEA) or a meaningful GeneRatio (top 10% of overlap, for ORA).
3. **A reasonable database choice**: the pathway is from a curated, well-cited collection (Hallmark, Reactome, KEGG, GO), not a one-off oddity from C7.
4. **Cross-method consistency**, ideally: both GSEA and ORA pick up the same pathway, or pathways from related biology cluster together.
5. **Biological plausibility**: the pathway is something that could reasonably be affected by the experimental perturbation. "Steroid response in dex-treated cells" is plausible. "Cardiac muscle development in cultured fibroblasts" is suspicious.

## Red flags that should pause a claim

- **Leading edge has < 5 genes.** Statistically a coincidence is hard to rule out.
- **Top pathway is GO:BP "regulation of biological process"** or similar — generic catch-all term whose enrichment usually means the gene list is noisy.
- **Single gene appears in many top pathways**. Probably driven by one outlier in your DE results, not real coordinated biology.
- **All significant pathways are from one large collection** (e.g., GO:BP), no overlap with Hallmark or Reactome — the broader collection is over-fitting.
- **Pathway name doesn't match the experiment**. Cancer pathway in non-cancer tissue, immune-cell-specific pathway in epithelial cells. Investigate before celebrating.
- **Conflicting NES signs across closely-related pathways**. E.g., one Reactome pathway says "interferon signaling activated" and another says "interferon signaling suppressed". Probably the database curation is inconsistent; check both.
- **Significant but the running-enrichment plot is flat.** The score is being pulled by a few outlier ranks. Fragile result.
- **The "right" pathway shows up only with one specific database / threshold / ranking combination.** A robust signal should appear under multiple plausible analysis choices.

## Validation checklist before reporting

For each headline pathway you intend to claim:

- [ ] The leading edge has ≥ 10 coherent genes.
- [ ] The running-enrichment plot (GSEA) shows a clear peak/valley, not a flat curve with one spike.
- [ ] The pathway makes biological sense given the perturbation.
- [ ] At least one specific leading-edge gene is well-known in the relevant biology (passes the "have you heard of it" test).
- [ ] The pathway appears (or has close relatives) in multiple databases.
- [ ] The result is stable when you slightly vary the analysis (different ranking metric, different size filter).
- [ ] You can summarize the finding in one biology-aware sentence without sounding speculative.

If you can't tick most of these, the result is exploratory and should be labeled as such in the manuscript.

## Comparing your result to external evidence

Before publishing a "novel pathway" claim:

1. **Search the literature.** Has this pathway been reported in similar experiments? If yes, you've replicated; if no, ask why.
2. **Check independent datasets.** GEO, Recount, GTEx may have a similar perturbation. Run the same enrichment on their data and check for agreement.
3. **Verify with a single-gene assay.** Pick 2-3 leading-edge genes and check qPCR / Western. A pathway claim with a wet-lab vignette is much stronger than computational alone.
4. **Run a complementary tool.** g:Profiler, Enrichr, or REACTOME's online tool. Cross-tool agreement strengthens the claim; disagreement reveals tool-specific artifacts.

## Reproducibility requirements

For each enrichment result you publish, archive:

- The DE results table that was the input.
- The ranking metric used.
- The exact gene-set library version (`packageVersion("msigdbr")`, MSigDB version number).
- The clusterProfiler version (`packageVersion("clusterProfiler")`).
- The complete parameter set (`min_size`, `max_size`, `nPermSimple`, `pAdjustMethod`).
- The random seed.
- The complete results table (CSV), not just the top-N filtered version.

A colleague should be able to reproduce your headline pathway list exactly from these inputs.

## Methods-section language

Acceptable claim types:

| Claim | Backing required |
|-------|------------------|
| "Pathway X was significantly enriched (padj = Y)" | A standard GSEA / ORA result table |
| "Genes in pathway X were coordinately upregulated" | GSEA NES > 0, p.adjust < 0.05, leading edge ≥ 10 genes |
| "Pathway X is activated upon treatment" | Above + biological plausibility |
| "Pathway X drives the phenotype" | The above + functional validation (knockout, inhibitor, etc.) |
| "Gene Y mediates the pathway X effect" | All of the above plus targeted experimental follow-up |

Most enrichment analyses justify rows 1-3. Rows 4-5 require wet-lab work.

## What not to claim

- "Pathway X was disregulated" — vague; specify direction (activated/suppressed).
- "Our data show that pathway X causes Y" — enrichment is associative, not causal.
- "All genes in pathway X were upregulated" — almost never true; check the leading edge.
- "The most significant pathway is the most important" — significance is not biological priority. Rank by NES magnitude, by leading-edge size, by biological context.
- "No pathways were enriched" — usually means weak signal, wrong database, or bad gene-ID mapping. State *what you tested* and what was *not significant*; don't claim absence of biology.
