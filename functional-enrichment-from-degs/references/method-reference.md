# Method Reference

Code patterns for GSEA, ORA, and the common variations. Adapt these to your specific analysis.

## Default workflow

```r
library(clusterProfiler)
library(msigdbr)

# Load DE results
de <- read.csv("de_results.csv")              # cols: gene, log2FC, padj, [stat]

# Ranked list for GSEA
df <- de
if (any(duplicated(df$gene))) {
  df <- df[order(df$padj), ]
  df <- df[!duplicated(df$gene), ]
}
if ("stat" %in% colnames(df)) {
  ranks <- df$stat
} else {
  ranks <- sign(df$log2FC) * -log10(pmax(df$padj, 1e-300))
}
names(ranks) <- df$gene
ranks <- sort(ranks, decreasing = TRUE)

# Get Hallmark gene sets
m <- msigdbr(species = "Homo sapiens", category = "H")
term2gene <- m[, c("gs_name", "gene_symbol")]

# Run GSEA
gsea <- GSEA(
  geneList     = ranks,
  TERM2GENE    = term2gene,
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  nPermSimple  = 10000,
  seed         = TRUE
)
head(as.data.frame(gsea))
```

## Variation 1: Multi-category MSigDB

```r
m_hallmark <- msigdbr(species = "Homo sapiens", category = "H")
m_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
m_kegg     <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")

combined <- rbind(
  m_hallmark[, c("gs_name", "gene_symbol")],
  m_reactome[, c("gs_name", "gene_symbol")],
  m_kegg    [, c("gs_name", "gene_symbol")]
)
combined <- unique(combined)
```

## Variation 2: ORA on up- and downregulated genes separately

```r
# Filter significant genes
sig_up   <- de$gene[de$padj <= 0.05 & de$log2FC >=  1]
sig_down <- de$gene[de$padj <= 0.05 & de$log2FC <= -1]
background <- de$gene  # universe = all tested genes

# Up-regulated ORA
ora_up <- enricher(
  gene          = sig_up,
  TERM2GENE     = term2gene,
  universe      = background,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# Down-regulated ORA
ora_down <- enricher(
  gene          = sig_down,
  TERM2GENE     = term2gene,
  universe      = background,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH"
)

# Always run both directions separately, never combined.
```

## Variation 3: GSEA against KEGG via clusterProfiler's built-in `gseKEGG`

If you have Entrez IDs and want KEGG specifically (academic only):

```r
# Convert symbols to Entrez
library(org.Hs.eg.db)
gene_df <- bitr(names(ranks), fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ranks_entrez <- ranks[gene_df$SYMBOL]
names(ranks_entrez) <- gene_df$ENTREZID
ranks_entrez <- sort(ranks_entrez, decreasing = TRUE)

gsea_kegg <- gseKEGG(
  geneList     = ranks_entrez,
  organism     = "hsa",
  minGSSize    = 15,
  maxGSSize    = 500,
  pvalueCutoff = 0.05
)
```

## Variation 4: GO:BP via `enrichGO` / `gseGO`

```r
library(org.Hs.eg.db)

# Convert symbols to ENTREZIDs
gene_df <- bitr(de$gene, fromType = "SYMBOL",
                toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
de2 <- merge(de, gene_df, by.x = "gene", by.y = "SYMBOL")

ranks2 <- de2$stat; names(ranks2) <- de2$ENTREZID
ranks2 <- sort(ranks2, decreasing = TRUE)

gsea_go <- gseGO(
  geneList     = ranks2,
  ont          = "BP",
  OrgDb        = "org.Hs.eg.db",
  keyType      = "ENTREZID",
  minGSSize    = 30,    # tighter for GO:BP
  maxGSSize    = 300,
  pvalueCutoff = 0.05
)
```

## Variation 5: Custom gene sets

Some questions require a gene set that isn't in MSigDB — a kinase substrate list, a specific Drosophila signature, a competitor's gene panel.

```r
# Two-column data.frame: term, gene
custom <- data.frame(
  term = c(rep("KINASE_X_SUBSTRATES", 25),
           rep("BIOMARKER_PANEL", 18)),
  gene = c(my_kinase_substrates, my_panel_genes)
)

gsea_custom <- GSEA(
  geneList     = ranks,
  TERM2GENE    = custom,
  minGSSize    = 5,     # smaller sets OK for custom hypothesis tests
  maxGSSize    = 5000,
  pvalueCutoff = 0.05
)
```

## Variation 6: Mouse data

```r
# Either use mouse MSigDB directly
m_mouse <- msigdbr(species = "Mus musculus", category = "H")
term2gene_mouse <- m_mouse[, c("gs_name", "gene_symbol")]

# Or map mouse symbols to human and use human MSigDB
# (use biomaRt or homologene; sketch below)
library(homologene)
mapping <- homologene(de$gene, inTax = 10090, outTax = 9606)
de_human <- merge(de, mapping[, c("10090", "9606")],
                  by.x = "gene", by.y = "10090")
# proceed with de_human and human MSigDB...
```

## Variation 7: Pre-ranked file format for GSEA-Java

If a collaborator wants the GSEA-Java tool's `.rnk` format:

```r
out <- data.frame(gene = names(ranks), rank = ranks)
write.table(out, "ranked.rnk", sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = FALSE)
```

## Variation 8: Plotting customization

```r
library(enrichplot)

# Dotplot of top pathways
dotplot(gsea, showCategory = 15)

# Running enrichment for one pathway
gseaplot2(gsea, geneSetID = "HALLMARK_INFLAMMATORY_RESPONSE")

# Enrichment map (pathway-pathway similarity)
emapplot(pairwise_termsim(gsea), showCategory = 20)

# Gene-concept network
cnetplot(gsea, foldChange = ranks, showCategory = 5)

# Ridge plot of fold-change distributions per pathway
ridgeplot(gsea)
```

## Variation 9: Comparing two conditions

```r
gsea_A <- run_gsea(ranks_A, term2gene)
gsea_B <- run_gsea(ranks_B, term2gene)

# Merge result tables
merged <- merge(
  as.data.frame(gsea_A)[, c("ID", "NES", "p.adjust")],
  as.data.frame(gsea_B)[, c("ID", "NES", "p.adjust")],
  by = "ID", suffixes = c("_A", "_B"), all = TRUE
)

# Scatter the NES values to see direction agreement
plot(merged$NES_A, merged$NES_B, pch = 20)
abline(0, 1, col = "red")
```

## Variation 10: Methods paragraph template

```
"Functional enrichment analysis was performed using clusterProfiler
(v{VERSION}, Yu et al. 2012) against MSigDB ({VERSION}, Liberzon et al.
2015) via msigdbr (v{VERSION}). Gene set enrichment analysis (GSEA;
Subramanian et al. 2005) was performed on the full ranked gene list,
ranked by the DESeq2 Wald statistic, with {N_PERM} permutations. Over-
representation analysis (ORA) was performed separately on up- and
downregulated gene lists ({UP_N} and {DOWN_N} genes respectively;
DESeq2 padj ≤ 0.05, |log2FC| ≥ 1) against the {BACKGROUND_N} tested
genes as background. Gene sets from MSigDB {COLLECTIONS} were tested
with size restricted to 15-500 members. Benjamini-Hochberg-adjusted
p-values < 0.05 were considered significant."
```
