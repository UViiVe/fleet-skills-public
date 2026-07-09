# CellChat Visualization Reference

Six plot families cover almost every CellChat figure you will need. This guide explains when each is most informative and how to customize it without departing from the bundled scripts.

## Plot Selection at a Glance

| Plot | Best for | CellChat function |
|---|---|---|
| Circle plot | Global topology: who-talks-to-whom counts or strengths | `netVisual_circle()` |
| Chord diagram | Detailed flows; pathway-specific routing | `netVisual_chord_cell()` / `netVisual_chord_gene()` |
| Bubble plot | Comparing many L–R pairs across source/target combinations | `netVisual_bubble()` |
| Heatmap | Overall outgoing/incoming pattern by cell type, or pathway contribution | `netAnalysis_signalingRole_heatmap()` |
| Scatter | Identifying dominant senders vs receivers in one glance | `netAnalysis_signalingRole_scatter()` |
| Hierarchy | Single-pathway flow with source/target column layout | `netVisual_aggregate(..., layout = "hierarchy")` |

## Colour Conventions

CellChat's default palette is fine but limited. Override consistently across plots:

```r
cell_types <- levels(cellchat@idents)
colors     <- scPalette(length(cell_types))   # CellChat built-in
names(colors) <- cell_types

# Hand-picked colors keyed by cell type
colors <- c(
    "CD4 T"        = "#E41A1C",
    "CD14+ Mono"   = "#377EB8",
    "B"            = "#4DAF4A",
    "CD8 T"        = "#984EA3",
    "NK"           = "#FF7F00",
    "FCGR3A+ Mono" = "#A65628",
    "DC"           = "#F781BF",
    "Platelet"     = "#999999"
)

netVisual_circle(cellchat@net$count,  color.use = colors, ...)
netVisual_chord_cell(cellchat,        color.use = colors, ...)
```

Reuse the same `colors` named vector across every plot for a coherent figure panel.

## Circle Plot Recipes

```r
# Vanilla — counts.
netVisual_circle(cellchat@net$count, weight.scale = TRUE)

# Vertices sized by cell count (interpretable at a glance).
sizes <- as.numeric(table(cellchat@idents))
netVisual_circle(
    cellchat@net$count,
    vertex.weight = sizes,
    weight.scale  = TRUE,
    label.edge    = FALSE,
    title.name    = "Number of Interactions"
)

# Show only the top 20 strongest edges.
mat <- cellchat@net$weight
mat[mat < sort(mat, decreasing = TRUE)[20]] <- 0
netVisual_circle(mat, weight.scale = TRUE)

# Out-degree view from a specific cell type.
netVisual_circle(
    cellchat@net$weight,
    sources.use  = c("CD14+ Mono", "FCGR3A+ Mono"),
    weight.scale = TRUE,
    title.name   = "Monocyte Outgoing Communication"
)
```

## Chord Diagram Recipes

```r
# Full aggregated weights.
netVisual_chord_cell(cellchat, net = cellchat@net$weight)

# Pathway-specific routing (one pathway, all cells).
netVisual_chord_gene(
    cellchat,
    signaling  = "MHC-II",
    slot.name  = "net",
    title.name = "MHC-II Signaling"
)

# Restrict to monocyte → T cell axis.
netVisual_chord_cell(
    cellchat,
    net          = cellchat@net$weight,
    sources.use  = c("CD14+ Mono"),
    targets.use  = c("CD4 T", "CD8 T"),
    title.name   = "Monocyte → T cell Communication"
)
```

## Bubble Plot Recipes

```r
# All non-isolated interactions.
netVisual_bubble(cellchat, remove.isolate = TRUE)

# Source × target slice.
netVisual_bubble(
    cellchat,
    sources.use    = c("CD14+ Mono"),
    targets.use    = c("CD4 T", "CD8 T", "NK"),
    remove.isolate = TRUE
)

# Specific pathways only.
netVisual_bubble(
    cellchat,
    signaling      = c("MHC-II", "TNF", "CCL"),
    remove.isolate = TRUE
)

# Order by source.
netVisual_bubble(cellchat, remove.isolate = TRUE, sort.by.source = TRUE)
```

## Heatmap Recipes

```r
# Outgoing patterns across all pathways.
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")

# Incoming patterns.
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

# Restrict to selected pathways.
netAnalysis_signalingRole_heatmap(
    cellchat,
    signaling = c("MHC-II", "TNF", "CCL", "CXCL"),
    pattern   = "outgoing"
)

# Per-pathway centrality (sender, receiver, mediator, influencer).
netAnalysis_signalingRole_network(cellchat, signaling = cellchat@netP$pathways)
```

## Scatter Recipes

```r
# Single overview of senders vs receivers across all pathways.
netAnalysis_signalingRole_scatter(cellchat)

# Restrict to one signaling category (e.g. secreted).
netAnalysis_signalingRole_scatter(cellchat, signaling = c("TNF", "CCL", "CXCL"))
```

## Pathway-Specific Deep Dive

```r
pathway <- "MHC-II"

# Three layouts side-by-side for inspection.
par(mfrow = c(1, 3))
netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")
netVisual_aggregate(cellchat, signaling = pathway, layout = "hierarchy")

# Which L–R pairs carry this pathway?
netAnalysis_contribution(cellchat, signaling = pathway)

# Ligand/receptor expression panels.
plotGeneExpression(cellchat, signaling = pathway)
```

## Multi-Condition Plots

After merging CellChat objects from two conditions:

```r
compareInteractions(merged, show.legend = TRUE, group = c(1, 2))
netVisual_diffInteraction(merged, weight.scale = TRUE)
rankNet(merged, mode = "comparison", stacked = TRUE)
netAnalysis_signalingChanges_scatter(merged, idents.use = "CD14+ Mono")
```

## Assembling Publication Figures

```r
# ggplot-based outputs (bubble, scatter) compose with patchwork.
library(patchwork)
p1 <- netVisual_bubble(cellchat, remove.isolate = TRUE)
p2 <- netAnalysis_signalingRole_scatter(cellchat)
combined <- p1 + p2 + plot_layout(widths = c(2, 1))
ggsave("combined.png", combined, width = 16, height = 8, dpi = 300)

# Base R plots (circle, chord) write directly to a device.
png("combined_networks.png", width = 16, height = 8, units = "in", res = 300)
par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))
netVisual_circle(cellchat@net$count,  weight.scale = TRUE, title.name = "Interaction Count")
netVisual_circle(cellchat@net$weight, weight.scale = TRUE, title.name = "Interaction Strength")
dev.off()
```

## Citing the Visualizations

Every plot above comes from CellChat:

- Jin S, et al. **Inference and analysis of cell-cell communication using CellChat.** *Nature Communications* (2021) 12:1088.
- Jin S, Plikus MV, Nie Q. **CellChat for systematic analysis of cell-cell communication from single-cell and spatially resolved transcriptomics.** *Nature Protocols* (2024).
- Function-by-function reference: https://github.com/jinworks/CellChat
