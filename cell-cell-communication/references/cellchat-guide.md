# CellChat v2: Parameter and Workflow Reference

This guide collects the day-to-day choices that come up when running CellChat: which database subset to use, how to tune communication probability, how to compare conditions, and how to handle awkward datasets (rare cell types, very large samples). It assumes you have already worked through the four-step workflow in `SKILL.md`.

## Ligand–Receptor Database

CellChat ships three categories of annotated interactions in `CellChatDB.human` and `CellChatDB.mouse`:

| Category | Signal type | Representative examples |
|---|---|---|
| Secreted Signaling | Diffusible ligands acting on distant cells | TNF, IL-6, CCL/CXCL chemokines, growth factors |
| ECM-Receptor | Extracellular matrix components engaging cell-surface receptors | Collagens, laminins, fibronectins |
| Cell-Cell Contact | Membrane-anchored ligand/receptor pairs needing physical adjacency | MHC-I/II, Notch–Jagged/Delta, CD80–CD28 |

Restrict to one category with `subsetDB()`:

```r
db_secreted <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
cellchat@DB <- db_secreted
```

The pipeline's `run_cellchat_analysis(..., db_category = "Secreted Signaling")` does the same in one call.

## Communication Probability

`computeCommunProb()` controls the central statistic — the mean expression of each ligand–receptor pair within and between cell groups. CellChat exposes three estimators:

| `type` | Behavior | When to use |
|---|---|---|
| `triMean` (default) | Conservative trimmed average; favors strong, consistent signals | Default for publication-quality figures |
| `truncatedMean` (`trim = 0.1`) | More permissive; surfaces weaker patterns | Large datasets or exploratory passes |
| `median` | Most conservative | When you need maximum specificity (e.g. clinical claims) |

Useful companion arguments:

- `population.size = TRUE` — weights probabilities by cell type abundance. Recommended when group sizes differ by an order of magnitude or more.
- `raw.use = TRUE` — runs on raw counts. Generally leave `FALSE`; CellChat is designed around normalized expression.

After computing probabilities, `filterCommunication(cellchat, min.cells = 10)` removes interactions involving cell groups with fewer than `min.cells` cells. Lower it to retain rare populations.

## Multi-Condition Comparison

CellChat's killer feature is differential signaling: which pathways are gained/lost between conditions.

```r
# 1. Run separately per condition (uses the standard workflow each time).
cellchat_ctrl <- run_cellchat_analysis(seurat_ctrl,  species = "human")
cellchat_trt  <- run_cellchat_analysis(seurat_trt,   species = "human")

# 2. Lift to a shared cell-type label space if the two conditions
#    have non-identical level sets.
shared <- union(levels(cellchat_ctrl@idents), levels(cellchat_trt@idents))
cellchat_ctrl <- liftCellChat(cellchat_ctrl, shared)
cellchat_trt  <- liftCellChat(cellchat_trt,  shared)

# 3. Merge.
merged <- mergeCellChat(
    list(Control = cellchat_ctrl, Treated = cellchat_trt),
    add.names = c("Control", "Treated")
)

# 4. Compare.
compareInteractions(merged, show.legend = TRUE)
netVisual_diffInteraction(merged, weight.scale = TRUE)
rankNet(merged, mode = "comparison", stacked = TRUE)
```

If your conditions are tagged inside one Seurat object:

```r
seurat_list   <- SplitObject(seurat_obj, split.by = "condition")
cellchat_ctrl <- run_cellchat_analysis(seurat_list[["control"]])
cellchat_trt  <- run_cellchat_analysis(seurat_list[["treated"]])
```

## Rare Cell Types

Cell types with fewer than ~10 cells produce unstable communication probabilities and can break the centrality step. Three remedies, in order of preference:

1. **Merge into a parent category** before `createCellChat()`:
   ```r
   seurat_obj$celltype_broad <- as.character(seurat_obj$celltype)
   seurat_obj$celltype_broad[seurat_obj$celltype %in% c("pDC", "cDC1", "cDC2")] <- "Dendritic"
   ```
2. **Raise the threshold**: `filterCommunication(cellchat, min.cells = 20)`.
3. **Drop the cell type outright**:
   ```r
   keep <- !seurat_obj$celltype %in% c("Platelets")
   seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[keep])
   ```

## Large Datasets (> 50k cells)

CellChat aggregates by cell group, so the raw cell count is less critical than total RAM usage during overexpression testing and the centrality decomposition. For very large datasets:

```r
# Downsample to 5,000 cells per type — group means are usually stable above ~500.
seurat_sub <- subset(seurat_obj, downsample = 5000)

# Population-size weighting recovers some signal lost by downsampling.
cellchat <- computeCommunProb(cellchat, type = "triMean", population.size = TRUE)

# Free memory between steps if R is close to the limit.
gc()
```

## Pathway Deep Dive

Once the main pipeline finishes, drill into individual pathways:

```r
cellchat@netP$pathways                              # list active pathways
pathway <- "MHC-II"

# Three layouts for the same pathway.
netVisual_aggregate(cellchat, signaling = pathway, layout = "hierarchy")
netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
netVisual_aggregate(cellchat, signaling = pathway, layout = "chord")

# Which L-R pairs drive this pathway?
netAnalysis_contribution(cellchat, signaling = pathway)

# Expression of the underlying ligands/receptors.
plotGeneExpression(cellchat, signaling = pathway)
```

## Extracting Results Programmatically

```r
df_all   <- subsetCommunication(cellchat)
df_mono  <- subsetCommunication(cellchat, sources.use = "CD14+ Mono")
df_tcell <- subsetCommunication(cellchat, targets.use = "CD8 T")
df_tnf   <- subsetCommunication(cellchat, signaling   = "TNF")
df_pw    <- subsetCommunication(cellchat, slot.name   = "netP")  # pathway-level

# Aggregate matrices already on the object:
cellchat@net$count    # cell type × cell type interaction counts
cellchat@net$weight   # cell type × cell type weighted strength
```

## CellChat Object Anatomy

For reference when you need to access internals directly:

| Slot | Contents |
|---|---|
| `cellchat@data.signaling` | Expression restricted to signaling genes |
| `cellchat@idents` | Cell type labels as a factor |
| `cellchat@DB` | The ligand–receptor database used |
| `cellchat@LR$LRsig` | Overexpressed L–R pairs |
| `cellchat@net$prob` | 3D array: source × target × L–R pair |
| `cellchat@net$pval` | Permutation p-values for each interaction |
| `cellchat@net$count` | Aggregated interaction counts |
| `cellchat@net$weight` | Aggregated interaction weights |
| `cellchat@netP$pathways` | Active signaling pathways |
| `cellchat@netP$prob` | 3D array: source × target × pathway |
| `cellchat@netP$centr` | Per-pathway centrality (out/in/flowbet/info) |

## Primary References

- Jin S, et al. **Inference and analysis of cell-cell communication using CellChat.** *Nature Communications* 12, 1088 (2021). doi:10.1038/s41467-021-21246-9
- Jin S, Plikus MV, Nie Q. **CellChat for systematic analysis of cell-cell communication from single-cell and spatially resolved transcriptomics.** *Nature Protocols* (2024). doi:10.1038/s41596-024-01045-4
- Tutorials: https://github.com/jinworks/CellChat/tree/main/tutorial
