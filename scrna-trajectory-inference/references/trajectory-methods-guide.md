# Trajectory Inference Methods

A working summary of the main trajectory methods, with notes on when each is the right choice.

## PAGA (partition-based graph abstraction)

PAGA represents a dataset as a graph where nodes are clusters and edges are cluster-to-cluster connectivities derived from the neighbor graph. It reconciles discrete clustering with continuous trajectory inference: clusters that are connected on the graph are likely to be transitionally adjacent in the data.

**Use when:**

- You want a quick topological overview before committing to a specific trajectory method.
- The biology has clear cluster structure (most differentiating systems).
- You're going to layer pseudotime or velocity on top.

**Parameters:**

- `n_neighbors`: 30 is a good default. Higher → smoother graph. Lower → more local.
- `groups`: the cluster column to abstract over.
- `threshold`: connectivity threshold for edge drawing in plots (default 0.05). Lower → more edges shown.

**Outputs:**

- `adata.uns["paga"]["connectivities"]` — cluster × cluster connectivity matrix.
- `adata.uns["paga"]["pos"]` — 2D positions for plotting.

**Reading the graph:**

- Strong edges between clusters X and Y indicate they're transitionally connected.
- Use the graph topology to verify expected lineage structure before pseudotime.
- An edge between cell types you don't expect (e.g., T cells ↔ erythroid) suggests doublets or batch effects.

## Diffusion pseudotime (DPT)

DPT orders cells by random-walk distance from a chosen root cell along the diffusion map. It works best on branching trajectories rooted at a known starting point.

**Use when:**

- You have a biologically clear root (e.g., HSC for hematopoiesis, ductal cells for endocrine pancreas).
- The trajectory has a tree-like topology.
- You want a fast, scanpy-native pseudotime estimate.

**Parameters:**

- `n_branchings`: for explicit branch detection. Usually leave at 0 and rely on the diffusion map for topology.

**Output:**

- `adata.obs["dpt_pseudotime"]` — pseudotime ∈ [0, 1].

**Pitfalls:**

- Wrong root cell inverts pseudotime. The trajectory gene direction reverses.
- Disconnected components produce infinity pseudotimes for cells unreachable from the root. Raise `n_neighbors` or check for batch effects.
- For long trajectories with multiple distinct fates, DPT can compress one fate against another. Inspect with the PAGA graph.

## Monocle3

Monocle3 builds a graph over cells and finds a principal-curve / principal-graph through it. Better than PAGA + DPT for very smooth continuous trajectories with no clear cluster structure (e.g., spermatogenesis, organoid differentiation).

R-only. Convert via `SeuratDisk::Convert("data.h5ad", dest = "h5Seurat")` or `zellkonverter`. This skillpack doesn't include a Monocle3 wrapper; use the official tutorial if needed.

**Parameters:**

- `learn_graph(use_partition = TRUE)`: handle disconnected partitions.
- `order_cells(root_cells = ...)`: anchor pseudotime to a specific cell.

## Slingshot

Slingshot (Bioconductor) fits principal curves through a reduced representation and a graph derived from clustering. Stable lineage assignments, good for branching but not too complex trajectories. R-only.

## Method choice

| Method      | Best for                                 | Pros                                          | Cons                                       |
| ----------- | ---------------------------------------- | --------------------------------------------- | ------------------------------------------ |
| PAGA + DPT  | Branching trajectories with clear root   | Fast, scanpy-native, topology-aware           | Sensitive to root choice                   |
| Monocle3    | Smooth continuous, no clear clusters     | Principal-graph, robust to noise              | R-only; slower                             |
| Slingshot   | Branching with moderate complexity       | Stable lineage assignments                    | R-only; struggles with deep branching      |
| scVelo + CellRank | Datasets with spliced/unspliced layers | Directionality from biology, not assumed root | Requires velocity-capable chemistry        |

For most production workflows in Python: PAGA + DPT + scVelo + CellRank. This is what `run_trajectory_analysis.py` does.

## Trajectory genes

Genes whose expression correlates with pseudotime. Computed by Spearman correlation against `dpt_pseudotime`, BH-corrected for multiple testing.

Filter for:

- `abs(correlation) > 0.3`.
- `FDR < 0.05`.
- Expression in ≥ 10% of cells along the trajectory (avoids ranking noise from rare-expression genes).

Direction:

- Positive correlation: increases along pseudotime (later in differentiation).
- Negative correlation: decreases along pseudotime (earlier).

For branching trajectories, compute per-branch correlations. Mixing branches dilutes the signal.

## Gene trends along pseudotime

GAM (generalized additive model) or LOESS fits of expression vs pseudotime. Smoothed trends are more interpretable than scatter plots.

Parameters:

- `n_bins`: 20–50 bins along pseudotime. More bins → smoother but noisier per-bin estimates.
- Smoothing: GAM with cubic-spline basis is the default.

Use these to identify gene activation/repression windows and to find genes that switch behavior at branch points.

## Cell-cycle confounding

Cycling cells appear at unexpected pseudotime positions because the cell cycle is its own circular trajectory. To handle:

- Filter cycling cells with `sc.tl.score_genes_cell_cycle()`.
- Regress out cell cycle scores during preprocessing.
- Or compute pseudotime on the non-cycling subset.

## Citations

1. Wolf FA, Hamey FK, Plass M, et al. *PAGA: graph abstraction reconciles clustering with trajectory inference.* Genome Biol 20:59 (2019).
2. Haghverdi L, Büttner M, Wolf FA, Buettner F, Theis FJ. *Diffusion pseudotime robustly reconstructs lineage branching.* Nat Methods 13:845–848 (2016).
3. Cao J, Spielmann M, Qiu X, et al. *The single-cell transcriptional landscape of mammalian organogenesis (Monocle3).* Nature 566:496–502 (2019).
4. Street K, Risso D, Fletcher RB, et al. *Slingshot: cell lineage and pseudotime inference for single-cell transcriptomics.* BMC Genomics 19:477 (2018).
5. Saelens W, Cannoodt R, Todorov H, Saeys Y. *A comparison of single-cell trajectory inference methods.* Nat Biotechnol 37:547–554 (2019).
