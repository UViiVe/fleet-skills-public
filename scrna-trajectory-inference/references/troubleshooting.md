# Trajectory Inference Troubleshooting

## Root cell selection

### `KeyError: root cell type not found`

```python
print(adata.obs[cluster_key].unique())
```

The string must match exactly — Python is case-sensitive, watch for trailing spaces.

### Pseudotime looks inverted

You picked the wrong root cell. The trajectory still works, but `dpt_pseudotime` runs the other direction. Verify with a known marker:

```python
sc.pl.scatter(adata, x="dpt_pseudotime", y="EARLY_MARKER")
```

If early marker peaks at pseudotime = 1 instead of 0, swap the root.

### `IndexError: invalid root cell index`

Usually a category dtype issue. Force categorical:

```python
adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")
```

## Neighbor graph

### Pseudotime has infinity values

Disconnected components. Cells unreachable from the root get `inf`. Raise `n_neighbors`:

```python
sc.pp.neighbors(adata, n_neighbors=50, n_pcs=30)
sc.tl.diffmap(adata, n_comps=10)
sc.tl.dpt(adata)
```

If that doesn't fix it, you have a real batch effect that integration didn't fully correct. Run `scvi`/`Harmony` first.

### PAGA graph shows no connectivity (all edges below threshold)

Either the clustering is too granular (clusters too small to overlap) or the neighbor graph is too sparse. Try:

- Coarsen the clustering: drop resolution from 0.8 to 0.4.
- Raise `n_neighbors` for the graph (try 50).
- Lower the visualization threshold: `sc.pl.paga(adata, threshold=0.01)`.

## Diffusion map

### `X_diffmap` is all zeros / degenerate

Cells too similar in low-dim space. Recompute PCA with more variable genes or a different `n_pcs`. Or upstream: re-run HVG selection with more features.

### Strange branches in diffusion map

Diffusion maps can introduce spurious branches when there are outlier cells. Filter outliers first (high QC variance, low neighbor density).

## scVelo

### `ModuleNotFoundError: No module named 'scvelo'`

```bash
pip install scvelo
```

The skill auto-detects and skips velocity if scVelo isn't available; install only if you actually have spliced/unspliced layers.

### `AttributeError: 'AnnData' object has no attribute 'layers'` (spliced/unspliced)

Your data wasn't quantified with intronic reads. Check:

```python
print(list(adata.layers.keys()))
```

If `spliced` / `unspliced` aren't there, re-run alignment with STARsolo `--soloFeatures Gene Velocyto` or kallisto|bustools `lamanno`.

### Dynamical model fails to converge

The script falls back to stochastic mode and logs the choice. Common causes:

- Median unspliced fraction < 0.05 — chemistry doesn't support velocity.
- Too few cells per cluster (< 50) — kinetic parameters can't be estimated.
- Too few variable features — bump `n_top_genes` to 3000.

### Velocity arrows go the wrong way

Phase portrait check:

```python
scv.pl.scatter(adata, var_names=["KNOWN_INDUCED_GENE"])
```

If the gene shows induction dynamics (unspliced > spliced in early cells) but arrows point away from those cells, the fit is wrong. Try:

- Stochastic mode (`scv.tl.velocity(adata, mode="stochastic")`).
- More HVGs.
- Tighter QC on the upstream object (drop ambient soup).

### `velocity_pseudotime` and `dpt_pseudotime` disagree

Common when the root cell choice is biased. Velocity-based pseudotime infers direction from data; DPT requires you to provide it. They should agree on cluster *ordering* even if absolute values differ.

## CellRank

### `ModuleNotFoundError: No module named 'cellrank'`

```bash
pip install cellrank
```

### Terminal states don't match expected mature cell types

Auto-detection via Schur decomposition can pick transient states. Override manually:

```python
g.set_terminal_states({
    "alpha": adata.obs_names[adata.obs["cell_type"] == "Alpha"],
    "beta":  adata.obs_names[adata.obs["cell_type"] == "Beta"],
})
g.compute_fate_probabilities()
```

### `MemoryError` in GPCCA

GPCCA scales poorly. Try macrostates:

```python
g.compute_macrostates(n_states=20)
g.predict_terminal_states()
```

Or subsample for the kernel computation.

### Driver genes look like batch effects

The kernel is using batch-correlated similarity. Re-integrate the data, or constrain the analysis to one batch.

## Pseudotime gene analysis

### `trajectory_genes.csv` is empty

Pseudotime range too narrow or all genes correlated weakly. Common cause: wrong root cell collapsing the trajectory.

```python
print(adata.obs["dpt_pseudotime"].describe())
```

Expect range covering most of [0, 1]. If it's [0, 0.1], the root is sitting in a cluster of cells that's far from the rest of the manifold.

### Gene trends are spiky / non-monotonic

Sparse pseudotime bins. Either:

- Increase `n_bins`.
- Restrict to one branch (compute trends per branch instead of full trajectory).
- Smooth with a higher-degree GAM.

## Cell-cycle confounding

### Cycling cells appear at unexpected pseudotime positions

The cell cycle is its own circular trajectory and confuses DPT. Score and regress:

```python
s_genes = ["MCM5", "PCNA", "TYMS", ...]   # canonical S-phase genes
g2m_genes = ["HMGB2", "CDK1", ...]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pp.regress_out(adata, ["S_score", "G2M_score"])
sc.pp.scale(adata)
sc.tl.pca(adata, n_comps=50)
# then re-run trajectory
```

Or restrict to non-cycling cells (`phase == "G1"`).

## Performance

### Trajectory analysis is slow

- `n_pcs=30` is enough for most datasets; using 100 is wasteful.
- Subsample to 50k cells for prototyping if you have > 100k.
- Set `sc.settings.n_jobs = max(1, os.cpu_count() - 2)`.

### `MemoryError` on large data

- Reduce `n_dcs` (diffusion components): 10 instead of 30.
- Subset to HVGs before trajectory: `adata = adata[:, adata.var.highly_variable].copy()`.
- Use `backed` mode for the upstream AnnData.

## When to ask for help

- scVelo Discourse: <https://discourse.scverse.org/c/scvelo>
- CellRank GitHub: <https://github.com/theislab/cellrank/issues>
- scanpy GitHub: <https://github.com/scverse/scanpy/issues>

Include scanpy/scvelo/cellrank versions, a minimal reproducible example, and the full traceback.
