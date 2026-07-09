# CellRank for fate mapping

CellRank uses Markov chains on the cell-similarity graph to estimate fate probabilities. Each cell gets a probability of ending up in each terminal state. The framework integrates multiple sources of directionality — RNA velocity, pseudotime, real-time labels — through kernels.

CellRank 2 is a major rewrite. This guide covers v2 (≥ 2.0).

## When to use CellRank

- You have a directionality signal (velocity, pseudotime, or experimental time points).
- You want fate probabilities, not just trajectory ordering.
- You want to identify driver genes per terminal fate.
- You have a complex branching trajectory where simple pseudotime collapses fates together.

## Kernels (directionality source)

| Kernel              | Source of directionality                              | Use when                                                  |
| ------------------- | ----------------------------------------------------- | --------------------------------------------------------- |
| `VelocityKernel`    | RNA velocity                                          | scVelo ran successfully                                   |
| `PseudotimeKernel`  | Diffusion pseudotime or any monotonic ordering        | No velocity; you have DPT or experimental timepoints     |
| `ConnectivityKernel`| Symmetric similarity graph                            | Backup; rarely the primary signal                         |
| `RealTimeKernel`    | Discrete time-resolved labels                         | Multi-timepoint datasets (day 0, day 3, day 7...)         |
| `CytoTRACEKernel`   | Differentiation potential from gene-count distribution | When no velocity, no clear root cell                     |

Combine kernels with weighted sums:

```python
from cellrank.kernels import VelocityKernel, ConnectivityKernel
vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined = 0.8 * vk + 0.2 * ck
```

The 80/20 velocity/connectivity mix is a common stabilizer when velocity is noisy.

## Terminal states

Terminal states are absorbing nodes of the Markov chain — cells that don't leave (in the model's view). Detected via Schur decomposition.

```python
from cellrank.estimators import GPCCA
g = GPCCA(combined)
g.fit(n_states=10, cluster_key="clusters")
g.predict_terminal_states()
```

`n_states` controls Schur decomposition granularity. Higher → more candidate states; choose by inspecting eigenvalue gap.

Auto-detection can be wrong. Always validate against expected biology:

```python
g.set_terminal_states({"alpha": [cells], "beta": [cells]})
```

## Fate probabilities

Once terminal states are set:

```python
g.compute_fate_probabilities()
g.plot_fate_probabilities(same_plot=False)
```

`adata.obsm["lineages_fwd"]` holds per-cell probabilities for each terminal fate. Rows sum to 1.

Use these to:

- Compute mean fate probability per cluster (composition of cluster X destined for fate Y).
- Identify decision-point cells (probabilities ~ uniform across fates) for downstream DE.
- Compare experimental conditions: do treated cells have different fate distributions?

## Driver genes

```python
drivers = g.compute_lineage_drivers(lineages="alpha")
```

Ranks genes by correlation with fate probability. Top drivers are candidate regulators of that fate.

## Initial states

Symmetric to terminal: cells most likely to be the source. Often the same set as your `root_cell_type`. Useful as a sanity check on the kernel.

## When CellRank breaks

- **Wrong terminal states:** check the PAGA graph against expected mature cell types. Use `g.set_terminal_states({...})` to override.
- **All cells assigned the same fate:** kernel too smoothed, or no real branching. Try a different kernel mix.
- **Drivers look like batch effects:** integration was incomplete. Re-integrate before CellRank.
- **Memory blow-ups:** GPCCA scales poorly past ~ 100k cells. Subsample or use the macrostate approximation.

## Macrostates (large datasets)

For > 100k cells, `g.compute_macrostates(n_states=20)` summarizes the chain at a coarse level before fate computation. Faster, less memory-heavy.

## Citations

1. Lange M, Bergen V, Klein M, et al. *CellRank for directed single-cell fate mapping.* Nat Methods 19:159–170 (2022).
2. Weiler P, Lange M, Klein M, Pe'er D, Theis F. *CellRank 2: unified fate mapping in multiview single-cell data.* Nat Methods 21:1196–1205 (2024).
3. Reuter B, Weber M, Fackeldey K, Röblitz S, Garcia ME. *Generalized Markov state modeling method for nonequilibrium biomolecular dynamics.* J Chem Phys 150:174103 (2019).
