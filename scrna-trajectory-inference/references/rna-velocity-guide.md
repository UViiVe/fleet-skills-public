# RNA Velocity (scVelo)

RNA velocity estimates the direction and rate of cell-state change from the ratio of unspliced to spliced mRNA. When unspliced (newly transcribed) counts of a gene exceed steady-state expectations, the gene is being induced — and the cell is moving toward a state with higher expression. The reverse implies repression.

`scVelo` implements this with dynamical and stochastic models. The dynamical model fits per-gene transcriptional kinetics (induction, repression, switching time) and gives the most biological signal when the chemistry supports it.

## When velocity is valid

**Works:**

- 10X Chromium v2 / v3 with intron-retained reads (default for `cellranger count`).
- Drop-seq with intronic reads.
- inDrop with intronic reads.
- Smart-seq2/3 (full-length capture).

**Doesn't work:**

- 5'-tag scRNA-seq without intronic coverage.
- Some CITE-seq protocols that drop intronic reads.
- Nuclei-only data (snRNA-seq) — unspliced fraction is enriched, scVelo assumptions break.

You can confirm by checking that `adata.layers["spliced"]` and `adata.layers["unspliced"]` both have non-trivial signal across genes (median unspliced fraction ~ 0.1–0.3 in cytoplasmic scRNA-seq).

## Generating spliced/unspliced layers

- **STARsolo** with `--soloFeatures Gene Velocyto`.
- **kallisto|bustools** with the `lamanno` workflow.
- **velocyto.py** on a CellRanger BAM (slower but works on existing data).

The skill expects the layers to already be present.

## Dynamical vs stochastic model

| Model       | When to use                                          | Behavior                                                  |
| ----------- | ---------------------------------------------------- | --------------------------------------------------------- |
| Dynamical   | Default for most data                                | Fits induction, repression, switching time per gene       |
| Stochastic  | Fallback when dynamical fails to converge            | Uses second moments; no explicit kinetic parameters       |
| Deterministic | Educational use, not recommended                   | Linear regression in spliced–unspliced space              |

`run_trajectory_analysis.py` tries dynamical first; on failure it falls back to stochastic and logs the choice. Don't override unless you have a reason.

## Preprocessing

scVelo has its own preprocessing (`scv.pp.filter_and_normalize`, `scv.pp.moments`). These computations live on top of whatever your AnnData already has — they don't modify `adata.X`.

```python
import scvelo as scv
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)
```

Use `n_top_genes=2000–3000`. Smaller numbers lose signal; larger numbers add noise.

## Outputs

- `adata.layers["velocity"]` — per-gene velocity per cell.
- `adata.obsm["velocity_umap"]` — velocity vectors projected onto UMAP.
- `adata.obs["velocity_pseudotime"]` — velocity-based pseudotime.
- `adata.obs["latent_time"]` (dynamical only) — gene-shared transcriptional time.
- `adata.var["fit_likelihood"]`, `fit_alpha`, `fit_beta`, `fit_gamma` (dynamical) — kinetic parameters.

## Visualization

- `scv.pl.velocity_embedding_stream(adata)` — flow streams on UMAP.
- `scv.pl.velocity_embedding_grid(adata)` — vector field on grid.
- `scv.pl.velocity_confidence(adata)` — per-cell confidence.
- `scv.pl.scatter(adata, gene)` — phase portrait (spliced vs unspliced).
- `scv.pl.heatmap(adata, var_names=top_genes, sortby="latent_time")` — gene trends.

## Top velocity genes

`scv.tl.rank_velocity_genes(adata, groupby="clusters")` ranks genes by how well their dynamics distinguish each cluster. Use as a candidate list for follow-up.

## Interpreting kinetic parameters

For the dynamical model:

- `fit_alpha`: transcription rate.
- `fit_beta`: splicing rate.
- `fit_gamma`: degradation rate.
- `fit_t_` / `fit_switching`: when the gene switches between induction and repression.

`fit_likelihood > 0.1` is the rough threshold for a well-fit gene. Lower-likelihood genes shouldn't drive interpretation.

## Failure modes

- **Velocity goes the wrong way.** Check the phase portrait of a known marker. If `scv.pl.scatter(adata, "MARKER")` shows the expected up-and-down dynamics but velocity arrows point the opposite direction, the kinetic fit failed — drop to stochastic mode.
- **Latent time is U-shaped.** Indicates branching not captured by the single latent time. Compute per-branch latent time instead.
- **Confidence is uniformly low.** Insufficient unspliced counts. Verify the chemistry generated intronic reads.
- **Pseudotime and velocity pseudotime disagree.** Check the PAGA graph — they should agree on cluster ordering even if exact values differ.

## Citations

1. La Manno G, Soldatov R, Zeisel A, et al. *RNA velocity of single cells.* Nature 560:494–498 (2018).
2. Bergen V, Lange M, Peidli S, Wolf FA, Theis FJ. *Generalizing RNA velocity to transient cell states through dynamical modeling.* Nat Biotechnol 38:1408–1414 (2020).
3. Soneson C, Srivastava A, Patro R, Stadler MB. *Preprocessing choices affect RNA velocity results for droplet scRNA-seq data.* PLoS Comput Biol 17:e1008585 (2021).
4. Bergen V, Soldatov RA, Kharchenko PV, Theis FJ. *RNA velocity — current challenges and future perspectives.* Mol Syst Biol 17:e10282 (2021).
