# CERES vs Chronos

DepMap shipped two distinct scoring algorithms over its history. CERES
(Meyers 2017) was the original, used through release 20Q4. **Chronos**
(Dempster 2021) replaced it from 21Q1 onward and is the current standard.

This document explains what each algorithm does, where they differ, and how
to reconcile data across the switch.

---

## What Both Algorithms Solve

The raw input to either is a matrix of guide-level log-fold-changes between
the post-screen and the pre-screen library — one number per (sgRNA, cell
line). This is noisy:

1. Multiple guides per gene must be **collapsed** to one gene-level number
2. **Copy-number** confounds it (high-copy amplicons have more cut sites and
   so produce stronger fitness defects unrelated to the gene's importance)
3. **Cell doubling time** scales the dynamic range (slow growers show
   smaller fold-changes overall)
4. **sgRNA efficacy** varies — a guide that doesn't cut produces no signal

Both algorithms model these jointly, but the modelling assumptions differ.

---

## CERES

**Reference**: Meyers RM, et al. (2017) *Nature Genetics* 49(12):1779–1784.

CERES is a constrained linear model. For each guide *i* in cell line *j*:

```
LFC[i,j] = (guide_efficacy[i]) × (gene_essentiality[gene_of_i, j])
           + (copy_number_effect[i,j])
           + noise
```

The constraint: gene-essentiality is rescaled so that the median of a
hand-curated essential reference set is `-1` and the median of a hand-curated
non-essential reference set is `0`. This anchoring is what allows pan-
essentials to behave the same across cell lines.

Strengths:
- Interpretable: each parameter has a clear meaning
- Copy-number correction is built in
- Anchored to a reference essential set

Weaknesses:
- Linear assumption is restrictive
- Each cell line is essentially independent — no population dynamics
- Saturation issues for very strong essentials (extreme LFC values get
  capped by the linear fit)
- Requires a curated essential set as a strong prior — biases against
  context-essentials that are not in the reference

---

## Chronos

**Reference**: Dempster JM, et al. (2021) *Genome Biology* 22:343.

Chronos models the *cell-population dynamics* of the screen directly. For a
cell line *j* with sgRNA *i* targeting gene *g*:

```
n_i,j(t) = n_i,j(0) × exp(growth_rate_j × t × (1 - guide_efficacy_i × gene_effect_g,j))
```

It then jointly infers `guide_efficacy_i`, `gene_effect_g,j`, and growth-
rate constants from the time-resolved guide counts (start, mid, end).

Strengths:
- Genuinely models population dynamics — the right physical model for what
  a fitness screen actually is
- Tighter dynamic range; better discrimination between mild and strong
  dependencies
- Less sensitive to copy-number confounding (corrected within the same
  framework)
- Publishes a **dependency probability** per (gene, line), not just a score
- Scale-anchored to pan-essentials similarly to CERES, but with stronger
  shrinkage

Weaknesses:
- More complex; harder to inspect what went wrong for a specific gene
- Probability calibration depends on the prior cell-line panel composition
- Needs the guide-level time series; cannot be applied to screens that
  only have endpoint counts

---

## What Differs Numerically

For the same (gene, cell line):
- Chronos values are **on the same overall scale** as CERES (`-1` ≈ median
  pan-essential), but the spread of scores is different
- Strong essentials: CERES might say `-1.8`, Chronos `-1.4` for the same line
- Mild effects: CERES is noisier; Chronos compresses these toward 0
- Sign convention is identical — negative = essential in both

The ranks tend to agree well at the top of the essentiality list. They
diverge for marginal cases — exactly the regime where "is this gene
selectively essential?" lives.

---

## Probability of Dependency

Only Chronos releases publish `CRISPRGeneDependency.csv`. It contains a
posterior probability per (gene, cell line) that the gene is a dependency.

How it works (informally): a normal-mixture prior is fit to the gene-effect
distribution (one mode at ≈ 0 for non-dependent genes, another at ≈ -1 for
dependent genes); the per-cell-line probability is the posterior assignment
to the "dependent" mode.

Recommendation: use the probability table instead of thresholding the raw
gene effect when you can. Threshold defaults:
- `prob ≥ 0.5` → call this a dependency (the standard)
- `prob ≥ 0.8` → high-confidence dependency

This is more robust than `gene_effect < -0.5` because it adapts the threshold
to the cell line's overall noise floor.

---

## Reconciling Across the CERES → Chronos Boundary

You cannot directly pool numeric scores from CERES and Chronos releases —
the rank order is largely preserved but the values shift. Strategies:

1. **Just don't pool.** Pick one release, document it, use it throughout.
2. **Rank-rescale.** Convert each release's gene-effect matrix to ranks
   within each gene (or each cell line), then operate on ranks. Loses
   magnitude information but is robust.
3. **Re-anchor to pan-essentials.** Compute the median pan-essential gene
   effect in each release and rescale so both are at `-1`. Imperfect
   because the spread of the rest of the distribution is also different.
4. **Use the probability table for Chronos and binarize CERES.** Threshold
   CERES at `-0.5` for dependency-yes/no, compare to Chronos `prob ≥ 0.5`.

The right answer is usually (1) for new analyses. (2) is good for figures
that span releases. (3) and (4) are last resorts.

---

## When Each Release Was Active

| Release | Algorithm |
|---|---|
| ≤ 20Q4 | CERES |
| 21Q1 onward | Chronos |

The 21Q1 release was the changeover. If you see a published analysis from
2021–2022 reporting "CERES scores", confirm by checking the release tag in
the methods — some 2021 papers used the 20Q4 release even after Chronos
shipped, because the analysis pipeline was already built around CERES.

---

## How To Tell What You've Got

Open `CRISPRGeneEffect.csv` and:
- Check the median of `RPS14`: if it sits at `-1.0 ± 0.2`, the data is
  scale-anchored (could be either algorithm)
- Check for `CRISPRGeneDependency.csv` in the same release: present →
  Chronos; absent → CERES
- Check the spread (e.g. `gene_effect.std()`): CERES typically a bit
  higher (~0.4 vs Chronos ~0.3)
- Check the release notes on the DepMap portal — the definitive answer

---

## Practical Code

```python
def detect_scoring_algorithm(release_dir: str) -> str:
    """Return 'CERES' or 'Chronos' based on files present in a release dir."""
    import os
    files = set(os.listdir(release_dir))
    if "CRISPRGeneDependency.csv" in files:
        return "Chronos"
    if "Achilles_gene_effect.csv" in files and "Achilles_gene_dependency.csv" not in files:
        return "CERES"
    # 21Q1 onwards uses CRISPRGeneEffect.csv + CRISPRGeneDependency.csv
    if "CRISPRGeneEffect.csv" in files and "CRISPRGeneDependency.csv" in files:
        return "Chronos"
    return "unknown"
```

Tag this on every output:

```python
print(f"Loaded {release} ({detect_scoring_algorithm(release_dir)})")
```
