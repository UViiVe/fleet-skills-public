# Selectivity Analysis

The interesting biology in DepMap is not "what is essential" — many genes
are essential — but "what is *selectively* essential in some context."
Selectivity is the property that makes a dependency actionable as a
therapeutic target: a gene whose loss kills tumour cells of type X but
spares everything else has the right cancer-specificity profile.

This document covers how to compute and interpret selectivity, with special
attention to the statistical pitfalls that show up at DepMap's sample size.

---

## Definitions

For each gene *g*:
- **Pan-essential**: median gene-effect < `-1` AND >= 80% of lines are
  dependent. These are the housekeeping fitness genes — translation,
  splicing, proteasome, basic transcription. Not actionable; filter out
  first.
- **Selectively essential**: a non-pan-essential gene with significantly
  stronger dependency in one cohort than in the background.
- **Non-essential**: no dependent lines.

Pan-essentials are the noise floor of "interesting selectivity" — you want
to identify them and remove them before any cohort-based ranking, because
otherwise they dominate every ranking.

---

## Pan-Essential Filter

```python
def pan_essential_genes(gene_effect, median_thr=-1.0, dep_frac_thr=0.80,
                        dep_score_thr=-0.5):
    median = gene_effect.median(axis=1)
    dep_frac = (gene_effect < dep_score_thr).mean(axis=1)
    return (median < median_thr) & (dep_frac >= dep_frac_thr)
```

Typical counts on a modern DepMap release:
- ~1,500–2,000 genes classified as pan-essential
- These overlap heavily with curated essential-gene lists (Hart, BAGEL,
  CEGv2)

If you are working with a "well-known cancer target" and it gets caught in
the pan-essential filter, that is informative — it means the target is
broadly essential in cancer cell lines and selectivity claims need careful
support. (See e.g. MYC: pan-essential under most thresholds.)

---

## Cohort-vs-Rest Selectivity Test

The basic test: given a cohort of cell lines and the rest of the panel,
ask whether the gene's effect distribution differs between the two.

### Wilcoxon rank-sum (default)

Non-parametric, fast, fits the typical DepMap n (5–500 per cohort):

```python
from scipy.stats import mannwhitneyu

def selectivity_score(gene_effect, subset_lines, background_lines=None):
    if background_lines is None:
        background_lines = gene_effect.columns.difference(subset_lines)
    sub_med = gene_effect[subset_lines].median(axis=1)
    bg_med = gene_effect[background_lines].median(axis=1)
    delta = sub_med - bg_med                # negative → cohort more dependent
    pvals = {}
    for gene in gene_effect.index:
        sub = gene_effect.loc[gene, subset_lines].dropna()
        bg  = gene_effect.loc[gene, background_lines].dropna()
        if len(sub) >= 3 and len(bg) >= 3:
            _, p = mannwhitneyu(sub, bg, alternative="less")
            pvals[gene] = p
    return pd.DataFrame({"delta": delta, "pval": pd.Series(pvals)}).dropna()
```

After computing, FDR-correct (`statsmodels.stats.multitest.multipletests`,
`method='fdr_bh'`) and join with median effect and dependency fraction.

### Permutation test (small cohorts)

When the cohort has <10 lines, the asymptotic Wilcoxon p-value is unreliable
— use a label-permutation test:

```python
def permutation_selectivity(gene_effect, gene, subset_lines, n_iter=10000, rng=None):
    rng = rng or np.random.default_rng(0)
    sub = gene_effect.loc[gene, subset_lines].dropna()
    all_lines = gene_effect.columns.dropna()
    observed = sub.median() - gene_effect.loc[gene, all_lines.difference(subset_lines)].median()
    n_sub = len(sub)
    count = 0
    for _ in range(n_iter):
        sample = rng.choice(all_lines, size=n_sub, replace=False)
        bg = all_lines.difference(sample)
        delta = gene_effect.loc[gene, sample].median() - gene_effect.loc[gene, bg].median()
        if delta <= observed:
            count += 1
    return observed, count / n_iter
```

10,000 iterations is enough for cohort sizes 3–10. Above that the asymptotic
test is fine.

---

## Effect-Size Thresholds

At DepMap's sample size, p-values get inflated. Always pair with effect-size
thresholds.

| Δ (median in cohort − median in background) | Interpretation |
|---|---|
| `< -0.5` | Strong cohort selectivity |
| `-0.5 to -0.3` | Moderate cohort selectivity |
| `-0.3 to -0.1` | Weak / suggestive |
| `> -0.1` | Not selective |

Recommend reporting hits with `delta < -0.3 AND FDR < 0.05`. Pure FDR
filtering brings in too many tiny effects.

---

## Confounders to Worry About

### Cell-line lineage compositions

If you compare AML lines vs the rest, every gene that is differentially
*expressed* between hematopoietic lines and solid lines will look
selectively essential. To control: define the background as "all
hematopoietic lines except AML" rather than "all lines."

### Expression confounding

A gene that is only expressed in 5 cell lines, all of which happen to be in
your cohort, will look selectively essential — but only because no other
line has it. Filter genes by minimum expression in both cohort and
background before testing:

```python
gene_effect_filt = gene_effect[
    (expression.loc[gene_effect.index, subset_lines].mean(axis=1) > 1.0) &
    (expression.loc[gene_effect.index, background_lines].mean(axis=1) > 0.1)
]
```

### Copy-number confounding

Chronos corrects for amplification effects, but high-amplitude amplifications
can still leak through. If a hit is a known oncogene that is amplified in the
cohort but not the background, the selectivity claim is fragile.

### Cell-line redundancy

Multiple sub-lines of the same parental (HeLa, A549, MCF7 derivatives) can
inflate counts. Most DepMap releases deduplicate, but it is worth checking
that your cohort doesn't have 5 different MCF7 derivatives counted as 5
independent observations.

---

## Reporting a Selective Essential

A defensible "this gene is selectively essential in cohort X" claim
requires:

1. The gene is **not in the pan-essential set** in the same release
2. **Δ < -0.3** between cohort median and background median
3. **FDR < 0.05** under Wilcoxon (or permutation for small cohorts)
4. The gene is **expressed in both cohort and background** at comparable
   levels (no expression confounding)
5. Result **replicates** in an orthogonal data layer where available — RNAi
   from DEMETER2, or Sanger Project Score CRISPR data
6. The cohort is **clinically meaningful** — a lineage, a subtype, a
   mutation status — not an ad-hoc cluster

A claim that satisfies these is publishable. One that satisfies only (3) is
a candidate worth investigating further but not a result.

---

## Worked Example: MLL1 in MLL-Rearranged AML

The textbook example. MLL-rearranged AML cell lines (~6 lines with a known
KMT2A fusion) are selectively dependent on the wild-type MLL1 allele
(`MLL1 (4297)` = `KMT2A`). The pattern:

- Median MLL1 gene effect in MLL-rearranged AML: ≈ `-1.5`
- Median MLL1 gene effect in non-rearranged AML: ≈ `-0.6`
- Median MLL1 gene effect across all lines: ≈ `-0.7`

So Δ vs all = `-0.8`, Δ vs lineage-matched = `-0.9`. Both pass the effect-
size threshold; lineage-matched is the cleaner comparison because it
controls for hematopoietic background.

Wilcoxon test against the lineage-matched background, n=6 vs n=~25, returns
a meaningful but not overwhelming p-value; permutation against the same
gives a similar number. Combined with the mechanistic prior (rearranged
allele requires WT for catalytic activity), this is a convincing call.

---

## Genome-Wide Selectivity Rankings

To produce a sortable list of selective essentials for a cohort:

```python
rank = (selectivity_df
        .query("not pan_essential")
        .query("delta < -0.1")
        .assign(score = lambda d: -d.delta * -np.log10(d.fdr.clip(lower=1e-20)))
        .sort_values("score", ascending=False))
```

The `score` is delta × log-significance, magnitude-weighted. Other reasonable
ranking metrics:
- Pure Δ (effect-size-only)
- Pure -log10(FDR) (significance-only)
- `cohort_dep_fraction - background_dep_fraction` (fraction-of-dependent
  lines difference)

There is no single right answer — pick what your downstream interpretation
requires.

---

## Outputs to Carry Forward

A typical selectivity output for downstream consumption:

```csv
gene,delta,pval,fdr,cohort_median,background_median,cohort_dep_frac,bg_dep_frac,is_pan_essential
KMT2A (4297),-0.92,1.2e-3,0.04,−1.51,−0.59,0.83,0.18,False
MEN1 (4221),-0.65,2.5e-3,0.06,−1.10,−0.45,0.67,0.12,False
...
```

This shape composes well with `lasso-biomarker-panel` (for predicting which
patients are likely to harbour the dependency) and with
`functional-enrichment-from-degs` (for pathway interpretation of the
selective set).
