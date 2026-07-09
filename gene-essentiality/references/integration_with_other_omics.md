# Integrating Essentiality With Other Omics

DepMap's value compounds when joined to other layers: expression
(predict who responds), mutation (find synthetic lethals), proteomics
(close the protein-level loop), drug response (validate computationally),
patient cohorts (translate to clinic). This document walks through the
common integrations, with code patterns and the pitfalls each one carries.

---

## The Universal Join

Every DepMap release ships matched omics for the cell-line panel. Join on
`ModelID`:

```python
import pandas as pd
gene_effect = pd.read_csv("CRISPRGeneEffect.csv", index_col=0)        # genes × cell lines
expression  = pd.read_csv("OmicsExpression…TPMLogp1.csv", index_col=0) # genes × cell lines
mutations   = pd.read_csv("OmicsSomaticMutationsMatrixHotspot.csv", index_col=0)
samples     = pd.read_csv("Model.csv", index_col=0)                    # cell lines × metadata

# Find the overlap of cell lines
shared = gene_effect.columns.intersection(expression.columns).intersection(samples.index)
gene_effect = gene_effect[shared]
expression  = expression[shared]
```

Most releases align cleanly. Where they don't, prefer using the latest
`Model.csv` mapping over the per-file headers, which can lag.

---

## Essentiality vs Expression

The most common integration. "Is the gene more essential in the cell lines
that express it more highly?" In Spearman rank terms:

```python
def expression_essentiality_corr(gene_effect, expression, gene):
    # Find row IDs in each matrix
    ge_id = next(i for i in gene_effect.index if i.startswith(f"{gene} ("))
    ex_id = next(i for i in expression.index if i.startswith(f"{gene} ("))
    # Intersect cell lines
    shared = gene_effect.columns.intersection(expression.columns)
    ge = gene_effect.loc[ge_id, shared].dropna()
    ex = expression.loc[ex_id, shared.intersection(ge.index)].dropna()
    common = ge.index.intersection(ex.index)
    # INVERT the gene effect for intuitive sign
    return stats.spearmanr(ex.loc[common], -ge.loc[common])
```

Interpretation after inversion:
- **positive rho** → genes essential in lines that express them (typical
  for lineage-specific TFs, e.g. PAX8 in ovarian)
- **negative rho** → genes essential only when *not* expressed (rare; often
  represents synthetic-lethal-via-paralog relationships)
- **rho ≈ 0** → no expression dependence

For a genome-wide screen: do this per gene, sort by |rho|, then look at
known interesting cases.

### Watch out for

- **Pan-essentials**: their essentiality doesn't vary with expression because
  it doesn't vary at all. They will show low |rho| and noise; filter them
  out first.
- **Low-expression genes**: if a gene is barely detectable in most lines,
  expression is essentially binary "expressed/not", and the rho is driven
  by a small fraction of expressing lines.
- **Sample size**: Spearman p-values get tiny at n=1100. Use rho thresholds
  (|rho| > 0.3 for "interesting") rather than p-value thresholds.

---

## Essentiality vs Mutation Status

"Is the gene more essential in cell lines carrying mutation X?" Wilcoxon
between mutant and wildtype groups:

```python
def mut_dependency(gene_effect, mutation_matrix, dep_gene, mut_gene):
    ge_id = next(i for i in gene_effect.index if i.startswith(f"{dep_gene} ("))
    mut_id = next(i for i in mutation_matrix.index if i.startswith(f"{mut_gene} ("))
    shared = gene_effect.columns.intersection(mutation_matrix.columns)
    mut_status = mutation_matrix.loc[mut_id, shared].fillna(0).astype(bool)
    ge = gene_effect.loc[ge_id, shared].dropna()
    mut = mut_status.reindex(ge.index)
    sub_mut = ge[mut]
    sub_wt = ge[~mut]
    if len(sub_mut) < 3 or len(sub_wt) < 3:
        return None
    delta = sub_mut.median() - sub_wt.median()
    _, p = stats.mannwhitneyu(sub_mut, sub_wt, alternative="less")
    return {"delta": delta, "pval": p, "n_mut": len(sub_mut), "n_wt": len(sub_wt)}
```

The classic example: KRAS-mutant lines depend on KRAS, RAS-effector pathway
genes (SHOC2, PHB, …), and Replication Protein A complex.

### Watch out for

- **Low-frequency mutations**: <5 mutant lines makes the test underpowered.
  Pool across hotspot positions or accept permutation testing.
- **Lineage confounding**: KRAS mutations cluster in lung, colon, pancreas.
  Any gene that is selectively essential in those lineages will look KRAS-
  mutation-associated. Stratify by lineage or use a mixed model.
- **Co-mutations**: a "TP53-dependent" hit may actually track with the
  copy-number changes that co-occur with TP53 mutations.

---

## Essentiality vs Proteomics

CCLE proteomics (Nusinow 2020) gives protein abundance per cell line for
~12,000 proteins. Use it the same way as expression:

```python
ge = gene_effect.loc[ge_id, shared]
pr = proteomics.loc[pr_id, shared]
common = ge.dropna().index.intersection(pr.dropna().index)
rho, p = stats.spearmanr(pr.loc[common], -ge.loc[common])
```

Differences from expression:
- **Smaller n** — proteomics covers ~378 lines, much less than 1,100
- **More biological signal in many cases** — protein abundance often
  predicts function better than mRNA
- **Less complete coverage** — only ~12,000 proteins quantified, no read-
  through to genes below the detection threshold
- **Plate-batch confounders** — TMT-based proteomics can have residual
  batch effects; the published file is batch-corrected but inspect with
  caution

Use proteomics as a **confirmatory** layer for expression-essentiality
correlations, not a primary screen.

---

## Essentiality vs Copy Number

`OmicsCNGene.csv` gives gene-level log2 copy number. Two distinct
questions:

1. **Are amplified genes more essential?** Spearman with copy number,
   inverted.
2. **Are deletions creating dependencies elsewhere?** Synthetic lethality
   — find genes whose essentiality is enriched in lines where a paralog
   is deleted.

For (1), check Chronos's built-in CN correction first — it removes most of
the spurious amplification-driven essentiality, so any residual signal is
genuine. For (2), pair-up paralogs explicitly:

```python
# Example: MAGOH / MAGOHB synthetic lethality
ge = gene_effect.loc["MAGOH (4116)"]
cn_paralog = copy_number.loc["MAGOHB (55110)"]
del_lines = cn_paralog < -1.0
print(ge[del_lines].median(), ge[~del_lines].median())
```

---

## Essentiality vs Drug Response

Drug response (`OmicsDrugResponseAuc.csv` for GDSC, or PRISM tables)
correlates with essentiality of the drug's target only loosely:

- The drug may engage off-targets — its sensitivity profile reflects those
  too
- Drug exposure may be transient, while genetic knockout is constitutive
- Some "lethal" knockouts aren't drug-tractable (no small molecule yet)

For a known drug-target relationship, expect Spearman |rho| ~0.3 at best
between essentiality of the target and AUC of the drug — anything stronger
suggests the screen captured something cleaner than typical.

```python
# Are BRAF-inhibitor-responsive lines BRAF-dependent?
braf_ge = gene_effect.loc["BRAF (673)"]
plx4720_auc = drug_auc.loc["PLX4720"]   # AUC: lower = more sensitive
shared = braf_ge.index.intersection(plx4720_auc.index)
rho, p = stats.spearmanr(braf_ge.loc[shared], plx4720_auc.loc[shared])
# Expect positive correlation: more essential (more negative) ↔ more sensitive (lower AUC)
```

---

## Essentiality vs Lineage / Tumour Type

The simplest stratification: group cell lines by `OncotreeLineage` or
`OncotreeSubtype` and ask which genes show lineage-specific dependency.

```python
lineage = samples.groupby("OncotreeLineage").apply(lambda d: d.index.tolist())
for name, lines in lineage.items():
    if len(lines) < 10:
        continue
    delta = gene_effect[lines].median(axis=1) - gene_effect.median(axis=1)
    top = delta.nsmallest(10)
    print(f"{name}: top selective essentials → {list(top.index)}")
```

Watch lineage size: <10 lines is unstable. The full ranking should
incorporate a Wilcoxon p-value per gene per lineage (see
`selectivity_analysis.md`).

---

## Translating to Patient Cohorts (TCGA / GTEx)

The endgame for many DepMap analyses is "which patients harbour the
dependency I just identified?" Workflow:

1. From DepMap: gene G is selectively essential in cell-line cohort X
   (defined by mutation, expression signature, etc.)
2. Build a classifier (`lasso-biomarker-panel`) on cell-line baseline
   features that predicts essentiality of G
3. Apply the classifier to TCGA patient tumours
4. Look at the predicted-dependent tumours' clinical features and
   outcomes

Caveats:
- Cell lines lose patient context: tumour microenvironment, treatment
  history, immune contexture are absent
- TCGA expression is RNA-seq; CCLE used a different protocol — apply
  ComBat or quantile-normalize before transferring
- Survival differences in the predicted-dependent cohort do not prove the
  dependency translates to in-vivo response — they only generate hypotheses

---

## Cross-Release Integration

Pooling data from multiple DepMap releases is risky:
- Cell-line panel composition changes
- Algorithm changes (CERES → Chronos in 21Q1)
- Annotation revisions (some lines re-classified by lineage between
  releases)

Safe patterns:
- Use one release end-to-end for an analysis
- For trend analysis across releases: rank-normalize within release before
  comparing

Unsafe patterns:
- Pooling absolute scores from CERES and Chronos
- Pooling "the same cell line" across releases without checking that the
  identifier still refers to the same lineage assignment

---

## What To Carry Forward

Per integration, a defensible output table:

```csv
dep_gene,feature_gene,feature_type,n_lines,rho_or_delta,fdr,direction
KMT2A,KMT2A,expression,1100,0.41,1e-15,positive
KRAS,KRAS,mutation,1100,-0.68,2e-9,mutant_more_essential
MAGOH,MAGOHB,copy_number,1100,-0.55,3e-7,paralog_loss_dependency
```

This shape composes with downstream skills (`multi-omics-integration`,
`lasso-biomarker-panel`, `pathway-enrichment`).
