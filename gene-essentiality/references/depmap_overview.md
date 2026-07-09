# DepMap Overview

The Cancer Dependency Map (DepMap) is a Broad Institute project that has
performed genome-wide CRISPR-Cas9 knockout screens across >1,100 cancer cell
lines (Project Achilles), plus a separate RNAi knockdown effort (DEMETER2).
The Wellcome Sanger Institute runs a parallel CRISPR effort (Project Score)
that DepMap integrates as of 2024.

The portal at https://depmap.org/portal/ exposes:
- Per-gene, per-cell-line **gene-effect scores** (the headline quantity)
- Per-gene, per-cell-line **dependency probabilities** (Chronos releases)
- Companion **omics layers** for each cell line: CCLE expression, mutations,
  copy number, methylation, proteomics (Nusinow 2020), metabolomics
- **Drug-response data** from PRISM (primary + secondary screens) and CTRP /
  GDSC links

This document covers what is in a release, how to fetch it programmatically,
and what release-to-release changes you should care about.

---

## Anatomy of a Release

Each release is tagged by year + quarter: `21Q1`, `22Q2`, `24Q4`, etc. The key
files (current naming, post-2023):

| File | Contents | Shape | Notes |
|---|---|---|---|
| `CRISPRGeneEffect.csv` | Gene effect, Chronos-scored | genes × cell lines | Row index: `SYMBOL (entrez_id)` |
| `CRISPRGeneDependency.csv` | Posterior probability of dependency | genes × cell lines | Chronos-only |
| `Model.csv` | Cell-line metadata | cell lines × features | OncotreeLineage, OncotreeSubtype, mutations summary, lineage |
| `OmicsExpressionProteinCodingGenesTPMLogp1.csv` | Bulk RNA-seq, log2(TPM+1) | genes × cell lines | Protein-coding only |
| `OmicsSomaticMutationsMatrixHotspot.csv` | Hotspot mutation calls | genes × cell lines, boolean | Use `OmicsSomaticMutations.csv` for full calls |
| `OmicsCNGene.csv` | Copy-number ratios | genes × cell lines | log2 |
| `Achilles_logfold_change.csv` | Raw guide-level fold-changes | guides × cell lines | Input to scoring |

Older releases (pre-22Q1) use ACH-XXXXXX style cell-line IDs and have a
slightly different file layout. The current schema standardised on `ModelID`
and an Oncotree-based lineage system.

---

## Programmatic Access

The portal exposes a download endpoint per release. Pattern:

```
https://depmap.org/portal/api/download/files/<RELEASE>/<FILENAME>
```

The official path also fronts a Figshare mirror; both are stable URLs.

```python
import pandas as pd
url = "https://depmap.org/portal/api/download/files/24Q2/CRISPRGeneEffect.csv"
ge = pd.read_csv(url, index_col=0)
```

The `load_depmap.py` helper in this skill caches under `~/.cache/depmap/<release>/`
and converts to parquet on first use — subsequent reads are ~100× faster.

---

## Quantities in the Gene-Effect Matrix

Cell lines are columns; genes are rows. The cell value is a real number.

- **0 ≈ no fitness effect**
- **-1 ≈ median effect of a known pan-essential** (RPL/RPS, POLR2A, etc.)
- **negative = essential** (knockout reduces fitness)
- **positive = pro-growth / suppressor** (knockout improves fitness; rare in
  cancer lines)

This sign convention is the single most common source of analysis errors.
Treat the data as essentiality requires negation; see
`sign_conventions.md` for the full discussion.

---

## DepMap Releases You Will Encounter

| Release window | Scoring | Notes |
|---|---|---|
| ≤ 20Q4 | CERES | Copy-number-correction model from Meyers 2017 |
| 21Q1 → 22Q4 | **Chronos** | Population-dynamics model; tighter dynamic range |
| 23Q1 → 24Q4 | Chronos | Schema rewritten (`Model.csv`, OmicsExpressionProteinCodingGenes…) |
| Sanger Project Score | CERES / BAGEL2 / Chronos | Separate panel; integrated into DepMap from 24Q1 |

**Recommendation**: use the latest release unless reproducing prior work. For
the latter, match release tag exactly and document the choice — DepMap
revises cell-line annotations and the panel on every release.

---

## What the Algorithm Does NOT Account For

Gene-effect scores are *cancer-cell-line* dependencies under standard 2D
culture. They will not capture:
- Tumour-microenvironment dependencies (stroma, immune, vasculature)
- Drug-resistance phenotypes that arise during treatment
- In-vivo growth dependencies that are absent in vitro
- Patient-derived organoid effects (only a small subset of DepMap lines are
  organoids)
- Synthetic-lethal combinations (DepMap is a single-knockout screen)

These caveats matter for therapeutic translation. A gene that is broadly
essential in 2D culture is often essential because of generic cell-fitness
machinery (translation, mitosis, splicing) and is not actionable; a gene that
is selectively essential in a context — especially one defined by a
clinically meaningful biomarker — is.

---

## Cell-Line Panel

The panel is ~1,100+ lines spanning the full range of solid and liquid
tumours. Coverage is unequal:
- Common cancers (lung, breast, colon) have hundreds of lines each
- Rare cancers may have 1–3 lines, making statistical claims unreliable
- Pediatric cancers are under-represented relative to adult cancers
- Female patients are under-represented in some lineages

Check `samples.OncotreeLineage.value_counts()` before claiming "lineage-
specific" anything. Anything below ~10 lines for the cohort calls for a
permutation test rather than asymptotic statistics.

---

## Integrated Data Layers

Every release ships matched omics for the same set of cell lines:

| Layer | File | Source | Notes |
|---|---|---|---|
| RNA-seq | `OmicsExpression*.csv` | CCLE | log2(TPM+1) |
| Mutations | `OmicsSomaticMutations*.csv` | CCLE | WES; hotspot subset for matrix file |
| Copy number | `OmicsCNGene.csv` | CCLE | gene-level log2 |
| Methylation | `OmicsAbsoluteCNGene.csv` (legacy) | Sanger | Probe-level still on legacy paths |
| Proteomics | `ProteinMatrixAveragedNormalizedAbundance.csv` | Nusinow 2020 | TMT, log2 ratios |
| Metabolomics | `Metabolomics.csv` | Li 2019 | Z-scored intensities |

Use the same `ModelID` join key across files.

---

## Drug Response Layers

- **PRISM** (Broad): two phases — primary (4,500 compounds × ~480 lines) and
  secondary / repurposing (1,448 compounds × 578 lines)
- **CTRP v2** (Broad): 481 compounds × 887 lines, 16-point dose curves
- **GDSC** (Sanger, integrated as `OmicsDrugResponseAuc.csv`)

Each ships AUC and IC50 per (compound, cell line). Drug response correlates
with essentiality of the drug target only if the cell line is sensitive AND
the drug's mechanism actually engages that target — frequently neither is
true.

---

## Practical Notes

- **First download is slow** (~500 MB), subsequent reads are cached
- **Always read the release notes** at https://depmap.org/portal/data_page/?tab=allData
  before pooling data from multiple releases
- **Use ModelID as the join key** in any modern release; `cell_line_name` is
  ambiguous (multiple lines share names across labs)
- **Gene identifiers are SYMBOL (entrez_id)** in the latest releases; index
  on the composite, not just the symbol, to avoid silent collisions
