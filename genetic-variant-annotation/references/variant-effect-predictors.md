# Variant effect predictors — SIFT, PolyPhen, AlphaMissense, REVEL, CADD

A practical guide to the missense-pathogenicity prediction landscape:
what each predictor does, when to prefer one over another, how to
calibrate their score ranges, and how to combine them sensibly. The goal
is to help you choose without falling into "more predictors = better"
mythology — most predictors share training data and their
non-independence undermines naive ensembling.

---

## The basic taxonomy

Missense predictors fall into three generations:

1. **First generation** (SIFT, PolyPhen-2) — sequence-conservation and
   protein-feature-based; trained on small curated datasets. Cheap to
   compute, available everywhere, but not state-of-the-art.

2. **Ensembl methods** (REVEL, MetaSVM, etc.) — combine many first-
   generation predictors as input features. Trained on ClinVar /
   HGMD. Sharper than any single predictor but inherit the training-set
   bias (ClinVar overrepresents well-studied genes).

3. **Deep-learning structure-aware predictors** (AlphaMissense, EVE,
   ESM-2-based methods) — leverage protein-language-model embeddings
   and/or AlphaFold structures. Best accuracy in recent benchmarks,
   especially on under-studied genes where ClinVar coverage is sparse.

CADD is in a separate category — it is a *general*-purpose deleteriousness
score that covers any variant type (coding + non-coding + splice +
indel), not just missense. It is the workhorse for first-pass
prioritisation across the whole VCF.

---

## SIFT (Sorting Intolerant From Tolerant)

SIFT scores a missense substitution by how conserved the affected
position is across orthologues. Score range: 0 → 1. Interpretation:

- `< 0.05` — predicted deleterious
- `≥ 0.05` — predicted tolerated

The convention is "low score = bad", which is the opposite of most other
predictors. Triple-check the direction before filtering.

Strengths: simple, fast, available for nearly every annotated coding
variant via VEP's default annotation set.

Limitations: ignores protein structure, indels, splicing; conservation
estimates are noisy in fast-evolving regions; calibrated against an old
training set.

---

## PolyPhen-2

Adds protein-structural features (Pfam domain, solvent accessibility,
PSIC scores) to the conservation signal. Score range: 0 → 1, with two
output flavours:

- `HumDiv` — trained on Mendelian disease; useful for damaging /
  benign discrimination on rare variants.
- `HumVar` — trained on common variants; better calibrated for
  prediction of pathogenicity vs benign polymorphism.

Interpretation thresholds:

| Score | HumDiv label | HumVar label |
|---|---|---|
| > 0.957 | probably damaging | probably damaging |
| > 0.453 (HumDiv) / > 0.447 (HumVar) | possibly damaging | possibly damaging |
| ≤ above | benign | benign |

Available for nearly every annotated coding variant via VEP. Use
together with SIFT — the two are often in agreement, and divergent calls
are a flag for manual review.

---

## AlphaMissense

DeepMind's 2023 predictor combining AlphaFold structural embeddings with
a protein-language-model backbone. Score range: 0 → 1. Two outputs:

- `AlphaMissense_score` — continuous score.
- `AlphaMissense_class` — categorical: `likely_pathogenic`,
  `likely_benign`, or `ambiguous`.

Thresholds (default from the publication):

- `score > 0.564` → likely pathogenic
- `score < 0.34` → likely benign
- otherwise → ambiguous

Strengths: best published accuracy on independent (non-ClinVar)
benchmarks; calibrated thresholds; covers ~71 M possible missense
substitutions across the human proteome.

Limitations: covers missense only (not stop-gain, frameshift, splice);
requires the AlphaMissense plugin and data file (~10 GB) for VEP;
licence is restricted for commercial use — verify before deploying in a
biotech setting.

**This is the recommended primary missense predictor for new pipelines.**

---

## REVEL

Ensemble method combining outputs of 13 individual predictors (SIFT,
PolyPhen-2, MutationAssessor, FATHMM, VEST, MetaLR, MetaSVM, and
others) using a random forest. Score range: 0 → 1.

Thresholds:

- `> 0.7` — strong pathogenic evidence (ACMG PP3 supporting; sometimes
  PS3 with caveats)
- `> 0.5` — moderate pathogenic evidence
- `< 0.3` — strong benign evidence

Available via the REVEL VEP plugin (data file ~50 MB compressed).

Strengths: well-calibrated, widely used in clinical literature, faster
than AlphaMissense for filtering.

Limitations: trained heavily on ClinVar — overestimates accuracy on
well-studied genes and underperforms on under-studied ones.

---

## CADD (Combined Annotation Dependent Depletion)

The "universal" deleteriousness score. Trains a logistic regression on
labelled "proxy-neutral" (simulated) vs "proxy-deleterious" (observed)
variants using hundreds of annotation features. Available for SNVs and
indels across the whole genome (coding + non-coding + splice).

Score outputs:

- `RawScore` — raw model output; not interpretable in isolation.
- `PHRED` — scaled `-10 * log10(rank/total)` so that PHRED=10 means
  top 10 %, PHRED=20 means top 1 %, etc. Use the PHRED score.

Thresholds:

| PHRED | Interpretation |
|---|---|
| ≥ 20 | Top 1 % of human variants; deleterious candidate |
| 15–20 | Top 3 % deleterious; supporting evidence |
| 10–15 | Top 10 %; weak |
| < 10 | Not flagged |

Strengths: universal coverage (any variant, any genomic region);
straightforward calibration via PHRED.

Limitations: not designed for missense in isolation — AlphaMissense or
REVEL beat CADD on missense alone. Use CADD for the first-pass
genome-wide filter; layer specialised predictors on top for missense.

---

## Splice predictors

For variants that look benign by missense prediction but might disrupt
splicing, the relevant predictors are different:

- **SpliceAI** — deep learning predictor of donor/acceptor gain or loss.
  Score range 0 → 1 per variant per nearby splice site. Thresholds:
  `> 0.8` strong, `0.5–0.8` moderate, `< 0.2` low evidence. Strongly
  recommended for any clinical pipeline.
- **MaxEntScan** — older information-theoretic donor/acceptor scoring.
  Useful as a corroborating signal but superseded by SpliceAI accuracy.
- **GeneSplicer** — neural network splice-site finder; legacy.

Splice variants must be considered any time you see a missense or
silent variant within ~20 bp of an exon boundary or within an intron.

---

## When to use which — recommended matrix

| Use case | Primary predictor | Secondary |
|---|---|---|
| Generic genome-wide first-pass | CADD PHRED > 20 | — |
| Missense in a clinical pipeline | AlphaMissense | REVEL as cross-check |
| Missense in a research pipeline (cost-sensitive) | REVEL | SIFT + PolyPhen via VEP `--everything` |
| Splice variants | SpliceAI | MaxEntScan |
| Under-studied gene | AlphaMissense | — (avoid REVEL alone) |
| Cancer somatic missense | REVEL or AlphaMissense | + COSMIC recurrence (separate axis) |
| Indels | CADD PHRED | + frame-shift consequence + LoF gene constraint |

---

## Ensembling — be careful

A naive "any predictor says pathogenic" or "all predictors agree" rule
sounds principled but is often flawed because most predictors share
training data and feature sets. REVEL already ensembles 13 underlying
predictors, so combining REVEL with SIFT and PolyPhen is mostly tautology.

A defensible approach:

1. Pick one **structure-aware deep learning** predictor (AlphaMissense or
   EVE) as the primary missense signal.
2. Pick one **ClinVar-trained ensemble** (REVEL) as the cross-check.
3. Use **CADD PHRED** for non-missense / non-splice and as a global
   sanity filter.
4. Use **SpliceAI** independently for any variant near a splice site.
5. **Disagreement is information.** When AlphaMissense flags
   pathogenic and REVEL flags benign (or vice versa), promote the
   variant to manual review — do not auto-collapse to a single answer.

---

## Calibrating against ClinVar

When deploying a new threshold for a clinical pipeline, calibrate against
ClinVar:

1. Pull all ClinVar Pathogenic + Likely Pathogenic variants in your
   genes of interest.
2. Pull ClinVar Benign + Likely Benign in the same genes.
3. Compute the predictor's true positive rate (TPR) and false positive
   rate (FPR) at candidate thresholds.
4. Pick the threshold that gives ≥ 90 % TPR at ≤ 10 % FPR for clinical
   reporting, or your group's local equivalent.

For AlphaMissense the publication's default thresholds usually need no
re-calibration in human genes. For REVEL the `> 0.7` threshold is the
clinical-grade default. For CADD the `PHRED > 20` threshold is the
genome-wide convention.

Re-calibrate annually as the underlying training data evolves; do not
freeze thresholds set in 2018.
