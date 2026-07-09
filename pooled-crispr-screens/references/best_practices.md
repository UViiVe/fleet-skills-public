# Pooled CRISPR Screen Best Practices

Guidance on designing, executing, and analyzing pooled CRISPR screens with single-cell RNA-seq readout. Most failures of a published screen trace back to one or two design decisions made before any cells are infected — this document tries to list those decisions explicitly.

---

## Library Design

### sgRNA coverage per gene

| Guides per gene | Use case | Trade-off |
|---|---|---|
| 1 | Pilot or genome-scale survey | Cannot distinguish on-target from off-target; one bad guide kills the gene call |
| 3–5 | **Standard** | Robust on-target inference; allows guide-level concordance testing |
| 6–10 | Saturation screens, targeted high-confidence libraries | Diminishing returns past ~6; cost grows linearly |

Three is the empirical floor for treating "majority-of-guides agree" as a hit criterion. Five gives you headroom to drop a guide for poor activity without re-running.

### Negative controls

- Aim for ≥100 non-targeting sgRNAs distributed throughout the library (not clustered in plate position).
- Match the GC content of the non-targeting set to the targeting set.
- Consider including a small (5–10) panel of safe-harbor-targeting sgRNAs (e.g. AAVS1) as a second control layer — they catch problems specific to having Cas9 bound to DNA at all.

### Positive controls

Reserve 5–10 well-characterized perturbations as positive controls. Choose them across phenotype classes:
- Essential genes (will deplete in CRISPRko, show stress signature in CRISPRi)
- A strong transcription factor with a known signature
- A pathway-defining gene whose downstream signature is known in your cell type

If positive controls do not validate, the rest of the screen cannot be trusted.

### Library size budget

| Scale | sgRNA count | Approximate gene count |
|---|---|---|
| Focused screen | 1k–5k | 10–500 |
| Pathway-scale | 5k–10k | 500–1500 |
| Subgenome | 20k–50k | 4000–10000 |
| Genome-wide | 50k–100k | 17000–20000 |

---

## MOI and Doublets

Target MOI is **0.3–0.5**. This puts most cells at 0 or 1 sgRNA and keeps the sgRNA-doublet rate manageable.

| MOI | Theoretical doublet rate | Typical observed | What it means downstream |
|---|---|---|---|
| 0.3 | ~4% | 5–8% | Clean singlets, lower cell yield per perturbation |
| 0.5 | ~7% | 8–12% | Balanced |
| 0.8 | ~13% | 15–20% | Need careful doublet filtering |
| 1.0 | ~18% | >18% | Treat as combinatorial screen or expect to discard many cells |

Titrate MOI on a small pilot (0.1, 0.3, 0.5, 0.8) and measure: sgRNA capture rate, computed doublet rate from feature-barcode UMI distributions, and post-selection viability.

---

## Cell Numbers

Minimum cells per perturbation: **50–100** after QC. To plan, work backwards:

```
target_cells = n_perturbations × cells_per_perturbation × safety_factor
```

`safety_factor` accounts for 30–50% drop-out across QC, doublet removal, and library skew. A 2x factor is conservative-safe; 1.5x is the floor.

Example: 1000-gene library × 50 cells/perturbation × 2 = **100,000 cells** to capture.

---

## Replicates

| Replicate type | What it tells you | Verdict |
|---|---|---|
| Technical (sequencing splits) | Sequencing-depth variance | Cheap but largely uninformative |
| Biological from same cell pool | Sample-handling variance | Acceptable minimum |
| Independent transductions of separate cell batches | Real reproducibility | **Aim for this** |

Two independent biological replicates are the published minimum. Three is the practical minimum if you want to claim discoveries rather than just describe candidates.

---

## Culture and Selection

- Allow 3–5 days post-transduction before puromycin (or equivalent) selection
- Select for 5–7 days; confirm complete killing of un-transduced controls
- Recover 2–3 days post-selection before phenotype induction or harvest
- Maintain ≥500 cells per sgRNA throughout culture; passage at consistent density and document doublings
- Strong-survival phenotypes dominate the population over time — limit culture to 7–14 days post-transduction or harvest early for lethal perturbations

---

## 10X Capture

- Cell concentration: 700–1200 cells/µL per 10X guidelines
- Viability >90% at capture; FACS-remove dead cells with DAPI/7-AAD if not
- Singlet gating by FSC-A vs FSC-H removes physical doublets before they reach the chip
- Feature Barcoding library typically needs 12–15 PCR cycles for amplification; verify with Bioanalyzer

---

## sgRNA Assignment Quality

Expected mapping rates by technology:

| Technology | Good | Acceptable | Investigate |
|---|---|---|---|
| 10X Feature Barcoding | >60% | 40–60% | <40% |
| Direct guide capture (Replogle 2020) | >70% | 50–70% | <50% |
| Separate sgRNA-seq library | >80% | 60–80% | <60% |

Singlet calling: a cell is a singlet if **>70% of sgRNA UMIs come from one guide**, and that guide has **≥2 UMIs**. Lower thresholds (60% / 1 UMI) are acceptable only when capture is poor and data is precious.

---

## Statistical Power

For DE between perturbed and control:
- 20 cells per group is the floor for a t-test
- 50+ cells is comfortable
- 100+ cells is needed to detect modest effect sizes (|log2FC| < 0.5)

Multiple-testing correction: Benjamini–Hochberg FDR for screen-wide; Bonferroni for small targeted screens; permutation-based FDR when the test's null distribution is suspect.

---

## Hit Validation

Tier the validation so that the most expensive experiments only happen for the most defensible candidates.

1. **Computational concordance**: ≥2 sgRNAs per gene produce the same direction and similar magnitude; target gene moves in the expected direction; phenotype is distinct from unrelated perturbations.
2. **Arrayed re-test**: re-introduce individual sgRNAs into bulk culture and re-measure by qPCR / Western / bulk RNA-seq.
3. **Orthogonal modality**: switch CRISPRi to CRISPRko (or vice versa); rescue with cDNA over-expression.
4. **Mechanistic / functional**: pathway perturbation, in vivo follow-up.

---

## Modality-Specific Notes

### CRISPRi

- Expected: 70–90% knockdown at the RNA level
- dCas9-KRAB or dCas9-KRAB-MeCP2 are the standard repressors
- Target the promoter-proximal region (TSS ± 200 bp)
- Phenotypes appear quickly (24–72 h) — good for essential genes

### CRISPRa

- Expected: 2–100x upregulation
- dCas9-VPR or dCas9-SAM
- Target ~TSS -200 to +100 bp
- Activation magnitude is locus-dependent and cell-type-dependent; include known-activatable positive controls

### CRISPRko

- Indels at early exons maximize frameshift probability
- Allow 7–10 days for protein turnover
- Essential genes will deplete — interpret missing perturbations as phenotype, not failure

---

## Common Pitfalls

| Pitfall | Symptom | Fix |
|---|---|---|
| Over-passaging | Loss of slow-growing perturbations | Harvest earlier; document doublings |
| Low viral titer | Sparse library representation | Re-make virus; re-titer |
| Poor cell viability at capture | High ambient RNA; sgRNA misassignment | Dead-cell removal pre-capture |
| Too-strict QC | True perturbation phenotypes (stressed cells) lost | Two-stage QC: strict on controls, relaxed on perturbed |
| Over-corrected batches | Biological signal squashed alongside batch | Compare pre/post on positive controls |
| Skipping target validation | Pathway analysis on non-functional sgRNAs | Always validate before enrichment |

---

## Recommended Reading

- Dixit et al. (2016) *Cell* — the original Perturb-seq paper
- Adamson et al. (2016) *Cell* — multiplexed scCRISPR
- Replogle et al. (2020) *Nat Biotech* — direct guide capture
- Norman et al. (2019) *Science* — genetic interaction manifolds
- Doench (2018) *Nat Rev Genet* — experimental design walkthrough
- Sanson et al. (2018) *Nat Commun* — optimized genome-wide libraries
