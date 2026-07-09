# QC Checklists — Pre-Run, Mid-Run, Post-Run

Concrete checks at each stage of a genomics experiment, organized so you can paste them into a project tracker.

---

## Power analysis QC

Before running the power calculation:

- [ ] CV or dispersion estimate matches the biological system.
  - Cell lines: 0.1–0.2.
  - Inbred mice: 0.2–0.3.
  - Outbred animals or primary cells: 0.3–0.4.
  - Human cohorts: 0.4–0.5.
  - Cross-check with [cv_tissue_database.csv](cv_tissue_database.csv) or your pilot.
- [ ] Pilot data used if available — pilot dispersion beats literature CV.
- [ ] Sequencing depth is realistic for the assay (not ambitious).
- [ ] Target effect size is biologically meaningful and detectable with your assay.

After running the power calculation:

- [ ] Power ≥ 0.80 at the proposed design (or ≥ 0.90 for grants).
- [ ] Sensitivity analyzed — power tested across a range of plausible effect sizes.
- [ ] Trade-off explored (depth vs. n) at a fixed budget.
- [ ] CV source, software, version, and parameters recorded for reproducibility.

---

## Sample size QC

- [ ] At least 3 biological replicates per condition (n = 2 is invalid).
- [ ] Buffer added for failures: 10–20% extra (~`ceiling(n × 1.15)`).
- [ ] Biological replicates prioritized over technical replicates.
- [ ] Total cost within budget given the chosen depth.
- [ ] Sample availability confirmed — recruitment, breeding, archive access.
- [ ] Sequencing facility capacity confirmed for the planned timeline.

For paired designs:

- [ ] Pilot data supports a real within-pair correlation (> 0.5).
- [ ] Failure buffer increased (20–30%) — losing one timepoint costs the whole pair.

---

## Batch design QC (the most critical step)

Pre-assignment:

- [ ] "Batch" defined clearly — is it a processing day, library prep plate, sequencing run, or all of the above?
- [ ] Maximum batch size determined by lab capacity.
- [ ] All variables to balance listed in priority order.

Confounding check (must pass):

- [ ] Every batch contains samples from every experimental condition.
- [ ] `check_confounding()` returns PASS (chi-square p > 0.05).
- [ ] Visual inspection of the contingency table — no empty cells in the batch × condition matrix.

Balance checks:

- [ ] Sex balanced across batches.
- [ ] Age, site, and other priority covariates balanced.
- [ ] Sample counts per batch within ~20% of each other (last batch may be smaller).
- [ ] Reference/control samples included in each batch if multi-batch processing.

Practical:

- [ ] Randomization applied within batches; processing order documented.
- [ ] Plate edge wells avoided when possible (edge effects).
- [ ] Metadata recording plan in place — lab knows what fields to fill in.

---

## Metadata recording — what to capture per sample

Mandatory:

- [ ] Sample ID.
- [ ] Batch ID.
- [ ] Processing date (YYYY-MM-DD).
- [ ] Operator name or initials.
- [ ] Reagent lot numbers — especially library prep kit.
- [ ] Instrument or equipment ID.
- [ ] Plate / chip ID and well position.

Strongly recommended:

- [ ] Library prep batch (if separate from extraction).
- [ ] Sequencing flowcell, lane, run date.
- [ ] Protocol deviations (any change from SOP).
- [ ] Environmental notes (room temperature, humidity if RNA stability matters).

The default failure mode: skip metadata "this once", batch effects emerge in analysis, no recovery path.

---

## Multiple-testing method QC

- [ ] Method matches the scale of testing (BH-FDR for genomics, Bonferroni for < 100 tests, IHW with informative covariate).
- [ ] Threshold pre-specified before seeing results.
- [ ] FDR vs. FWER chosen based on tolerance for false positives.
- [ ] For IHW: covariate is independent of the test statistic under the null (use mean expression, not fold-change).

Post-correction:

- [ ] Number of significant results is plausible — not 0, not 60% of tests.
- [ ] QQ plot shows well-calibrated p-values (no severe inflation or deflation).
- [ ] Top hits make biological sense — known biology, consistent direction.

---

## Grant proposal checklist

- [ ] Power analysis present with explicit CV source, effect size, α, target power.
- [ ] Sample size calculation shown (or referenced supplementary).
- [ ] Budget justification ties n × depth × cost to the funding ask.
- [ ] Batch handling described (must mention balanced batches and no confounding).
- [ ] Multiple-testing method pre-specified.
- [ ] Statistical analysis plan attached (or summarized).
- [ ] Power analysis methods cited (Hart et al. 2013 for RNA-seq; Vieth et al. 2017 for scRNA-seq; Schurch et al. 2016 for replicate guidance).
- [ ] Preliminary data supports the assumed effect size and variability.

---

## Red flags

Power analysis:

- ⚠️ Power < 0.60 — study severely underpowered.
- ⚠️ Required n > 20 per group for a "large" effect — re-check assumptions, especially CV.
- ⚠️ CV < 0.1 for any heterogeneous biological sample — unrealistic.
- ⚠️ Depth > 100 M reads for RNA-seq — diminishing returns; add samples instead.

Batch design:

- 🚨 Batch confounded with condition — STOP, do not proceed.
- ⚠️ All samples of one condition in the first batch — bad randomization.
- ⚠️ No reference samples in a multi-month processing timeline.
- ⚠️ Different operators consistently assigned to different conditions.

Sample size:

- ⚠️ n = 2 — statistically invalid.
- ⚠️ No failure buffer.
- ⚠️ Technical replicates being counted toward biological n.
- ⚠️ Paired data treated as independent in power calculation.

---

## Cross-references

- Design principles — [experimental_design_best_practices.md](experimental_design_best_practices.md).
- Power workflow — [power_analysis_guidelines.md](power_analysis_guidelines.md).
- Batch design — [batch_effect_mitigation.md](batch_effect_mitigation.md).
- Method selection — [multiple_testing_guide.md](multiple_testing_guide.md).
- Troubleshooting — [troubleshooting_guide.md](troubleshooting_guide.md).
