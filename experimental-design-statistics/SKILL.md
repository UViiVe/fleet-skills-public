---
name: experimental-design-statistics
id: "skill_0cad451229b04b71b97b33a9bc0638c2"
display-name: "Experimental Design & Statistics"
short-description: "Plan genomics experiments with power analysis, sample size estimation, batch-balanced layouts, and multiple-testing strategy. Use when designing a new RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methyl..."
category: methodology
keywords: [experimental design, statistics, power, sample size, t-test]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Plan genomics experiments with power analysis, sample size estimation, batch-balanced layouts, and multiple-testing strategy. Use when designing a new RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, or proteomics study before any samples are processed."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "Methodology guide (power analysis, sample size, multiple testing — BH-FDR, Cohen's d)"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Benjamini Y & Hochberg Y, J Royal Statistical Society B 1995, doi:10.1111/j.2517-6161.1995.tb02031.x; Cohen J, Statistical Power Analysis for the Behavioral Sciences, 2nd ed, 1988"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# Experimental Design and Statistical Planning

A pre-experiment workflow that turns a fuzzy plan ("we want to do RNA-seq on liver from control vs. drug-treated mice") into a defensible design with a chosen sample size, a balanced batch layout that survives a confounding check, and a written analysis plan that locks in the multiple-testing method before any data is collected.

The job is unglamorous but high-leverage: most failed genomics experiments fail at design, not at analysis. Too few replicates, a batch that's secretly the treatment, or an FDR threshold picked after the volcano plot — these are the patterns that turn a $40k experiment into a methods-section footnote.

---

## When to Use This Skill

✅ **Use this skill when:**
- Planning a new genomics experiment (any assay: bulk RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, proteomics).
- Justifying a sample size for a grant, IRB, or IACUC application.
- Translating a budget cap into a defensible design (depth vs. n trade-off).
- Generating a balanced batch layout when samples must be processed in groups.
- Locking in the multiple-testing correction method before seeing results.
- Pre-registering a study or writing a statistical analysis plan.

❌ **Do not use this skill for:**
- Running the analysis after samples are sequenced — use the assay-specific DE skills.
- Single fixed-n two-sample comparisons where standard t-test power suffices.
- Computational batch correction of already-collected data — design prevention is upstream of `ComBat`/`SVA`/`RUVSeq`.
- Post-hoc power analysis to "explain" a null result. This is statistically meaningless and a red flag in any review.

**Core principle:** Power calculations belong *before* the experiment. After data is collected, the only relevant numbers are effect sizes and confidence intervals.

---

## Quick Start

Typical end-to-end run for a 2-condition RNA-seq design with pilot data:

```r
source("scripts/load_example_data.R")
source("scripts/power_rnaseq.R")
source("scripts/sample_size_de.R")
source("scripts/batch_assignment.R")
source("scripts/batch_validation.R")
source("scripts/plot_power_curves.R")
source("scripts/export_design.R")

pilot   <- load_example_data()
req_n   <- samplesize_from_pilot(pilot$dds, fold_change = 1.5, power = 0.8)
batch   <- assign_samples_to_batches(pilot$metadata, batch_size = 8,
                                     balance_vars = c("condition"))
check_confounding(batch, "condition")
plot_power_vs_samplesize(cv = pilot$cv$median, fold_changes = c(1.5, 2, 3),
                         output_file = "out/power_vs_n")
export_complete_design(batch, design_params = list(
  assay = "RNA-seq", conditions = c("control", "treated"),
  n_per_group = req_n$required_n_per_group, power = 0.8,
  alpha = 0.05, effect_size = 1.5, multiple_testing = "BH-FDR"
), output_dir = "out")
```

**Runtime:** 30 s – 3 min on a laptop. `powsimR` simulations (if invoked) can take 10–30 min.

---

## Installation

| Package | Version | License | Commercial use | Install |
|---|---|---|---|---|
| RNASeqPower | ≥1.36 | GPL-2+ | yes | `BiocManager::install("RNASeqPower")` |
| ssizeRNA | ≥1.4 | GPL-2+ | yes | `BiocManager::install("ssizeRNA")` |
| DESeq2 | ≥1.40 | LGPL | yes | `BiocManager::install("DESeq2")` |
| powsimR | ≥1.2.3 | GPL-3 | yes | `BiocManager::install("powsimR")` |
| OSAT | ≥1.38 | Artistic-2.0 | yes | `BiocManager::install("OSAT")` |
| IHW | ≥1.28 | Artistic-2.0 | yes | `BiocManager::install("IHW")` |
| qvalue | ≥2.32 | LGPL | yes | `BiocManager::install("qvalue")` |
| sva | ≥3.48 | Artistic-2.0 | yes | `BiocManager::install("sva")` |
| pasilla | ≥1.18 | Artistic-2.0 | yes | `BiocManager::install("pasilla")` (example data) |
| ggplot2 | ≥3.4 | MIT | yes | `install.packages("ggplot2")` |
| ggprism | ≥1.0.4 | GPL-3 | yes | `install.packages("ggprism")` |
| jsonlite | ≥1.8 | MIT | yes | `install.packages("jsonlite")` |

Quick one-shot install:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("RNASeqPower","ssizeRNA","DESeq2","powsimR","OSAT","IHW","qvalue","sva","pasilla"))
install.packages(c("ggplot2","ggprism","jsonlite"))
```

Full install notes including OS-specific compilation issues are in [references/software_requirements.md](references/software_requirements.md).

---

## Inputs

**Required (gather from clarification questions):**
- Assay type — bulk RNA-seq, scRNA-seq, ATAC-seq, ChIP-seq, methylation, proteomics, or other.
- Number of conditions and planned n per group (or "unknown — calculate from power").
- Target effect size (fold-change to detect) and target power (0.80 standard, 0.90 for grants).
- α (typically 0.05) and intended multiple-testing method (or "need recommendation").

**Optional but high-value:**
- Pilot data — a `DESeqDataSet` RDS, a count matrix CSV, or a `SingleCellExperiment` RDS. Pilot dispersion estimates beat literature CV estimates every time.
- Covariates to balance across batches (sex, age, collection site, genetic background).
- Budget constraints (max samples, max sequencing depth, batch processing capacity).
- Existing metadata table with sample IDs and condition labels for batch assignment.

If no pilot data is available, the skill falls back to literature CV values in [references/cv_tissue_database.csv](references/cv_tissue_database.csv).

---

## Outputs

All written to a single output directory (default `design_results/`):

| File | Format | Downstream consumer |
|---|---|---|
| `power_analysis_results.csv` | CSV | Methods section, grant proposal |
| `sample_size_recommendation.txt` | text | Quick summary for stakeholders |
| `power_vs_n.png` / `.svg` | image | Figure for grant / pre-registration |
| `batch_layout_for_lab.csv` | CSV | Lab bench — sample-to-batch assignment |
| `batch_design_validation.txt` | text | Confounding test result, balance summary |
| `batch_design.png` / `.svg` | image | Visual layout check |
| `statistical_analysis_plan.md` | Markdown | Pre-registration, reviewer-facing plan |
| `lab_protocol_checklist.md` | Markdown | Bench checklist, metadata recording |
| `design_parameters.json` | JSON | Machine-readable design record |
| `batch_design.rds` | RDS | Load with `readRDS()` for downstream batch correction |
| `design_parameters.rds` | RDS | Load for analysis-stage validation |

The `.rds` files are the load-bearing artifacts: downstream DE skills read them to verify the analysis matches the pre-specified plan.

---

## Clarification Questions

### 1. ⚠️ CRITICAL: ASK THIS FIRST — Pilot data or literature estimates?

Do you have pilot data we can use for variability estimates?
- **Yes — RDS file (DESeqDataSet)** → preferred path, gives realistic dispersion.
- **Yes — count matrix (CSV/TSV)** → we'll build a DESeq2 object and extract dispersion.
- **No** → we'll use tissue-specific CV from [references/cv_tissue_database.csv](references/cv_tissue_database.csv). Tell me the tissue and organism.

This question is first because the answer changes which calculation path runs. A pilot-based estimate can shift the recommended n by 2–4× compared to a generic literature CV.

### 2. Assay type and experimental structure

- Assay: bulk RNA-seq / scRNA-seq / ATAC-seq / ChIP-seq / methylation / proteomics / other.
- Conditions: 2 (case-control), 3+ (multi-group), or factorial?
- Samples: independent, paired (before/after), or repeated measures?
- Covariates to balance: sex, age, site, genetic background, anything else?

### 3. Effect size and power targets

- Target fold-change: large (≥2×), moderate (1.5–2×), small (1.2–1.5×).
- Target power: 0.80 (standard), 0.90 (grants, clinical), or custom.
- α: 0.05 (standard), 0.01 (stringent), or custom.
- Multiple testing: BH-FDR (default), IHW (more powerful with covariate), Bonferroni (small candidate panels), or need guidance.

### 4. Practical constraints

- Budget: maximum samples or maximum sequencing reads available.
- Sample availability — limited? Include a 10–20% failure buffer?
- Batch structure: how many samples can the lab process simultaneously?
- Sequencing depth target: RNA-seq 15–30M reads, ATAC-seq 25–50M, scRNA-seq 50–100K reads/cell.

### 5. Primary objective

Power-only / sample-size-only / batch-design-only / multiple-testing guidance / complete design / budget optimization.

Full clarification rubric: [references/experimental_design_best_practices.md](references/experimental_design_best_practices.md).

---

## Standard Workflow

The scripts are designed to be sourced and called as functions, not copy-pasted. Use them rather than writing inline power or batch-assignment code — the OSAT optimization and SVG-fallback plotting are not trivial to reproduce.

### Step 1 — Load parameters

With pilot data:

```r
source("scripts/load_example_data.R")
pilot <- load_example_data()         # example data, or replace with your RDS
```

Without pilot data:

```r
cv_db <- load_cv_database()
# choose tissue row; cv_db$CV holds the median estimate
```

**Verify:** `"✓ Example pilot data loaded successfully!"` or the CV table prints.

### Step 2 — Calculate the design

Power analysis at a planned n:

```r
source("scripts/power_rnaseq.R")
power <- calc_power_rnaseq(depth = 20, n_per_group = 6,
                           cv = pilot$cv$median, fold_change = 2, alpha = 0.05)
```

Required n from pilot dispersion:

```r
source("scripts/sample_size_de.R")
req_n <- samplesize_from_pilot(pilot$dds, fold_change = 1.5, power = 0.8, fdr = 0.05)
```

Balanced batch layout (OSAT-optimized):

```r
source("scripts/batch_assignment.R")
batch <- assign_samples_to_batches(pilot$metadata, batch_size = 8,
                                   balance_vars = c("condition", "sex"))
```

For scRNA-seq sample sizes use `scripts/sample_size_scrna.R`; for ATAC-seq use `scripts/power_atacseq.R`.

**The cardinal rule:** batch must not be confounded with condition. The next step verifies this.

### Step 3 — Validate and visualize

```r
source("scripts/batch_validation.R")
check_confounding(batch, "condition")          # MUST return "PASS"
check_balance(batch, c("sex", "age_group"))    # chi-square p > 0.05 ideal

source("scripts/plot_power_curves.R")
plot_power_vs_samplesize(cv = pilot$cv$median, fold_changes = c(1.5, 2, 3),
                         output_file = "design_results/power_vs_n")
visualize_batch_design(batch, condition_var = "condition",
                       output_file = "design_results/batch_design")
```

Both PNG and SVG are written. If `svglite` is missing the SVG falls back gracefully to base R's `svg()` device — do **not** try to install `svglite` to "fix" this.

If `check_confounding()` returns FAIL, regenerate the batch layout with stricter constraints. Do not proceed.

### Step 4 — Multiple-testing strategy (lock in *now*)

```r
source("scripts/multiple_testing.R")
recommend_method(n_tests = 20000, effect_heterogeneity = "moderate", sample_size = 6)
```

Record the recommended method in the statistical analysis plan. **Switching methods after seeing results is a form of p-hacking.**

### Step 5 — Export everything

```r
source("scripts/export_design.R")
params <- list(
  assay = "RNA-seq", conditions = c("control", "treated"),
  n_per_group = 6, power = 0.85, alpha = 0.05,
  effect_size = 2, multiple_testing = "BH-FDR"
)
export_complete_design(batch, params, output_dir = "design_results")
```

**Verify:** `"=== Export Complete ==="`. The exported `.rds` files are required for downstream batch correction and analysis-time validation.

---

## When Scripts Fail

Follow this escalation hierarchy. Most failures are missing packages.

1. **Fix and retry (90% of cases).** Install the missing Bioconductor or CRAN package the error message names. Re-source the script and re-run.
2. **Modify the script (5%).** If a parameter range needs to change (e.g., the default `n_range = 2:15` is too narrow for your CV), edit the script file and document the change in your analysis plan.
3. **Use the script as reference (4%).** If a dependency truly cannot be installed (institution firewall blocking Bioconductor), read the script, adapt the algorithm, and cite the source script in your code.
4. **Write from scratch (1%).** Only if the above fail and you can articulate why in writing. Document in the analysis plan.

Never skip directly to inline code without trying the script first. The scripts encode design decisions (OSAT optimization, SVG fallback, NA handling in p-value adjustment) that are easy to get wrong.

---

## Decision Points

| Decision | Default | When to override |
|---|---|---|
| Pilot CV vs. literature CV | Pilot if available | If pilot n < 3 per group, blend with literature (use higher value) |
| Sample size vs. depth | Add samples first | Only add depth below 10M reads (RNA-seq) or 25M (ATAC-seq) |
| Power target | 0.80 | 0.90 for grants, clinical, expensive validation |
| FDR threshold | q < 0.05 | q < 0.01 for high-stakes, q < 0.10 for exploratory pilots |
| Multiple-testing method | BH-FDR | IHW if mean expression varies widely; Bonferroni if < 100 tests |
| Batch size | Lab capacity | Multiple of n_conditions × 2 keeps balance clean |
| Buffer for failures | 10% | 20–30% for paired designs, precious samples, or long timelines |

Full rationale and references: [references/experimental_design_best_practices.md](references/experimental_design_best_practices.md) and [references/power_analysis_guidelines.md](references/power_analysis_guidelines.md).

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| Power < 0.80 at max feasible n | Effect size too small or CV too high | Revise target fold-change upward; verify CV; consider paired design; or accept lower power and document |
| `check_confounding()` returns FAIL | Batch correlates with condition | Regenerate with `assign_samples_to_batches`, ensure `condition` is in `balance_vars`, verify batch_size allows balance |
| Required n exceeds sample availability | Pilot dispersion higher than expected | Investigate sources of variation (outliers, hidden batches); use paired design; revise effect-size target |
| `OSAT::osat` fails or hangs | Too many balance variables for batch size | Reduce balance variables to condition + 1–2 highest-priority covariates; increase batch size |
| CV estimate from pilot wildly variable | Low-count genes inflate per-gene CV | Filter genes with mean count < 10 before CV calculation; use median CV, not mean |
| `ssizeRNA_vary` returns n = 2 | Very large effect or near-zero dispersion | n = 2 is statistically invalid; floor at n = 3, prefer n ≥ 4 |
| SVG export warning printed | `svglite` not installed | Ignored — the script falls back to base R `svg()`. Do not install `svglite`. PNG always saved |
| IHW returns same as BH | Covariate uninformative under null | Try mean expression instead of variance; if no informative covariate, use BH |
| Multi-group power feels too low | Pairwise comparisons share α budget | Increase n by √k for k groups, or focus power on the primary contrast |
| `powsimR` install takes 30 min | Many heavyweight dependencies | This is normal. Run install once, then cache |
| Post-hoc power question | User wants to "salvage" a null result | Decline. Post-hoc power is not informative. Suggest sensitivity analysis instead |
| Batch processing delayed by weeks | Operational slip | Record dates as a covariate; add cross-batch reference samples for remaining batches |

Deeper troubleshooting: [references/troubleshooting_guide.md](references/troubleshooting_guide.md).

---

## Best Practices

- **Biological replicates beat technical replicates.** Modern sequencing technical noise is small. Don't sequence the same library twice and call it n = 2.
- **Pilot studies pay for themselves.** A 6-sample pilot can rescue a 40-sample main study by exposing high dispersion before the budget is spent.
- **Document the design before sequencing, not after.** The `statistical_analysis_plan.md` artifact is the pre-registration. Lock it in.
- **Randomize, don't alphabetize.** Sample order, plate positions, processing day — all randomized within batches.
- **Record reagent lots.** Lot-to-lot variation is a real batch effect. Capture it in metadata so it can be modeled later.
- **Match the multiple-testing method to the experiment scale, not your preference.** GWAS: Bonferroni at 5×10⁻⁸. RNA-seq: BH-FDR. Candidate panel of 20 genes: Bonferroni.

---

## Suggested Next Steps

After completing the design and collecting samples, downstream:

1. **bulk-rnaseq-counts-to-de-deseq2** — DE analysis using the pre-specified design formula and FDR method.
2. **bulk-rnaseq-differential-expression** — broader DE comparisons including edgeR and limma-voom.
3. **functional-enrichment-from-degs** — pathway analysis once DE genes are identified, using the analysis plan's significance threshold.
4. **pathway-enrichment** — ORA or GSEA with the pre-specified background gene set.
5. **bulk-omics-clustering** — exploratory clustering / PCA to verify batch correction worked.
6. **data-analysis-best-practices** — data validation checks before running analysis.

Each downstream skill can load `design_parameters.rds` to verify the running analysis matches the plan.

---

## Related Skills

- **Upstream:** none — this is typically the first artifact in a project.
- **Companion:** [data-analysis-best-practices](../data-analysis-best-practices/SKILL.md) for the QC/validation steps that follow data delivery.
- **Downstream DE:** the assay-specific DE skills listed in Suggested Next Steps.
- **Computational batch correction (if design prevention falls short):** look for ComBat/SVA/RUVSeq workflows.

---

## References

**Bundled with this skill:**
- [references/experimental_design_best_practices.md](references/experimental_design_best_practices.md) — design principles, pitfalls, budget optimization.
- [references/power_analysis_guidelines.md](references/power_analysis_guidelines.md) — power workflow, pilot vs. literature, assay-specific guidance.
- [references/batch_effect_mitigation.md](references/batch_effect_mitigation.md) — preventing confounding, balanced layouts.
- [references/multiple_testing_guide.md](references/multiple_testing_guide.md) — choosing BH/Bonferroni/IHW/q-value.
- [references/qc_guidelines.md](references/qc_guidelines.md) — pre-run and post-run QC checklists.
- [references/troubleshooting_guide.md](references/troubleshooting_guide.md) — common failure modes and fixes.
- [references/software_requirements.md](references/software_requirements.md) — install instructions, license summary.
- [references/cv_tissue_database.csv](references/cv_tissue_database.csv) — tissue/assay-specific CV estimates from published data.

**Scripts:**
- `scripts/load_example_data.R` — pilot data loader, CV database loader.
- `scripts/power_rnaseq.R` — RNA-seq power and sample size via RNASeqPower.
- `scripts/power_atacseq.R` — ATAC-seq power via ssizeRNA.
- `scripts/power_pilot_based.R` — pilot-dispersion-driven power, PROPER wrapper.
- `scripts/sample_size_de.R` — DE sample size with FDR control (ssizeRNA).
- `scripts/sample_size_scrna.R` — scRNA-seq sample size and cells-per-sample heuristics.
- `scripts/batch_assignment.R` — OSAT-optimized batch layouts, plate position assignment.
- `scripts/batch_validation.R` — confounding test, balance test, visual layout.
- `scripts/multiple_testing.R` — BH/Bonferroni/qvalue/IHW with method recommender.
- `scripts/plot_power_curves.R` — power vs n, power vs FC, sample-size heatmap (PNG+SVG).
- `scripts/export_design.R` — bundled export of analysis plan, lab protocol, RDS, JSON.

**Primary citations:**
- Hart SN, Therneau TM, Zhang Y, Poland GA, Kocher J-P. *Calculating sample size estimates for RNA sequencing data.* J Comput Biol 2013;20(12):970–978.
- Schurch NJ et al. *How many biological replicates are needed in an RNA-seq experiment and which differential expression tool should you use?* RNA 2016;22(6):839–851.
- Leek JT, Scharpf RB, Bravo HC, et al. *Tackling the widespread and critical impact of batch effects in high-throughput data.* Nat Rev Genet 2010;11(10):733–739.
- Benjamini Y, Hochberg Y. *Controlling the false discovery rate.* J R Stat Soc B 1995;57(1):289–300.
- Ignatiadis N, Klaus B, Zaugg JB, Huber W. *Data-driven hypothesis weighting increases detection power.* Nat Methods 2016;13(7):577–580.
- Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biol 2014;15(12):550.
- Yasrebi H. *OSAT: A tool for sample-to-batch allocations in genomics experiments.* BMC Genomics 2012;13:689.
