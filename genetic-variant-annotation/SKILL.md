---
name: genetic-variant-annotation
id: "skill_12f31bcbd8e2477daaef943493c52137"
display-name: "Genetic Variant Annotation"
short-description: "Annotate a VCF with functional consequences, clinical significance, population frequencies, and pathogenicity scores using Ensembl VEP or SNPEff. Use when you have called variants and need gene imp..."
category: genetics
keywords: [variant, SNP, VCF, annotation, GWAS]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Annotate a VCF with functional consequences, clinical significance, population frequencies, and pathogenicity scores using Ensembl VEP or SNPEff. Use when you have called variants and need gene impact, ClinVar links, or CADD/REVEL prioritisation."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "Ensembl VEP (Variant Effect Predictor) / SNPEff"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "McLaren W et al., Genome Biology 2016, doi:10.1186/s13059-016-0974-4"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
  - version: 1.0.0
    date: 2026-05-20
    event: initial-import
    summary: "Initial public import."
status: active
promoted: true
---

# Genetic Variant Annotation

Take a called VCF (germline, somatic, or population) and add the functional context: gene consequence, ClinVar pathogenicity, population frequencies, conservation and pathogenicity scores (CADD, REVEL, SIFT, PolyPhen), and ACMG-style classification. Two annotation engines are wrapped under a common interface — **Ensembl VEP** (most comprehensive clinical and regulatory annotation, large cache) and **SNPEff** (fast, lightweight, 38,000+ genome databases). The skill automatically picks the right tool from organism, use case, and compute budget.

Outputs include the annotated VCF, a tidy variants table, gene-level burden summaries, six standard plots, an Excel report, and a pickled analysis object that plugs into downstream pathway-enrichment, protein-structure, and literature-mining skills.

## When to Use This Skill

✅ **Use when:**
- You have **VCF files** from variant calling (GATK, bcftools, FreeBayes, DeepVariant, …)
- You want **functional consequence** (missense, stop_gained, splice_acceptor, …) per variant
- You need **clinical interpretation**: ClinVar significance, ACMG classification, gnomAD frequencies
- You want **pathogenicity prioritisation** (CADD, REVEL, SIFT, PolyPhen) for downstream filtering
- You are working with **human, mouse, or any of the 38,000+ genomes** SNPEff supports
- You want **gene-level aggregation** (variant burden per gene) for pathway analysis

❌ **Do not use when:**
- You are still **calling variants** (use GATK HaplotypeCaller, DeepVariant, bcftools call first)
- You need **structural variant annotation** (limited support; use AnnotSV or SVAFotate)
- You need **non-coding regulatory fine-mapping** (use functional annotation skills with cell-type-specific tracks)
- The variants are **CNVs or large indels** > 50 bp (annotation is approximate)
- You want **trans-eQTL** or **regulatory landscape** rather than per-variant consequence — see [gwas-to-function-twas](../gwas-to-function-twas/)

## Quick Start

ClinVar pathogenic BRCA1/BRCA2 demo (~10 variants):

```python
from pathlib import Path
import sys
sys.path.insert(0, "scripts")

from load_example_data import load_clinvar_pathogenic_sample
from validate_vcf import validate_vcf, vcf_summary_stats
from select_tool import select_annotation_tool
from run_snpeff import run_snpeff, check_snpeff_installation
from parse_snpeff_output import parse_snpeff_vcf
from annotate_genes import annotate_gene_summary
from filter_variants import filter_by_consequence, filter_by_frequency
from plot_variant_distribution import (
    plot_consequence_distribution, plot_impact_by_chromosome,
    plot_pathogenicity_scores, plot_allele_frequency, plot_gene_burden,
)
from export_results import export_all

data = load_clinvar_pathogenic_sample()
input_vcf = data["vcf_path"]

assert validate_vcf(input_vcf)["is_valid"]
print(vcf_summary_stats(input_vcf))

tool = select_annotation_tool(
    organism="human", use_case="clinical",
    resources="standard",
    annotation_priorities=["clinical_significance", "pathogenicity_predictions"],
)
# either VEP or SNPEff path is selected by select_annotation_tool()
annotated = run_snpeff(input_vcf, "annotated.snpeff.vcf.gz",
                      genome="GRCh38.105", canon=True, threads=4)
df = parse_snpeff_vcf("annotated.snpeff.vcf.gz")
genes = annotate_gene_summary(df)

high_impact = filter_by_consequence(df, consequence_levels=["HIGH", "MODERATE"])
rare = filter_by_frequency(df, max_af=0.01)

output_dir = Path("results")
output_dir.mkdir(exist_ok=True)
plot_consequence_distribution(df, output_dir / "consequence_distribution.png")
plot_impact_by_chromosome(df,   output_dir / "impact_by_chromosome.png")
plot_pathogenicity_scores(df,   output_dir / "pathogenicity_scores.png")
plot_allele_frequency(df,       output_dir / "allele_frequency.png")
plot_gene_burden(genes,         output_file=output_dir / "gene_burden.png")

export_all(df=df, gene_df=genes, original_vcf=input_vcf,
           output_dir=output_dir, tool_name=tool)
```

Expected: ~6 HIGH-impact variants, ~8 pathogenic ClinVar classifications, runtime ≤ 5 min once the cache / database is installed.

## Installation

Two annotation engines; pick one. The selection helper `select_annotation_tool()` chooses by organism + use case.

| Tool | Best for | Install | Cache size |
|---|---|---|---|
| Ensembl VEP | Human clinical, comprehensive annotation | `conda install -c bioconda ensembl-vep` | 15–20 GB / species |
| SNPEff | Non-model organisms, quick analysis | `conda install -c bioconda snpeff` | 2–5 GB / genome |

Python dependencies:

```bash
pip install "pandas>=2.0" "numpy>=1.24" pysam cyvcf2 \
    matplotlib openpyxl
```

Minimum versions: Python 3.9, VEP ≥ 110, SNPEff ≥ 5.1.

**One-line auto-install via `install_tools.py`** if conda is available:

```python
from install_tools import install_vep, install_snpeff
install_snpeff()           # ~2 min
install_vep()              # ~5 min, then download cache separately
```

**Cache / database setup:**

```bash
vep_install -a c -s homo_sapiens -y GRCh38           # VEP, ~30–60 min, ~15 GB
snpEff download GRCh38.105                            # SNPEff, ~10 min, ~3 GB
```

Detailed setup, plugin lists, and offline-cluster recipes: `references/installation_guide.md`.

**Licences:** VEP is Apache 2.0 (commercial use permitted); SNPEff is LGPL-3 (commercial use permitted).

## Inputs

**Required:**

- **VCF file** (uncompressed or `.gz`/`.bgz`), VCF 4.x specification, coordinate-sorted, indexed when bgzipped.
- **Genome assembly declaration**: `GRCh38` / `GRCh37` for human; `GRCm39` / `GRCm38` for mouse; or any of the 38,000+ SNPEff IDs (e.g. `GRCz11.105` for zebrafish).

**Optional but recommended:**

- **Reference FASTA** for full validation (reference-allele check against the VCF `REF` column).
- **BED file** to restrict annotation to specific regions (capture intervals, ROI panels).
- **Custom gene lists** for filtering or prioritisation (one HGNC symbol per line).

**VCF schema expectations:**

- `CHROM` matches the contig naming of the chosen genome (`chr1` vs `1` — convert with `bcftools annotate --rename-chrs` if mismatched).
- `REF` and `ALT` use the IUPAC ACGT alphabet (no IUPAC ambiguity codes; multi-allelic split is fine).
- `QUAL`, `FILTER`, `INFO` are present — empty values OK.

## Outputs

| File | Format | Contents |
|---|---|---|
| `analysis_object.pkl` | pickle | Complete analysis object for downstream skills (`pickle.load`) |
| `annotated.vep.vcf.gz` / `annotated.snpeff.vcf.gz` | VCF | Annotated VCF (annotations in INFO/CSQ or INFO/ANN) |
| `all_variants.csv` | CSV | Per-variant tidy table |
| `high_impact_variants.csv` | CSV | HIGH-impact subset |
| `moderate_impact_variants.csv` | CSV | HIGH + MODERATE subset |
| `rare_variants.csv` | CSV | Variants with AF < threshold (default 0.01) |
| `filtered_high_impact.vcf.gz` | VCF | Filtered VCF (preserves format for downstream tools) |
| `gene_summary.csv` | CSV | Per-gene burden, top impact, ACMG counts |
| `consequence_distribution.png/svg` | image | Variant types histogram |
| `impact_by_chromosome.png/svg` | image | Impact severity by chromosome |
| `pathogenicity_scores.png/svg` | image | CADD / REVEL / SIFT score distributions |
| `allele_frequency.png/svg` | image | Population-frequency distribution |
| `gene_burden.png/svg` | image | Top genes by variant burden |
| `summary_report.xlsx` | Excel | Multi-sheet report |

Downstream consumers of `analysis_object.pkl`: pathway-enrichment skills, protein-structure mapping, literature-mining workflows.

## Clarification Questions

⚠️ **CRITICAL: ASK THIS FIRST**

**Q1. Where is the VCF, what organism is it from, and what genome build did variant calling use?**
The most catastrophic failure mode in this skill is a build mismatch (everything annotates as `intergenic_variant`). Open the VCF header and find the `##reference=` line, or ask the user explicitly. Confirm `chr` vs no-`chr` chromosome prefix at the same time.

Then:

**Q2.** Use case — clinical diagnostics, population genetics, cancer somatic, or research prioritisation? Drives tool selection (VEP for clinical / regulatory; SNPEff for non-model organisms or quick scans) and the plugin / annotation set.

**Q3.** Which annotation priorities? Multiple-select: clinical significance (ClinVar / OMIM), pathogenicity (CADD / REVEL / SIFT / PolyPhen), population frequency (gnomAD / 1000G), regulatory impact (ENCODE), basic consequences only.

**Q4.** Compute budget — high-performance (32+ GB RAM, full VEP cache), standard (16–32 GB, either tool), or limited (< 16 GB, SNPEff)?

**Q5.** Want filtering applied during annotation, or annotate first and filter separately? The latter is more flexible; the former is faster.

**Q6.** Output requirements: variant prioritisation list, gene-level summary for pathway enrichment, Excel clinical report, or all three?

## Standard Workflow

The workflow is a 4-step chain. Use the scripts as functions — do not inline the underlying CLI invocations.

1. **Load and validate.** `validate_vcf(input_vcf)` → checks header, contig definitions, REF allele agreement (if reference FASTA provided). `vcf_summary_stats(input_vcf)` → total variants, SNV / indel counts, transitions/transversions ratio.
2. **Annotate.** `select_annotation_tool(...)` decides VEP or SNPEff; the matching `check_*` / `install_*` helpers guarantee the tool is present; `run_vep()` or `run_snpeff()` produces the annotated VCF; `parse_vep_vcf()` / `parse_snpeff_vcf()` flattens annotations into a DataFrame.
3. **Filter, summarise, and plot.** `filter_by_consequence()`, `filter_by_frequency()`, `prioritize_variants()`, `annotate_gene_summary()`, plus the five plotters in `plot_variant_distribution.py`.
4. **Export.** `export_all()` writes every CSV, every plot, the Excel workbook, the analysis pickle, and a filtered VCF.

`prioritize_variants.py` adds an ACMG-style classification (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign) when the relevant evidence fields are populated.

## When Scripts Fail

1. **Fix and retry (~90%)** — install missing Python package, run `install_tools.py` to add a missing engine, point to the correct cache directory.
2. **Modify the script (~5%)** — to add a new VEP plugin or change the SNPEff genome ID, edit the small param block at the top of `run_vep.py` / `run_snpeff.py`. Document the change.
3. **Use as reference (~4%)** — to add a custom annotation source (e.g. an in-house ClinGen export), read `parse_vep_output.py` for the CSQ-field parsing pattern.
4. **Write from scratch (~1%)** — only for niche structural-variant pipelines that this skill explicitly does not cover.

⚠️ **Do not** fall back to VEP's API mode if the local binary is missing — it is 10–100× slower and yields fewer annotations. Use `install_tools.py` instead.

## Decision Points

| Decision | Default | Override when |
|---|---|---|
| Annotation engine | VEP for human clinical, SNPEff for non-model / limited resources | `select_annotation_tool()` returns the recommendation; override with explicit `tool=` argument |
| Human assembly | GRCh38 | Legacy clinical pipelines: GRCh37, but document choice |
| VEP cache mode | Local cache (`--cache`) | Offline cluster: pre-stage cache; never fall back to API |
| SNPEff genome | GRCh38.105 | Specific Ensembl release for reproducibility |
| Filter consequence threshold | HIGH + MODERATE | Research: only HIGH. Population studies: all consequences |
| Frequency threshold | AF < 0.01 (rare) | Common-variant analyses: drop the filter |
| CADD / REVEL plugin | enabled when available | Skip if not needed — saves cache space |
| Plot DPI | 300 | Publication: 600 |
| Excel report | always | Skip for very large VCFs (> 10 M variants) — file size explodes |

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| `vep: command not found` | Tool not installed or not on PATH | `python install_tools.py --tool vep` (conda) |
| `snpEff: command not found` | SNPEff not on PATH | `python install_tools.py --tool snpeff` |
| All variants annotated as `intergenic_variant` | Genome build or contig-naming mismatch | Check `##reference=`; convert chromosome names with `bcftools annotate --rename-chrs` |
| "Cannot find genome database" | Cache or database missing | `vep_install -a c -s homo_sapiens -y GRCh38` or `snpEff download GRCh38.105` |
| No CADD / REVEL scores | Plugins not installed | Install via `vep_install -a p -g CADD,REVEL`; ensure plugin data files are downloaded |
| Memory error mid-annotation | Large VCF or default buffer too high | Add `--buffer_size 1000` to VEP; split VCF by chromosome with `bcftools view -r` |
| "Invalid VCF format" | Header truncation, missing tabs | Re-validate with `bcftools view -h`; regenerate VCF if header is corrupted |
| `Reference allele mismatch` | VCF REF disagrees with reference genome | Likely contig-naming or build mismatch; verify with `samtools faidx` |
| Output VCF much smaller than input | Annotation filter dropped non-coding variants | Use `parse_*_vcf()` without the `--no_intergenic` SNPEff flag if you want everything kept |
| `parse_vep_vcf` produces NaN columns | Plugin output absent from VCF | Verify with `bcftools view -h | grep CSQ` to check which fields are in the CSQ description |
| Pathogenicity score histogram empty | No plugins ran | Re-run annotation with `--plugin CADD,...` |
| ACMG classifier always returns VUS | Required evidence fields (PVS1 / PS1 / PM1 / PP3) missing | Confirm annotation included ClinVar + AF + CADD + splice score plugins |
| Excel report fails on huge VCF | openpyxl row limit (1,048,576) | Use the CSV outputs; skip the Excel sheet |

## Best Practices

- **Validate before annotating.** A 5-second `validate_vcf()` check catches the 80 % of failures that would otherwise eat 30 minutes of compute.
- **Match the contig naming convention.** `chr1` (Ensembl + UCSC modern) vs `1` (1000 Genomes legacy) — pick one and convert with `bcftools annotate --rename-chrs` if needed.
- **Co-install the most important plugins** the first time: CADD, REVEL, dbNSFP, gnomAD. Once installed they cost nothing per run.
- **Avoid VEP API mode.** It is documented as a fallback but is unfit for production — 10–100× slower and partial annotations.
- **Aggregate per gene** for any downstream pathway analysis. The `annotate_gene_summary()` output is the input format every enrichment tool expects.
- **Save the pickled analysis object.** It is the only convenient way to reload all annotations across processes.
- **Document the tool, version, genome build, and plugins used** in any clinical report — reviewers and regulators will ask.

## Suggested Next Steps

- **Pathway enrichment** on the gene burden table → use a functional-enrichment skill on the gene list ordered by impact or burden.
- **Protein-structure mapping** of HIGH-impact missense variants on AlphaFold or PDB structures, including PTM-aware visualisation.
- **Literature mining** of the top-prioritised genes for prior clinical evidence.
- **PRS or MR follow-up** — see [polygenic-risk-score-prs-catalog](../polygenic-risk-score-prs-catalog/) or [mendelian-randomization-twosamplemr](../mendelian-randomization-twosamplemr/) for population-level analyses.
- **Drug-target hypothesis assembly** → forward the prioritised gene list to [genetic-target-hypothesis](../genetic-target-hypothesis/) for multi-evidence prioritisation.

## Related Skills

- [gwas-to-function-twas](../gwas-to-function-twas/) — Gene-level effects from GWAS via TWAS + COLOC + MR.
- [mendelian-randomization-twosamplemr](../mendelian-randomization-twosamplemr/) — Causal inference using genetic instruments.
- [polygenic-risk-score-prs-catalog](../polygenic-risk-score-prs-catalog/) — Polygenic-score application to genotype data.
- [genetic-target-hypothesis](../genetic-target-hypothesis/) — Multi-evidence drug-target prioritisation pipeline.

## References

- **[Ensembl VEP]**: McLaren W, Gil L, Hunt SE, Riat HS, Ritchie GR, Thormann A, Flicek P, Cunningham F. "The Ensembl Variant Effect Predictor." *Genome Biology* 2016;17:122. doi:10.1186/s13059-016-0974-4 — [link](https://doi.org/10.1186/s13059-016-0974-4)
- **[SNPEff]**: Cingolani P, Platts A, Wang LL, et al. "A program for annotating and predicting the effects of single nucleotide polymorphisms — SnpEff." *Fly* 2012;6(2):80–92. doi:10.4161/fly.19695
- **[ACMG/AMP variant classification]**: Richards S, Aziz N, Bale S, et al. "Standards and guidelines for the interpretation of sequence variants." *Genetics in Medicine* 2015;17:405–424. PubMed:25741868
- **[CADD]**: Rentzsch P, Witten D, Cooper GM, Shendure J, Kircher M. "CADD: predicting the deleteriousness of variants throughout the human genome." *Nucleic Acids Research* 2019;47:D886–D894.
- **[REVEL]**: Ioannidis NM, Rothstein JH, Pejaver V, et al. "REVEL: an ensemble method for predicting the pathogenicity of missense variants." *Am J Hum Genet* 2016;99:877–885.
- **[dbNSFP]**: Liu X, Li C, Mou C, Dong Y, Tu Y. "dbNSFP v4." *Genome Medicine* 2020;12:103.
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)

### Online resources

- Ensembl VEP: <https://www.ensembl.org/info/docs/tools/vep/>
- SNPEff: <https://pcingola.github.io/SnpEff/>
- ClinVar: <https://www.ncbi.nlm.nih.gov/clinvar/>
- gnomAD: <https://gnomad.broadinstitute.org/>
- CADD: <https://cadd.gs.washington.edu/>
- REVEL: <https://sites.google.com/site/revelgenomics/>
- dbNSFP: <https://sites.google.com/site/jpopgen/dbNSFP>