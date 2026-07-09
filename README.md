# Fleet Skills Public Library

A curated collection of open-source agentic-science skills for literature review,
bioinformatics, statistics, reporting, and scientific quality control.

Most repository content is MIT licensed. Some skills rely on third-party tools
with their own licenses and citation requirements; see `NOTICE.md` and each
skill's references before redistributing bundled outputs or dependencies.

---

## What is a Skill?

A **skill** is a reusable, documented workflow that an AI agent can follow to complete a
well-defined scientific or analytical task. Each skill contains:

- **`SKILL.md`** - The skill itself: when to use it, inputs/outputs, step-by-step workflow,
  common issues, and suggested next steps
- **`references/`** - Supporting docs, domain guides, decision trees
- **`scripts/`** - Executable Python/R/shell code used by the skill

Skills follow the [SKILLSLIB schema](skill-creator/) for consistency and discoverability.

---

## Skills in This Library

### Literature & APIs

| Skill | Description |
|---|---|
| [arxiv-fetch](arxiv-fetch/) | Fetch and parse arXiv preprints by ID or search query |
| [biorxiv-fetch](biorxiv-fetch/) | Fetch bioRxiv/medRxiv preprints |
| [ncbi-api](ncbi-api/) | Search PubMed, fetch abstracts, retrieve PMC full-text |
| [semantic-scholar-api](semantic-scholar-api/) | Search Semantic Scholar, build citation networks |

### Transcriptomics

| Skill | Description |
|---|---|
| [bulk-rnaseq-counts-to-de-deseq2](bulk-rnaseq-counts-to-de-deseq2/) | End-to-end DESeq2 differential expression from raw counts |
| [cell-cell-communication](cell-cell-communication/) | Ligand-receptor cell communication analysis (CellChat/LIANA) |
| [scrna-trajectory-inference](scrna-trajectory-inference/) | Pseudotime / trajectory inference from scRNA-seq data |
| [scrnaseq-scanpy-core-analysis](scrnaseq-scanpy-core-analysis/) | Core scRNA-seq pipeline in Python (Scanpy) |
| [scrnaseq-seurat-core-analysis](scrnaseq-seurat-core-analysis/) | Core scRNA-seq pipeline in R (Seurat) |
| [splicecraft](splicecraft/) | Alternative splicing and isoform analysis |

### Functional Genomics

| Skill | Description |
|---|---|
| [functional-enrichment-from-degs](functional-enrichment-from-degs/) | GO/KEGG/Reactome enrichment from differential expression results |
| [gene-essentiality](gene-essentiality/) | Gene essentiality scoring from CRISPR screen data |
| [pooled-crispr-screens](pooled-crispr-screens/) | MAGeCK-based pooled CRISPR screen analysis |

### Genetics

| Skill | Description |
|---|---|
| [genetic-variant-annotation](genetic-variant-annotation/) | Annotate VCF variants with functional and clinical context |

### Statistics & Methodology

| Skill | Description |
|---|---|
| [data-analysis-best-practices](data-analysis-best-practices/) | QC, normalization, and statistical best practices reference |
| [ensemble-analyze](ensemble-analyze/) | Multi-model ensemble analysis and consensus scoring |
| [experimental-design-statistics](experimental-design-statistics/) | Power analysis, sample size, and experimental design |
| [horizon](horizon/) | Large-scale literature benchmark and retrieval pipeline |
| [lasso-biomarker-panel](lasso-biomarker-panel/) | LASSO-based biomarker panel selection from omics data |

### Clinical

| Skill | Description |
|---|---|
| [peer-review-gate](peer-review-gate/) | Structured peer-review quality gate for scientific outputs |
| [verified-science-review](verified-science-review/) | High-confidence systematic literature review with citation verification |

### Reporting

| Skill | Description |
|---|---|
| [docx-generation](docx-generation/) | Generate structured DOCX reports from agent outputs |
| [pdf-report-generation](pdf-report-generation/) | Generate PDF reports (LaTeX/Pandoc pipeline) |

### Wet Lab

| Skill | Description |
|---|---|
| [pcr-primer-design](pcr-primer-design/) | PCR primer design with Primer3 and validation |

### Meta / Skill Creation

| Skill | Description |
|---|---|
| [skill-creator](skill-creator/) | Meta-skill for designing and packaging new skills |

---

## Usage

### Standalone

Each `SKILL.md` is a self-contained document. Read it, follow the steps, run the scripts.
No agent framework is required for offline use.

### Agent Frameworks

Copy a skill directory into your framework's skill/plugin location and point the
agent at the `SKILL.md` file. Keep generated outputs outside the repository, or
in paths ignored by `.gitignore`.

---

## Contributing

Skills that meet the [SKILLSLIB schema](skill-creator/) and pass validation are welcome.

1. Fork this repo
2. Add your skill under a new directory
3. Run the validator: `python3 skill-creator/scripts/validate_skill.py --skill-dir <your-skill>/`
4. Submit a PR

---

## License

Repository-authored content is MIT licensed unless otherwise noted in individual
`SKILL.md` frontmatter. Third-party tools, datasets, APIs, and citations retain
their own licenses and terms.

---

Maintained by Fleet Skills contributors.
