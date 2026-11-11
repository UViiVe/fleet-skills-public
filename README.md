# Fleet Skills Public Library

A curated collection of **36** open-source agentic-science skills for literature review, bioinformatics, statistics, reporting, and scientific quality control.

Most repository content is MIT licensed. Some skills rely on third-party tools with their own licenses and citation requirements; see `NOTICE.md` and each skill's references before redistributing bundled outputs or dependencies.

---

## What is a Skill?

A **skill** is a reusable, documented workflow that an AI agent can follow to complete a well-defined scientific or analytical task. Each skill contains:

- **`SKILL.md`** - The skill itself: when to use it, inputs/outputs, step-by-step workflow, common issues, and suggested next steps
- **`references/`** - Supporting docs, domain guides, decision trees
- **`scripts/`** - Executable Python/R/shell code used by the skill

Skills follow the [SKILLSLIB schema](skill-creator/) for consistency and discoverability.

---

## Skills in This Library

### Regulatory Genomics

| Skill | Description |
|---|---|
| [alphagenome2](alphagenome2/) | Bring-your-own-key AlphaGenome 2 noncommercial research workflow: condition metadata, one bounded reference prediction, and explicit privacy/interpretation gates. |

### Methodology

| Skill | Description |
|---|---|
| [abstract-lookup](abstract-lookup/) | Maximizes abstract coverage for research papers using a waterfall approach across multiple APIs (PubMed → Europe PMC → Semantic Scholar → Crossref). Use when building... |
| [arxiv-fetch](arxiv-fetch/) | Library card skill - programmatic access to arXiv preprints. Both metadata (Atom API) and PDFs are fully open: no Cloudflare, no auth, no browser required. Direct HTTP... |
| [arxiv-search](arxiv-search/) | Search arXiv for academic papers using the official arXiv API. Use when: (1) looking up recent preprints or papers by topic, (2) fetching paper metadata (title, author... |
| [biorxiv-fetch](biorxiv-fetch/) | Library card skill - programmatic access to bioRxiv/medRxiv preprint METADATA (title, authors, abstract, DOI, date) via the official REST API, AND a VERIFIED browser-a... |
| [biorxiv-search](biorxiv-search/) | Search bioRxiv and medRxiv preprints via the official biorxiv.org API. No API key required. Use when looking for recent life sciences preprints on a topic — genomics,... |
| [chemrxiv-fetch](chemrxiv-fetch/) | Library card skill — programmatic access to ChemRxiv chemistry preprints. Metadata via Crossref API (fully open). PDFs are behind Cloudflare: headless browser clears t... |
| [data-analysis-best-practices](data-analysis-best-practices/) | Pre-analysis validation, missing-data handling, multiple-testing correction, batch-effect checks, reproducibility, and seed/version pinning for any data analysis. Use... |
| [docling-parse](docling-parse/) | Parse PDFs into structured markdown, JSON, or HTML using IBM Docling — handles tables, figures, multi-column layouts, and scientific papers without any API key. Use wh... |
| [ensemble-analyze](ensemble-analyze/) | Run a scientific question or dataset analysis through a configurable panel of frontier LLMs (Gemini, GPT-5.5, Claude Opus, Qwen local) in parallel, synthesize a consen... |
| [europepmc-search](europepmc-search/) | Search Europe PMC for preprints and published papers across 30+ sources in one API call. No API key required. Covers bioRxiv, medRxiv, Research Square, Preprints.org,... |
| [experimental-design-statistics](experimental-design-statistics/) | Plan genomics experiments with power analysis, sample size estimation, batch-balanced layouts, and multiple-testing strategy. Use when designing a new RNA-seq, scRNA-s... |
| [ncbi-api](ncbi-api/) | Search PubMed, fetch abstracts, retrieve full-text links, and download open-access papers via NCBI E-utilities and PMC APIs. Use when searching biomedical literature b... |
| [openalex-search](openalex-search/) | Search OpenAlex for scholarly works including preprints across all major servers and repositories. No API key required. 450M+ works indexed across arXiv, bioRxiv, medR... |
| [paperclip](paperclip/) | Search and read from Paperclip — a virtual filesystem of 8M+ full-text biomedical papers (PubMed Central, bioRxiv, medRxiv, arXiv). Use for literature search, paper re... |
| [peer-review-gate](peer-review-gate/) | Require independent review before scaled, expensive, irreversible, or high-risk scientific work. Use when a second-review checkpoint should happen before execution. |
| [pipette-skillgraph](pipette-skillgraph/) | Query the Pipette.bio SkillGraph MCP server to discover bioinformatics workflows, pipeline transitions, and tool documentation across 78 skills and 483 pipeline edges.... |
| [semantic-scholar-api](semantic-scholar-api/) | Search academic papers, fetch citations/references, get author profiles, and retrieve paper details via the Semantic Scholar Academic Graph API. Use when searching cro... |
| [skill-creator](skill-creator/) | Draft, validate, and package a reusable agentic-science skill with public-safe metadata, references, scripts, and validation checks. Use when designing a new workflow. |
| [summarize-pdf](summarize-pdf/) | Summarize PDF reusable workflow. Use when this workflow is needed. |
| [verified-science-review](verified-science-review/) | High-confidence systematic literature review pipeline for biomedical research. Multi-pass architecture: retrieve → verify citations → synthesize → audit. Every PMID ve... |

### Transcriptomics

| Skill | Description |
|---|---|
| [bulk-rnaseq-counts-to-de-deseq2](bulk-rnaseq-counts-to-de-deseq2/) | End-to-end DESeq2 pipeline that turns a raw count matrix into a ranked, shrunken, QC-checked DE result table. Use when you have integer counts and biological replicate... |
| [cell-cell-communication](cell-cell-communication/) | Infer ligand-receptor signaling networks from annotated single-cell RNA-seq data and produce publication-ready communication plots. Use when you have a Seurat object w... |
| [scrna-trajectory-inference](scrna-trajectory-inference/) | Reconstruct differentiation trajectories from scRNA-seq using PAGA, diffusion pseudotime, scVelo RNA velocity, and CellRank fate probabilities. Use when ordering cells... |
| [scrnaseq-scanpy-core-analysis](scrnaseq-scanpy-core-analysis/) | End-to-end scRNA-seq analysis in Python using Scanpy and the scverse stack — QC, normalization, integration, clustering, annotation, pseudobulk DE. Use when analyzing... |
| [scrnaseq-seurat-core-analysis](scrnaseq-seurat-core-analysis/) | End-to-end scRNA-seq analysis in R using Seurat v5 — QC, SCTransform, Harmony / CCA / RPCA integration, clustering, annotation, pseudobulk DE. Use when working in R/Bi... |
| [spatial-transcriptomics](spatial-transcriptomics/) | Run a complete 10x Visium pipeline — QC, clustering, spatially variable genes, neighborhood enrichment, co-occurrence — on Space Ranger output, H5AD, or a built-in tis... |

### Functional Genomics

| Skill | Description |
|---|---|
| [functional-enrichment-from-degs](functional-enrichment-from-degs/) | Pathway / gene-set enrichment on differential-expression results using clusterProfiler. GSEA (rank-based, all genes) is the default; ORA (cutoff-based, significant gen... |
| [gene-essentiality](gene-essentiality/) | Query, interpret, and integrate DepMap (Cancer Dependency Map) CRISPR knockout essentiality scores — including CERES, Chronos, gene-effect inversion, pan-essential vs... |
| [pooled-crispr-screens](pooled-crispr-screens/) | Tiered analysis of pooled CRISPR screens with single-cell readout: fast t-test screening, on-target validation, and batch-corrected DE. Use when you have Perturb-seq,... |

### Genetics

| Skill | Description |
|---|---|
| [genetic-variant-annotation](genetic-variant-annotation/) | Annotate a VCF with functional consequences, clinical significance, population frequencies, and pathogenicity scores using Ensembl VEP or SNPEff. Use when you have cal... |

### Clinical

| Skill | Description |
|---|---|
| [lasso-biomarker-panel](lasso-biomarker-panel/) | Select a minimal (5–15 feature) biomarker panel from high-dimensional omics with nested-CV penalized logistic regression and stability selection; validate on an indepe... |

### Reporting

| Skill | Description |
|---|---|
| [docx-generation](docx-generation/) | Generate professional, scientifically-styled Word (.docx) reports with python-docx — title block, tables, embedded figures, callouts, headers/footers, and consistent t... |
| [pdf-report-generation](pdf-report-generation/) | Generate professional, scientifically-styled PDF reports with ReportLab Platypus — title block, embedded figures, styled tables, callout boxes, headers, footers with p... |

### Wet Lab

| Skill | Description |
|---|---|
| [pcr-primer-design](pcr-primer-design/) | Design and validate primers for PCR, qPCR, TaqMan probes, multiplex, and sequencing assays with thermodynamic QC and MIQE-compliant reporting. Use when you need oligos... |
| [splicecraft](splicecraft/) | Programmatic plasmid design, cloning simulation, and primer design via the SpliceCraft local API. Use when construct work should stay local. |

---

## Count Verification

This repository currently contains **35** skill directories with `SKILL.md` files.

```bash
find . -mindepth 2 -maxdepth 2 -name SKILL.md | wc -l
```
