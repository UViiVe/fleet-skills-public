# Template Catalog

JSON spec schemas for the report types the builder handles out of the box,
plus three worked examples and notes on extending the catalog with a new
template.

← Back to [SKILL.md](../SKILL.md)
Related: [Style Guide](style_guide.md)

---

## Spec Schema Overview

A report spec is a JSON object with one mandatory `sections` array. Every
field except `title` and `sections` is optional.

```json
{
  "title": "string",
  "subtitle": "string (optional)",
  "attribution": "string (optional)",
  "header_text": "string (default: title)",
  "accent_color": "hex string (default: D4A04A)",
  "page_size": "letter | a4 (default: letter)",
  "sections": [
    {
      "heading": "string",
      "level": 1,
      "paragraphs": ["string", "string"],
      "table": {
        "headers": ["string"],
        "rows": [["any"]],
        "col_widths": [1.5, 1.0, 1.5],
        "bold_first_column": true
      },
      "figure": {
        "path": "path/to/image.png",
        "caption": "string",
        "width_in": 5.0
      },
      "callout": {
        "bold_prefix": "Key finding:",
        "text": "string"
      },
      "bullets": ["string"],
      "numbered": ["string"],
      "divider": false,
      "page_break_before": false
    }
  ]
}
```

A section may contain any combination of `paragraphs`, `table`, `figure`,
`callout`, `bullets`, `numbered`. They render in the order written above.

Use `level: 1` for top-level headings (H1) and `level: 2` for
subheadings (H2). The builder chooses font and spacing automatically.

---

## Template 1 — Scientific Analysis Report

For RNA-seq, proteomics, or any analysis with tables and a few figures.

```json
{
  "title": "Differential Expression Analysis",
  "subtitle": "RNA-seq | Treated vs Control | HeLa cell line",
  "attribution": "Generated 2026-05-20",
  "sections": [
    {
      "heading": "Executive Summary",
      "paragraphs": [
        "We compared HeLa cells treated with compound X (n = 3) against vehicle control (n = 3). 421 genes are significantly differentially expressed at FDR < 0.05; 187 are upregulated, 234 downregulated.",
        "The top hits cluster in MAPK and apoptosis pathways, consistent with the proposed mechanism of action."
      ]
    },
    {
      "heading": "Methods",
      "paragraphs": [
        "Raw reads were aligned to GRCh38 with STAR 2.7. Counts were summarized at the gene level with featureCounts. Differential expression was tested with DESeq2 using a Wald test on a model with treatment as the only covariate."
      ],
      "bullets": [
        "Aligner: STAR 2.7.10b (default parameters)",
        "Quantifier: featureCounts (Subread 2.0.6)",
        "DE method: DESeq2 1.42"
      ]
    },
    {
      "heading": "Top Differentially Expressed Genes",
      "level": 1,
      "table": {
        "headers": ["Gene", "log2FC", "FDR"],
        "rows": [
          ["TP53", "2.31", "3.36e-06"],
          ["MYC", "1.94", "8.71e-05"],
          ["BAX", "1.62", "4.10e-04"],
          ["BCL2", "-1.81", "5.20e-04"]
        ],
        "col_widths": [1.5, 1.2, 1.3],
        "bold_first_column": true
      }
    },
    {
      "heading": "Figure 1 — Volcano plot",
      "figure": {
        "path": "figures/volcano.png",
        "caption": "Volcano plot of all genes. Red: FDR < 0.05 and |log2FC| > 1. Top 10 most significant labeled."
      }
    },
    {
      "heading": "Discussion",
      "paragraphs": [
        "The increase in TP53 expression alongside elevated BAX and reduced BCL2 fits a classical apoptosis induction signature."
      ],
      "callout": {
        "bold_prefix": "Key finding:",
        "text": "Apoptotic markers BAX (up) and BCL2 (down) suggest compound X induces programmed cell death."
      }
    }
  ]
}
```

---

## Template 2 — Literature Review / Background

For background sections, target-disease overviews, or grant introductions.
Heavy on prose, light on figures.

```json
{
  "title": "HER2-Positive Breast Cancer: Background",
  "subtitle": "Target overview",
  "attribution": "Compiled 2026-05-20",
  "sections": [
    {
      "heading": "Disease Context",
      "paragraphs": [
        "HER2-positive breast cancer accounts for approximately 15–20% of invasive breast cancers. The HER2 (ERBB2) gene amplification drives constitutive receptor signaling through the MAPK and PI3K-AKT pathways.",
        "Despite targeted therapies (trastuzumab, pertuzumab, T-DM1), resistance develops in a substantial fraction of patients within 12–24 months of treatment initiation."
      ]
    },
    {
      "heading": "Existing Therapies",
      "table": {
        "headers": ["Therapy", "Modality", "First-line approval"],
        "rows": [
          ["Trastuzumab", "Monoclonal antibody", "1998"],
          ["Pertuzumab", "Monoclonal antibody", "2012"],
          ["T-DM1", "Antibody-drug conjugate", "2013"],
          ["T-DXd", "Antibody-drug conjugate", "2019"]
        ],
        "col_widths": [2.0, 2.2, 1.5]
      }
    },
    {
      "heading": "Unmet Need",
      "paragraphs": [
        "Patients who progress on T-DXd have few remaining options. There is a clear opportunity for novel formats (bispecifics, T-cell engagers) that can address resistance mechanisms."
      ]
    }
  ]
}
```

---

## Template 3 — Methods Protocol

For a procedural document — wet-lab protocol, computational pipeline,
or SOP. Heavy on numbered steps.

```json
{
  "title": "qPCR Standard Operating Procedure",
  "subtitle": "Quantification with SYBR Green I",
  "attribution": "Version 2.3 — 2026-05-20",
  "sections": [
    {
      "heading": "Materials",
      "bullets": [
        "SYBR Green Master Mix (Applied Biosystems, A25742)",
        "Forward and reverse primers (10 µM each in TE)",
        "cDNA template (5 ng/µL working stock)",
        "96-well qPCR plates and optical seals"
      ]
    },
    {
      "heading": "Reaction Setup",
      "numbered": [
        "Thaw all reagents on ice and vortex briefly to mix.",
        "Prepare a master mix for n + 1 reactions (n = number of samples).",
        "Aliquot 18 µL master mix into each well.",
        "Add 2 µL of cDNA template to each well.",
        "Seal the plate, centrifuge briefly to settle the reactions, and load into the instrument."
      ]
    },
    {
      "heading": "Cycling Program",
      "paragraphs": [
        "95°C for 10 min initial denaturation, followed by 40 cycles of 95°C for 15 s and 60°C for 60 s. Add a melt curve from 60 to 95°C in 0.5°C increments at the end."
      ]
    },
    {
      "heading": "Acceptance Criteria",
      "callout": {
        "bold_prefix": "Run accepted when:",
        "text": "Standard curve R² > 0.98, efficiency 90–110%, NTC Cq > 35, and the melt curve shows a single peak."
      }
    }
  ]
}
```

---

## Adding a New Template

To add a fourth template (say, "Clinical Case Report"), do three things:

1. **Add a section ordering convention.** Settle on the order of subsections
   the template implies (e.g., Patient → History → Findings → Treatment →
   Outcome).
2. **Document the schema below.** Show the section structure, mention any
   non-standard widgets (sidebar tables, multi-column lists).
3. **Add an example spec** to this file with realistic content.

The builder doesn't need code changes — it always honours the spec schema
above.

---

## Spec Validation Rules

The builder rejects specs that fail these:

- `title` is non-empty
- `sections` is a non-empty array
- Every section has either `heading` or content (`paragraphs`, `table`,
  `figure`, `callout`, `bullets`, `numbered`)
- Image paths referenced by `figure.path` exist on disk
- Table rows are all the same length as `headers`
- `col_widths` sum to ≤ 6.5 in (US Letter) or ≤ 6.27 in (A4)
- `accent_color` is a 6-char hex string
- `level` is 1 or 2

Validation failures are raised before any file is written.

---

## Style Inheritance

A spec inherits everything from the style guide:
- Arial throughout
- Body color BODY (`2C2A26`)
- Section headings HEADING (`111111`)
- Captions and attribution MUTED (`8A8378`)
- Tables, dividers, header underline ACCENT (`D4A04A` by default)

Override per-spec with `accent_color`. Per-element overrides (a single
red heading) require dropping to python-docx; the spec schema does not
support them.

---

Last updated: 2026-05-20.
