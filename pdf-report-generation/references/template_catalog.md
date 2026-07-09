# Template Catalog

JSON spec schemas for the report types the builder handles, with three
worked examples and notes on adding a new template.

The spec schema is intentionally the same as docx-generation so a single
spec can render to either format.

← Back to [SKILL.md](../SKILL.md)
Related: [Style Guide](style_guide.md)

---

## Spec Schema Overview

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
      "paragraphs": ["string"],
      "bullets": ["string"],
      "numbered": ["string"],
      "table": {
        "headers": ["string"],
        "rows": [["any"]],
        "col_widths_in": [1.5, 1.0, 1.5],
        "bold_first_column": true
      },
      "figure": {
        "path": "path/to/image.png",
        "caption": "string",
        "width_in": 5.5
      },
      "callout": {"text": "string"},
      "divider": false,
      "page_break_before": false
    }
  ]
}
```

A section can mix any of `paragraphs`, `bullets`, `numbered`, `table`,
`figure`, `callout`. They render in that order.

Set `level: 1` for top-level headings (H1) or `level: 2` for subheadings
(H2).

---

## Template 1 — Scientific Analysis Report

```json
{
  "title": "Differential Expression Analysis",
  "subtitle": "RNA-seq | Treated vs Control | HeLa",
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
        "Raw reads were aligned to GRCh38 with STAR 2.7. Counts were summarized at the gene level with featureCounts. Differential expression was tested with DESeq2 (Wald test, FDR cutoff 0.05)."
      ],
      "bullets": [
        "Aligner: STAR 2.7.10b",
        "Quantifier: featureCounts 2.0.6",
        "DE method: DESeq2 1.42"
      ]
    },
    {
      "heading": "Top Differentially Expressed Genes",
      "table": {
        "headers": ["Gene", "log2FC", "FDR"],
        "rows": [
          ["TP53", "2.31", "3.36e-06"],
          ["MYC", "1.94", "8.71e-05"],
          ["BAX", "1.62", "4.10e-04"],
          ["BCL2", "-1.81", "5.20e-04"]
        ],
        "col_widths_in": [1.8, 1.4, 1.5],
        "bold_first_column": true
      }
    },
    {
      "heading": "Figure 1 — Volcano plot",
      "figure": {
        "path": "figures/volcano.png",
        "caption": "Volcano plot of all genes. Red: FDR < 0.05 and |log2FC| > 1."
      }
    },
    {
      "heading": "Discussion",
      "paragraphs": [
        "TP53 upregulation alongside BAX (up) and BCL2 (down) fits a classical apoptosis induction signature."
      ],
      "callout": {
        "text": "Apoptotic markers BAX (up) and BCL2 (down) indicate compound X induces programmed cell death via the TP53 axis."
      }
    }
  ]
}
```

---

## Template 2 — Literature Review / Target Overview

```json
{
  "title": "HER2-Positive Breast Cancer: Target Overview",
  "subtitle": "Background and unmet need",
  "attribution": "Compiled 2026-05-20",
  "sections": [
    {
      "heading": "Disease Context",
      "paragraphs": [
        "HER2-positive breast cancer accounts for 15–20% of invasive breast cancers. HER2 (ERBB2) gene amplification drives constitutive receptor signaling through MAPK and PI3K-AKT.",
        "Despite the success of trastuzumab and downstream therapies, resistance develops in a meaningful fraction of patients within 12–24 months of treatment initiation."
      ]
    },
    {
      "heading": "Approved Therapies",
      "table": {
        "headers": ["Therapy", "Modality", "First-line approval"],
        "rows": [
          ["Trastuzumab", "Monoclonal antibody", "1998"],
          ["Pertuzumab", "Monoclonal antibody", "2012"],
          ["T-DM1", "Antibody-drug conjugate", "2013"],
          ["T-DXd", "Antibody-drug conjugate", "2019"]
        ],
        "col_widths_in": [2.0, 2.2, 1.5]
      }
    },
    {
      "heading": "Unmet Need",
      "paragraphs": [
        "Patients who progress on T-DXd have limited remaining options. There is a clear opportunity for novel formats (bispecifics, T-cell engagers) addressing the resistance mechanism."
      ],
      "callout": {
        "text": "Late-line HER2+ patients (post-T-DXd progression) are the highest-need subset; novel formats targeting CD3-recruited cytotoxicity show the most promise."
      }
    }
  ]
}
```

---

## Template 3 — Methods / SOP

```json
{
  "title": "qPCR Standard Operating Procedure",
  "subtitle": "SYBR Green I quantification",
  "attribution": "Version 2.3 — 2026-05-20",
  "sections": [
    {
      "heading": "Materials",
      "bullets": [
        "SYBR Green Master Mix (Applied Biosystems A25742)",
        "Forward and reverse primers (10 µM each in TE)",
        "cDNA template (5 ng/µL working stock)",
        "96-well qPCR plates and optical seals"
      ]
    },
    {
      "heading": "Reaction Setup",
      "numbered": [
        "Thaw all reagents on ice and vortex briefly to mix.",
        "Prepare a master mix for n + 1 reactions.",
        "Aliquot 18 µL master mix into each well.",
        "Add 2 µL cDNA template per well.",
        "Seal the plate, centrifuge briefly, and load into the instrument."
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
        "text": "Standard curve R² > 0.98, efficiency 90–110%, NTC Cq > 35, single melt-curve peak."
      }
    }
  ]
}
```

---

## Adding a New Template

To document a new report type (say, "Clinical Case Report"):

1. **Decide the section sequence** that template will use. Settle on
   conventions for what each section contains.
2. **Document the schema** as a worked example below.
3. **No code change is needed** — the builder dispatches on the spec
   schema, not on a `type` field.

---

## Spec Validation Rules

The builder rejects specs that fail:

- `title` is non-empty
- `sections` is a non-empty array
- Every section has either `heading` or at least one content block
- Image paths referenced by `figure.path` exist on disk
- Table rows match `headers` length
- `col_widths_in` sum to ≤ 6.83 in (Letter) or ≤ 6.6 in (A4)
- `accent_color` is a 6-char hex string
- `level` is 1 or 2

Validation failures are raised before any file is written.

---

## Inheritance from Style Guide

Every spec inherits:
- Helvetica throughout
- Body color BODY_COLOR (`2C2A26`)
- Section headings HEADING_COLOR (`111111`)
- Captions and attribution MUTED_COLOR (`8A8378`)
- Tables, dividers, header underline ACCENT (`D4A04A` by default)

Override `accent_color` per spec. Per-element overrides (a single red
heading) require dropping to ReportLab; the spec schema does not support
them.

---

Last updated: 2026-05-20.
