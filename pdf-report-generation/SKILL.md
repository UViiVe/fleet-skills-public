---
name: pdf-report-generation
id: "skill_6e01fd0d1d034cdb9157d0d686906517"
display-name: "PDF Report Generation"
short-description: "Generate professional, scientifically-styled PDF reports with ReportLab Platypus — title block, embedded figures, styled tables, callout boxes, headers, footers with page numbers, and section divid..."
category: reporting
keywords: [pdf, report, latex, generation, document]
last-updated: "2026-07-09"
version: 1.0.0
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Generate professional, scientifically-styled PDF reports with ReportLab Platypus — title block, embedded figures, styled tables, callout boxes, headers, footers with page numbers, and section dividers. Use when the user explicitly asks for a PDF, a polished report, or a self-contained deliverable to share outside the platform."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
status: active
promoted: true
---

# PDF Report Generation

Builds .pdf reports with ReportLab Platypus: US Letter pages, Helvetica
throughout, accent-color section dividers and table headers, alternating
row shading, centered figures with italic captions, header with title and
accent underline, footer with page numbers.

`scripts/build_report.py` turns a JSON spec into a finished PDF.
`scripts/validate_output.py` opens the PDF with pypdf and runs structural
QC before declaring it shippable.

## When to Use This Skill

✅ Use when:
- The user explicitly asks for a PDF, "report", or "polished deliverable"
- The output goes to PIs, collaborators, grant reviewers, or anyone
  outside the platform
- The deliverable is final — not meant to be edited
- The user wants multi-page output with figures and tables that pages
  consistently

❌ Avoid when:
- The user will edit it — use docx-generation
- The user will present it — use pptx-generation
- They just want a markdown summary in chat
- They need to extract or merge existing PDFs — use pypdf directly

## Quick Start

```python
from scripts.build_report import build_report

build_report(
    spec_path="examples/de_analysis_spec.json",
    output_path="report_de_analysis.pdf",
)
```

A 5–8 page report with title block, executive summary, methods, results
table, and an embedded figure renders in ~2 seconds.

The spec schema is identical to docx-generation, so the same content
renders cleanly to both formats. Minimal example:

```json
{
  "title": "Differential Expression Analysis",
  "subtitle": "RNA-seq | Treated vs Control",
  "attribution": "Generated 2026-05-20",
  "sections": [
    {"heading": "Executive Summary",
     "paragraphs": ["421 genes are significantly upregulated at FDR < 0.05."]},
    {"heading": "Results",
     "table": {"headers": ["Gene", "log2FC", "FDR"],
               "rows": [["TP53", "2.31", "3.36e-06"]]}}
  ]
}
```

## Installation

| Package | Version | License | Commercial use | Install |
|---------|---------|---------|----------------|---------|
| reportlab | ≥ 4.0 | BSD-3 | yes | `pip install reportlab` |
| pypdf | ≥ 4.0 | BSD-3 | yes | `pip install pypdf` |
| Pillow | ≥ 10.0 | HPND | yes | `pip install Pillow` (image dimension checks) |

```bash
pip install reportlab pypdf Pillow
```

## Inputs

**Required**
- A spec — dict or JSON file with `title` and non-empty `sections` array
- Image files referenced by `figure.path` must exist before build

**Optional**
- `subtitle`, `attribution` — title block content
- `header_text` — header text (default: title)
- `accent_color` — hex used for section dividers, table headers, and the
  header underline (default `D4A04A`)
- `page_size` — `letter` (default) or `a4`
- Per-table `col_widths` in inches, `bold_first_column` flag

**Style canon** (all reporting skills): US Letter (8.5 × 11 in), 60-pt
side margins, 52-pt top/bottom margins, Helvetica body, no Unicode
superscripts (write `3.36e-06`, not 3.36 × 10⁻⁶), figures sized to ≥ 80%
of content width (typically 5.5 in wide), raster images embedded at 300
DPI.

## Outputs

**Primary**
- A `.pdf` file at `output_path`

**Validation artifacts**
- Structured QC summary (page count, text extractability, image count,
  table count, font usage, image DPI)

Downstream consumers: any PDF reader, journal submission portals,
electronic lab notebooks, archival systems.

## Clarification Questions

### Q1. ⚠️ CRITICAL: ASK THIS FIRST — Is PDF the right format?

If the user will edit it: pick docx-generation. If they'll present it:
pptx-generation. If they need to share an immutable final deliverable
or archive a result: PDF is the right call.

Ambiguous requests like "make a report" usually want PDF when the
recipient is external, docx when the recipient is internal.

### Q2. What sections and content?

Default scientific report:
1. Title block
2. Executive Summary
3. Methods (input, tools, parameters)
4. Results (tables and figures with captions)
5. Discussion
6. Appendix (optional)

Confirm with the user which to include.

### Q3. Page size — Letter or A4?

US recipients expect Letter; European recipients expect A4. Default to
Letter unless told otherwise.

### Q4. Branding — accent color?

Default gold (`D4A04A`) is neutral. If the organization has a specific
brand color, swap it in by passing `accent_color` in the spec.

### Q5. Any figures, tables, or callouts to embed?

Get source PNG paths, table data, and callout text. For tables with
many rows, show the top N and summarize the rest — the full table
should go in a CSV alongside the PDF.

## Standard Workflow

1. **Compose the spec** — dict or JSON file matching
   `references/template_catalog.md`
2. **Build** with `scripts/build_report.py` — runs styled paragraphs,
   tables, figures, dividers, and adds the canvas-callback header/footer
3. **Validate** with `scripts/validate_output.py` — page count, text
   extractability, font consistency
4. **Open in a reader** to verify layout

The build script always wraps figure + caption in `KeepTogether` to
prevent orphan captions across pages.

## When Scripts Fail

1. **`KeyError: 'title'`** — spec is missing required keys. Compare against
   `references/template_catalog.md`.
2. **`FileNotFoundError` on figure** — verify the image exists at the
   given path; relative paths resolve from the working directory.
3. **Table overflows page width** — sum of `col_widths` exceeds the
   ~492-pt usable content width (Letter, 60-pt side margins). Tighten
   widths or reduce columns.
4. **Black boxes where superscripts should be** — Unicode superscripts in
   plain text. Use `<super>` XML tags inside Paragraph strings, or write
   plain `e-notation` (`3.36e-06`).
5. **Page numbers stuck at 1** — canvas callback wasn't passed to
   `doc.build`. Confirm `onFirstPage` and `onLaterPages` are set.
6. **Figure split from its caption** — wrap in `KeepTogether` (the
   helper does this).

## Decision Points

| Choice | Default | When to change |
|--------|---------|----------------|
| Page size | Letter | A4 for European recipients |
| Margins | 60/52 pt | 50/40 pt if content is dense |
| Body font size | 10.5 pt | 11 pt for older readers |
| Accent color | `D4A04A` | Match organization brand |
| Figure width | 5.5 in | Smaller for portraits, larger for landscape |
| Header text | Report title | Project name when relevant |
| Footer page-number style | `Page X` | `Page X of Y` requires a two-pass build |
| Repeat header on long tables | Yes | Always — `repeatRows=1` |

## Common Issues

| Issue | Likely cause | Fix |
|-------|-------------|-----|
| PDF under 5 KB | Empty story | Re-run with non-empty `sections` |
| First page text extraction empty | Image-only content | Add at least one Paragraph |
| Table grid invisible | Missing `GRID` in TableStyle | Use the `styled_table` helper |
| Header text overflowing | Title too long | Truncate via `header_text` in spec |
| Footer overlaps content | Bottom margin too small | Restore `bottomMargin=52` |
| Page numbers stuck | Callback not registered | Pass both `onFirstPage` and `onLaterPages` |
| Orphan caption | Manual flowable insertion | Use the `add_figure` helper which uses `KeepTogether` |
| Wrong font in some paragraphs | Style not applied | Use the `Body`, `Caption`, `SectionHead` styles consistently |
| Image looks pixelated | Embedded at too-large dimensions | Limit width to 5.5 in; resize source if huge |
| File size > 10 MB | Many high-DPI raster images | Resample sources to actual display size at 300 DPI |
| `KeepTogether overflowed frame` | Figure + caption taller than page | Reduce figure dimensions |
| Unicode µ missing | Helvetica doesn't include it | Use `<font name="Symbol">m</font>` or stick to `ug/mL` |

## Best Practices

- **Use Helvetica.** ReportLab embeds it correctly on every system. Don't
  rely on system fonts.
- **Set `colWidths` on every table.** Without them, ReportLab guesses
  poorly and tables look uneven.
- **Always center figures, images, and tables** by setting `hAlign =
  "CENTER"`. The helpers do this for you.
- **`KeepTogether` figure + caption.** A caption marooned on the next
  page is the most common pdf layout mistake.
- **Validate before shipping.** Run `validate_output.py` and surface the
  summary to the user.
- **300-DPI raster images.** Lower DPI looks fine on screen and terrible
  when printed.

## Suggested Next Steps

- **Archive the PDF** alongside the analysis it describes.
- **Send for review.** Recipient can annotate with their PDF reader of
  choice.
- **Build a parallel editable version** with docx-generation when the
  recipient also wants to comment in-line.
- **For an oral version** — use pptx-generation with the same spec.

## Related Skills

- **docx-generation** — editable Word output, same style canon
- **pptx-generation** — slide deck, same style canon
- Any analysis skill that produces tables and figures — PDF is typically
  the final deliverable

## References

### Local documentation
- [references/style_guide.md](references/style_guide.md) — page geometry,
  typography, color tokens, ReportLab style sheets, canvas-callback
  patterns
- [references/template_catalog.md](references/template_catalog.md) — JSON
  spec schema, three worked examples, instructions for extending the
  templates

### Tools and tutorials
- ReportLab User Guide: https://www.reportlab.com/docs/reportlab-userguide.pdf
- ReportLab Platypus reference: https://www.reportlab.com/docs/reportlab-reference.pdf
- pypdf documentation: https://pypdf.readthedocs.io/
- PDF/A archival standards: https://www.iso.org/standard/57229.html