---
name: docx-generation
id: "skill_77a6487db1004cd1bcf59b222490ac0c"
display-name: "DOCX Report Generation"
short-description: "Generate professional, scientifically-styled Word (.docx) reports with python-docx — title block, tables, embedded figures, callouts, headers/footers, and consistent typography. Use when the user e..."
category: reporting
keywords: [docx, word, report, document, generation]
last-updated: "2026-07-09"
version: 1.0.0
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Generate professional, scientifically-styled Word (.docx) reports with python-docx — title block, tables, embedded figures, callouts, headers/footers, and consistent typography. Use when the user explicitly asks for an editable Word document or a deliverable that will go through review cycles."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
status: active
promoted: true
---

# Word Document Generation

Builds .docx reports with python-docx using a consistent scientific style:
US Letter page, Arial 11pt body, accent-color section heads, gold table
headers with alternating row shading, centered figures with italic
captions, and an automatic page footer.

The accompanying scripts in `scripts/` build a finished report from a
JSON spec, and validate the output before the user sees it.

## When to Use This Skill

✅ Use when:
- The user explicitly asks for a Word document, `.docx`, or "editable report"
- Collaborators need to comment, track changes, or edit the deliverable
- The output goes into Word-driven workflows (journal manuscripts, grant submissions, internal reviews)
- The user wants to add their own content or figures after generation

❌ Avoid when:
- The output is final and should not be edited — use **pdf-report-generation** instead
- The user just needs a brief markdown reply in chat
- The deliverable is a slide deck — use **pptx-generation**
- The data is purely tabular — a CSV/XLSX is more useful

## Quick Start

```python
from scripts.build_report import build_report

build_report(
    spec_path="examples/de_analysis_spec.json",
    output_path="document_de_analysis.docx",
)
```

A 4–6 page report with title block, executive summary, methods, a results
table, an embedded figure, and a footer with page numbers renders in
under 2 seconds.

A minimal JSON spec looks like this:

```json
{
  "title": "Differential Expression Analysis",
  "subtitle": "RNA-seq | Treated vs Control | HeLa",
  "attribution": "Generated 2026-05-20",
  "sections": [
    {"heading": "Executive Summary",
     "paragraphs": ["421 genes are significantly upregulated at FDR < 0.05."]},
    {"heading": "Results",
     "table": {
       "headers": ["Gene", "log2FC", "FDR"],
       "rows": [["TP53", "2.31", "3.36e-06"],
                ["MYC",  "1.94", "8.71e-05"]]
     }},
    {"heading": "Figure 1",
     "figure": {"path": "figures/volcano.png",
                "caption": "Volcano plot of all genes."}}
  ]
}
```

## Installation

| Package | Version | License | Commercial use | Install |
|---------|---------|---------|----------------|---------|
| python-docx | ≥ 1.1 | MIT | yes | `pip install python-docx` |
| Pillow | ≥ 10.0 | HPND | yes | `pip install Pillow` (for image dimension checks) |
| matplotlib | ≥ 3.6 | matplotlib | yes | optional, for figure generation |

```bash
pip install python-docx Pillow
```

## Inputs

**Required**
- A report specification — Python dict or JSON file with at minimum a `title`
  and a list of `sections`
- For embedded figures: image files referenced by relative or absolute path

**Optional**
- `subtitle` — accent-color subtitle below the main title
- `attribution` — italic muted line under the subtitle
- `header_text` — what appears in the running header (defaults to title)
- `accent_color` — hex color used for table headers, section dividers, and
  the header underline (default `D4A04A`)
- `body_color` — hex color for body text (default `2C2A26`)
- Tables can specify `col_widths` (in inches) and `bold_first_column`

**Style canon** (all reporting skills): US Letter page (8.5 × 11 in),
1-inch margins, Arial body, no Unicode superscripts (use `3.36e-06`
rather than 3.36 × 10⁻⁶), figures sized to ≥ 80% of the text area
(default 5.0 in), 300 DPI raster images embedded.

## Outputs

**Primary**
- A `.docx` file at `output_path`

**Optional validation artifacts**
- `<output>.pdf` — LibreOffice-rendered PDF for visual inspection
- `<output>.validation.json` — structured QC report listing paragraph
  count, table count, embedded image count, and any flagged issues

Downstream consumers: Word, Office 365, Google Docs, journal submission
portals, anything that accepts OOXML.

## Clarification Questions

### Q1. ⚠️ CRITICAL: ASK THIS FIRST — Is Word really the right format?

If the user just needs a document, the format determines the workflow:

- **Word** — recipient will edit, comment, or track changes
- **PDF** — recipient should not edit; final deliverable (use pdf-report-generation)
- **Markdown** — recipient will view in chat or a wiki
- **Slides** — recipient will present; use pptx-generation

Don't generate .docx unless the user explicitly named one of those terms or
the workflow obviously needs editability.

### Q2. What sections and content?

Default scientific report structure:
1. Title block
2. Executive Summary (2-3 paragraphs)
3. Methods (input data, tools, parameters)
4. Results (tables and figures with captions)
5. Discussion / Interpretation
6. Appendix (optional — full lists, parameter tables, supplementary figures)

Confirm with the user which to include and what data goes into each.

### Q3. Branding — accent color, attribution line?

The default accent (`D4A04A`, gold) is neutral. If the user has an
organization color, swap it in. The attribution defaults to the date; ask
if a person, department, or version string should appear.

### Q4. Any figures or tables to embed?

Get the source PNG paths and the table data (header + rows). For tables
with many rows, show the top N and summarize the rest — full data should
go to CSV.

### Q5. Where should the file land?

A relative path in the current working directory is fine. The script
verifies the parent directory exists.

## Standard Workflow

1. **Compose the spec** — Python dict or JSON file matching the schema in
   `references/template_catalog.md`
2. **Build** with `scripts/build_report.py` — runs all styling helpers and
   writes the .docx
3. **Validate** with `scripts/validate_output.py` — reopens the file,
   counts paragraphs/tables/images, checks file size and content density
4. **Optional visual check** — convert to PDF with LibreOffice and open in
   a reader to verify layout

The build script is intentionally idempotent: re-running with the same
spec produces a byte-identical file (modulo embedded timestamps where
Word always changes the digest).

## When Scripts Fail

1. **`KeyError: 'title'`** — the spec is missing required keys. Compare against
   `references/template_catalog.md`.
2. **`FileNotFoundError` on figure** — verify the image path exists and is
   readable; relative paths resolve from the script's working directory.
3. **Tables render with no borders** — the OOXML borders block needs to be
   applied. The helper does this; if you're calling python-docx directly,
   see `references/style_guide.md`.
4. **Fonts not applied** — python-docx does not inherit `font.name` reliably
   into table cells. The helper sets it on every run; if you bypass the
   helper, set `run.font.name = "Arial"` explicitly on every run.
5. **Image too wide for page** — drop `image.width` below 5.5 in for US
   Letter with 1-inch margins.
6. **Unicode rendering glitches** — Arial supports Greek (α, β), math (±, ×,
   °, ≤, ≥), and most accented characters. Avoid combining diacritics; use
   precomposed characters.

## Decision Points

| Choice | Default | When to change |
|--------|---------|----------------|
| Page size | US Letter (8.5 × 11 in) | A4 (8.27 × 11.69 in) for European recipients |
| Margins | 1 in all sides | 0.75 in only when content is very dense |
| Body font | Arial 11 pt | Bump to 12 pt for grant-style readability |
| Accent color | Gold `D4A04A` | Match recipient's brand color |
| Figure width | 5.0 in | 6.0 in for landscape figures; 4.0 in for portraits |
| Tables alignment | Center | Left when the table is narrow and inline |
| Header text | Report title | Project / version when the recipient cares about that |

## Common Issues

| Issue | Likely cause | Fix |
|-------|-------------|-----|
| Document under 5 KB after save | Empty paragraphs only | Re-run with a non-empty spec |
| Tables look uneven | `col_widths` not set | Provide widths in inches; sum should be ≤ 6.5 |
| Black background in table header | Used `w:val="solid"` instead of `clear` | Use the `styled_table` helper |
| Image not centered | Default left-align on paragraph | Set `paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER` |
| Font ignored in some runs | python-docx doesn't inherit | Use the `add_body` helper which sets `font.name` per run |
| Hyperlink not blue | Missing run-level style | Use `add_hyperlink` helper; it sets `Hyperlink` style |
| Footer page number frozen at 1 | Hardcoded "1" instead of PAGE field | Use the field code XML, see style guide |
| Bullet doesn't indent | Used `"• "` text prefix | Use `style="List Bullet"` |
| Heading shows in Calibri | python-docx default heading style | Override via the `add_section_heading` helper |
| Validation reports zero paragraphs | Wrote to a different path | Confirm `output_path` matches what validator opens |
| LibreOffice PDF conversion fails | LibreOffice not installed | Skip the visual check or install it: `brew install libreoffice` |
| Cells look cramped | Missing `tcMar` padding | Apply the `set_cell_padding` helper to every cell |

## Best Practices

- **One helper module, called everywhere.** All styling decisions live in
  `scripts/build_report.py` so the look-and-feel is consistent.
- **Always validate before shipping.** Run `validate_output.py` and
  surface its summary to the user; treat any flagged issue as blocking.
- **No Unicode superscripts in scientific notation.** Write `3.36e-06`, not
  `3.36 × 10⁻⁶` — superscripts often render incorrectly in non-Arial fonts
  on the recipient's machine.
- **Embed figures at 300 DPI.** Lower DPI looks blurry when the recipient
  zooms in. The `add_figure` helper does not resample but does verify DPI
  via Pillow.
- **Tables top N + summary.** A 500-row table buries the lede; show top 20
  and a brief total-row.

## Suggested Next Steps

- **Send for review.** The .docx is ready to share with collaborators.
- **Convert to PDF when finalized.** Either via Word (Save as PDF) or
  programmatically with LibreOffice (`libreoffice --headless --convert-to
  pdf`).
- **Use pdf-report-generation** for a parallel PDF deliverable when
  publishing or archiving.
- **For slide-deck output** — use pptx-generation with the same spec
  schema (the section/table/figure layout maps closely).

## Related Skills

- **pdf-report-generation** — same style canon, output is uneditable PDF
- **pptx-generation** — same content model, slide layout instead of pages
- Any analysis skill that produces tables and figures — the .docx is
  typically the last stop in the pipeline

## References

### Local documentation
- [references/style_guide.md](references/style_guide.md) — typography, color tokens,
  spacing, OOXML snippets for borders/shading/page numbers
- [references/template_catalog.md](references/template_catalog.md) — JSON schemas
  for the report spec, three worked examples, and notes on extending the templates

### Tools and tutorials
- python-docx documentation: https://python-docx.readthedocs.io/
- Office Open XML reference (ECMA-376): https://www.ecma-international.org/publications/standards/Ecma-376.htm
- LibreOffice headless conversion: https://help.libreoffice.org/latest/en-US/text/shared/guide/start_parameters.html