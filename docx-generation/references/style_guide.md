# Word Document Style Guide

The canonical typography, spacing, color, and OOXML-fragment reference for
docx output. Every helper in `scripts/build_report.py` follows this guide.

← Back to [SKILL.md](../SKILL.md)
Related: [Template Catalog](template_catalog.md)

---

## Page Setup

US Letter is the default. A4 is supported for European recipients.

| Property | Value |
|----------|-------|
| Page width | 8.5 in (US Letter) or 8.27 in (A4) |
| Page height | 11 in (US Letter) or 11.69 in (A4) |
| Top margin | 1 in |
| Bottom margin | 1 in |
| Left margin | 1 in |
| Right margin | 1 in |
| Orientation | Portrait |

```python
from docx.shared import Inches

section = doc.sections[0]
section.page_width = Inches(8.5)
section.page_height = Inches(11)
section.top_margin = Inches(1)
section.bottom_margin = Inches(1)
section.left_margin = Inches(1)
section.right_margin = Inches(1)
```

## Color Tokens

Use hex strings in OOXML attributes (`"D4A04A"`) and `RGBColor` constants
in python-docx (`RGBColor(0xD4, 0xA0, 0x4A)`).

| Token | Hex | Use |
|-------|-----|-----|
| HEADING | `111111` | section headings |
| BODY | `2C2A26` | body text |
| ACCENT | `D4A04A` | table headers, dividers, header underline |
| MUTED | `8A8378` | captions, attribution, footnotes |
| TABLE_HEADER_FILL | `D4A04A` | table header background |
| TABLE_ALT_FILL | `F9F7F3` | alternating row shading |
| TABLE_BORDER | `D5CFC5` | table grid lines |
| LINK | `0563C1` | hyperlinks |
| WHITE | `FFFFFF` | table header text |

To match a different accent (organization branding), change `ACCENT`,
`TABLE_HEADER_FILL`, and the header-underline color simultaneously.

## Typography

Arial for everything — universally available, renders predictably across
Windows, macOS, Linux, and the web.

| Element | Font | Size (pt) | Weight | Color |
|---------|------|-----------|--------|-------|
| Report title | Arial | 26 | Bold | HEADING |
| Subtitle | Arial | 11 | Normal | ACCENT |
| Attribution | Arial | 10 | Italic | MUTED |
| H1 (section) | Arial | 18 | Bold | HEADING |
| H2 (subsection) | Arial | 14 | Bold | HEADING |
| Body | Arial | 11 | Normal | BODY |
| Caption | Arial | 9 | Italic | MUTED |
| Table header | Arial | 10 | Bold | WHITE |
| Table cell | Arial | 10 | Normal | BODY |
| Footer | Arial | 8 | Normal | MUTED |

```python
FONT = "Arial"
```

## Spacing

| Element | Space before | Space after |
|---------|--------------|-------------|
| Title | 20 pt | 4 pt |
| Subtitle | 0 | 2 pt |
| Attribution | 0 | 16 pt |
| H1 | 24 pt | 10 pt |
| H2 | 16 pt | 8 pt |
| Body paragraph | 0 | 8 pt |
| List item | 0 | 4 pt |
| Figure | 12 pt | 4 pt |
| Caption | 4 pt | 14 pt |
| Section divider | 4 pt | 10 pt |

Line spacing on body paragraphs: 15 pt (set as
`paragraph_format.line_spacing = Pt(15)`).

## Header and Footer

The header carries the report title in muted text with a thin accent
underline. The footer centers the word "Page" followed by an auto-updating
page number.

### Header underline

```python
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

def apply_header_underline(header_paragraph, color="D4A04A", size="6"):
    pPr = header_paragraph._element.get_or_add_pPr()
    pBdr = OxmlElement("w:pBdr")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), size)
    bottom.set(qn("w:color"), color)
    bottom.set(qn("w:space"), "4")
    pBdr.append(bottom)
    pPr.append(pBdr)
```

### Footer with page number

```python
def add_page_number(footer_paragraph):
    from docx.shared import Pt, RGBColor

    run = footer_paragraph.add_run("Page ")
    run.font.name = "Arial"
    run.font.size = Pt(8)
    run.font.color.rgb = RGBColor(0x8A, 0x83, 0x78)

    # PAGE field code
    begin = OxmlElement("w:fldChar")
    begin.set(qn("w:fldCharType"), "begin")
    instr = OxmlElement("w:instrText")
    instr.set(qn("xml:space"), "preserve")
    instr.text = " PAGE "
    end = OxmlElement("w:fldChar")
    end.set(qn("w:fldCharType"), "end")

    for el in (begin, instr, end):
        r = footer_paragraph.add_run()
        r._element.append(el)
```

## Section Headings

```python
def add_section_heading(doc, text, level=1):
    from docx.shared import Pt, RGBColor
    heading = doc.add_heading(text, level=level)
    for run in heading.runs:
        run.font.name = "Arial"
        run.font.color.rgb = RGBColor(0x11, 0x11, 0x11)
    space = {1: (24, 10), 2: (16, 8)}.get(level, (10, 6))
    heading.paragraph_format.space_before = Pt(space[0])
    heading.paragraph_format.space_after = Pt(space[1])
    return heading
```

## Tables

Tables use accent-colored header fill, white bold header text, alternating
warm-off-white rows, and thin warm-gray borders. Cell padding is required —
without it, the cells look cramped.

```python
def styled_table(doc, headers, rows, col_widths=None, bold_first_column=True):
    from docx.shared import Pt, Inches, RGBColor
    from docx.enum.table import WD_TABLE_ALIGNMENT

    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER

    if col_widths:
        for i, w in enumerate(col_widths):
            for row in table.rows:
                row.cells[i].width = Inches(w)

    # Header
    for i, label in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(label)
        run.font.name = "Arial"
        run.font.size = Pt(10)
        run.bold = True
        run.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        _apply_shading(cell, "D4A04A")

    # Body
    for ri, row_data in enumerate(rows):
        for ci, value in enumerate(row_data):
            cell = table.rows[ri + 1].cells[ci]
            cell.text = ""
            run = cell.paragraphs[0].add_run(str(value))
            run.font.name = "Arial"
            run.font.size = Pt(10)
            run.font.color.rgb = RGBColor(0x2C, 0x2A, 0x26)
            if bold_first_column and ci == 0:
                run.bold = True
            if ri % 2 == 1:
                _apply_shading(cell, "F9F7F3")

    _apply_borders(table, color="D5CFC5", size="4")
    for row in table.rows:
        for cell in row.cells:
            _apply_cell_padding(cell, top=80, bottom=80, left=120, right=120)
    return table
```

The OOXML helpers (`_apply_shading`, `_apply_borders`,
`_apply_cell_padding`) live in `scripts/build_report.py`.

## Figures

Embed at the source resolution; size to about 5 inches wide for a US-
Letter page with 1-inch margins (that's roughly 80% of the 6.5-inch text
column). Use PNG, not JPG, for figures with sharp edges.

```python
def add_figure(doc, image_path, caption, width_in=5.0):
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    from docx.shared import Inches, Pt, RGBColor

    img_para = doc.add_paragraph()
    img_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    img_para.add_run().add_picture(image_path, width=Inches(width_in))
    img_para.paragraph_format.space_before = Pt(12)
    img_para.paragraph_format.space_after = Pt(4)

    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = cap.add_run(caption)
    run.font.name = "Arial"
    run.font.size = Pt(9)
    run.italic = True
    run.font.color.rgb = RGBColor(0x8A, 0x83, 0x78)
    cap.paragraph_format.space_after = Pt(14)
```

## Hyperlinks

```python
def add_hyperlink(paragraph, url, text):
    from docx.shared import RGBColor

    part = paragraph.part
    r_id = part.relate_to(
        url,
        "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
        is_external=True,
    )
    link = OxmlElement("w:hyperlink")
    link.set(qn("r:id"), r_id)

    r = OxmlElement("w:r")
    rPr = OxmlElement("w:rPr")
    style = OxmlElement("w:rStyle")
    style.set(qn("w:val"), "Hyperlink")
    color = OxmlElement("w:color")
    color.set(qn("w:val"), "0563C1")
    underline = OxmlElement("w:u")
    underline.set(qn("w:val"), "single")
    rPr.append(style)
    rPr.append(color)
    rPr.append(underline)
    r.append(rPr)

    t = OxmlElement("w:t")
    t.text = text
    r.append(t)
    link.append(r)
    paragraph._element.append(link)
```

## Callout Boxes

Single-cell tables with an accent-colored left border. Use for the most
important finding on the page.

```python
def callout_box(doc, body_text, bold_prefix=None):
    table = doc.add_table(rows=1, cols=1)
    cell = table.rows[0].cells[0]
    cell.text = ""
    p = cell.paragraphs[0]
    if bold_prefix:
        b = p.add_run(bold_prefix + " ")
        b.font.name = "Arial"
        b.font.size = Pt(10)
        b.bold = True
        b.font.color.rgb = RGBColor(0x2C, 0x2A, 0x26)
    t = p.add_run(body_text)
    t.font.name = "Arial"
    t.font.size = Pt(10)
    t.font.color.rgb = RGBColor(0x2C, 0x2A, 0x26)
    _apply_shading(cell, "FAF9F3")
    _apply_border_per_edge(
        cell,
        {"start": ("D4A04A", "18"),
         "top":   ("D5CFC5", "4"),
         "bottom":("D5CFC5", "4"),
         "end":   ("D5CFC5", "4")},
    )
    _apply_cell_padding(cell, top=120, bottom=120, left=160, right=160)
    return table
```

## Section Dividers

A thin accent-colored bottom border on an empty paragraph:

```python
def add_divider(doc):
    para = doc.add_paragraph()
    pPr = para._element.get_or_add_pPr()
    pBdr = OxmlElement("w:pBdr")
    bottom = OxmlElement("w:bottom")
    bottom.set(qn("w:val"), "single")
    bottom.set(qn("w:sz"), "6")
    bottom.set(qn("w:color"), "D4A04A")
    bottom.set(qn("w:space"), "1")
    pBdr.append(bottom)
    pPr.append(pBdr)
    para.paragraph_format.space_before = Pt(4)
    para.paragraph_format.space_after = Pt(10)
```

## Page Breaks

```python
from docx.enum.text import WD_BREAK

para = doc.add_paragraph()
para.add_run().add_break(WD_BREAK.PAGE)
```

## Bullet and Numbered Lists

Use the built-in styles, never raw Unicode bullets:

```python
doc.add_paragraph("First", style="List Bullet")
doc.add_paragraph("Second", style="List Bullet")

doc.add_paragraph("Step 1", style="List Number")
doc.add_paragraph("Step 2", style="List Number")
```

Then apply font on the runs:

```python
for run in para.runs:
    run.font.name = "Arial"
    run.font.size = Pt(11)
    run.font.color.rgb = RGBColor(0x2C, 0x2A, 0x26)
```

## Numbers and Scientific Notation

- p-values, FDR, q-values: write as `3.36e-06`. **Do not** use Unicode
  superscripts (3.36 × 10⁻⁶) — they don't render reliably in every viewer.
- Confidence intervals: `(95% CI 1.12 — 1.45)`
- Percentages: `21.4%` not `21.4 %`
- Units: `5.2 µg/mL` (Arial includes µ; non-breaking space between number
  and unit)
- Greek letters: precomposed (α, β, γ, σ²); avoid combining diacritics

## Validation Rules

After saving, the validation helper checks:

| Check | Threshold |
|-------|-----------|
| File size | > 5 KB |
| Paragraph count | ≥ 5 |
| Total text length | ≥ 100 characters |
| Every table has borders | yes |
| Every figure paragraph is centered | yes |
| All runs have `font.name` set | yes |
| Image DPI | ≥ 300 (if Pillow available) |

The skill will refuse to declare success if any of these fail.

---

Last updated: 2026-05-20.
