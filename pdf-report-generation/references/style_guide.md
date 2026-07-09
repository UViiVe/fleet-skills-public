# PDF Style Guide

The canonical page geometry, typography, color, and ReportLab idioms for
the .pdf output. The helpers in `scripts/build_report.py` follow these
rules.

← Back to [SKILL.md](../SKILL.md)
Related: [Template Catalog](template_catalog.md)

---

## Page Geometry

| Property | Letter | A4 |
|----------|--------|----|
| Page size | 8.5 × 11 in (612 × 792 pt) | 8.27 × 11.69 in (595 × 842 pt) |
| Left margin | 60 pt | 60 pt |
| Right margin | 60 pt | 60 pt |
| Top margin | 52 pt | 52 pt |
| Bottom margin | 52 pt | 52 pt |
| Usable content width | 492 pt (~6.83 in) | 475 pt (~6.6 in) |
| Usable content height | 688 pt | 738 pt |

```python
from reportlab.lib.pagesizes import letter, A4
from reportlab.platypus import SimpleDocTemplate

doc = SimpleDocTemplate(
    output_path,
    pagesize=letter,
    leftMargin=60, rightMargin=60,
    topMargin=52, bottomMargin=52,
    title="Report Title",
    author="",
)
```

Set `title` in the `SimpleDocTemplate` constructor so the PDF metadata
matches what's in the visible title block.

## Color Tokens

```python
from reportlab.lib.colors import HexColor

HEADING_COLOR = HexColor("#111111")
BODY_COLOR    = HexColor("#2C2A26")
MUTED_COLOR   = HexColor("#8A8378")
ACCENT        = HexColor("#D4A04A")
ACCENT_LIGHT  = HexColor("#E8CC8A")
WARM_BG       = HexColor("#FAF9F3")
ALT_ROW       = HexColor("#F9F7F3")
BORDER        = HexColor("#D5CFC5")
LINK_COLOR    = HexColor("#0563C1")
WHITE         = HexColor("#FFFFFF")
```

For charts (when rendered via ReportLab Graphics), cycle through:

```python
CHART_PALETTE = [
    HexColor("#D4A04A"),  # accent
    HexColor("#0279EE"),  # blue
    HexColor("#75A025"),  # green
    HexColor("#FF9400"),  # orange
    HexColor("#FD9BED"),  # pink
    HexColor("#111111"),  # heading
]
```

To match an organization's brand, change `ACCENT` (and the first entry
of `CHART_PALETTE`).

## Typography

Use ReportLab's built-in Helvetica fonts — they embed cleanly across
viewers without font-substitution surprises.

| Element | Font | Size | Color | Notes |
|---------|------|------|-------|-------|
| Report title | Helvetica-Bold | 26 pt | HEADING | leading 32 |
| Subtitle | Helvetica | 11 pt | ACCENT | |
| Attribution | Helvetica-Oblique | 10 pt | MUTED | |
| Section head (H1) | Helvetica-Bold | 18 pt | HEADING | space before 24 pt, after 10 pt |
| Subhead (H2) | Helvetica-Bold | 13 pt | HEADING | space before 16 pt, after 8 pt |
| Body | Helvetica | 10.5 pt | BODY | justified, leading 15 |
| Caption | Helvetica-Oblique | 9 pt | MUTED | centered |
| Table header | Helvetica-Bold | 10 pt | WHITE | |
| Table cell | Helvetica | 10 pt | BODY | |
| Footer | Helvetica | 8 pt | MUTED | centered |
| Callout | Helvetica | 10.5 pt | BODY | inside warm-bg box |

Use `<sub>` and `<super>` XML tags for sub/superscripts inside Paragraph
strings. **Never** use Unicode superscripts (²³⁴, ¹⁰⁻⁶) — they render as
black boxes. Better still, write scientific notation as `3.36e-06`.

## Style Sheet

```python
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER

styles = getSampleStyleSheet()

styles.add(ParagraphStyle(
    name="ReportTitle", fontName="Helvetica-Bold",
    fontSize=26, leading=32,
    textColor=HEADING_COLOR, spaceBefore=0, spaceAfter=6,
))

styles.add(ParagraphStyle(
    name="Subtitle", fontName="Helvetica", fontSize=11,
    textColor=ACCENT, spaceAfter=4,
))

styles.add(ParagraphStyle(
    name="Attribution", fontName="Helvetica-Oblique", fontSize=10,
    textColor=MUTED_COLOR, spaceAfter=8,
))

styles.add(ParagraphStyle(
    name="SectionHead", fontName="Helvetica-Bold", fontSize=18,
    textColor=HEADING_COLOR, spaceBefore=24, spaceAfter=10,
))

styles.add(ParagraphStyle(
    name="SubHead", fontName="Helvetica-Bold", fontSize=13,
    textColor=HEADING_COLOR, spaceBefore=16, spaceAfter=8,
))

styles.add(ParagraphStyle(
    name="Body", fontName="Helvetica", fontSize=10.5,
    textColor=BODY_COLOR, alignment=TA_JUSTIFY,
    spaceAfter=8, leading=15,
))

styles.add(ParagraphStyle(
    name="Caption", fontName="Helvetica-Oblique", fontSize=9,
    textColor=MUTED_COLOR, alignment=TA_CENTER, spaceAfter=14,
))

styles.add(ParagraphStyle(
    name="TableHeader", fontName="Helvetica-Bold", fontSize=10,
    textColor=WHITE, alignment=TA_CENTER,
))

styles.add(ParagraphStyle(
    name="TableCell", fontName="Helvetica", fontSize=10,
    textColor=BODY_COLOR,
))

styles.add(ParagraphStyle(
    name="Callout", fontName="Helvetica", fontSize=10.5,
    textColor=BODY_COLOR, leading=15,
))
```

## Header and Footer

ReportLab fits the header and footer through a canvas callback. The same
callback runs on every page; pass it as both `onFirstPage` and
`onLaterPages`.

```python
def make_page_callback(header_text, accent_color=ACCENT):
    def page_callback(canvas, doc):
        canvas.saveState()
        w, h = doc.pagesize

        # Header
        canvas.setFont("Helvetica", 9)
        canvas.setFillColor(MUTED_COLOR)
        canvas.drawString(60, h - 40, header_text)
        canvas.setStrokeColor(accent_color)
        canvas.setLineWidth(1)
        canvas.line(60, h - 48, w - 60, h - 48)

        # Footer
        canvas.setStrokeColor(BORDER)
        canvas.setLineWidth(0.75)
        canvas.line(60, 40, w - 60, 40)
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(MUTED_COLOR)
        canvas.drawCentredString(w / 2, 26, f"Page {doc.page}")

        canvas.restoreState()
    return page_callback

doc.build(story, onFirstPage=make_page_callback("Report Title"),
          onLaterPages=make_page_callback("Report Title"))
```

Always wrap drawing code in `saveState`/`restoreState` so styles don't
leak into Platypus content rendered on the same page.

## Title Block

Reports open directly with content — no full-bleed cover page. The title
block is a flowable on page 1:

```python
from reportlab.platypus import Paragraph, Spacer

story.append(Spacer(1, 40))
story.append(Paragraph("Report Title", styles["ReportTitle"]))
story.append(Paragraph("Subtitle text", styles["Subtitle"]))
story.append(Spacer(1, 6))
story.append(Paragraph(f"<i>Generated 2026-05-20</i>", styles["Attribution"]))
story.append(Spacer(1, 24))
```

## Tables

```python
from reportlab.platypus import Table, TableStyle

def styled_table(headers, rows, col_widths, accent_color=ACCENT,
                 bold_first_column=True):
    header_paras = [Paragraph(f"<b>{h}</b>", styles["TableHeader"])
                    for h in headers]
    body_paras = [
        [Paragraph(str(c), styles["TableCell"]) for c in row]
        for row in rows
    ]
    if bold_first_column:
        for row in body_paras:
            row[0] = Paragraph(
                f"<b>{row[0].text}</b>", styles["TableCell"]
            )
    data = [header_paras, *body_paras]

    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.hAlign = "CENTER"
    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), accent_color),
        ("TEXTCOLOR", (0, 0), (-1, 0), WHITE),
        ("ALIGN", (0, 0), (-1, 0), "CENTER"),
        ("GRID", (0, 0), (-1, -1), 0.5, BORDER),
        ("BOX", (0, 0), (-1, -1), 0.75, BORDER),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ]
    for ri in range(2, len(data), 2):
        style_cmds.append(("BACKGROUND", (0, ri), (-1, ri), ALT_ROW))
    t.setStyle(TableStyle(style_cmds))
    return t
```

Column widths sum to ≤ 492 pt (Letter content width). `repeatRows=1`
makes the header re-render on continuation pages.

## Figures

Always wrap figure + caption in `KeepTogether` so the caption doesn't
strand on the next page.

```python
from reportlab.platypus import Image, KeepTogether, Spacer

def add_figure(image_path, caption, width_in=5.5):
    img = Image(image_path, width=width_in * 72, height=width_in * 72 * 0.6)
    img.hAlign = "CENTER"
    cap = Paragraph(f"<i>{caption}</i>", styles["Caption"])
    return KeepTogether([img, Spacer(1, 4), cap])
```

For accurate aspect ratio, query the image dimensions before fixing
height:

```python
from PIL import Image as PILImage

def add_figure(image_path, caption, width_in=5.5):
    with PILImage.open(image_path) as im:
        w_px, h_px = im.size
    aspect = h_px / w_px
    width_pt = width_in * 72
    img = Image(image_path, width=width_pt, height=width_pt * aspect)
    img.hAlign = "CENTER"
    cap = Paragraph(f"<i>{caption}</i>", styles["Caption"])
    return KeepTogether([img, Spacer(1, 4), cap])
```

## Callout Boxes

A warm-off-white single-cell table with an accent-colored left border:

```python
def callout_box(text, accent_color=ACCENT, width_in=6.0):
    cell = Paragraph(text, styles["Callout"])
    t = Table([[cell]], colWidths=[width_in * 72])
    t.hAlign = "CENTER"
    t.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), WARM_BG),
        ("BOX", (0, 0), (-1, -1), 0.5, BORDER),
        ("LINEBEFORE", (0, 0), (0, -1), 3, accent_color),
        ("TOPPADDING", (0, 0), (-1, -1), 12),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 12),
        ("LEFTPADDING", (0, 0), (-1, -1), 14),
        ("RIGHTPADDING", (0, 0), (-1, -1), 14),
    ]))
    return t
```

## Section Dividers

```python
from reportlab.platypus import HRFlowable

def divider(accent_color=ACCENT, width_in=6.8):
    return HRFlowable(
        width=width_in * 72, thickness=1, color=accent_color,
        spaceBefore=4, spaceAfter=10,
    )
```

## Bullet and Numbered Lists

```python
from reportlab.platypus import ListFlowable, ListItem

def bullet_list(items, style):
    return ListFlowable(
        [ListItem(Paragraph(t, style)) for t in items],
        bulletType="bullet",
        leftIndent=18,
    )

def numbered_list(items, style):
    return ListFlowable(
        [ListItem(Paragraph(t, style)) for t in items],
        bulletType="1",
        leftIndent=22,
    )
```

## Notation Conventions

- p-values, FDR: `3.36e-06`. **Do not** use Unicode superscripts.
- Confidence intervals: `(95% CI 1.12 — 1.45)`
- Percentages: `21.4%`
- Units with non-breaking space: `5.2&nbsp;µg/mL` (ReportLab respects
  `&nbsp;`)
- Greek: precomposed (α, β, γ) — Helvetica supports them via the
  standard Adobe glyph set

## KeepTogether and Page Breaks

```python
from reportlab.platypus import PageBreak

# Force a new page (e.g., between Methods and Results)
story.append(PageBreak())

# Bind a heading to the first paragraph that follows it so the heading
# never orphans at the bottom of a page
story.append(KeepTogether([heading_para, first_body_para]))
```

## Validation Rules

| Check | Threshold |
|-------|-----------|
| File size | > 5 KB |
| Page count | ≥ 1 |
| Extractable text on page 1 | non-empty |
| Total extractable text | > 200 chars |
| Embedded image count | matches expected (when spec known) |
| Image DPI | ≥ 300 (when Pillow available) |
| Title metadata | matches visible title |

---

Last updated: 2026-05-20.
