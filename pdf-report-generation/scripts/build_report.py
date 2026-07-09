"""Build a styled PDF report from a JSON spec using ReportLab Platypus.

US Letter (or A4) page, Helvetica throughout, accent-color section
dividers and table headers, alternating row shading, centered figures
with italic captions, header with title and accent underline, footer
with page numbers.

Inputs (CLI):
    --spec PATH
    --out PATH
Outputs:
    A .pdf at --out.

Example invocation
------------------
    python -m scripts.build_report --spec examples/de_analysis.json --out report.pdf
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional, Sequence

try:
    from reportlab.lib.colors import HexColor
    from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
    from reportlab.lib.pagesizes import A4, letter
    from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
    from reportlab.platypus import (
        HRFlowable, Image, KeepTogether, ListFlowable, ListItem,
        PageBreak, Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle,
    )
except ImportError:
    sys.exit("Install with: pip install reportlab")

try:
    from PIL import Image as PILImage
    _HAS_PIL = True
except ImportError:
    _HAS_PIL = False


HEADING_COLOR = HexColor("#111111")
BODY_COLOR = HexColor("#2C2A26")
MUTED_COLOR = HexColor("#8A8378")
DEFAULT_ACCENT = HexColor("#D4A04A")
ALT_ROW = HexColor("#F9F7F3")
BORDER = HexColor("#D5CFC5")
WARM_BG = HexColor("#FAF9F3")
WHITE = HexColor("#FFFFFF")

PAGE_SIZES = {"letter": letter, "a4": A4}
LEFT_MARGIN = 60
RIGHT_MARGIN = 60
TOP_MARGIN = 52
BOTTOM_MARGIN = 52


def build_report(
    spec: Optional[Dict[str, Any]] = None,
    spec_path: Optional[str] = None,
    output_path: str = "report.pdf",
) -> str:
    if spec is None and spec_path is None:
        raise ValueError("Provide spec or spec_path.")
    if spec is None:
        with open(spec_path) as f:
            spec = json.load(f)
    validate_spec(spec)

    accent_hex = spec.get("accent_color", "D4A04A")
    accent_color = HexColor(f"#{accent_hex}")
    pagesize = PAGE_SIZES[spec.get("page_size", "letter")]
    header_text = spec.get("header_text", spec["title"])

    styles = _make_styles(accent_color)
    story = _build_story(spec, styles, accent_color)

    parent = os.path.dirname(os.path.abspath(output_path))
    if parent:
        os.makedirs(parent, exist_ok=True)

    doc = SimpleDocTemplate(
        output_path,
        pagesize=pagesize,
        leftMargin=LEFT_MARGIN, rightMargin=RIGHT_MARGIN,
        topMargin=TOP_MARGIN, bottomMargin=BOTTOM_MARGIN,
        title=spec["title"], author="",
    )
    callback = _page_callback(header_text, accent_color)
    doc.build(story, onFirstPage=callback, onLaterPages=callback)

    print(f"✓ pdf written: {output_path} "
          f"({os.path.getsize(output_path)} bytes)")
    return output_path


def validate_spec(spec: Dict[str, Any]) -> None:
    if not isinstance(spec, dict):
        raise TypeError("Spec must be a dict.")
    if not spec.get("title"):
        raise ValueError("Spec is missing `title`.")
    sections = spec.get("sections")
    if not isinstance(sections, list) or not sections:
        raise ValueError("Spec is missing a non-empty `sections` array.")
    for i, s in enumerate(sections):
        if not isinstance(s, dict):
            raise ValueError(f"Section {i} is not a dict.")
        if not any(
            s.get(k)
            for k in ("heading", "paragraphs", "table", "figure",
                      "callout", "bullets", "numbered")
        ):
            raise ValueError(f"Section {i} has no content.")
        if "table" in s and s["table"]:
            t = s["table"]
            if not t.get("headers") or not isinstance(t.get("rows"), list):
                raise ValueError(f"Section {i} table malformed.")
            ncols = len(t["headers"])
            for r in t["rows"]:
                if len(r) != ncols:
                    raise ValueError(f"Section {i} table row mismatch.")
        if "figure" in s and s["figure"]:
            path = s["figure"].get("path")
            if not path or not os.path.exists(path):
                raise FileNotFoundError(f"Figure {path!r} not found (section {i}).")

    accent = spec.get("accent_color", "D4A04A")
    if not (isinstance(accent, str) and len(accent) == 6):
        raise ValueError(f"accent_color {accent!r} must be a 6-char hex string.")
    page_size = spec.get("page_size", "letter")
    if page_size.lower() not in PAGE_SIZES:
        raise ValueError(f"page_size must be one of {list(PAGE_SIZES)}.")


def _make_styles(accent_color) -> Dict[str, ParagraphStyle]:
    base = getSampleStyleSheet()
    out = {}
    out["ReportTitle"] = ParagraphStyle(
        name="ReportTitle", fontName="Helvetica-Bold", fontSize=26, leading=32,
        textColor=HEADING_COLOR, spaceBefore=0, spaceAfter=6,
    )
    out["Subtitle"] = ParagraphStyle(
        name="Subtitle", fontName="Helvetica", fontSize=11,
        textColor=accent_color, spaceAfter=4,
    )
    out["Attribution"] = ParagraphStyle(
        name="Attribution", fontName="Helvetica-Oblique", fontSize=10,
        textColor=MUTED_COLOR, spaceAfter=8,
    )
    out["SectionHead"] = ParagraphStyle(
        name="SectionHead", fontName="Helvetica-Bold", fontSize=18,
        textColor=HEADING_COLOR, spaceBefore=24, spaceAfter=10,
    )
    out["SubHead"] = ParagraphStyle(
        name="SubHead", fontName="Helvetica-Bold", fontSize=13,
        textColor=HEADING_COLOR, spaceBefore=16, spaceAfter=8,
    )
    out["Body"] = ParagraphStyle(
        name="Body", fontName="Helvetica", fontSize=10.5,
        textColor=BODY_COLOR, alignment=TA_JUSTIFY,
        spaceAfter=8, leading=15,
    )
    out["Caption"] = ParagraphStyle(
        name="Caption", fontName="Helvetica-Oblique", fontSize=9,
        textColor=MUTED_COLOR, alignment=TA_CENTER, spaceAfter=14,
    )
    out["TableHeader"] = ParagraphStyle(
        name="TableHeader", fontName="Helvetica-Bold", fontSize=10,
        textColor=WHITE, alignment=TA_CENTER,
    )
    out["TableCell"] = ParagraphStyle(
        name="TableCell", fontName="Helvetica", fontSize=10,
        textColor=BODY_COLOR,
    )
    out["TableCellBold"] = ParagraphStyle(
        name="TableCellBold", fontName="Helvetica-Bold", fontSize=10,
        textColor=BODY_COLOR,
    )
    out["Callout"] = ParagraphStyle(
        name="Callout", fontName="Helvetica", fontSize=10.5,
        textColor=BODY_COLOR, leading=15,
    )
    out["ListItem"] = ParagraphStyle(
        name="ListItem", fontName="Helvetica", fontSize=10.5,
        textColor=BODY_COLOR, leading=14,
    )
    return out


def _build_story(spec: Dict[str, Any], styles: Dict[str, ParagraphStyle],
                 accent_color) -> List[Any]:
    story: List[Any] = []
    story.append(Spacer(1, 40))
    story.append(Paragraph(spec["title"], styles["ReportTitle"]))
    if spec.get("subtitle"):
        story.append(Paragraph(spec["subtitle"], styles["Subtitle"]))
    if spec.get("attribution"):
        story.append(Spacer(1, 4))
        story.append(
            Paragraph(f"<i>{spec['attribution']}</i>", styles["Attribution"])
        )
    story.append(Spacer(1, 16))

    for section in spec["sections"]:
        if section.get("page_break_before"):
            story.append(PageBreak())

        if section.get("heading"):
            level = section.get("level", 1)
            head_style = styles["SectionHead"] if level == 1 else styles["SubHead"]
            story.append(Paragraph(section["heading"], head_style))

        for para in section.get("paragraphs", []) or []:
            story.append(Paragraph(para, styles["Body"]))

        if section.get("bullets"):
            story.append(_bullet_list(section["bullets"], styles["ListItem"]))
        if section.get("numbered"):
            story.append(_numbered_list(section["numbered"], styles["ListItem"]))

        if section.get("table"):
            story.append(_styled_table(
                headers=section["table"]["headers"],
                rows=section["table"]["rows"],
                col_widths_in=section["table"].get("col_widths_in"),
                bold_first_column=section["table"].get("bold_first_column", True),
                styles=styles, accent_color=accent_color,
            ))
            story.append(Spacer(1, 12))

        if section.get("figure"):
            story.append(_add_figure(
                image_path=section["figure"]["path"],
                caption=section["figure"]["caption"],
                width_in=section["figure"].get("width_in", 5.5),
                styles=styles,
            ))

        if section.get("callout"):
            story.append(_callout_box(
                text=section["callout"]["text"],
                styles=styles,
                accent_color=accent_color,
            ))
            story.append(Spacer(1, 12))

        if section.get("divider"):
            story.append(_divider(accent_color))

    return story


def _bullet_list(items: Sequence[str], style: ParagraphStyle) -> ListFlowable:
    return ListFlowable(
        [ListItem(Paragraph(t, style)) for t in items],
        bulletType="bullet",
        leftIndent=18,
    )


def _numbered_list(items: Sequence[str], style: ParagraphStyle) -> ListFlowable:
    return ListFlowable(
        [ListItem(Paragraph(t, style)) for t in items],
        bulletType="1",
        leftIndent=22,
    )


def _styled_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
    col_widths_in: Optional[Sequence[float]],
    bold_first_column: bool,
    styles: Dict[str, ParagraphStyle],
    accent_color,
) -> Table:
    header_paras = [
        Paragraph(f"<b>{h}</b>", styles["TableHeader"]) for h in headers
    ]
    body_rows = []
    for row in rows:
        cells = []
        for ci, value in enumerate(row):
            text = str(value)
            style = styles["TableCellBold"] if (bold_first_column and ci == 0) else styles["TableCell"]
            cells.append(Paragraph(text, style))
        body_rows.append(cells)
    data = [header_paras, *body_rows]

    if col_widths_in:
        col_widths = [w * 72 for w in col_widths_in]
    else:
        content_width = 612 - LEFT_MARGIN - RIGHT_MARGIN
        col_widths = [content_width / len(headers)] * len(headers)

    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.hAlign = "CENTER"

    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), accent_color),
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


def _add_figure(image_path: str, caption: str, width_in: float,
                styles: Dict[str, ParagraphStyle]) -> KeepTogether:
    width_pt = width_in * 72
    if _HAS_PIL:
        with PILImage.open(image_path) as im:
            w_px, h_px = im.size
        aspect = h_px / w_px
        height_pt = width_pt * aspect
    else:
        height_pt = width_pt * 0.6
    img = Image(image_path, width=width_pt, height=height_pt)
    img.hAlign = "CENTER"
    cap = Paragraph(f"<i>{caption}</i>", styles["Caption"])
    return KeepTogether([img, Spacer(1, 4), cap])


def _callout_box(text: str, styles: Dict[str, ParagraphStyle],
                 accent_color) -> Table:
    cell = Paragraph(text, styles["Callout"])
    t = Table([[cell]], colWidths=[6.0 * 72])
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


def _divider(accent_color) -> HRFlowable:
    return HRFlowable(
        width=6.8 * 72, thickness=1, color=accent_color,
        spaceBefore=4, spaceAfter=10,
    )


def _page_callback(header_text: str, accent_color):
    def cb(canvas, doc):
        canvas.saveState()
        w, h = doc.pagesize
        # Header
        canvas.setFont("Helvetica", 9)
        canvas.setFillColor(MUTED_COLOR)
        canvas.drawString(LEFT_MARGIN, h - 40, header_text)
        canvas.setStrokeColor(accent_color)
        canvas.setLineWidth(1)
        canvas.line(LEFT_MARGIN, h - 48, w - RIGHT_MARGIN, h - 48)
        # Footer
        canvas.setStrokeColor(BORDER)
        canvas.setLineWidth(0.75)
        canvas.line(LEFT_MARGIN, 40, w - RIGHT_MARGIN, 40)
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(MUTED_COLOR)
        canvas.drawCentredString(w / 2, 26, f"Page {doc.page}")
        canvas.restoreState()
    return cb


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--spec", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()
    build_report(spec_path=args.spec, output_path=args.out)
    print(f"\n✓ PDF report built from {args.spec} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
