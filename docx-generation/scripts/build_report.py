"""Build a styled Word (.docx) report from a JSON spec.

Implements the style guide: US Letter page, Arial body, accent-color
section dividers and table headers, alternating row shading, centered
figures with italic captions, footer with auto-updating page number.

Inputs (CLI):
    --spec PATH (JSON or .py with `spec = ...`)
    --out PATH
Outputs:
    A .docx file at --out.

Example invocation
------------------
    python -m scripts.build_report --spec examples/de_analysis.json --out report.docx
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any, Dict, List, Optional, Sequence

try:
    from docx import Document
    from docx.enum.section import WD_ORIENT
    from docx.enum.table import WD_TABLE_ALIGNMENT
    from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_BREAK
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    from docx.shared import Inches, Pt, RGBColor
except ImportError:
    sys.exit("Install with: pip install python-docx")


FONT = "Arial"
BODY_COLOR = RGBColor(0x2C, 0x2A, 0x26)
HEADING_COLOR = RGBColor(0x11, 0x11, 0x11)
MUTED_COLOR = RGBColor(0x8A, 0x83, 0x78)
LINK_COLOR_HEX = "0563C1"
WHITE = RGBColor(0xFF, 0xFF, 0xFF)


def build_report(
    spec: Optional[Dict[str, Any]] = None,
    spec_path: Optional[str] = None,
    output_path: str = "report.docx",
) -> str:
    """Render `spec` (or read from `spec_path`) into a .docx at `output_path`."""
    if spec is None and spec_path is None:
        raise ValueError("Provide spec or spec_path.")
    if spec is None:
        with open(spec_path) as f:
            spec = json.load(f)
    validate_spec(spec)

    doc = Document()
    _apply_page_setup(doc, spec.get("page_size", "letter"))
    _apply_header_footer(
        doc,
        header_text=spec.get("header_text", spec["title"]),
        accent_color=spec.get("accent_color", "D4A04A"),
    )

    _add_title_block(doc, spec)
    accent = spec.get("accent_color", "D4A04A")

    for section in spec["sections"]:
        if section.get("page_break_before"):
            doc.add_paragraph().add_run().add_break(WD_BREAK.PAGE)

        if "heading" in section and section["heading"]:
            _add_section_heading(doc, section["heading"], level=section.get("level", 1))

        for paragraph in section.get("paragraphs", []) or []:
            _add_body(doc, paragraph)

        for bullet in section.get("bullets", []) or []:
            _add_list_item(doc, bullet, style="List Bullet")

        for item in section.get("numbered", []) or []:
            _add_list_item(doc, item, style="List Number")

        if "table" in section and section["table"]:
            _styled_table(
                doc,
                headers=section["table"]["headers"],
                rows=section["table"]["rows"],
                col_widths=section["table"].get("col_widths"),
                bold_first_column=section["table"].get("bold_first_column", True),
                header_fill=accent,
            )

        if "figure" in section and section["figure"]:
            _add_figure(
                doc,
                image_path=section["figure"]["path"],
                caption=section["figure"]["caption"],
                width_in=section["figure"].get("width_in", 5.0),
            )

        if "callout" in section and section["callout"]:
            _callout_box(
                doc,
                body_text=section["callout"]["text"],
                bold_prefix=section["callout"].get("bold_prefix"),
                accent=accent,
            )

        if section.get("divider"):
            _add_divider(doc, accent_color=accent)

    parent = os.path.dirname(os.path.abspath(output_path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    doc.save(output_path)
    print(f"✓ docx written: {output_path} "
          f"({os.path.getsize(output_path)} bytes)")
    return output_path


def validate_spec(spec: Dict[str, Any]) -> None:
    """Reject malformed specs before building."""
    if not isinstance(spec, dict):
        raise TypeError("Spec must be a dict.")
    if not spec.get("title"):
        raise ValueError("Spec is missing `title`.")
    sections = spec.get("sections")
    if not isinstance(sections, list) or not sections:
        raise ValueError("Spec is missing a non-empty `sections` array.")

    for i, section in enumerate(sections):
        if not isinstance(section, dict):
            raise ValueError(f"Section {i} is not a dict.")
        if not any(
            section.get(k)
            for k in ("heading", "paragraphs", "table", "figure",
                      "callout", "bullets", "numbered")
        ):
            raise ValueError(f"Section {i} has no content.")
        if "table" in section and section["table"]:
            t = section["table"]
            if not t.get("headers") or not isinstance(t.get("rows"), list):
                raise ValueError(f"Section {i} table is malformed.")
            ncols = len(t["headers"])
            for r in t["rows"]:
                if len(r) != ncols:
                    raise ValueError(
                        f"Section {i} table row {r!r} length != {ncols}."
                    )
        if "figure" in section and section["figure"]:
            path = section["figure"].get("path")
            if not path or not os.path.exists(path):
                raise FileNotFoundError(f"Figure {path!r} not found (section {i}).")

    accent = spec.get("accent_color", "D4A04A")
    if not (isinstance(accent, str) and len(accent) == 6):
        raise ValueError(f"accent_color {accent!r} must be a 6-char hex string.")


# ---- helpers --------------------------------------------------------------


def _apply_page_setup(doc, page_size: str) -> None:
    s = doc.sections[0]
    if page_size.lower() == "a4":
        s.page_width = Inches(8.27)
        s.page_height = Inches(11.69)
    else:
        s.page_width = Inches(8.5)
        s.page_height = Inches(11.0)
    s.top_margin = Inches(1)
    s.bottom_margin = Inches(1)
    s.left_margin = Inches(1)
    s.right_margin = Inches(1)

    style = doc.styles["Normal"]
    style.font.name = FONT
    style.font.size = Pt(11)
    style.font.color.rgb = BODY_COLOR


def _apply_header_footer(doc, header_text: str, accent_color: str) -> None:
    section = doc.sections[0]
    header = section.header
    header.is_linked_to_previous = False
    hp = header.paragraphs[0]
    hp.text = ""
    run = hp.add_run(header_text)
    run.font.name = FONT
    run.font.size = Pt(9)
    run.font.color.rgb = MUTED_COLOR
    _paragraph_bottom_border(hp, color=accent_color, size="6")

    footer = section.footer
    footer.is_linked_to_previous = False
    fp = footer.paragraphs[0]
    fp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    _paragraph_top_border(fp, color="D5CFC5", size="4")
    label = fp.add_run("Page ")
    label.font.name = FONT
    label.font.size = Pt(8)
    label.font.color.rgb = MUTED_COLOR
    _add_page_field(fp)


def _add_page_field(paragraph) -> None:
    for tag, attr in (("begin", None), ("instr", " PAGE "), ("end", None)):
        run = paragraph.add_run()
        if tag == "instr":
            instr = OxmlElement("w:instrText")
            instr.set(qn("xml:space"), "preserve")
            instr.text = attr
            run._element.append(instr)
        else:
            ch = OxmlElement("w:fldChar")
            ch.set(qn("w:fldCharType"), tag)
            run._element.append(ch)


def _paragraph_bottom_border(paragraph, color: str, size: str) -> None:
    pPr = paragraph._element.get_or_add_pPr()
    pBdr = OxmlElement("w:pBdr")
    el = OxmlElement("w:bottom")
    el.set(qn("w:val"), "single")
    el.set(qn("w:sz"), size)
    el.set(qn("w:color"), color)
    el.set(qn("w:space"), "4")
    pBdr.append(el)
    pPr.append(pBdr)


def _paragraph_top_border(paragraph, color: str, size: str) -> None:
    pPr = paragraph._element.get_or_add_pPr()
    pBdr = OxmlElement("w:pBdr")
    el = OxmlElement("w:top")
    el.set(qn("w:val"), "single")
    el.set(qn("w:sz"), size)
    el.set(qn("w:color"), color)
    el.set(qn("w:space"), "4")
    pBdr.append(el)
    pPr.append(pBdr)


def _add_title_block(doc, spec: Dict[str, Any]) -> None:
    accent = spec.get("accent_color", "D4A04A")
    r, g, b = int(accent[0:2], 16), int(accent[2:4], 16), int(accent[4:6], 16)

    title = doc.add_paragraph()
    t = title.add_run(spec["title"])
    t.font.name = FONT
    t.font.size = Pt(26)
    t.bold = True
    t.font.color.rgb = HEADING_COLOR
    title.paragraph_format.space_before = Pt(20)
    title.paragraph_format.space_after = Pt(4)

    if spec.get("subtitle"):
        sub = doc.add_paragraph()
        s = sub.add_run(spec["subtitle"])
        s.font.name = FONT
        s.font.size = Pt(11)
        s.font.color.rgb = RGBColor(r, g, b)
        sub.paragraph_format.space_after = Pt(2)

    if spec.get("attribution"):
        attr = doc.add_paragraph()
        a = attr.add_run(spec["attribution"])
        a.font.name = FONT
        a.font.size = Pt(10)
        a.italic = True
        a.font.color.rgb = MUTED_COLOR
        attr.paragraph_format.space_after = Pt(16)


def _add_section_heading(doc, text: str, level: int) -> None:
    heading = doc.add_heading(text, level=level)
    for run in heading.runs:
        run.font.name = FONT
        run.font.color.rgb = HEADING_COLOR
    space_before, space_after = (24, 10) if level == 1 else (16, 8)
    heading.paragraph_format.space_before = Pt(space_before)
    heading.paragraph_format.space_after = Pt(space_after)


def _add_body(doc, text: str) -> None:
    para = doc.add_paragraph()
    run = para.add_run(text)
    run.font.name = FONT
    run.font.size = Pt(11)
    run.font.color.rgb = BODY_COLOR
    para.paragraph_format.space_after = Pt(8)
    para.paragraph_format.line_spacing = Pt(15)


def _add_list_item(doc, text: str, style: str) -> None:
    para = doc.add_paragraph(style=style)
    run = para.add_run(text)
    run.font.name = FONT
    run.font.size = Pt(11)
    run.font.color.rgb = BODY_COLOR


def _styled_table(
    doc,
    headers: Sequence[str],
    rows: Sequence[Sequence[Any]],
    col_widths: Optional[Sequence[float]],
    bold_first_column: bool,
    header_fill: str,
) -> None:
    table = doc.add_table(rows=1 + len(rows), cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER

    if col_widths:
        for i, w in enumerate(col_widths):
            for row in table.rows:
                row.cells[i].width = Inches(w)

    for i, label in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(label)
        run.font.name = FONT
        run.font.size = Pt(10)
        run.bold = True
        run.font.color.rgb = WHITE
        _apply_shading(cell, header_fill)

    for ri, row_data in enumerate(rows):
        for ci, value in enumerate(row_data):
            cell = table.rows[ri + 1].cells[ci]
            cell.text = ""
            run = cell.paragraphs[0].add_run(str(value))
            run.font.name = FONT
            run.font.size = Pt(10)
            run.font.color.rgb = BODY_COLOR
            if bold_first_column and ci == 0:
                run.bold = True
            if ri % 2 == 1:
                _apply_shading(cell, "F9F7F3")

    _apply_borders(table, color="D5CFC5", size="4")
    for row in table.rows:
        for cell in row.cells:
            _apply_cell_padding(cell, top=80, bottom=80, left=120, right=120)


def _apply_shading(cell, color_hex: str) -> None:
    shading = OxmlElement("w:shd")
    shading.set(qn("w:fill"), color_hex)
    shading.set(qn("w:val"), "clear")
    cell._element.get_or_add_tcPr().append(shading)


def _apply_borders(table, color: str, size: str) -> None:
    tbl = table._tbl
    tblPr = tbl.tblPr if tbl.tblPr is not None else OxmlElement("w:tblPr")
    borders = OxmlElement("w:tblBorders")
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), size)
        el.set(qn("w:color"), color)
        el.set(qn("w:space"), "0")
        borders.append(el)
    tblPr.append(borders)


def _apply_cell_padding(cell, top: int, bottom: int, left: int, right: int) -> None:
    tcPr = cell._element.get_or_add_tcPr()
    tcMar = OxmlElement("w:tcMar")
    for edge, val in (("top", top), ("bottom", bottom),
                      ("start", left), ("end", right)):
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:w"), str(val))
        el.set(qn("w:type"), "dxa")
        tcMar.append(el)
    tcPr.append(tcMar)


def _add_figure(doc, image_path: str, caption: str, width_in: float) -> None:
    img_para = doc.add_paragraph()
    img_para.alignment = WD_ALIGN_PARAGRAPH.CENTER
    img_para.add_run().add_picture(image_path, width=Inches(width_in))
    img_para.paragraph_format.space_before = Pt(12)
    img_para.paragraph_format.space_after = Pt(4)

    cap = doc.add_paragraph()
    cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = cap.add_run(caption)
    run.font.name = FONT
    run.font.size = Pt(9)
    run.italic = True
    run.font.color.rgb = MUTED_COLOR
    cap.paragraph_format.space_after = Pt(14)


def _callout_box(doc, body_text: str, bold_prefix: Optional[str], accent: str) -> None:
    table = doc.add_table(rows=1, cols=1)
    cell = table.rows[0].cells[0]
    cell.text = ""
    p = cell.paragraphs[0]
    if bold_prefix:
        b = p.add_run(bold_prefix + " ")
        b.font.name = FONT
        b.font.size = Pt(10)
        b.bold = True
        b.font.color.rgb = BODY_COLOR
    t = p.add_run(body_text)
    t.font.name = FONT
    t.font.size = Pt(10)
    t.font.color.rgb = BODY_COLOR
    _apply_shading(cell, "FAF9F3")

    tcBorders = OxmlElement("w:tcBorders")
    edges = {
        "start": (accent, "18"),
        "top": ("D5CFC5", "4"),
        "bottom": ("D5CFC5", "4"),
        "end": ("D5CFC5", "4"),
    }
    for edge, (color, sz) in edges.items():
        el = OxmlElement(f"w:{edge}")
        el.set(qn("w:val"), "single")
        el.set(qn("w:sz"), sz)
        el.set(qn("w:color"), color)
        el.set(qn("w:space"), "0")
        tcBorders.append(el)
    cell._element.get_or_add_tcPr().append(tcBorders)
    _apply_cell_padding(cell, top=120, bottom=120, left=160, right=160)


def _add_divider(doc, accent_color: str) -> None:
    para = doc.add_paragraph()
    _paragraph_bottom_border(para, color=accent_color, size="6")
    para.paragraph_format.space_before = Pt(4)
    para.paragraph_format.space_after = Pt(10)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--spec", required=True, help="Path to JSON spec")
    p.add_argument("--out", required=True, help="Output .docx path")
    args = p.parse_args()
    build_report(spec_path=args.spec, output_path=args.out)
    print(f"\n✓ Word report built from {args.spec} -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
