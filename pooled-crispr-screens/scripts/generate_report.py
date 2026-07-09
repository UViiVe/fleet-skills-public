"""PDF report generator for pooled CRISPR-screen results.

Produces a four-section PDF (Introduction, Methods, Results with figures,
Conclusions). Requires reportlab; falls back gracefully if absent.

Example:
    python generate_report.py --in adata_processed.h5ad \\
        --summary perturbation_summary.csv --out-dir results/
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from typing import Optional

try:
    from reportlab.lib import colors
    from reportlab.lib.pagesizes import letter
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import inch
    from reportlab.lib.enums import TA_CENTER
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
        Image, PageBreak, HRFlowable,
    )
    HAS_REPORTLAB = True
except ImportError:
    HAS_REPORTLAB = False


_PRIMARY = "#1B4F72"
_ACCENT  = "#E74C3C"
_DARK    = "#2C3E50"
_LIGHT   = "#F2F3F4"
_MED     = "#BDC3C7"


def _styles():
    s = getSampleStyleSheet()
    s.add(ParagraphStyle("RTitle",   parent=s["Title"],   fontSize=24, textColor=colors.HexColor(_PRIMARY), spaceAfter=6,  alignment=TA_CENTER, fontName="Helvetica-Bold"))
    s.add(ParagraphStyle("RSub",     parent=s["Normal"],  fontSize=12, textColor=colors.HexColor(_DARK),    spaceAfter=20, alignment=TA_CENTER))
    s.add(ParagraphStyle("RSec",     parent=s["Heading1"], fontSize=16, textColor=colors.HexColor(_PRIMARY), spaceBefore=18, spaceAfter=8,  fontName="Helvetica-Bold"))
    s.add(ParagraphStyle("RSub2",    parent=s["Heading2"], fontSize=13, textColor=colors.HexColor(_DARK),    spaceBefore=12, spaceAfter=6,  fontName="Helvetica-Bold"))
    s.add(ParagraphStyle("RBody",    parent=s["Normal"],  fontSize=10, textColor=colors.HexColor(_DARK),    spaceAfter=6, leading=14))
    s.add(ParagraphStyle("RStat",    parent=s["Normal"],  fontSize=28, textColor=colors.HexColor(_ACCENT),  alignment=TA_CENTER, fontName="Helvetica-Bold"))
    s.add(ParagraphStyle("RStatLbl", parent=s["Normal"],  fontSize=9,  textColor=colors.HexColor(_DARK),    alignment=TA_CENTER, spaceAfter=8))
    s.add(ParagraphStyle("RCap",     parent=s["Normal"],  fontSize=9,  textColor=colors.HexColor(_DARK),    alignment=TA_CENTER, spaceAfter=12, fontName="Helvetica-Oblique"))
    s.add(ParagraphStyle("RFoot",    parent=s["Normal"],  fontSize=8,  textColor=colors.HexColor(_MED),     alignment=TA_CENTER))
    return s


def _embed_figure(elements, path, caption, styles, max_width_in=6.0):
    if not os.path.exists(path):
        elements.append(Paragraph(
            f"<i>[Figure not available: {os.path.basename(path)}]</i>",
            styles["RBody"],
        ))
        return
    img = Image(path)
    aspect = img.imageHeight / img.imageWidth
    w = min(max_width_in * inch, 6.5 * inch)
    h = w * aspect
    if h > 5.0 * inch:
        h = 5.0 * inch
        w = h / aspect
    img.drawWidth, img.drawHeight = w, h
    elements.append(img)
    elements.append(Paragraph(caption, styles["RCap"]))


def _table(header, rows, col_widths=None):
    data = [header] + rows
    cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_PRIMARY)),
        ("TEXTCOLOR",  (0, 0), (-1, 0), colors.white),
        ("FONTNAME",   (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",   (0, 0), (-1, 0), 9),
        ("FONTNAME",   (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE",   (0, 1), (-1, -1), 9),
        ("GRID",       (0, 0), (-1, -1), 0.5, colors.HexColor(_MED)),
        ("VALIGN",     (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
    ]
    for i in range(2, len(data), 2):
        cmds.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor(_LIGHT)))
    t = Table(data, colWidths=col_widths, repeatRows=1)
    t.setStyle(TableStyle(cmds))
    return t


def _find_figure(output_dir: str, filename: str) -> str:
    for sub in ("", "figures", "plots"):
        candidate = os.path.join(output_dir, sub, filename) if sub else os.path.join(output_dir, filename)
        if os.path.exists(candidate):
            return candidate
    return os.path.join(output_dir, filename)


def generate_report(
    adata=None,
    perturbation_summary=None,
    de_results=None,
    validation_df=None,
    output_dir: str = "results",
    output_file: Optional[str] = None,
    parameters: Optional[dict] = None,
):
    """Render the PDF; returns the output path or None if reportlab is missing."""
    if not HAS_REPORTLAB:
        print("  reportlab not installed — PDF skipped (text report still available)")
        return None
    if output_file is None:
        output_file = os.path.join(output_dir, "crispr_screen_report.pdf")

    styles = _styles()
    el = []

    n_cells = adata.n_obs if adata is not None else 0
    n_genes = adata.n_vars if adata is not None else 0
    n_pert = n_hits = 0
    hit_rate = 0.0
    if perturbation_summary is not None and len(perturbation_summary) > 0:
        n_pert = len(perturbation_summary)
        if "is_hit" in perturbation_summary.columns:
            n_hits = int(perturbation_summary["is_hit"].sum())
            hit_rate = n_hits / n_pert * 100
    mean_de = 0.0
    if de_results and perturbation_summary is not None and "n_de_genes" in perturbation_summary.columns:
        mean_de = float(perturbation_summary["n_de_genes"].mean())

    # Title
    el += [
        Spacer(1, 0.8 * inch),
        Paragraph("Pooled CRISPR Screen<br/>Analysis Report", styles["RTitle"]),
        Paragraph("Perturb-seq / CROP-seq tiered workflow", styles["RSub"]),
        Paragraph(datetime.now().strftime("%B %d, %Y"), styles["RSub"]),
        Spacer(1, 0.3 * inch),
        HRFlowable(width="80%", thickness=1, color=colors.HexColor(_PRIMARY), spaceBefore=10, spaceAfter=20),
    ]

    # Stat row 1
    stat_rows = [
        [Paragraph(f"<b>{n_cells:,}</b>",         styles["RStat"]),
         Paragraph(f"<b>{n_pert}</b>",            styles["RStat"]),
         Paragraph(f"<b>{n_hits}</b>",            styles["RStat"])],
        [Paragraph("Cells analysed",              styles["RStatLbl"]),
         Paragraph("Perturbations tested",        styles["RStatLbl"]),
         Paragraph("Validated hits",              styles["RStatLbl"])],
    ]
    t = Table(stat_rows, colWidths=[2.2 * inch] * 3)
    t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ALIGN", (0, 0), (-1, -1), "CENTER")]))
    el.append(t)
    el.append(Spacer(1, 0.3 * inch))

    # Stat row 2
    stat_rows2 = [
        [Paragraph(f"<b>{hit_rate:.1f}%</b>",     styles["RStat"]),
         Paragraph(f"<b>{mean_de:.0f}</b>",       styles["RStat"]),
         Paragraph(f"<b>{n_genes:,}</b>",         styles["RStat"])],
        [Paragraph("Hit rate",                    styles["RStatLbl"]),
         Paragraph("Mean DE genes / hit",         styles["RStatLbl"]),
         Paragraph("Genes measured",              styles["RStatLbl"])],
    ]
    t2 = Table(stat_rows2, colWidths=[2.2 * inch] * 3)
    t2.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("ALIGN", (0, 0), (-1, -1), "CENTER")]))
    el.append(t2)

    # Introduction
    el.append(PageBreak())
    el.append(Paragraph("1. Introduction", styles["RSec"]))
    el.append(Paragraph(
        "Pooled CRISPR screens with single-cell RNA-seq readout (Perturb-seq, CROP-seq, "
        "ECCITE-seq) systematically interrogate gene function at scale. Each cell carries "
        "a unique sgRNA; the resulting transcriptome links a perturbation to its "
        "downstream effects.",
        styles["RBody"],
    ))
    el.append(Paragraph(
        "This report summarises a tiered analysis: a fast t-test screen across every "
        "perturbation, a target-gene validation step that filters out non-functional "
        "guides, and (where applicable) a batch-corrected NB GLM tier for the surviving "
        "hits.",
        styles["RBody"],
    ))

    # Methods
    el.append(Spacer(1, 0.2 * inch))
    el.append(Paragraph("2. Methods", styles["RSec"]))
    el.append(Paragraph("2.1 Pipeline", styles["RSub2"]))
    el.append(Paragraph(
        "<b>Tier 1 (screening):</b> Welch's t-test of each perturbation against the "
        "non-targeting control on log-normalized counts.",
        styles["RBody"],
    ))
    el.append(Paragraph(
        "<b>Tier 2 (validation):</b> the target gene's own log2FC must move in the "
        "expected direction (down for CRISPRi/ko, up for CRISPRa) at the expected "
        "magnitude, with p < 0.05.",
        styles["RBody"],
    ))
    el.append(Paragraph(
        "<b>Tier 3 (refinement):</b> glmGamPoi negative-binomial GLM with batch as a "
        "fixed effect on the top validated hits.",
        styles["RBody"],
    ))
    if parameters:
        el.append(Paragraph("2.2 Parameters", styles["RSub2"]))
        el.append(_table(["Parameter", "Value"],
                          [[str(k), str(v)] for k, v in parameters.items()],
                          col_widths=[2.5 * inch, 4.0 * inch]))

    # Results
    el.append(PageBreak())
    el.append(Paragraph("3. Results", styles["RSec"]))
    el.append(Paragraph("3.1 Summary", styles["RSub2"]))
    if n_pert == 0:
        el.append(Paragraph("<b>No perturbation results available.</b>", styles["RBody"]))
    else:
        el.append(Paragraph(
            f"Screening tested <b>{n_pert} perturbations</b> across <b>{n_cells:,} cells</b>, "
            f"identifying <b>{n_hits} validated hits</b> ({hit_rate:.1f}% hit rate).",
            styles["RBody"],
        ))

        el.append(Paragraph("3.2 Top hits", styles["RSub2"]))
        if "is_hit" in perturbation_summary.columns:
            hits = perturbation_summary[perturbation_summary["is_hit"]]
        else:
            hits = perturbation_summary
        sort_col = "outlier_fraction" if "outlier_fraction" in hits.columns else hits.columns[0]
        top = hits.nlargest(10, sort_col) if len(hits) else hits.head(10)
        if len(top):
            rows = []
            for i, (_, r) in enumerate(top.iterrows(), 1):
                g = str(r.get("gene", r.get("perturbation", "N/A")))
                frac = f"{r['outlier_fraction']:.1%}" if "outlier_fraction" in r.index else "N/A"
                nde = str(int(r["n_de_genes"])) if "n_de_genes" in r.index else "N/A"
                rows.append([str(i), g, frac, nde])
            el.append(_table(
                ["Rank", "Perturbation", "Outlier frac.", "DE genes"],
                rows,
                col_widths=[0.5 * inch, 2.5 * inch, 1.2 * inch, 1.2 * inch],
            ))
            el.append(Spacer(1, 0.15 * inch))

        if validation_df is not None and len(validation_df) > 0:
            el.append(Paragraph("3.3 Target validation", styles["RSub2"]))
            n_ok = int(validation_df["target_affected"].sum()) if "target_affected" in validation_df.columns else 0
            n_v = len(validation_df)
            el.append(Paragraph(
                f"<b>{n_ok}/{n_v}</b> perturbations passed on-target validation "
                f"({n_ok/n_v*100:.1f}%).",
                styles["RBody"],
            ))

        el.append(Paragraph("3.4 Visualisations", styles["RSub2"]))
        _embed_figure(el, _find_figure(output_dir, "hit_summary.png"),
                      "<b>Figure 1.</b> Screen overview: top perturbations by outlier fraction (left); "
                      "DE-gene count vs outlier-cell count (right).",
                      styles)
        if len(top):
            g_top = str(top.iloc[0].get("gene", top.iloc[0].get("perturbation", "")))
            if g_top:
                _embed_figure(el, _find_figure(output_dir, f"volcano_{g_top}.png"),
                              f"<b>Figure 2.</b> Volcano plot for top hit {g_top}. "
                              "Red indicates FDR-significant DE genes.",
                              styles)

    # Conclusions
    el.append(PageBreak())
    el.append(Paragraph("4. Conclusions", styles["RSec"]))
    if n_hits > 0:
        el.append(Paragraph(
            f"This screen identified <b>{n_hits} validated hits</b> out of {n_pert} "
            f"perturbations ({hit_rate:.1f}% hit rate). Each hit showed both significant "
            "transcriptional changes and confirmed on-target effects.",
            styles["RBody"],
        ))
    else:
        el.append(Paragraph(
            "No validated hits were identified. Possible causes: weak perturbation effects, "
            "insufficient cells per perturbation, or overly stringent thresholds.",
            styles["RBody"],
        ))

    el.append(Paragraph("4.1 Caveats", styles["RSub2"]))
    el.append(Paragraph(
        "• sgRNA efficiency varies — incomplete knockdown or activation produces false negatives.<br/>"
        "• High-MOI screens contain cells with multiple guides; interpret with care.<br/>"
        "• Paralog or feedback buffering can mask perturbation effects at the RNA level.<br/>"
        "• Hit calling depends on chained thresholds (min DE genes, log2FC, FDR, outlier "
        "fraction); record them in the run manifest.",
        styles["RBody"],
    ))

    el.append(Paragraph("4.2 Suggested next steps", styles["RSub2"]))
    el.append(Paragraph(
        "• Pathway / GSEA enrichment on per-hit DE gene lists.<br/>"
        "• Gene-regulatory-network inference from the perturbed-cell subset.<br/>"
        "• Cluster perturbations by transcriptional signature.<br/>"
        "• Validate top hits orthogonally (qPCR, Western, arrayed re-test).",
        styles["RBody"],
    ))

    el += [
        Spacer(1, 0.5 * inch),
        HRFlowable(width="60%", thickness=0.5, color=colors.HexColor(_MED), spaceBefore=8, spaceAfter=8),
        Paragraph(
            f"Generated by pooled-crispr-screens skill | {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            styles["RFoot"],
        ),
    ]

    doc = SimpleDocTemplate(
        output_file, pagesize=letter,
        leftMargin=0.75 * inch, rightMargin=0.75 * inch,
        topMargin=0.75 * inch, bottomMargin=0.75 * inch,
        title="Pooled CRISPR Screen Analysis Report",
        author="pooled-crispr-screens",
    )
    doc.build(el)
    print(f"  → {os.path.basename(output_file)}")
    return output_file


def main() -> None:
    if not HAS_REPORTLAB:
        sys.exit("Install with: pip install reportlab")
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--out-dir", default="results/")
    args = parser.parse_args()

    try:
        import scanpy as sc
        import pandas as pd
    except ImportError:
        sys.exit("Install with: pip install scanpy pandas")

    adata = sc.read_h5ad(args.in_path)
    summary = pd.read_csv(args.summary)
    out_path = generate_report(adata=adata, perturbation_summary=summary, output_dir=args.out_dir)
    if out_path:
        print(f"✓ PDF written: {out_path}")


if __name__ == "__main__":
    main()
