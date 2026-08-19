#!/usr/bin/env python3
"""
parse_pdf.py — Docling PDF Parser (single file)

Parse a PDF (local path or URL) into markdown, JSON, or HTML using IBM Docling.
Returns a metadata dict with word_count, table_count, figure_count, page_count,
title, and parse_time_s.

Usage:
    python parse_pdf.py --input paper.pdf --output paper.md --format markdown
    python parse_pdf.py --input https://arxiv.org/pdf/2408.09869 --format json --output paper.json
    python parse_pdf.py --input scan.pdf --output scan.md --use-ocr

Part of the docling-parse fleet skill.
Ekker Lab | Zevo 🦓 | 2026-05-27
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Dict, Optional

# Force CPU device on Apple Silicon — Docling's layout model uses float64 which
# MPS (Metal) does not support. CPU fallback is reliable and fast enough for fleet use.
os.environ["DOCLING_DEVICE"] = "cpu"
try:
    import torch
    if torch.backends.mps.is_available():
        torch.set_default_device("cpu")
except ImportError:
    pass


def parse_pdf(
    input_path: str,
    output_format: str = "markdown",
    output_path: Optional[str] = None,
    export_figures: bool = False,
    export_tables: bool = False,
    use_ocr: bool = False,
    fast_mode: bool = False,
    table_detection_threshold: float = 0.5,
) -> Dict[str, Any]:
    """
    Parse a PDF (local path or URL) using IBM Docling.

    Args:
        input_path:               Local file path or URL to PDF.
        output_format:            "markdown", "json", or "html". Default: "markdown".
        output_path:              Where to write output. If None, returns content string only.
        export_figures:           Extract figure images to figures/ directory.
        export_tables:            Extract tables as CSV files to tables/ directory.
        use_ocr:                  Enable OCR for scanned PDFs (requires docling[ocr]).
        fast_mode:                Text-only extraction; skips layout analysis. Faster but less structured.
        table_detection_threshold: Confidence threshold for table detection (0.0–1.0).

    Returns:
        Dict with keys: word_count, table_count, figure_count, page_count, title,
                        parse_time_s, output_path (if written), content (str, always).
    """
    # --- Import check ---
    try:
        from docling.document_converter import DocumentConverter
        from docling.datamodel.pipeline_options import PdfPipelineOptions
        from docling.datamodel.base_models import InputFormat
        from docling.document_converter import PdfFormatOption
    except ImportError as exc:
        print(
            "✗ Docling not installed. Run: pip install docling",
            file=sys.stderr,
        )
        raise ImportError("Docling is required. Install with: pip install docling") from exc

    start_time = time.time()

    # --- Configure pipeline ---
    pipeline_options = PdfPipelineOptions()
    pipeline_options.do_ocr = use_ocr

    # Table detection threshold (PdfPipelineOptions may support this in newer versions)
    # Fall back gracefully if attribute not present
    if hasattr(pipeline_options, "table_structure_options"):
        if hasattr(pipeline_options.table_structure_options, "do_cell_matching"):
            pipeline_options.table_structure_options.do_cell_matching = (
                table_detection_threshold >= 0.5
            )

    if fast_mode:
        # Disable expensive layout models for fast text extraction
        pipeline_options.do_table_structure = False
        if hasattr(pipeline_options, "do_figure_classification"):
            pipeline_options.do_figure_classification = False

    # --- Build converter ---
    converter = DocumentConverter(
        format_options={
            InputFormat.PDF: PdfFormatOption(pipeline_options=pipeline_options)
        }
    )

    # --- Convert ---
    try:
        result = converter.convert(input_path)
    except Exception as exc:
        print(f"✗ Conversion failed for {input_path}: {exc}", file=sys.stderr)
        raise

    doc = result.document
    parse_time_s = round(time.time() - start_time, 2)

    # --- Export content ---
    output_format = output_format.lower()
    if output_format == "markdown":
        content = doc.export_to_markdown()
    elif output_format == "json":
        content = json.dumps(doc.export_to_dict(), indent=2, ensure_ascii=False)
    elif output_format == "html":
        content = doc.export_to_html()
    else:
        raise ValueError(f"Unknown output_format: {output_format!r}. Use 'markdown', 'json', or 'html'.")

    # --- Write output ---
    if output_path:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        Path(output_path).write_text(content, encoding="utf-8")

    # --- Gather metadata ---
    # Page count
    page_count = len(result.document.pages) if hasattr(result.document, "pages") else 0

    # Table count
    table_count = 0
    try:
        table_count = len(list(doc.tables))
    except Exception:
        pass

    # Figure count
    figure_count = 0
    try:
        figure_count = len(list(doc.pictures))
    except Exception:
        pass

    # Word count (from markdown export — reliable baseline)
    words = content.split() if output_format != "json" else doc.export_to_markdown().split()
    word_count = len(words)

    # Title — try metadata, fall back to filename
    title = None
    try:
        meta = doc.origin
        if meta and hasattr(meta, "filename") and meta.filename:
            title = Path(meta.filename).stem
    except Exception:
        pass
    if not title and input_path:
        title = Path(input_path).stem if not input_path.startswith("http") else input_path.split("/")[-1]

    # --- Optional: export figures ---
    if export_figures:
        figures_dir = Path(output_path).parent / "figures" if output_path else Path("figures")
        figures_dir.mkdir(parents=True, exist_ok=True)
        try:
            for i, picture in enumerate(doc.pictures):
                if hasattr(picture, "image") and picture.image:
                    img_path = figures_dir / f"figure_{i + 1}.png"
                    if hasattr(picture.image, "save"):
                        picture.image.pil_image.save(str(img_path))
        except Exception as exc:
            print(f"  ⚠ Figure export partially failed: {exc}", file=sys.stderr)

    # --- Optional: export tables as CSV ---
    if export_tables:
        try:
            import pandas as pd
            tables_dir = Path(output_path).parent / "tables" if output_path else Path("tables")
            tables_dir.mkdir(parents=True, exist_ok=True)
            for i, table in enumerate(doc.tables):
                try:
                    df = table.export_to_dataframe()
                    csv_path = tables_dir / f"table_{i + 1}.csv"
                    df.to_csv(str(csv_path), index=False)
                except Exception as exc:
                    print(f"  ⚠ Table {i+1} export failed: {exc}", file=sys.stderr)
        except ImportError:
            print("  ⚠ pandas not installed — table CSV export skipped.", file=sys.stderr)

    metadata = {
        "word_count": word_count,
        "table_count": table_count,
        "figure_count": figure_count,
        "page_count": page_count,
        "title": title,
        "parse_time_s": parse_time_s,
        "output_path": output_path,
        "content": content,
    }

    # --- Summary line ---
    out_label = output_path or "(string only)"
    print(
        f"✓ PDF parsed: {page_count} pages, {table_count} tables → {out_label} "
        f"({parse_time_s}s)"
    )

    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Parse a PDF into markdown, JSON, or HTML using IBM Docling.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Local path or URL to input PDF."
    )
    parser.add_argument(
        "--output", "-o", default=None,
        help="Output file path. If omitted, content is printed to stdout."
    )
    parser.add_argument(
        "--format", "-f", default="markdown",
        choices=["markdown", "json", "html"],
        help="Output format (default: markdown)."
    )
    parser.add_argument(
        "--use-ocr", action="store_true",
        help="Enable OCR for scanned PDFs (requires: pip install docling[ocr])."
    )
    parser.add_argument(
        "--fast-mode", action="store_true",
        help="Fast text-only extraction — skips layout analysis. Faster but less structured."
    )
    parser.add_argument(
        "--export-figures", action="store_true",
        help="Extract figures to figures/ directory next to output file."
    )
    parser.add_argument(
        "--export-tables", action="store_true",
        help="Extract tables as CSV files to tables/ directory next to output file."
    )
    parser.add_argument(
        "--table-threshold", type=float, default=0.5,
        help="Table detection confidence threshold (0.0–1.0, default: 0.5)."
    )

    args = parser.parse_args()

    meta = parse_pdf(
        input_path=args.input,
        output_format=args.format,
        output_path=args.output,
        export_figures=args.export_figures,
        export_tables=args.export_tables,
        use_ocr=args.use_ocr,
        fast_mode=args.fast_mode,
        table_detection_threshold=args.table_threshold,
    )

    if not args.output:
        print(meta["content"])

    # Print metadata summary to stderr so content can be piped cleanly
    print(
        f"  word_count={meta['word_count']} | figure_count={meta['figure_count']} | "
        f"title={meta['title']!r}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
