#!/usr/bin/env python3
"""
parse_batch.py — Docling Batch PDF Parser

Process a directory of PDFs and write one output file per PDF.
Supports markdown, JSON, and HTML output formats.

Usage:
    python parse_batch.py --input-dir ./pdfs/ --output-dir ./parsed/ --format markdown
    python parse_batch.py --input-dir ./pdfs/ --output-dir ./json/ --format json --use-ocr
    python parse_batch.py --input-dir ./pdfs/ --output-dir ./out/ --batch-size 5

Only .pdf files in the input directory are processed.
Files that fail are logged and counted but do not stop the batch.

Part of the docling-parse fleet skill.
Ekker Lab | Zevo 🦓 | 2026-05-27
"""

import argparse
import json
import sys
import time
from pathlib import Path
from typing import List


def parse_batch(
    input_dir: str,
    output_dir: str,
    output_format: str = "markdown",
    use_ocr: bool = False,
    fast_mode: bool = False,
    export_tables: bool = False,
    export_figures: bool = False,
    batch_size: int = 0,
    table_detection_threshold: float = 0.5,
) -> dict:
    """
    Parse all PDFs in a directory using IBM Docling.

    Args:
        input_dir:    Directory containing .pdf files.
        output_dir:   Directory to write output files.
        output_format: "markdown", "json", or "html".
        use_ocr:      Enable OCR for scanned PDFs.
        fast_mode:    Text-only extraction, skips layout analysis.
        export_tables: Extract tables as CSV files alongside each output.
        export_figures: Extract figures alongside each output.
        batch_size:   Max PDFs to process per batch (0 = all at once).
        table_detection_threshold: Table detection confidence threshold.

    Returns:
        Dict with n_success, n_failed, failed_files, total_time_s, output_dir.
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

    # --- Discover PDFs ---
    input_path = Path(input_dir)
    if not input_path.exists():
        raise FileNotFoundError(f"Input directory not found: {input_dir}")

    pdf_files: List[Path] = sorted(input_path.glob("*.pdf"))
    if not pdf_files:
        print(f"⚠ No .pdf files found in {input_dir}", file=sys.stderr)
        return {
            "n_success": 0, "n_failed": 0, "failed_files": [],
            "total_time_s": 0, "output_dir": output_dir,
        }

    print(f"Found {len(pdf_files)} PDF(s) in {input_dir}")

    # --- Prepare output directory ---
    out_path = Path(output_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # --- Configure Docling once (reuse converter for all files) ---
    pipeline_options = PdfPipelineOptions()
    pipeline_options.do_ocr = use_ocr
    if fast_mode:
        pipeline_options.do_table_structure = False

    converter = DocumentConverter(
        format_options={
            InputFormat.PDF: PdfFormatOption(pipeline_options=pipeline_options)
        }
    )

    # --- Extension map ---
    ext_map = {"markdown": ".md", "json": ".json", "html": ".html"}
    ext = ext_map.get(output_format, ".md")

    # --- Process ---
    n_success = 0
    n_failed = 0
    failed_files: List[str] = []
    start_total = time.time()

    # Apply batch_size chunking
    batches = [pdf_files]
    if batch_size > 0:
        batches = [
            pdf_files[i : i + batch_size]
            for i in range(0, len(pdf_files), batch_size)
        ]

    for batch_idx, batch in enumerate(batches):
        if len(batches) > 1:
            print(
                f"  Batch {batch_idx + 1}/{len(batches)}: {len(batch)} files"
            )

        for pdf_path in batch:
            stem = pdf_path.stem
            out_file = out_path / f"{stem}{ext}"
            t0 = time.time()

            try:
                result = converter.convert(str(pdf_path))
                doc = result.document

                # Export
                fmt = output_format.lower()
                if fmt == "markdown":
                    content = doc.export_to_markdown()
                elif fmt == "json":
                    content = json.dumps(
                        doc.export_to_dict(), indent=2, ensure_ascii=False
                    )
                elif fmt == "html":
                    content = doc.export_to_html()
                else:
                    content = doc.export_to_markdown()

                out_file.write_text(content, encoding="utf-8")

                # Optional: tables
                if export_tables:
                    try:
                        import pandas as pd
                        tables_dir = out_path / f"{stem}_tables"
                        tables_dir.mkdir(exist_ok=True)
                        for i, tbl in enumerate(doc.tables):
                            df = tbl.export_to_dataframe()
                            df.to_csv(
                                str(tables_dir / f"table_{i+1}.csv"), index=False
                            )
                    except Exception as te:
                        print(
                            f"  ⚠ Table export failed for {pdf_path.name}: {te}",
                            file=sys.stderr,
                        )

                # Optional: figures
                if export_figures:
                    try:
                        figs_dir = out_path / f"{stem}_figures"
                        figs_dir.mkdir(exist_ok=True)
                        for i, pic in enumerate(doc.pictures):
                            if hasattr(pic, "image") and pic.image:
                                img_p = figs_dir / f"figure_{i+1}.png"
                                if hasattr(pic.image, "pil_image"):
                                    pic.image.pil_image.save(str(img_p))
                    except Exception as fe:
                        print(
                            f"  ⚠ Figure export failed for {pdf_path.name}: {fe}",
                            file=sys.stderr,
                        )

                elapsed = round(time.time() - t0, 2)
                # Table count for progress line
                try:
                    tbl_count = len(list(doc.tables))
                except Exception:
                    tbl_count = 0
                page_count = len(result.document.pages) if hasattr(result.document, "pages") else "?"
                print(
                    f"  ✓ {pdf_path.name} → {out_file.name} "
                    f"({page_count}p, {tbl_count}t, {elapsed}s)"
                )
                n_success += 1

            except Exception as exc:
                elapsed = round(time.time() - t0, 2)
                print(
                    f"  ✗ FAILED: {pdf_path.name} — {exc} ({elapsed}s)",
                    file=sys.stderr,
                )
                failed_files.append(str(pdf_path))
                n_failed += 1

    total_time_s = round(time.time() - start_total, 2)
    print(
        f"✓ Batch complete: {n_success} PDFs parsed, {n_failed} failed → {output_dir} "
        f"({total_time_s}s total)"
    )

    return {
        "n_success": n_success,
        "n_failed": n_failed,
        "failed_files": failed_files,
        "total_time_s": total_time_s,
        "output_dir": output_dir,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Batch parse all PDFs in a directory using IBM Docling.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input-dir", "-i", required=True,
        help="Directory containing .pdf files to process."
    )
    parser.add_argument(
        "--output-dir", "-o", required=True,
        help="Directory to write output files."
    )
    parser.add_argument(
        "--format", "-f", default="markdown",
        choices=["markdown", "json", "html"],
        help="Output format for all files (default: markdown)."
    )
    parser.add_argument(
        "--use-ocr", action="store_true",
        help="Enable OCR for scanned PDFs (requires: pip install docling[ocr])."
    )
    parser.add_argument(
        "--fast-mode", action="store_true",
        help="Text-only extraction — skips layout analysis."
    )
    parser.add_argument(
        "--export-tables", action="store_true",
        help="Extract tables as CSV files to <stem>_tables/ next to output."
    )
    parser.add_argument(
        "--export-figures", action="store_true",
        help="Extract figures to <stem>_figures/ next to output."
    )
    parser.add_argument(
        "--batch-size", type=int, default=0,
        help="Process in batches of N PDFs (0 = all at once, default: 0)."
    )
    parser.add_argument(
        "--table-threshold", type=float, default=0.5,
        help="Table detection confidence threshold (default: 0.5)."
    )

    args = parser.parse_args()

    results = parse_batch(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        output_format=args.format,
        use_ocr=args.use_ocr,
        fast_mode=args.fast_mode,
        export_tables=args.export_tables,
        export_figures=args.export_figures,
        batch_size=args.batch_size,
        table_detection_threshold=args.table_threshold,
    )

    if results["n_failed"] > 0:
        print("\nFailed files:", file=sys.stderr)
        for f in results["failed_files"]:
            print(f"  - {f}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
