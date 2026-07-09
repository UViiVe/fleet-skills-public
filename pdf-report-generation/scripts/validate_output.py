"""Reopen a generated .pdf and check structural QC.

Verifies file size, page count, text extractability, embedded image
count and DPI (when Pillow is available).

Inputs (CLI):
    --pdf PATH
Outputs:
    JSON summary; exits 0 on pass, 1 if any check fails.

Example invocation
------------------
    python -m scripts.validate_output --pdf report.pdf --out report.validation.json
"""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
import zipfile
from typing import Any, Dict, List

try:
    from pypdf import PdfReader
except ImportError:
    sys.exit("Install with: pip install pypdf")

try:
    from PIL import Image
    _HAS_PIL = True
except ImportError:
    _HAS_PIL = False


CHECKS = {
    "min_file_size_bytes": 5_000,
    "min_pages": 1,
    "min_text_length": 200,
    "min_first_page_text": 10,
    "min_image_dpi": 300,
}


def validate_pdf(path: str) -> Dict[str, Any]:
    failures: List[str] = []
    warnings: List[str] = []

    if not os.path.exists(path):
        return {"path": path, "passes": False,
                "failures": [f"File does not exist: {path}"], "warnings": []}

    file_size = os.path.getsize(path)
    if file_size < CHECKS["min_file_size_bytes"]:
        failures.append(
            f"File size {file_size} B < {CHECKS['min_file_size_bytes']} B"
        )

    try:
        reader = PdfReader(path)
    except Exception as e:
        return {"path": path, "passes": False,
                "failures": [f"Could not open PDF: {e}"], "warnings": []}

    pages = list(reader.pages)
    if len(pages) < CHECKS["min_pages"]:
        failures.append(f"Only {len(pages)} pages")

    full_text = ""
    page_summaries: List[Dict[str, Any]] = []
    for i, page in enumerate(pages, 1):
        try:
            text = page.extract_text() or ""
        except Exception:
            text = ""
        full_text += text + "\n"
        page_summaries.append({"page": i, "text_length": len(text.strip())})

    if pages and len(pages[0].extract_text() or "") < CHECKS["min_first_page_text"]:
        failures.append("First page has effectively no extractable text")

    if len(full_text.strip()) < CHECKS["min_text_length"]:
        failures.append(
            f"Total extracted text {len(full_text.strip())} < "
            f"{CHECKS['min_text_length']} chars"
        )

    image_records = _scan_images(path)
    if _HAS_PIL:
        for img in image_records:
            if img.get("dpi") is not None and img["dpi"] < CHECKS["min_image_dpi"]:
                warnings.append(
                    f"Image {img['name']} DPI={img['dpi']} "
                    f"< {CHECKS['min_image_dpi']}"
                )
    else:
        warnings.append(
            "Pillow not installed — skipping image DPI check. "
            "`pip install Pillow` to enable."
        )

    return {
        "path": path,
        "file_size_bytes": file_size,
        "page_count": len(pages),
        "total_text_length": len(full_text.strip()),
        "pages": page_summaries,
        "image_count": len(image_records),
        "images": image_records,
        "passes": not failures,
        "failures": failures,
        "warnings": warnings,
        "thresholds": CHECKS,
    }


def _scan_images(pdf_path: str) -> List[Dict[str, Any]]:
    """Extract embedded images via the underlying object stream.

    pypdf exposes image objects per-page through `page.images`. Each image's
    DPI isn't recorded in the PDF itself, but Pillow can infer dimensions.
    """
    records: List[Dict[str, Any]] = []
    try:
        reader = PdfReader(pdf_path)
    except Exception:
        return records

    for i, page in enumerate(reader.pages, 1):
        try:
            images = page.images
        except Exception:
            images = []
        for img_obj in images:
            info: Dict[str, Any] = {
                "name": img_obj.name,
                "page": i,
            }
            if _HAS_PIL:
                try:
                    im = Image.open(io.BytesIO(img_obj.data))
                    info["width_px"] = im.width
                    info["height_px"] = im.height
                    dpi = im.info.get("dpi")
                    info["dpi"] = round(min(dpi)) if dpi else None
                    info["format"] = im.format
                except Exception as e:  # pragma: no cover
                    info["error"] = str(e)
            else:
                info["dpi"] = None
            records.append(info)
    return records


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--pdf", required=True)
    p.add_argument("--out", help="Optional JSON output path")
    args = p.parse_args()
    report = validate_pdf(args.pdf)
    text = json.dumps(report, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(text)
        print(f"✓ Validation report written: {args.out}")
    else:
        print(text)
    if report["passes"]:
        print(f"\n✓ Validation passed: {args.pdf}")
        return 0
    print(f"\n✗ Validation FAILED: {args.pdf}")
    for fail in report["failures"]:
        print(f"  - {fail}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
