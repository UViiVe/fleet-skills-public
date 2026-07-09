"""Open a generated .docx and run structural checks.

Verifies file size, paragraph and table counts, total text length, image
embedding, and (when Pillow is available) embedded image DPI. Emits a
structured JSON summary so callers can flag failures programmatically.

Inputs (CLI):
    --docx PATH
Outputs:
    Prints summary, exit code 0 on pass, 1 if any check fails.

Example invocation
------------------
    python -m scripts.validate_output --docx report.docx --out report.validation.json
"""

from __future__ import annotations

import argparse
import io
import json
import os
import sys
from typing import Any, Dict, List

try:
    from docx import Document
except ImportError:
    sys.exit("Install with: pip install python-docx")

try:
    from PIL import Image
    _HAS_PIL = True
except ImportError:
    _HAS_PIL = False


CHECKS = {
    "min_file_size_bytes": 5000,
    "min_paragraphs": 5,
    "min_text_length": 100,
    "min_image_dpi": 300,
}


def validate_docx(path: str) -> Dict[str, Any]:
    """Run every check; return a structured report."""
    failures: List[str] = []
    warnings: List[str] = []

    if not os.path.exists(path):
        return {
            "path": path,
            "passes": False,
            "failures": [f"File does not exist: {path}"],
            "warnings": [],
        }

    file_size = os.path.getsize(path)
    if file_size < CHECKS["min_file_size_bytes"]:
        failures.append(
            f"File size {file_size} B < {CHECKS['min_file_size_bytes']} B "
            "(suggests empty document)"
        )

    doc = Document(path)
    paragraph_count = len(doc.paragraphs)
    table_count = len(doc.tables)
    all_text = "\n".join(p.text for p in doc.paragraphs).strip()
    text_length = len(all_text)

    if paragraph_count < CHECKS["min_paragraphs"]:
        failures.append(
            f"Only {paragraph_count} paragraphs (< {CHECKS['min_paragraphs']})"
        )
    if text_length < CHECKS["min_text_length"]:
        failures.append(
            f"Only {text_length} characters of text (< {CHECKS['min_text_length']})"
        )

    image_records = _scan_images(path)
    if image_records:
        if _HAS_PIL:
            for img in image_records:
                if img["dpi"] is not None and img["dpi"] < CHECKS["min_image_dpi"]:
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
        "paragraph_count": paragraph_count,
        "table_count": table_count,
        "text_length": text_length,
        "image_count": len(image_records),
        "images": image_records,
        "passes": not failures,
        "failures": failures,
        "warnings": warnings,
        "thresholds": CHECKS,
    }


def _scan_images(docx_path: str) -> List[Dict[str, Any]]:
    """List embedded images and (when Pillow available) their dimensions/DPI."""
    import zipfile

    records: List[Dict[str, Any]] = []
    with zipfile.ZipFile(docx_path) as z:
        for name in z.namelist():
            if name.startswith("word/media/"):
                info: Dict[str, Any] = {"name": os.path.basename(name)}
                if _HAS_PIL:
                    with z.open(name) as f:
                        try:
                            img = Image.open(io.BytesIO(f.read()))
                            info["width_px"] = img.width
                            info["height_px"] = img.height
                            dpi = img.info.get("dpi")
                            info["dpi"] = round(min(dpi)) if dpi else None
                            info["format"] = img.format
                        except Exception as e:  # pragma: no cover
                            info["error"] = str(e)
                else:
                    info["dpi"] = None
                records.append(info)
    return records


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--docx", required=True)
    p.add_argument("--out", help="Optional JSON output path")
    args = p.parse_args()

    report = validate_docx(args.docx)
    text = json.dumps(report, indent=2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(text)
        print(f"✓ Validation report written: {args.out}")
    else:
        print(text)

    if report["passes"]:
        print(f"\n✓ Validation passed: {args.docx}")
        return 0
    print(f"\n✗ Validation FAILED: {args.docx}")
    for fail in report["failures"]:
        print(f"  - {fail}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
