#!/usr/bin/env python3
"""
validate_output.py — Docling Output Quality Checker

Validates a parsed markdown or JSON file for minimum quality thresholds:
  - Minimum word count (>100)
  - At least one section heading
  - Table count reported
  - Basic encoding sanity check

Usage:
    python validate_output.py --input paper.md
    python validate_output.py --input paper.json --min-words 200

Exit codes:
    0 — All checks pass
    1 — One or more checks failed

Part of the docling-parse fleet skill.
Ekker Lab | Zevo 🦓 | 2026-05-27
"""

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple


def validate_markdown(content: str, min_words: int = 100) -> List[Tuple[str, bool, str]]:
    """
    Validate a parsed markdown string.

    Returns a list of (check_name, passed, detail) tuples.
    """
    checks = []

    # 1. Word count
    words = content.split()
    word_count = len(words)
    wc_pass = word_count >= min_words
    checks.append((
        "Word count",
        wc_pass,
        f"{word_count} words (minimum: {min_words})"
    ))

    # 2. Section headings (# heading or ## heading)
    headings = re.findall(r"^#{1,6}\s+.+", content, re.MULTILINE)
    heading_pass = len(headings) > 0
    checks.append((
        "Section headings",
        heading_pass,
        f"{len(headings)} heading(s) found"
    ))

    # 3. Table count (Markdown pipe tables)
    table_rows = re.findall(r"^\|.+\|$", content, re.MULTILINE)
    # Rough table count: count separator rows |---|---|
    table_separators = re.findall(r"^\|[-:| ]+\|$", content, re.MULTILINE)
    table_count = len(table_separators)
    checks.append((
        "Table count",
        True,  # Tables are informational, not a failure condition
        f"{table_count} table(s) detected"
    ))

    # 4. Encoding sanity (no replacement characters)
    replacement_chars = content.count("\ufffd")
    encoding_pass = replacement_chars < 10
    checks.append((
        "Encoding sanity",
        encoding_pass,
        f"{replacement_chars} Unicode replacement character(s) found"
    ))

    # 5. Not empty / not trivially short
    nonempty_pass = len(content.strip()) > 50
    checks.append((
        "Non-empty content",
        nonempty_pass,
        f"{len(content.strip())} characters in output"
    ))

    return checks


def validate_json(content: str, min_words: int = 100) -> List[Tuple[str, bool, str]]:
    """
    Validate a Docling lossless JSON output.

    Returns a list of (check_name, passed, detail) tuples.
    """
    checks = []

    # 1. Valid JSON
    try:
        data = json.loads(content)
        json_pass = True
        json_detail = "Valid JSON"
    except json.JSONDecodeError as e:
        checks.append(("JSON validity", False, f"JSON parse error: {e}"))
        return checks
    checks.append(("JSON validity", json_pass, json_detail))

    # 2. Expected Docling schema
    schema_name = data.get("schema_name", "")
    schema_pass = "Docling" in schema_name or "docling" in schema_name
    checks.append((
        "Docling schema",
        schema_pass,
        f"schema_name: {schema_name!r}"
    ))

    # 3. Word count (from texts array)
    texts = data.get("texts", [])
    all_text = " ".join(t.get("text", "") for t in texts if isinstance(t, dict))
    word_count = len(all_text.split())
    wc_pass = word_count >= min_words
    checks.append((
        "Word count",
        wc_pass,
        f"{word_count} words across {len(texts)} text elements (minimum: {min_words})"
    ))

    # 4. Section headings
    headers = [
        t for t in texts
        if isinstance(t, dict) and t.get("label", "") in (
            "section_header", "title", "chapter_title"
        )
    ]
    heading_pass = len(headers) > 0
    checks.append((
        "Section headings",
        heading_pass,
        f"{len(headers)} header element(s) found"
    ))

    # 5. Table count (informational)
    tables = data.get("tables", [])
    checks.append((
        "Table count",
        True,
        f"{len(tables)} table(s) in JSON"
    ))

    # 6. Page count
    pages = data.get("pages", {})
    page_count = len(pages)
    checks.append((
        "Page count",
        page_count > 0,
        f"{page_count} page(s) found"
    ))

    return checks


def validate_output(
    input_path: str,
    min_words: int = 100,
) -> Dict:
    """
    Validate a parsed output file (markdown or JSON).

    Args:
        input_path: Path to .md or .json file.
        min_words:  Minimum word count to pass (default: 100).

    Returns:
        Dict with keys: passed (bool), checks (list), summary (str).
    """
    path = Path(input_path)
    if not path.exists():
        print(f"✗ File not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    content = path.read_text(encoding="utf-8", errors="replace")
    suffix = path.suffix.lower()

    # Dispatch to format-specific validator
    if suffix in (".md", ".markdown", ".txt"):
        checks = validate_markdown(content, min_words=min_words)
        fmt = "markdown"
    elif suffix == ".json":
        checks = validate_json(content, min_words=min_words)
        fmt = "json"
    else:
        # Try markdown validation for unknown types
        checks = validate_markdown(content, min_words=min_words)
        fmt = "unknown"

    # --- Print QC summary ---
    print(f"\n{'='*55}")
    print(f"  QC Report: {path.name} ({fmt})")
    print(f"{'='*55}")

    all_pass = True
    for check_name, passed, detail in checks:
        icon = "✓" if passed else "✗"
        status = "PASS" if passed else "FAIL"
        print(f"  {icon} [{status}]  {check_name:25s} {detail}")
        if not passed:
            all_pass = False

    print(f"{'='*55}")
    overall = "PASS ✓" if all_pass else "FAIL ✗"
    print(f"  Overall: {overall}")
    print(f"{'='*55}\n")

    # Extract table count for summary
    table_count = 0
    for name, _, detail in checks:
        if "Table count" in name or "table" in name.lower():
            nums = re.findall(r"\d+", detail)
            if nums:
                table_count = int(nums[0])
                break

    summary = (
        f"{'PASS' if all_pass else 'FAIL'} | "
        + " | ".join(
            f"{name}: {'PASS' if p else 'FAIL'}"
            for name, p, _ in checks
        )
    )

    return {
        "passed": all_pass,
        "checks": checks,
        "table_count": table_count,
        "summary": summary,
        "input_path": input_path,
        "format": fmt,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Validate a Docling-parsed markdown or JSON output file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--input", "-i", required=True,
        help="Path to parsed .md or .json file to validate."
    )
    parser.add_argument(
        "--min-words", type=int, default=100,
        help="Minimum word count to pass (default: 100)."
    )

    args = parser.parse_args()

    result = validate_output(
        input_path=args.input,
        min_words=args.min_words,
    )

    sys.exit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
