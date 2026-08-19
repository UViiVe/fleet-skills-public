---
id: skill_1b5dbc71962ec5e7ce6f0e69697dcacd
name: docling-parse
display-name: Docling Parse
short-description: Parse PDFs into structured markdown, JSON, or HTML using IBM Docling
  — handles tables, figures, multi-column layouts, and scientific papers without any
  API key. Use when ingesting PDFs for RAG pipelines, structured data extraction,
  or document analysis.
category: methodology
keywords: docling, PDF, document parsing, markdown, tables, structured extraction
version: 1.0.1
last-updated: '2026-08-19'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
provenance:
  wraps: IBM Docling
  upstream_source: https://github.com/docling-project/docling
  citation: Auer et al., arXiv:2408.09869, 2024
changelog:
- version: 1.0.1
  date: '2026-07-07'
  event: citation-provenance
  summary: Added provenance block and References section. Vienna public subset.
description: Parse PDFs into structured markdown, JSON, or HTML using IBM Docling
  — handles tables, figures, multi-column layouts, and scientific papers without any
  API key. Use when ingesting PDFs for RAG pipelines, structured data extraction,
  or document analysis.
---

# Docling PDF Parser

Parse PDFs into structured markdown, JSON, or HTML using IBM Docling.
No API key required. Fully self-hosted. Works on Apple Silicon and x86 systems.

---

## When to Use

✅ **Use when:**
- Ingesting PDFs for RAG pipelines or document ingestion pipelines
- Extracting tables from scientific papers
- Converting complex multi-column layouts to clean markdown
- No API key available or desired (fully self-hosted)
- Apple Silicon or CPU-only environments (Docling runs without GPU)
- Parsing PDFs with figures, formulas, or complex layouts

❌ **Avoid when:**
- Simple single-column text-only PDFs (pypdf is faster)
- Scanned PDFs with no text layer and no GPU available (needs OCR model)
- Real-time parsing in latency-sensitive pipelines (Docling is thorough but not instant)
- PDF is password-protected without pre-decryption (use pikepdf first)

---

## Installation

If Docling is already installed in a venv, activate it before running scripts. Otherwise, install fresh:

```bash
python3.13 -m venv ~/.venvs/docling
source ~/.venvs/docling/bin/activate
pip install docling

# Or call Python directly:
~/.venvs/docling/bin/python3 scripts/parse_pdf.py --input paper.pdf --output paper.md

# Optional — OCR support for scanned PDFs (requires model download ~200MB on first use):
pip install docling[ocr]
```

**⚠️ macOS Homebrew note:** Python 3.14 (default `pip`) has a broken pyexpat on macOS 15+. Use Python 3.13 explicitly. Create a venv — do not install into system Python.

Python 3.10+ required. Works on macOS (arm64/x86_64), Linux, Windows.
First run downloads layout model (~200MB) — normal behavior.

---

## Quick Start

```python
from scripts.parse_pdf import parse_pdf

result = parse_pdf(
    input_path="paper.pdf",
    output_format="markdown",   # or "json" or "html"
    output_path="paper.md"
)
print(result["table_count"], "tables extracted")
```

Or via CLI:

```bash
python scripts/parse_pdf.py --input paper.pdf --output paper.md --format markdown

# Batch:
python scripts/parse_batch.py --input-dir ./pdfs/ --output-dir ./parsed/ --format markdown

# Validate output quality:
python scripts/validate_output.py --input paper.md
```

---

## Inputs

| Parameter | Type | Description |
|---|---|---|
| `input_path` | str | Local file path or URL to PDF |
| `output_format` | str | `"markdown"` \| `"json"` \| `"html"` (default: `"markdown"`) |
| `output_path` | str | Where to write output (optional — returns string if omitted) |
| `export_figures` | bool | Extract figure images to `figures/` directory (default: `False`) |
| `export_tables` | bool | Extract tables as separate CSV files to `tables/` directory (default: `False`) |
| `use_ocr` | bool | Enable OCR for scanned PDFs (default: `False`) |
| `fast_mode` | bool | Text-only extraction, skips layout analysis (default: `False`) |
| `table_detection_threshold` | float | Confidence threshold for table detection (default: `0.5`) |

---

## Outputs

| Output | Description |
|---|---|
| Primary file | `.md` / `.json` / `.html` written to `output_path` |
| `figures/` directory | Extracted figure images (if `export_figures=True`) |
| `tables/` directory | CSV files per table (if `export_tables=True`) |
| Metadata dict | `{word_count, table_count, figure_count, page_count, title, parse_time_s}` |

---

## Compute Notes

Runs on any modern CPU (Apple Silicon or x86). For scanned PDFs requiring OCR, a GPU is recommended but not required.

---

## Core API Usage

```python
from docling.document_converter import DocumentConverter
from docling.datamodel.pipeline_options import PdfPipelineOptions

# Basic conversion
converter = DocumentConverter()
result = converter.convert("paper.pdf")
doc = result.document

# Export options
markdown = doc.export_to_markdown()
html = doc.export_to_html()
json_data = doc.export_to_dict()         # for JSON serialization
text = doc.export_to_text()

# With OCR enabled
pipeline_options = PdfPipelineOptions()
pipeline_options.do_ocr = True
converter = DocumentConverter(
    format_options={InputFormat.PDF: PdfFormatOption(pipeline_options=pipeline_options)}
)
```

---

## Common Issues

| Issue | Cause | Fix |
|---|---|---|
| `ImportError: docling not found` | Not installed | `pip install docling` |
| Blank output on scanned PDF | No text layer, OCR not enabled | `pip install docling[ocr]` and pass `use_ocr=True` |
| Tables missing from output | Table detection threshold too low | Pass `table_detection_threshold=0.5` (or lower) |
| Very slow on large PDFs | Full layout analysis | Set `fast_mode=True` for text-only extraction |
| Figure extraction fails | No `figures/` dir | Script creates it automatically — check write permissions |
| Memory error on large PDF | Default batch size too large | Split PDF or reduce `batch_size` parameter |
| Unicode errors in output | Encoding issue in source PDF | Pass `encoding="utf-8-sig"` |
| Wrong column order in tables | Complex spanning cells | Known Docling limitation for merged cells — flag with [VERIFY] |
| Output file empty but no error | PDF is password-protected | Decrypt first with `pikepdf` |
| Title not extracted | No metadata in PDF | Falls back to filename; set title manually in metadata |
| Slow first run | Model download on first use | Docling downloads layout model on first call (~200MB) — normal |
| `parse_batch.py` skips files | Non-PDF files in directory | Only `.pdf` files are processed; others silently skipped |
| `AttributeError: export_to_dict` | Older Docling version | Upgrade: `pip install --upgrade docling` |
| High RAM usage on batch | Many large PDFs in memory | Process in smaller batches; use `parse_batch.py` with `--batch-size` |
| JSON output very large | Lossless JSON includes coordinates | Use `export_to_markdown()` for lighter output |
| OCR model not found | First OCR run, no internet | Pre-download with `docling --download-models` |

---

## Suggested Next Steps

- **paperclip**: deposit parsed document to library
- **pdf-report-generation**: generate PDF from structured data (reverse direction)
- **abstract-lookup**: combine with parsed paper for citation enrichment

---

## References

- **IBM Docling**: Auer C, et al. "Docling Technical Report." *arXiv*. 2024. arXiv:2408.09869. https://arxiv.org/abs/2408.09869
- **GitHub source**: https://github.com/docling-project/docling
- **Documentation**: https://docling-project.github.io/docling/
- See `references/` directory for detailed local reference docs
