# Docling Overview

## What Is Docling?

Docling is an open-source Python library developed by IBM Research (AI for Knowledge team, IBM Research Zurich) and now hosted under the LF AI & Data Foundation. It converts diverse document formats into a unified "DoclingDocument" representation, then exports to multiple AI-friendly formats.

**Key strength:** Advanced PDF understanding — including page layout detection, reading order, table structure extraction, figure classification, formula recognition, and OCR support — without requiring any external API or GPU.

**Tagline:** "Get your documents ready for gen AI"

---

## Architecture

### Pipeline Overview

```
Input (PDF/DOCX/PPTX/...)
        │
        ▼
 DocumentConverter
        │
   ┌────┴────┐
   │ Layout  │  ← Heron model (default, faster) or older layout model
   │ Model   │
   └────┬────┘
        │
   ┌────┴────────┐
   │ Table       │  ← TableFormer model: extracts table structure
   │ Structure   │
   └────┬────────┘
        │
   ┌────┴──────┐
   │ OCR       │  ← Optional: EasyOCR / Tesseract for scanned pages
   │ (optional)│
   └────┬──────┘
        │
        ▼
 DoclingDocument (unified representation)
        │
   ┌────┴──────────────────────────────────────┐
   │              Export formats               │
   │  Markdown │ JSON │ HTML │ DocTags │ Text  │
   └───────────────────────────────────────────┘
```

### Core Components

1. **DocumentConverter** — Entry point. Accepts local paths or URLs. Orchestrates the pipeline.
2. **PdfPipelineOptions** — Configuration object for PDF parsing (OCR, table detection, figure extraction, layout options).
3. **DoclingDocument** — Unified document object. All exports derive from this.
4. **Format-specific parsers** — PDF, DOCX, PPTX, XLSX, HTML, LaTeX, plain text, audio (ASR), images.

### Model Stack (PDF pipeline)

| Model | Purpose | Notes |
|---|---|---|
| **Heron** | Page layout detection | Default since ~2025; faster than original |
| **TableFormer** | Table structure recognition | IBM Research model; handles complex tables |
| **EasyOCR / Tesseract** | OCR for scanned text | Optional; requires `docling[ocr]` |
| **GraniteDocling** (258M) | VLM-based parsing | Alternative pipeline via `--pipeline vlm` |

---

## Installation

```bash
# Core install
pip install docling

# With OCR support
pip install docling[ocr]

# With all extras
pip install "docling[all]"
```

**Requirements:**
- Python 3.10+
- macOS (arm64/x86_64), Linux (x86_64/arm64), Windows
- No GPU required for standard parsing
- GPU (CUDA/MPS) speeds up OCR and VLM pipelines

**First run:** Downloads Heron layout model and TableFormer (~200–400MB total). Cached in `~/.cache/docling/` or equivalent.

---

## Key Features

### Document Understanding
- **Multi-column layouts** — correctly identifies and linearizes reading order
- **Tables** — extracts structure including spanning cells (with known limitations on complex merges)
- **Figures** — detects and can extract figure images
- **Formulas** — identifies mathematical expressions
- **Code blocks** — preserves code formatting
- **Headers/footers** — identifies and can exclude

### Input Formats
- PDF (primary strength)
- DOCX, PPTX, XLSX
- HTML, Markdown, AsciiDoc, LaTeX, plain text
- Images (PNG, JPEG, TIFF, BMP, WEBP) — with OCR
- Audio (MP3, WAV) — ASR transcription
- WebVTT (transcript files)
- USPTO patents, JATS articles, XBRL financial reports

### Output Formats
- **Markdown** — human-readable, LLM-friendly
- **JSON (lossless)** — full fidelity including coordinates, bounding boxes, component metadata
- **HTML** — structured web representation
- **DocTags** — structured markup for AI; preserves document semantics
- **Plain text** — raw text without structure
- **WebVTT** — time-coded text (for transcript formats)

---

## Performance Characteristics

| Scenario | Typical Speed | Notes |
|---|---|---|
| Simple text PDF (10 pages) | 2–5s | Fast mode or standard |
| Complex scientific paper (20 pages) | 10–30s | Full layout + table extraction |
| Large PDF (100 pages) | 60–180s | Full pipeline |
| Scanned PDF + OCR (10 pages) | 30–60s | Depends on OCR model and hardware |
| Batch (100 PDFs) | Variable | Memory-bounded; use `parse_batch.py` |

*Times on Apple Silicon M1/M2 without GPU acceleration.*

---

## Integrations

- **LangChain** — `DoclingLoader` with `ExportType.DOC_CHUNKS` or `ExportType.MARKDOWN`
- **LlamaIndex** — native document loader
- **Haystack** — document converter component
- **CrewAI** — tool integration
- **MCP server** — connect any agent via Model Context Protocol

---

## Fleet Usage Notes

- All fleet Apple Silicon nodes (Zevo, Atlas, Uvy, SARAH, Pip, Galen) can run Docling without GPU
- GPU nodes (Buster, Sparq) recommended only for OCR-heavy workloads
- Docling is local-only — no data leaves the machine (important for sensitive lab data)
- Model cache is per-user: `~/.cache/docling/` — shared if multiple agents run as same OS user

---

## Citation

```bibtex
@techreport{Docling,
  author = {Deep Search Team},
  month = {8},
  title = {Docling Technical Report},
  url = {https://arxiv.org/abs/2408.09869},
  eprint = {2408.09869},
  doi = {10.48550/arXiv.2408.09869},
  version = {1.0.0},
  year = {2024}
}
```

---

*Reference compiled by Zevo 🦓 | 2026-05-27*
