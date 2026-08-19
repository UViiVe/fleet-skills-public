# Docling Output Formats

Detailed reference for all export formats supported by Docling's DoclingDocument.

---

## Overview

All Docling exports derive from the `DoclingDocument` object returned by `converter.convert(source).document`.

```python
from docling.document_converter import DocumentConverter

converter = DocumentConverter()
result = converter.convert("paper.pdf")
doc = result.document   # DoclingDocument — all exports from here
```

---

## Markdown

**Method:** `doc.export_to_markdown()`
**File extension:** `.md`

### Characteristics
- Human-readable, LLM-friendly
- Preserves headings, lists, tables (as Markdown tables), code blocks
- Reading order is linearized correctly from multi-column layouts
- Figures represented as `![Figure N](...)` placeholders (or omitted if not extracted)
- Tables rendered as GFM (GitHub Flavored Markdown) pipe tables

### Output Example

```markdown
## Abstract

This paper presents a novel approach to document parsing...

## Methods

### Data Collection

| Sample | Count | Notes |
|--------|-------|-------|
| Train  | 1000  | Filtered |
| Test   | 200   | Held out |

### Model Architecture

```python
model = TransformerEncoder(...)
```
```

### API Options

```python
# Basic
md = doc.export_to_markdown()

# Strict mode (only confirmed text, no inferred content)
md = doc.export_to_markdown(strict_text=True)

# With image references
md = doc.export_to_markdown(image_mode=ImageRefMode.EMBEDDED)
```

### Use When
- Feeding into LLMs or RAG pipelines
- Human review of extracted content
- Library deposit (paperclip skill)
- RDD curation input

---

## JSON (Lossless)

**Method:** `doc.export_to_dict()` → then `json.dumps()`
**File extension:** `.json`

### Characteristics
- Full fidelity — preserves ALL information from the parse
- Includes bounding boxes, page coordinates, confidence scores
- Tables stored as structured objects with row/column metadata
- Figures stored with location data
- Larger file size than markdown (10–50x)

### Top-Level Structure

```json
{
  "schema_name": "DoclingDocument",
  "version": "1.3.0",
  "name": "paper",
  "origin": {
    "mimetype": "application/pdf",
    "filename": "paper.pdf"
  },
  "furniture": { ... },
  "body": {
    "children": [
      {
        "cref": "#/texts/0"
      },
      ...
    ]
  },
  "texts": [
    {
      "self_ref": "#/texts/0",
      "label": "section_header",
      "text": "Abstract",
      "prov": [
        {
          "page_no": 1,
          "bbox": {"l": 72.0, "t": 742.0, "r": 540.0, "b": 728.0},
          "charspan": [0, 8]
        }
      ]
    }
  ],
  "tables": [
    {
      "self_ref": "#/tables/0",
      "data": {
        "table_cells": [
          {"row_span": 1, "col_span": 1, "start_row_offset_idx": 0,
           "end_row_offset_idx": 1, "start_col_offset_idx": 0,
           "end_col_offset_idx": 1, "text": "Sample", "column_header": true}
        ],
        "num_rows": 3,
        "num_cols": 3
      }
    }
  ],
  "pictures": [ ... ],
  "pages": {
    "1": {
      "size": {"width": 612.0, "height": 792.0},
      "image": null
    }
  }
}
```

### Key Fields

| Field | Description |
|---|---|
| `texts` | All text elements with labels and bounding boxes |
| `tables` | Structured table data with cell-level detail |
| `pictures` | Figure metadata and location |
| `pages` | Page dimensions, images (if extracted) |
| `body.children` | Document tree — reading order via `cref` references |
| `label` | Element type: `section_header`, `paragraph`, `list_item`, `code`, `formula`, etc. |
| `prov[].bbox` | Bounding box: `l` (left), `t` (top), `r` (right), `b` (bottom) in points |
| `prov[].page_no` | 1-indexed page number |

### Use When
- Building vector stores with coordinate-level provenance
- UI rendering with document highlighting
- Advanced table extraction pipelines
- Archival with full fidelity

---

## HTML

**Method:** `doc.export_to_html()`
**File extension:** `.html`

### Characteristics
- Hierarchical structure with semantic HTML5 elements
- Tables as `<table>` elements
- Headings as `<h1>–<h6>`
- Figures with `<figure>` and `<figcaption>`
- Inline styles for layout preservation

### Use When
- Web rendering of parsed documents
- Feed into HTML parsers for further processing
- Document viewers

---

## DocTags

**Method:** `doc.export_to_doctags()`
**File extension:** `.dt` (custom)

### Characteristics
- Docling's native structured markup format
- Designed for AI model training and fine-tuning
- Encodes document semantics with explicit tags
- Used in GraniteDocling VLM training

### Format Example

```
<document>
<section>
<title>Abstract</title>
<paragraph>This paper presents...</paragraph>
</section>
<table>
<row><cell>Sample</cell><cell>Count</cell></row>
<row><cell>Train</cell><cell>1000</cell></row>
</table>
</document>
```

### Use When
- Training or fine-tuning document understanding models
- When strict semantic structure preservation matters more than readability

---

## Plain Text

**Method:** `doc.export_to_text()`
**File extension:** `.txt`

### Characteristics
- Stripped of all formatting
- Just raw text, newline-separated
- Fastest export
- No table structure preserved

### Use When
- Classical NLP pipelines (tokenization, NER, etc.)
- Quick text extraction without structure
- Keyword search indexing

---

## Table Export (CSV)

Tables can be exported individually as CSV via the `parse_pdf.py` script with `export_tables=True`.

### Method

```python
# Iterate tables from the document
for i, table in enumerate(doc.tables):
    df = table.export_to_dataframe()   # pandas DataFrame
    df.to_csv(f"tables/table_{i}.csv", index=False)
```

### Notes
- Requires pandas
- Merged/spanning cells may be flattened incorrectly (known limitation)
- Column headers detected automatically when `column_header=True` in cell metadata

---

## Choosing the Right Format

| Use Case | Recommended Format |
|---|---|
| LLM / RAG pipeline | Markdown |
| Library deposit | Markdown |
| Vector store with provenance | JSON (lossless) |
| Table extraction | JSON + CSV export |
| Web rendering | HTML |
| Classical NLP | Plain text |
| Model training | DocTags |
| Human review | Markdown or HTML |

---

*Reference compiled by Zevo 🦓 | 2026-05-27*
