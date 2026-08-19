---
id: skill_1cda0f22c69211ebdbfbadcd38531a10
name: arxiv-search
display-name: arXiv Search
short-description: 'Search arXiv for academic papers using the official arXiv API.
  Use when: (1) looking up recent preprints or papers by topic, (2) fetching paper
  metadata (title, authors, abstract, PDF URL, DOI) for any domain, (3) scoping searches
  to specific arXiv categories (e.g. cs.AI, q-bio.GN, q-bio.BM, physics.bio-ph), (4)
  seeding literature reviews or research scouting runs. No API key required. Works
  on any machine with Python 3 and the arxiv package installed.'
category: methodology
keywords: arxiv, preprint, literature search, metadata, academic papers, API
version: 1.0.1
last-updated: '2026-08-19'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
provenance:
  wraps: arXiv API (arxiv Python package, lukasschwab/arxiv.py)
  citation: 'arXiv is maintained by Cornell University. No primary paper. Python client:
    https://github.com/lukasschwab/arxiv.py | API: https://arxiv.org/help/api'
changelog:
- version: 1.0.1
  date: '2026-07-07'
  event: citation-provenance
  summary: Added provenance block and References section. Vienna public subset.
- version: 1.0.0
  date: '2026-07-07'
  event: schema-compliance
  summary: 'Backfilled public schema: version, license, status, author. De-agentified
    for public release. Vienna public subset.'
description: 'Search arXiv for academic papers using the official arXiv API. Use when:
  (1) looking up recent preprints or papers by topic, (2) fetching paper metadata
  (title, authors, abstract, PDF URL, DOI) for any domain, (3) scoping searches to
  specific arXiv categories (e.g. cs.AI, q-bio.GN, q-bio.BM, physics.bio-ph), (4)
  seeding literature reviews or research scouting runs. No API key required. Works
  on any machine with Python 3 and the arxiv package installed.

  '
---

# arXiv Search

Query the arXiv API for academic papers. Free, no key required.

## Quick Start

```bash
python3 skills/arxiv-search/scripts/arxiv_search.py "VEGF-A peptide inhibitor" --max 10
```

Returns JSON array with: `arxiv_id`, `title`, `authors` (up to 6), `abstract` (500 chars), `published`, `updated`, `pdf_url`, `categories`, `journal_ref`, `doi`.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--max N` | 10 | Results (max 50) |
| `--sort` | `relevance` | `relevance` / `date` / `updated` |
| `--category` | none | arXiv category filter (e.g. `cs.AI`, `q-bio.GN`) |

## Category Reference

See `references/arxiv-categories.md` for full list. Key categories for research use:

- `cs.AI` — AI/ML general
- `cs.LG` — Machine learning
- `cs.CL` — NLP/language models
- `q-bio.GN` — Genomics
- `q-bio.BM` — Biomolecules
- `q-bio.QM` — Quantitative methods in biology
- `physics.bio-ph` — Biological physics

## Installation

```bash
pip3 install arxiv
```

Verify: `python3 -c "import arxiv; print('OK')"`

## Usage

Call the script and parse the JSON output:

```bash
python3 skills/arxiv-search/scripts/arxiv_search.py \
  "zebrafish CHCHD2 mitochondria" --max 5 --sort date
```

For programmatic use within Python, import directly:

```python
import sys
sys.path.insert(0, 'skills/arxiv-search/scripts')
from arxiv_search import search
results = search("sparse attention transformer", max_results=10, category="cs.LG")
```

## References

- **arXiv**: Maintained by Cornell University. No primary paper. API documentation: https://arxiv.org/help/api
- **arxiv Python package**: lukasschwab/arxiv.py — [https://github.com/lukasschwab/arxiv.py](https://github.com/lukasschwab/arxiv.py)
