---
id: skill_e829fe594ec35d84869daa5331db3de8
name: abstract-lookup
display-name: Abstract Lookup
short-description: Maximizes abstract coverage for research papers using a waterfall
  approach across multiple APIs (PubMed → Europe PMC → Semantic Scholar → Crossref).
  Use when building publication databases, analyzing scholar profiles, or needing
  >95% abstract coverage for biomedical literature. Significantly improves on single-source
  approaches (OpenAlex ~70% → waterfall ~95%+).
category: methodology
keywords: literature, abstracts, PubMed, Europe PMC, Semantic Scholar, Crossref, metadata,
  biomedical literature
version: 1.0.1
last-updated: '2026-08-19'
status: active
license: MIT
author:
  agent: UV-2001
  human: UV-1001
provenance:
  wraps: PubMed E-utilities / Europe PMC / Semantic Scholar / Crossref APIs (waterfall)
  citation: 'Sayers EW et al. Database resources of the National Center for Biotechnology
    Information. Nucleic Acids Research. 2022;50:D20–D26. doi:10.1093/nar/gkab1112
    | Europe PMC Consortium. Europe PMC: a full-text literature database for the life
    sciences and platform for innovation. Nucleic Acids Research. 2015;43:D1042–D1048.
    doi:10.1093/nar/gku1061 | Kinney RM et al. The Semantic Scholar Open Data Platform.
    arXiv:2301.10140'
changelog:
- version: 1.0.1
  date: '2026-07-07'
  event: citation-provenance
  summary: Added provenance block and References section. Vienna public subset.
description: Maximizes abstract coverage for research papers using a waterfall approach
  across multiple APIs (PubMed → Europe PMC → Semantic Scholar → Crossref). Use when
  building publication databases, analyzing scholar profiles, or needing >95% abstract
  coverage for biomedical literature. Significantly improves on single-source approaches
  (OpenAlex ~70% → waterfall ~95%+).
---

# Abstract Lookup — Waterfall Retrieval for Research Papers

This skill provides high-coverage abstract retrieval by trying multiple academic APIs in waterfall order until an abstract is found.

## When to Use

- Building publication databases for lab members or collaborators
- Analyzing scholar profiles with complete abstract coverage
- Literature review requiring full-text availability information
- Improving abstract coverage from incomplete sources (e.g., OpenAlex)
- Any task requiring >95% abstract success rate for biomedical papers

## How It Works

The waterfall approach tries sources in order of speed and biomedical relevance:

1. **PubMed E-utilities** — Free, fast, excellent biomedical coverage (~36M citations)
2. **Europe PMC** — Free, broader coverage including European journals (~40M records)
3. **Semantic Scholar** — Free (rate-limited), AI-enhanced, cross-disciplinary (200M+ papers)
4. **Crossref** — Free, DOI-based, fills remaining gaps (150M+ records)

Each source is tried only if previous sources failed, minimizing API calls.

## Usage

### Basic Usage

The core script is `scripts/fetch_abstracts.py`. It expects a JSON file with paper metadata:

```bash
python3 scripts/fetch_abstracts.py papers.json output_with_abstracts.json
```

**Input format** (minimal):
```json
[
  {
    "title": "Paper Title",
    "doi": "10.1234/journal.2024.12345",
    "pmid": "38123456",
    "year": 2024,
    "authors": "Smith J, Doe A"
  }
]
```

**Output format** (adds `abstract` field):
```json
[
  {
    "title": "Paper Title",
    "doi": "10.1234/journal.2024.12345",
    "pmid": "38123456",
    "year": 2024,
    "authors": "Smith J, Doe A",
    "abstract": "Full abstract text retrieved from best available source..."
  }
]
```

### Field Requirements

**Required:**
- `title` — Paper title (used for Semantic Scholar search if no DOI/PMID)

**Recommended:**
- `doi` — Digital Object Identifier (enables Semantic Scholar, Crossref)
- `pmid` — PubMed ID (enables PubMed, Europe PMC)

**Optional:**
- `year`, `authors`, `journal` — Preserved in output but not used for retrieval

### Example Workflow

```python
# 1. Load papers from existing source (e.g., OpenAlex)
papers = load_papers_from_openalex("Stephen Ekker")

# 2. Export to JSON
with open('ekker_papers.json', 'w') as f:
    json.dump(papers, f)

# 3. Run waterfall abstract fetcher
!python3 scripts/fetch_abstracts.py ekker_papers.json ekker_papers_with_abstracts.json

# 4. Check stats (printed to stderr)
# Success rate: 289/304 (95.1%)
```

## Performance

**Typical coverage improvement:**
- OpenAlex alone: 68-79% (biomedical researchers)
- With waterfall: 95%+ (biomedical literature)

**Rate limits:**
- PubMed: No key required, 3 requests/second recommended
- Europe PMC: No key required, generous limits
- Semantic Scholar: 100 requests/5 minutes (script includes delays)
- Crossref: No key required, polite usage recommended

**Speed:**
- ~0.5-1 second per paper (with rate limiting)
- 300 papers ≈ 5-10 minutes

## Output Stats

The script prints statistics showing which sources provided abstracts:

```
=== Abstract Retrieval Stats ===
Total papers: 304
  PubMed: 187
  Europe PMC: 45
  Semantic Scholar: 52
  Crossref: 15

Success rate: 289/304 (95.1%)
```

## Limitations

- **Not magic** — Some papers genuinely don't have abstracts (editorials, author responses, etc.)
- **Rate limits** — Semantic Scholar is slowest due to rate limiting
- **Network dependent** — Requires internet access
- **No API keys** — Using free tiers; commercial users may need keys for higher limits

## API Documentation

See `references/api_docs.md` for detailed API specifications and authentication options.

## References

- **PubMed / NCBI**: Sayers EW, Bolton EE, Brister JR, et al. "Database resources of the National Center for Biotechnology Information." *Nucleic Acids Research* 2022;50:D20–D26. doi:10.1093/nar/gkab1112 — [link](https://doi.org/10.1093/nar/gkab1112)
- **Europe PMC**: Europe PMC Consortium. "Europe PMC: a full-text literature database for the life sciences and platform for innovation." *Nucleic Acids Research* 2015;43:D1042–D1048. doi:10.1093/nar/gku1061 — [link](https://doi.org/10.1093/nar/gku1061)
- **Semantic Scholar**: Kinney RM, Anastasiades G, Authur R, et al. "The Semantic Scholar Open Data Platform." arXiv:2301.10140 — [link](https://arxiv.org/abs/2301.10140)
