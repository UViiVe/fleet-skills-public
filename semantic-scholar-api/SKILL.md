---
name: semantic-scholar-api
id: "skill_44a0758dcf404f0484a72d953e656f2a"
display-name: "Semantic Scholar API"
short-description: "Search academic papers, fetch citations/references, get author profiles, and retrieve paper details via the Semantic Scholar Academic Graph API. Use when searching cross-disciplinary literature, fi..."
category: methodology
keywords: [semantic scholar, citations, api, literature, papers]
last-updated: "2026-07-07"
version: "1.0.1"
status: active
license: MIT
description: "Search academic papers, fetch citations/references, get author profiles, and retrieve paper details via the Semantic Scholar Academic Graph API. Use when searching cross-disciplinary literature, finding citing/cited papers, building citation networks, getting author h-index and publication lists, or looking up paper metadata by DOI/ArXiv/PMID. Covers all fields (CS, biology, medicine, physics, etc.). Handles API key config, pagination, and bulk operations. Triggers on: Semantic Scholar, find papers citing, citation network, author profile, h-index, papers by author, related papers, recommendations, search all fields."
author: "Fleet Skills contributors"
provenance:
  wraps: "Semantic Scholar Academic Graph API (Allen Institute for AI)"
  api_docs: "https://api.semanticscholar.org/api-docs/"
  citation: "Kinney RM, et al. \"The Semantic Scholar Open Data Platform.\" arXiv. 2023. arXiv:2301.10140"
changelog:
  - version: "1.0.1"
    date: "2026-07-07"
    event: citation-provenance
    summary: "Backfilled public schema, added provenance block and References section."
  - version: "1.0.0"
    date: "2026-06-01"
    event: authored
    summary: "Initial skill."
promoted: true
---

# Semantic Scholar API

Base URL: `https://api.semanticscholar.org/graph/v1`

## API Key Setup

Unauthenticated: 100 requests/5 min (shared pool). With API key: 1 request/second dedicated.

Check if key configured:
```bash
test -n "$SEMANTIC_SCHOLAR_API_KEY" && echo "SEMANTIC_SCHOLAR_API_KEY is set"
```

Request a free key at: https://www.semanticscholar.org/product/api#api-key-form

Set key:
```bash
export SEMANTIC_SCHOLAR_API_KEY="YOUR_KEY_HERE"
```

Pass in header: `x-api-key: $SEMANTIC_SCHOLAR_API_KEY`

## Core Workflows

### 1. Search Papers
```
GET /paper/search?query=QUERY&fields=paperId,title,abstract,year,authors,citationCount,externalIds,openAccessPdf&limit=10
```
`query` supports natural language or field-specific terms.
`limit` max: 100. Use `offset` for pagination.

**Key fields to request:**
- `paperId` - S2 internal ID (use for follow-up calls)
- `externalIds` - contains `DOI`, `PubMed`, `ArXiv`, `PubMedCentral`
- `openAccessPdf` - `{url, status}` if available
- `citationCount`, `influentialCitationCount`
- `tldr` - AI-generated one-sentence summary

### 2. Get Paper Details by ID
```
GET /paper/{paper_id}?fields=paperId,title,abstract,year,authors,citations,references,openAccessPdf,externalIds,tldr
```
`paper_id` can be:
- S2 ID: `649def34f8be52c8b66281af98ae884c09aef38b`
- DOI: `DOI:10.1038/nature12160`
- PubMed ID: `PMID:23851394`
- ArXiv ID: `ARXIV:2106.01345`
- PMC ID: `PMC:3338380`

### 3. Get Citations (papers that cite this paper)
```
GET /paper/{paper_id}/citations?fields=paperId,title,year,authors,citationCount&limit=100
```
Returns `data[].citingPaper`. Paginate with `offset`.

### 4. Get References (papers cited by this paper)
```
GET /paper/{paper_id}/references?fields=paperId,title,year,authors&limit=100
```
Returns `data[].citedPaper`.

### 5. Author Search & Profile
```
# Search by name
GET /author/search?query=Jane Smith&fields=authorId,name,hIndex,citationCount,paperCount,affiliations

# Get author details
GET /author/{author_id}?fields=authorId,name,hIndex,citationCount,paperCount,affiliations,papers

# Get all papers by author
GET /author/{author_id}/papers?fields=paperId,title,year,citationCount,externalIds&limit=100
```

### 6. Paper Recommendations
```
# Papers similar to one paper
GET /recommendations/v1/papers/forpaper/{paper_id}?fields=paperId,title,year,citationCount&limit=10

# Papers similar to a set of papers (POST)
POST /recommendations/v1/papers/
Body: {"positivePaperIds": ["id1","id2"], "negativePaperIds": []}
```

### 7. Bulk Paper Fetch (POST)
```
POST /paper/batch
Body: {"ids": ["PMID:12345", "DOI:10.1038/...", "ARXIV:2106.01345"]}
Query: ?fields=paperId,title,abstract,externalIds,year
```
Batch limit: 500 IDs per request. Mix ID types freely.

## Rate Limiting

Without key: `time.sleep(3)` between requests.
With key: `time.sleep(1)` between requests.
On 429: check `Retry-After` header, back off accordingly.

## Practical Tips

- **Cross-link to PubMed**: Use `externalIds.PubMed` to get PMID for NCBI follow-up
- **Free PDF check**: `openAccessPdf.url` gives direct PDF link when available
- **Influence score**: `influentialCitationCount` is better than raw `citationCount` for impact
- **TLDR**: `tldr.text` gives a one-sentence AI summary - great for rapid screening

## Reference

See `references/fields.md` for complete field listings and pagination patterns.

---

## References

- **Semantic Scholar**: Kinney RM, et al. "The Semantic Scholar Open Data Platform." *arXiv*. 2023. arXiv:2301.10140
- **API documentation**: https://api.semanticscholar.org
- **API key request**: https://www.semanticscholar.org/product/api#api-key-form