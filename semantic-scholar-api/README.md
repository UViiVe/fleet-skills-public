# semantic-scholar-api Skill

Search academic papers, build citation networks, get author profiles, and retrieve paper metadata via the Semantic Scholar Academic Graph API. Works across all fields - biology, medicine, CS, physics, and more.

## What It Does

- Search papers by natural language query across all academic disciplines
- Get full citation and reference lists for any paper
- Look up author profiles: h-index, paper count, citation count, affiliations
- Fetch paper details by DOI, PubMed ID, ArXiv ID, or PMC ID
- Get paper recommendations based on one or more seed papers
- Bulk fetch up to 500 papers in a single request with mixed ID types

## API Key Setup

**Without a key:** 100 requests/5 min (shared public pool)  
**With a key:** 1 request/second dedicated (free, ~1 day approval time)

Request a free Semantic Scholar API key at: https://www.semanticscholar.org/product/api#api-key-form

Set it in your shell environment:
```bash
export SEMANTIC_SCHOLAR_API_KEY="YOUR_KEY_HERE"
```

Pass in API calls as header: `x-api-key: YOUR_KEY`

## Files

- `SKILL.md` - Full workflow instructions loaded when skill triggers
- `references/fields.md` - Complete paper/author field listings, pagination patterns, bulk fetch Python code

## Quick Reference

| Task | Endpoint |
|---|---|
| Search papers | `GET /paper/search?query=QUERY&fields=...` |
| Get paper by ID | `GET /paper/{paperId\|DOI:...\|PMID:...}` |
| Get citations | `GET /paper/{id}/citations` |
| Get references | `GET /paper/{id}/references` |
| Search authors | `GET /author/search?query=NAME` |
| Author papers | `GET /author/{authorId}/papers` |
| Bulk fetch | `POST /paper/batch` with `{"ids": [...]}` |

**Supported ID prefixes:** `PMID:`, `DOI:`, `ARXIV:`, `PMC:`, `ACL:`, `MAG:`

## Key Fields

- `externalIds` - get DOI, PubMed, ArXiv, PMC IDs for cross-linking
- `openAccessPdf.url` - direct PDF link when open access
- `tldr.text` - AI one-sentence summary for rapid screening
- `influentialCitationCount` - better impact metric than raw citation count

## Rate Limits

- Without key: `time.sleep(3)` between requests
- With key: `time.sleep(1)` between requests
- On 429: check `Retry-After` header

---

*Added to the skill library: March 23, 2026*
