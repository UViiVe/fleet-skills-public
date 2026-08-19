# API Documentation Reference

## PubMed E-utilities

**Base URL:** `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

**Key endpoints:**
- `esearch.fcgi` — Search for PMIDs
- `efetch.fcgi` — Fetch full records (including abstracts)

**Usage (no API key required for low volume):**
```bash
# Fetch abstract by PMID
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=38123456&retmode=xml&rettype=abstract"
```

**Rate limits:** 3 requests/second without key, 10/second with key

**Documentation:** https://www.ncbi.nlm.nih.gov/books/NBK25501/

---

## Europe PMC REST API

**Base URL:** `https://www.ebi.ac.uk/europepmc/webservices/rest/`

**Search endpoint:**
```bash
curl "https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=ext_id:38123456&format=json&resulttype=core"
```

**Query formats:**
- By PMID: `ext_id:38123456 src:med`
- By DOI: `doi:10.1234/journal.2024.123`
- By title: `title:"Exact Paper Title"`

**Rate limits:** Generous, no key required

**Documentation:** https://europepmc.org/RestfulWebService

---

## Semantic Scholar API

**Base URL:** `https://api.semanticscholar.org/graph/v1/`

**Endpoints:**
- Paper by DOI: `/paper/DOI:{doi}?fields=abstract`
- Paper search: `/paper/search?query={title}&fields=abstract&limit=1`

**Example:**
```bash
curl "https://api.semanticscholar.org/graph/v1/paper/DOI:10.1038/nature12345?fields=abstract"
```

**Rate limits:** 
- Free tier: 100 requests per 5 minutes
- Recommended: Add 200ms delay between requests

**Documentation:** https://api.semanticscholar.org/api-docs/graph

---

## Crossref API

**Base URL:** `https://api.crossref.org/works/`

**Usage:**
```bash
curl "https://api.crossref.org/works/10.1038/nature12345"
```

**Returns:** Full metadata including abstract when available (hit-or-miss)

**Rate limits:** 
- Polite pool: 50 requests/second (with User-Agent header)
- No API key required

**Documentation:** https://www.crossref.org/documentation/retrieve-metadata/rest-api/

---

## Authentication (Optional)

### PubMed/NCBI
Register for API key at: https://www.ncbi.nlm.nih.gov/account/settings/

Add to requests: `&api_key=YOUR_KEY`

### Semantic Scholar
Contact for API key at: https://www.semanticscholar.org/product/api

Higher rate limits available for academic/research use.

### Best Practice
Always include:
- Descriptive User-Agent header
- Email in requests for API support contact
- Rate limiting/delays in your code
