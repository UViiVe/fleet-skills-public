---
name: ncbi-api
id: "skill_7c8e238ec0144b008ad4123ab9497a8c"
display-name: "NCBI API"
short-description: "Search PubMed, fetch abstracts, retrieve full-text links, and download open-access papers via NCBI E-utilities and PMC APIs. Use when searching biomedical literature by keyword/MeSH/author/PMID, fe..."
category: methodology
keywords: [ncbi, pubmed, api, abstract, literature]
last-updated: "2026-07-07"
version: "1.0.1"
status: active
license: MIT
description: "Search PubMed, fetch abstracts, retrieve full-text links, and download open-access papers via NCBI E-utilities and PMC APIs. Use when searching biomedical literature by keyword/MeSH/author/PMID, fetching abstracts in bulk, checking open-access availability, or downloading PMC full-text PDFs. Covers PubMed, PubMed Central (PMC), Gene, Nucleotide, Protein, and SRA databases. Handles rate limits, API key config, and batch fetching. Triggers on: search PubMed, fetch abstract, PMID lookup, PMC download, NCBI search, find papers on, get full text."
author: "Fleet Skills contributors"
provenance:
  wraps: "NCBI E-utilities (PubMed, PubMed Central, Entrez)"
  api_docs: "https://www.ncbi.nlm.nih.gov/books/NBK25497/"
  citation: "Sayers EW, et al. \"Database resources of the National Center for Biotechnology Information.\" Nucleic Acids Research. 2022;50:D20\u2013D26. doi:10.1093/nar/gkab1112"
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

# NCBI API

NCBI E-utilities base URL: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

## API Key Setup

Unauthenticated: 3 requests/second. With API key: 10 requests/second.

Check if key is configured:
```bash
test -n "$NCBI_API_KEY" && echo "NCBI_API_KEY is set"
```

Set API key (get one free at https://www.ncbi.nlm.nih.gov/account/):
```bash
export NCBI_API_KEY="YOUR_KEY_HERE"
```

In API calls, append `&api_key=$NCBI_API_KEY` when set.

## Core Workflow

### 1. Search PubMed (esearch)
```
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi
  ?db=pubmed
  &term=QUERY
  &retmax=20
  &retmode=json
  &usehistory=y
```
Returns: `esearchresult.idlist` (array of PMIDs), `webenv`, `query_key` for large batch fetching.

**Query syntax:**
- Keywords: `CRISPR cancer`
- MeSH: `neoplasms[MeSH Terms]`
- Author: `Jane Smith[Author]`
- Date range: `2020:2026[pdat]`
- Combine: `TP53[Title] AND neoplasms[MeSH Terms]`

### 2. Fetch Abstracts (efetch)
```
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
  ?db=pubmed
  &id=PMID1,PMID2,PMID3
  &rettype=abstract
  &retmode=text
```
For JSON/XML structured data use `&retmode=xml` with `&rettype=medline`.

Batch limit: 200 PMIDs per request (use `&query_key` + `&WebEnv` for larger sets).

### 3. Fetch Summaries (esummary)
```
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi
  ?db=pubmed
  &id=PMID1,PMID2
  &retmode=json
```
Returns title, authors, journal, pubdate, doi, pmcid per record.

### 4. Check PMC Full-Text Availability
```
GET https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi
  ?dbfrom=pubmed
  &db=pmc
  &id=PMID
  &retmode=json
```
If `linksets[0].linksetdbs` contains `pubmed_pmc`, the paper is in PMC. Extract the PMCID.

### 5. Download PMC Full-Text PDF
** Never use PMC web URL for PDFs** - returns JS redirect stub.

Use PMC FTP:
```python
# Step 1: Get FTP path from OA file list
import requests
r = requests.get(f"https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC{pmcid}&format=pdf")
# Parse XML: <link format="pdf" href="ftp://..."/>

# Step 2: Convert ftp:// to https://ftp.ncbi.nlm.nih.gov/...
ftp_path = link_href.replace("ftp://ftp.ncbi.nlm.nih.gov", "https://ftp.ncbi.nlm.nih.gov")
pdf = requests.get(ftp_path)
# Validate: len(pdf.content) > 10000 and pdf.content[:4] == b'%PDF'
```

## Rate Limiting

Without API key: add `time.sleep(0.34)` between requests (3/s limit).
With API key: add `time.sleep(0.1)` between requests (10/s limit).
On 429: back off 5 seconds and retry once.

## Other NCBI Databases

Change `&db=` parameter:
- `gene` - gene records (use `&rettype=gene_table`)
- `nucleotide` - DNA/RNA sequences
- `protein` - protein sequences
- `sra` - sequencing datasets
- `pmc` - search full-text PMC articles directly

## Reference

See `references/endpoints.md` for full endpoint reference, XML field mappings, and bulk fetch patterns.

---

## References

- **NCBI E-utilities**: Sayers EW, et al. "Database resources of the National Center for Biotechnology Information." *Nucleic Acids Research*. 2022;50:D20–D26. doi:10.1093/nar/gkab1112
- **API documentation**: https://www.ncbi.nlm.nih.gov/books/NBK25497/
- **API key registration**: https://www.ncbi.nlm.nih.gov/account/
