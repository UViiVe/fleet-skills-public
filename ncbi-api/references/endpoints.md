# NCBI E-utilities - Endpoint Reference

## Base URL
`https://eutils.ncbi.nlm.nih.gov/entrez/eutils/`

## Endpoints

| Endpoint | Purpose | Key params |
|---|---|---|
| `esearch.fcgi` | Search a database, get IDs | `db`, `term`, `retmax`, `usehistory` |
| `efetch.fcgi` | Fetch full records by ID | `db`, `id`, `rettype`, `retmode` |
| `esummary.fcgi` | Fetch doc summaries by ID | `db`, `id`, `retmode=json` |
| `elink.fcgi` | Find related records across DBs | `dbfrom`, `db`, `id` |
| `einfo.fcgi` | Database field/index info | `db` |
| `egquery.fcgi` | Search all NCBI DBs at once | `term` |

## PubMed `rettype` options (efetch)
- `abstract` - plain-text abstract + citation
- `medline` - MEDLINE format (pipe-delimited)
- `xml` - full PubMed XML (most complete)
- `uilist` - just a list of PMIDs

## PubMed XML Field Mappings
```
MedlineCitation/PMID          → PMID
MedlineCitation/Article/ArticleTitle → Title
MedlineCitation/Article/Abstract/AbstractText → Abstract
MedlineCitation/Article/AuthorList/Author → Authors
MedlineCitation/Article/Journal/Title → Journal
MedlineCitation/Article/Journal/JournalIssue/PubDate → PubDate
PubmedData/ArticleIdList/ArticleId[@IdType="doi"] → DOI
PubmedData/ArticleIdList/ArticleId[@IdType="pmc"] → PMCID
```

## Large Batch Fetch (>200 records)
```python
import requests, time

# Step 1: esearch with usehistory=y
r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={
    "db": "pubmed", "term": "your query", "retmax": 0,
    "usehistory": "y", "retmode": "json"
})
data = r.json()
webenv = data["esearchresult"]["webenv"]
query_key = data["esearchresult"]["querykey"]
total = int(data["esearchresult"]["count"])

# Step 2: efetch in batches of 200
for start in range(0, total, 200):
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi", params={
        "db": "pubmed", "query_key": query_key, "WebEnv": webenv,
        "retstart": start, "retmax": 200, "rettype": "xml", "retmode": "xml"
    })
    # process r.content
    time.sleep(0.1)  # with API key; use 0.34 without
```

## PMC OA API (PDF/full-text download)
```
GET https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id=PMC{pmcid}&format=pdf
```
Response XML:
```xml
<OA>
  <records>
    <record id="PMC1234567" citation="...">
      <link format="pdf" href="ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa/pdf/...pdf"/>
    </record>
  </records>
</OA>
```
Convert `ftp://ftp.ncbi.nlm.nih.gov` → `https://ftp.ncbi.nlm.nih.gov` for HTTP download.

## Common Query Examples
```
# All papers by Jane Smith since 2020
Jane Smith[Author] AND 2020:2026[pdat]

# CRISPR cancer models
CRISPR[Title/Abstract] AND neoplasms[MeSH Terms]

# TP53 tumor suppressor
TP53[Gene/Protein Name] AND neoplasms[MeSH Terms]

# Papers with free full text
your query AND free full text[filter]

# Clinical trials only
your query AND clinical trial[pt]
```

## Error Codes
- `429` - rate limit hit, back off 5s
- `414` - query too long, split into smaller batches
- Empty `idlist` - no results, check query syntax
