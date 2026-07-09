# Semantic Scholar API — Fields & Pagination Reference

## Paper Fields (complete list)

| Field | Type | Description |
|---|---|---|
| `paperId` | string | S2 unique ID |
| `externalIds` | object | DOI, PubMed, ArXiv, PubMedCentral, CorpusId, ACL, DBLP, MAG |
| `url` | string | S2 paper page URL |
| `title` | string | Paper title |
| `abstract` | string | Full abstract |
| `venue` | string | Journal/conference name |
| `publicationVenue` | object | `{id, name, type, issn, url}` |
| `year` | int | Publication year |
| `referenceCount` | int | Number of references |
| `citationCount` | int | Total citations |
| `influentialCitationCount` | int | Highly-influential citations (S2 metric) |
| `isOpenAccess` | bool | Whether OA version exists |
| `openAccessPdf` | object | `{url, status}` — direct PDF URL if available |
| `fieldsOfStudy` | array | e.g. `["Biology", "Medicine"]` |
| `s2FieldsOfStudy` | array | More granular S2 categorization |
| `publicationTypes` | array | e.g. `["JournalArticle", "Review"]` |
| `publicationDate` | string | ISO date |
| `journal` | object | `{name, volume, pages}` |
| `authors` | array | `[{authorId, name}]` — add nested fields for more |
| `citations` | array | Papers citing this (nested Paper objects) |
| `references` | array | Papers cited by this (nested Paper objects) |
| `embedding` | object | S2 semantic embedding vector |
| `tldr` | object | `{model, text}` — AI one-sentence summary |

## Author Fields

| Field | Type | Description |
|---|---|---|
| `authorId` | string | S2 unique author ID |
| `name` | string | Author name |
| `affiliations` | array | Institution strings |
| `homepage` | string | Author webpage |
| `paperCount` | int | Total papers |
| `citationCount` | int | Total citations |
| `hIndex` | int | h-index |
| `papers` | array | Nested Paper objects |

## Nested Field Syntax

Request nested fields with dot notation:
```
?fields=authors.name,authors.affiliations,citations.title,citations.year
```

## Pagination

All list endpoints support `limit` and `offset`:
```
GET /paper/{id}/citations?limit=100&offset=0
GET /paper/{id}/citations?limit=100&offset=100   # page 2
```

Response includes:
```json
{
  "offset": 0,
  "next": 100,       // next offset (absent when no more results)
  "data": [...]
}
```

Loop until `next` is absent or `len(data) < limit`.

## Bulk Fetch Pattern (Python)
```python
import requests, time

API_KEY = "your_key"  # or None
HEADERS = {"x-api-key": API_KEY} if API_KEY else {}
BASE = "https://api.semanticscholar.org/graph/v1"

def batch_fetch(ids, fields="paperId,title,abstract,year,externalIds"):
    """Fetch up to 500 papers by mixed IDs."""
    r = requests.post(
        f"{BASE}/paper/batch",
        params={"fields": fields},
        json={"ids": ids},
        headers=HEADERS
    )
    r.raise_for_status()
    time.sleep(1 if API_KEY else 3)
    return r.json()

def get_all_citations(paper_id, fields="paperId,title,year,citationCount"):
    """Fetch all citations with pagination."""
    results, offset = [], 0
    while True:
        r = requests.get(
            f"{BASE}/paper/{paper_id}/citations",
            params={"fields": fields, "limit": 100, "offset": offset},
            headers=HEADERS
        )
        r.raise_for_status()
        data = r.json()
        results.extend([item["citingPaper"] for item in data["data"]])
        time.sleep(1 if API_KEY else 3)
        if "next" not in data:
            break
        offset = data["next"]
    return results
```

## ID Type Prefixes for Cross-DB Lookup
```
PMID:12345678          → PubMed ID
DOI:10.1038/nature12160 → DOI
ARXIV:2106.01345       → ArXiv
PMC:3338380            → PubMed Central
ACL:P19-1042           → ACL Anthology
MAG:2963403868         → Microsoft Academic Graph
```

## Common Search Query Examples
```
# Natural language
CRISPR cancer models

# Specific genes
TP53 tumor suppressor

# Find reviews
mitochondrial dynamics review

# Limit by year (use field filter in results, not query syntax)
# → filter results where year >= 2020 after fetching
```
