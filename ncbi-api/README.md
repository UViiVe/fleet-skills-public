# ncbi-api Skill

Search PubMed, fetch abstracts, retrieve full-text links, and download open-access papers via NCBI E-utilities and PMC APIs.

## What It Does

- Search PubMed by keyword, MeSH term, author, or date range
- Fetch abstracts and metadata for lists of PMIDs
- Check PMC full-text availability and download open-access PDFs
- Access all NCBI databases: Gene, Nucleotide, Protein, SRA, PMC
- Handle batch fetching (up to 200 PMIDs per request, unlimited via WebEnv)

## API Key Setup

**Without a key:** 3 requests/second (public rate limit)  
**With a key:** 10 requests/second (free, takes ~2 minutes to set up)

Get a free NCBI API key at: https://www.ncbi.nlm.nih.gov/account/

Set it in your shell environment:
```bash
export NCBI_API_KEY="YOUR_KEY_HERE"
```

## Files

- `SKILL.md` - Full workflow instructions loaded when skill triggers
- `references/endpoints.md` - Complete endpoint reference, XML field mappings, bulk fetch patterns, query syntax examples

## Quick Reference

| Task | Endpoint |
|---|---|
| Search PubMed | `esearch.fcgi?db=pubmed&term=QUERY` |
| Fetch abstracts | `efetch.fcgi?db=pubmed&id=PMID1,PMID2&rettype=abstract` |
| Fetch summaries (JSON) | `esummary.fcgi?db=pubmed&id=PMID1,PMID2&retmode=json` |
| Check PMC availability | `elink.fcgi?dbfrom=pubmed&db=pmc&id=PMID` |
| Download PMC PDF | PMC OA API → FTP path (see SKILL.md) |

 **PDF rule:** Always use the PMC FTP path (`https://ftp.ncbi.nlm.nih.gov/...`), never the PMC web URL - the web URL returns a JS redirect stub.

## Rate Limits

- Without key: `time.sleep(0.34)` between requests
- With key: `time.sleep(0.1)` between requests
- On 429: back off 5 seconds, retry once

---

*Added to the skill library: March 23, 2026*
