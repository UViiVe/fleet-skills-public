---
id: skill_d49e991b8083918b010c4d9bb68f113f
name: biorxiv-search
display-name: bioRxiv Search
short-description: Search bioRxiv and medRxiv preprints via the official biorxiv.org
  API. No API key required. Use when looking for recent life sciences preprints on
  a topic — genomics, structural biology, biotech, clinical.
category: methodology
keywords: biorxiv, medrxiv, preprint, literature search, life sciences, API
version: 1.0.1
last-updated: '2026-08-19'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
provenance:
  wraps: bioRxiv / medRxiv API (api.biorxiv.org)
  citation: Abdill RJ, Blekhman R. Tracking the popularity and outcomes of all bioRxiv
    preprints. eLife. 2019;8:e45133. doi:10.7554/eLife.45133
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
description: Search bioRxiv and medRxiv preprints via the official biorxiv.org API.
  No API key required. Use when looking for recent life sciences preprints on a topic
  — genomics, structural biology, biotech, clinical.
---

# biorxiv-search Skill

Search bioRxiv and medRxiv preprints via the official biorxiv.org API.
No API key required. No Python package required beyond `requests` (stdlib-adjacent).

## When to use
- User asks about recent bioRxiv/medRxiv preprints on a topic
- Periodic literature scouting on life sciences topics
- Checking for new structural biology, genomics, or biotech preprints

## Script
`skills/biorxiv-search/scripts/biorxiv_search.py`

## Usage

```bash
# Search bioRxiv for VEGF-A peptide papers in last 30 days
python3 skills/biorxiv-search/scripts/biorxiv_search.py "VEGF-A peptide" --server biorxiv --days 30 --max 10

# Search medRxiv for clinical trial results
python3 skills/biorxiv-search/scripts/biorxiv_search.py "rosacea treatment" --server medrxiv --days 60 --max 5

# JSON output for programmatic use
python3 skills/biorxiv-search/scripts/biorxiv_search.py "zebrafish CRISPR" --json
```

## Parameters
- `query` — space-separated terms; ALL must appear in title or abstract
- `--server` — `biorxiv` (default) or `medrxiv`
- `--days` — how many days back to search (default: 30)
- `--max` — max results to return (default: 10)
- `--json` — output raw JSON

## Output fields
- `doi`, `title`, `authors`, `date`, `category`, `abstract` (truncated to 400 chars), `url`, `server`

## Notes
- Search is client-side filtered: fetches all papers in date range, then filters by query terms
- For very broad date ranges or very common terms, use `--days` to limit fetch size
- Rate limit: ~3 req/sec (built-in 0.4s delay between pages)
- API docs: https://api.biorxiv.org

## References

- **bioRxiv**: Abdill RJ, Blekhman R. "Tracking the popularity and outcomes of all bioRxiv preprints." *eLife* 2019;8:e45133. doi:10.7554/eLife.45133 — [link](https://doi.org/10.7554/eLife.45133)
- API documentation: https://api.biorxiv.org
