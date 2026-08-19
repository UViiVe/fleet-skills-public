---
id: skill_e2110fdc0665fef64f9d2f7428741d21
name: europepmc-search
display-name: Europe PMC Search
short-description: Search Europe PMC for preprints and published papers across 30+
  sources in one API call. No API key required. Covers bioRxiv, medRxiv, Research
  Square, Preprints.org, PubMed, PMC, and more. Use when this workflow is needed.
category: methodology
keywords: Europe PMC, PubMed, preprints, literature search, biomedical papers, API
version: 1.0.1
last-updated: '2026-08-19'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
provenance:
  wraps: Europe PMC REST API (European Bioinformatics Institute)
  api_docs: https://europepmc.org/RestfulWebService
  citation: 'Europe PMC Consortium. "Europe PMC: a full-text literature database for
    the life sciences and platform for innovation." Nucleic Acids Research. 2015;43:D1042–D1048.
    doi:10.1093/nar/gku1061'
changelog:
- version: 1.0.1
  date: '2026-07-07'
  event: citation-provenance
  summary: Backfilled public schema, added provenance block and References section.
    Vienna public subset.
- version: 1.0.0
  date: '2026-06-01'
  event: authored
  summary: Initial skill.
description: Search Europe PMC for preprints and published papers across 30+ sources
  in one API call. No API key required. Covers bioRxiv, medRxiv, Research Square,
  Preprints.org, PubMed, PMC, and more.
---

# europepmc-search

Search Europe PMC for preprints and published papers across 30+ sources in one API call.
No API key required. Covers bioRxiv, medRxiv, Research Square, Preprints.org, PubMed, PMC, and more.

## When to use

- Finding preprints across ALL servers (not just bioRxiv/medRxiv) in one query
- Cross-server preprint monitoring for a topic (CF, rare disease, CRISPR, etc.)
- Full-text search including abstracts (stronger than bioRxiv API)
- Finding open-access papers on a topic
- Searching PubMed alongside preprints in one pass

## Script

`scripts/europepmc_search.py`

## Usage

```bash
# Preprints only, last 30 days
python3 ~/.openclaw/workspace/skills/europepmc-search/scripts/europepmc_search.py "cystic fibrosis CFTR" --preprints-only --days 30

# Summary stats only
python3 ~/.openclaw/workspace/skills/europepmc-search/scripts/europepmc_search.py "cystic fibrosis" --preprints-only --days 60 --summary

# All sources (published + preprints), last 14 days
python3 ~/.openclaw/workspace/skills/europepmc-search/scripts/europepmc_search.py "VEGF-A peptide" --days 14 --max 20

# Open access only
python3 ~/.openclaw/workspace/skills/europepmc-search/scripts/europepmc_search.py "zebrafish CRISPR" --open-access --max 10

# Specific source (PPR = preprints, MED = PubMed, PMC = full text)
python3 ~/.openclaw/workspace/skills/europepmc-search/scripts/europepmc_search.py "CHCHD2 mitochondria" --source MED --max 10

# List available source codes
python3 ~/.openclaw/workspace/skills/europepmc-search/scripts/europepmc_search.py --list-sources
```

## Arguments

| Argument | Description | Default |
|---|---|---|
| `keywords` | Search terms (Europe PMC query language) | — |
| `--days` | Search last N days | none (all time) |
| `--from` | Start date `YYYY-MM-DD` | — |
| `--to` | End date `YYYY-MM-DD` | — |
| `--max` | Max results (default: 25) | 25 |
| `--source` | Source code: `MED`, `PMC`, `PPR` | all |
| `--preprints-only` | Filter to preprints only (SRC:PPR) | off |
| `--open-access` | Open access papers only | off |
| `--pmid` | Fetch specific paper by PubMed ID | — |
| `--summary` | Print counts only, no full records | off |
| `--list-sources` | List source codes | — |

## Output

JSON with `total_hits_in_epmc` (full database count), `returned` (actual results), `by_source` (breakdown), and `results` array with: `title`, `authors`, `date`, `source`, `journal`, `doi`, `pmid`, `url`, `abstract` (500 char), `open_access`, `cited_by`.

## Source Codes

- `PPR` — Preprints (bioRxiv, medRxiv, Research Square, Preprints.org, etc.)
- `MED` — PubMed/MEDLINE (published, indexed)
- `PMC` — PubMed Central (full-text open access)
- `PAT` — Patents
- `CTX` — ClinicalTrials.gov

## Notes

- API: `https://www.ebi.ac.uk/europepmc/webservices/rest/search`
- Free, no key required; uses cursor-based pagination
- `sort` parameter not compatible with cursor pagination — results ordered by relevance
- No pip dependencies — uses `urllib` only
- Europe PMC query language supports field tags: `AUTH:`, `TITLE:`, `ABSTRACT:`, `JOURNAL:`, `OPEN_ACCESS:Y`

## Installed on

- Atlas ✅ (no pip deps needed)
- Copy the helper script from this skill directory when installing in another agent environment.

---

## References

- **Europe PMC**: Europe PMC Consortium. "Europe PMC: a full-text literature database for the life sciences and platform for innovation." *Nucleic Acids Research*. 2015;43:D1042–D1048. doi:10.1093/nar/gku1061
- **API documentation**: https://europepmc.org/RestfulWebService
- **API base URL**: https://www.ebi.ac.uk/europepmc/webservices/rest/search
