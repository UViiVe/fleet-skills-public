---
name: horizon
id: "skill_662c71fd86bd43d88c7dba2378937c47"
short-description: "Set up a configurable preprint radar for bioRxiv, medRxiv, arXiv, and ChemRxiv. Use when monitoring new papers against a public-safe relevance profile."
category: methodology
keywords: [horizon, literature, benchmark, paper, retrieval]
display-name: "Horizon Preprint Radar"
description: "Set up a configurable preprint radar for bioRxiv, medRxiv, arXiv, and ChemRxiv. Use when a project needs recurring literature monitoring, relevance-filtered preprint search, or a daily literature briefing."
version: "1.0.1"
author: "Fleet Skills contributors"
created: "2026-06-06"
last-updated: "2026-07-09"
tier: 1
changelog:
  - version: "1.0.1"
    date: "2026-07-09"
    event: public-cleanup
    summary: "Replaced internal deployment details and live profiles with generic public setup guidance."
  - version: "1.0.0"
    date: "2026-06-17"
    event: authored
    summary: "Initial preprint radar skill."
status: active
promoted: true
license: MIT
---

# Horizon Preprint Radar

Horizon is a pattern for recurring preprint monitoring. It fetches new records
from public preprint sources, scores them against a relevance profile, stores
kept records in a searchable index, and keeps a lightweight archive of seen
DOIs so the process can resume cleanly.

This public skill intentionally uses placeholders for servers, databases, and
profiles. Keep real team profiles and infrastructure runbooks in a private
repository.

---

## When to Use This Skill

Use when:

- A project needs daily or weekly awareness of new papers.
- A team wants a repeatable relevance profile instead of ad-hoc keyword alerts.
- You need a searchable local index of papers that passed a relevance threshold.
- You want to track what was discarded so future profile changes can backfill.

Do not use when:

- You only need a one-time literature search.
- The project does not have a stable topic profile yet.
- You cannot store API results or abstracts under the relevant source terms.

---

## Inputs

Required:

- `sources`: preprint sources to poll, such as `biorxiv`, `medrxiv`, `arxiv`, or `chemrxiv`.
- `profile`: a JSON relevance profile based on public-safe research interests.
- `state`: a JSON file that stores the most recent cursor/date for each source.

Optional:

- `keep_threshold`: score at or above which a record is kept without review.
- `triage_band`: score range for manual or LLM-assisted triage.
- `archive_dir`: location for monthly JSONL archives of seen records.
- `index_path`: local vector or keyword index for kept papers.

---

## Outputs

```text
<horizon-root>/
  state.json
  profiles/
    TEMPLATE.json
    <profile-name>.json
  archive/
    YYYY-MM.jsonl
  logs/
    YYYY-MM-DD.log
  index/
    kept-papers.index
```

---

## Profile Schema

Start from `assets/profiles/TEMPLATE.json` and replace the placeholder topics
with public-safe interests.

```json
{
  "profile_name": "example-profile",
  "version": 1,
  "updated": "YYYY-MM-DD",
  "keep_threshold": 0.72,
  "triage_band": [0.60, 0.72],
  "discard_below": 0.60,
  "interest_seeds": [
    "single-cell transcriptomics batch correction",
    "cell signaling disease models",
    "clinical variant interpretation standards"
  ],
  "hard_keep_keywords": [
    "single-cell", "mitochondria", "ClinVar"
  ],
  "hard_discard_keywords": [],
  "sources": ["biorxiv", "medrxiv", "arxiv", "chemrxiv"],
  "arxiv_categories": ["q-bio", "cs.AI", "cs.LG", "stat.ML"],
  "auto_refine": false
}
```

Keep private product names, unpublished target lists, customer names, and
internal project codenames out of public profiles.

---

## Standard Workflow

### Step 1 - Create a profile

```bash
cp assets/profiles/TEMPLATE.json profiles/example-profile.json
```

Edit `interest_seeds`, thresholds, and source lists.

### Step 2 - Run a dry run

```bash
python3 scripts/ingest.py \
  --sources biorxiv \
  --days 7 \
  --profile profiles/example-profile.json \
  --dry-run
```

Review kept and discarded examples before enabling scheduled runs.

### Step 3 - Run an ingest

```bash
python3 scripts/ingest.py \
  --sources all \
  --profile profiles/example-profile.json \
  --state state.json \
  --archive archive/ \
  --index index/kept-papers.index
```

### Step 4 - Schedule recurring runs

Use cron, launchd, systemd, or your scheduler of choice. Keep logs and generated
archives outside the git-tracked skill directory unless they are small synthetic
examples.

```bash
0 3 * * * /usr/bin/env bash /path/to/horizon/daily_ingest.sh
```

---

## Search

```bash
python3 scripts/search.py \
  --query "cell signaling disease model" \
  --index index/kept-papers.index \
  --k 10
```

Search output should include score, date, source, title, DOI, and an abstract
snippet.

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| Too many papers kept | Threshold too low or seed terms too broad | Raise `keep_threshold`; split broad interests into narrower phrases |
| Relevant papers discarded | Seed profile too narrow | Add positive seed examples and hard-keep keywords |
| Duplicate records | Same DOI appears across source/version feeds | Deduplicate by DOI and source version |
| Schedule reruns same window | State file not written or not mounted | Check `state.json` permissions and scheduler working directory |
| API rate limits | Polling too frequently | Add source-specific delay and retry handling |
| Index grows unexpectedly | Archives or logs stored inside index path | Separate `archive/`, `logs/`, and `index/` directories |
| Private terms in profile | Real project profile copied into public repo | Replace with synthetic or generic examples before commit |

---

## References

- bioRxiv API: `https://api.biorxiv.org/details/{server}/{date1}/{date2}/{cursor}/json`
- arXiv OAI-PMH: `https://oaipmh.arxiv.org/oai?verb=ListRecords&...`
- ChemRxiv public API: `https://chemrxiv.org/engage/chemrxiv/public-api/v1/items`
