---
id: skill_d8ac13f95359d2a45256d312676193b3
name: paperclip
display-name: Paperclip
short-description: Search and read from Paperclip — a virtual filesystem of 8M+ full-text
  biomedical papers (PubMed Central, bioRxiv, medRxiv, arXiv). Use for literature
  search, paper reading, figure analysis, corpus-wide queries, and synthesizing results
  across many papers. Auth via PAPERCLIP_API_KEY env var. One API key per account.
  Use when this workflow is needed.
category: methodology
keywords: paperclip, literature search, full text, biomedical papers, MCP, BYOK
version: 1.2.1
last-updated: '2026-08-19'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
provenance:
  wraps: Paperclip research library (paperclip.gxl.ai)
  api_docs: https://paperclip.gxl.ai
  citation: Paperclip is an external service — no primary paper. See https://paperclip.gxl.ai.
    BYOK.
changelog:
- version: 1.2.1
  date: '2026-07-07'
  event: citation-provenance
  summary: Added provenance block and References section. Vienna public subset.
- version: 1.2.0
  date: '2026-07-07'
  event: schema-compliance
  summary: 'De-agentified for public release: BYOK pattern, removed fleet key references
    and deployment table. Vienna public subset.'
- version: 1.1.0
  date: '2026-05-10'
  event: update
  summary: Added paper layout, recommended workflow, citation format, bash pipes,
    abstracts source.
- version: 1.0.0
  date: '2026-05-10'
  event: authored
  summary: Initial skill.
description: 'Search and read from Paperclip — a virtual filesystem of 8M+ full-text
  biomedical papers (PubMed Central, bioRxiv, medRxiv, arXiv). Use for literature
  search, paper reading, figure analysis, corpus-wide queries, and synthesizing results
  across many papers. Auth via PAPERCLIP_API_KEY env var. One API key per account.

  '
---

# Paperclip Skill

## What it is
Paperclip is a research library giving access to 8M+ biomedical papers across PubMed Central,
bioRxiv, medRxiv, and arXiv. It exposes an MCP endpoint at `https://paperclip.gxl.ai/mcp`.

## Authentication
- **Env var:** `PAPERCLIP_API_KEY` (must be set in your agent's environment config)
- **Auth method:** `x-api-key` header on all MCP requests
- **Obtain your key at:** https://paperclip.gxl.ai
- **Set PAPERCLIP_API_KEY** in your environment or in your agent's config file
- **One API key per account**
- If your key stops working, regenerate it at https://paperclip.gxl.ai

## How to call

### Via exec (recommended)
```bash
curl -s \
  -H "x-api-key: $PAPERCLIP_API_KEY" \
  -H "Content-Type: application/json" \
  -X POST \
  -d "{\"jsonrpc\":\"2.0\",\"method\":\"tools/call\",\"params\":{\"name\":\"paperclip\",\"arguments\":{\"command\":\"COMMAND_HERE\"}},\"id\":1}" \
  "https://paperclip.gxl.ai/mcp" \
  | python3 -c "import sys,json; print(json.load(sys.stdin)['result']['content'][0]['text'])"
```

Replace `COMMAND_HERE` with any Paperclip command.

## Command Reference

### Search & Discovery
| Command | Description |
|---------|-------------|
| `search "QUERY" -n N` | Semantic + keyword search, N results (default 10) |
| `search "QUERY" -s SOURCE` | Filter by source: `biorxiv`, `pubmed`, `arxiv`, `medrxiv` |
| `search "QUERY" --since 2023` | Filter by year |
| `search "QUERY" --author "Name"` | Filter by author |
| `grep PATTERN /papers/` | Regex search across full corpus |
| `lookup doi 10.1234/...` | Find paper by DOI |
| `lookup pmid 12345678` | Find paper by PubMed ID |
| `sql "SELECT title, year FROM documents WHERE ..."` | SQL query (200-row limit) |

### Paper Layout

Each paper lives at `/papers/<id>/` with this structure:

| Path | Contents |
|------|----------|
| `meta.json` | Metadata: title, authors, abstract, DOI, year, journal |
| `content.lines` | Full text, **line-numbered** — prefer this over full.txt |
| `sections/` | Named section files (Introduction, Methods, Results, etc.) |
| `figures/` | Figure files (use with `ask-image`) |
| `supplements/` | Supplementary material |

ID prefixes: `PMC` (PubMed Central), `bio_` (bioRxiv), `med_` (medRxiv), `arx_` (arXiv)

`/.gxl/` is a writable scratch space for intermediate results.
`/papers/` is read-only.

### Reading Papers
| Command | Description |
|---------|-------------|
| `cat /papers/bio_ID/meta.json` | Paper metadata (title, authors, abstract, DOI) |
| `cat /papers/bio_ID/meta.json \| jq .title` | Extract specific field |
| `head -40 /papers/bio_ID/content.lines` | First 40 lines — prefer over full.txt |
| `cat /papers/bio_ID/sections/Results` | Read a specific named section |
| `grep -i "term" /papers/bio_ID/content.lines` | Search within a paper |
| `scan /papers/bio_ID/content.lines "term1" "term2"` | Multi-pattern search within paper |
| `ask-image /papers/bio_ID/figures/fig1.png "q"` | Analyze a figure with vision |

### Analysis (Power Features)
| Command | Description |
|---------|-------------|
| `map --from SEARCH_ID "question"` | Ask same question across all search results (~15s, 20 papers parallel) |
| `reduce --from SEARCH_ID "question"` | Synthesize map results — strategies: `summarize`, `table`, `themes` |
| `ask-image /papers/bio_ID/fig1.png "question"` | Analyze figure with vision |
| `results SEARCH_ID` | View saved results from a previous search/map |

## Recommended Workflow

The server's own recommended pattern:

1. **Find:** `search "topic"` → get search ID from results
2. **One paper:** `head` / `grep` / `scan` on `/papers/<id>/content.lines`
3. **Many papers:** `search -n 10` → `map --from ID "question"` → synthesize directly
4. **Stats:** `sql "SELECT ..."`

**After `map`, synthesize directly — don't re-read papers individually.**

## Typical Workflows

### Quick literature search
```
search "VEGF-A peptide inhibitor UV redness" -n 10
```

### Deep sweep — question across 20 papers
```
search "VEGF-A angiogenesis peptide" -n 20   # get search ID from results
map --from s_XXXX "Does this paper describe a short peptide binding VEGF-A directly?"
reduce --from s_XXXX "table"
```

### Read a specific paper
```
head -40 /papers/bio_XXXX/content.lines
cat /papers/bio_XXXX/sections/Methods
grep -i "IC50" /papers/bio_XXXX/content.lines
```

### Pipes and bash
```
bash 'search "protein design" | grep -i "diffusion"'
```

### Find by author
```
search "protein design" --author "David Baker" -n 5
```

### Broader abstract-only search (larger corpus, no full text)
```
search "QUERY" -s abstracts
```

### Patent + literature combined
Pair with `epo-ops` and `uspto-api` skills for IP landscape work:
1. Paperclip for prior art literature
2. EPO OPS for patent family search
3. USPTO for US prosecution history

## Citations

Cite inline as **[1]**, **[2]**, etc. End responses with a references block:

```
--------
REFERENCES
[1] Authors. "Title." *Journal* vol, pages (year). doi:XX
    https://citations.gxl.ai/papers/<doc_id>#L<n>
```

Never expose raw `doc_id` in prose. Cite every factual claim.

## First Law Note
Paperclip returns real paper metadata and text. Do not paraphrase results in a way that
changes meaning. Quote titles and abstracts verbatim. Cite paper IDs (bio_XXXX) and DOIs.
Never claim a paper says something it doesn't — read the full text before asserting conclusions.
Use `head`/`grep`/`scan` on `content.lines` before making claims about paper content.

## Version History
- **v1.2.1** (2026-07-07) — Added provenance block and References section. Vienna public subset.
- **v1.2.0** (2026-07-07) — De-agentified for public release: BYOK pattern, removed fleet key references and deployment table.
- **v1.1.0** (2026-05-10) — Added paper layout (content.lines, sections/, figures/, supplements/, /.gxl/ scratch), recommended workflow (search→map→synthesize), citation format, bash pipes, abstracts source. Corrected full.txt → content.lines.
- **v1.0.0** (2026-05-10) — Initial skill build.

---

## References

- **Paperclip platform**: External service — no primary paper. See https://paperclip.gxl.ai. BYOK.
- **API endpoint**: https://paperclip.gxl.ai/mcp
