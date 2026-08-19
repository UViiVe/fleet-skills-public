# Paperclip Quick Reference

## One-liner to call any command

```bash
paperclip() {
  curl -s \
    -H "x-api-key: $PAPERCLIP_API_KEY" \
    -H "Content-Type: application/json" \
    -X POST \
    -d "{\"jsonrpc\":\"2.0\",\"method\":\"tools/call\",\"params\":{\"name\":\"paperclip\",\"arguments\":{\"command\":\"$*\"}},\"id\":1}" \
    "https://paperclip.gxl.ai/mcp" \
    | python3 -c "import sys,json; print(json.load(sys.stdin)['result']['content'][0]['text'])"
}
```

> Paste that function in your shell session, then call `paperclip search "VEGF" -n 5`.

---

## Most-used commands at a glance

| What you want | Command |
|---|---|
| Quick search | `search "QUERY" -n 10` |
| Source-filtered search | `search "QUERY" -s biorxiv` |
| Find by PMID | `lookup pmid 12345678` |
| Find by DOI | `lookup doi 10.1101/2024.01.01.000000` |
| Read abstract + title | `cat /papers/bio_ID/meta.json` |
| Read full text | `cat /papers/bio_ID/full.txt` |
| Search inside a paper | `scan /papers/bio_ID/full.txt "term1" "term2"` |
| Regex across corpus | `grep "VEGF.{0,20}peptide" /papers/` |
| SQL query | `sql "SELECT title, year FROM documents WHERE year > 2022 LIMIT 20"` |
| Analyze a figure | `ask-image /papers/bio_ID/fig1.png "What does panel B show?"` |

---

## Power workflow: sweep 20 papers for one question (~15s)

```bash
# Step 1 — search and note the search ID (printed as s_XXXX in results)
search "VEGF-A short peptide inhibitor" -n 20

# Step 2 — map one question across all 20 results in parallel
map --from s_XXXX "Does this paper report a peptide shorter than 20 AA that binds VEGF-A?"

# Step 3 — synthesize into a table
reduce --from s_XXXX "table"

# (Optional) view saved results later
results s_XXXX
```

---

## Source codes
| Code | Database |
|---|---|
| `pubmed` | PubMed Central (full-text OA) |
| `biorxiv` | bioRxiv preprints |
| `medrxiv` | medRxiv preprints |
| `arxiv` | arXiv (cross-disciplinary) |

---

## Worked example — UViiVe literature sweep

```bash
# Find papers on UV-induced VEGF and vascular redness
search "ultraviolet VEGF-A redness erythema vascular" -n 15

# Map across results
map --from s_XXXX "Does this paper report a mechanism linking UV exposure to VEGF-A upregulation in skin?"

# Synthesize
reduce --from s_XXXX "summarize"
```

---

## Combine with IP tools

1. `paperclip search "VEGF peptide" -n 20` → find prior art literature
2. `epo-ops` skill → patent family search
3. `uspto-api` skill → US prosecution history
4. Triangulate: anything claimed in both literature AND patents needs special handling

---

*See SKILL.md for full auth setup and command reference.*
