---
id: skill_507f3118a19d5d5f84ef80323b06f804
name: chemrxiv-fetch
display-name: ChemRxiv Paper Fetcher
short-description: 'Library card skill — programmatic access to ChemRxiv chemistry
  preprints. Metadata via Crossref API (fully open). PDFs are behind Cloudflare: headless
  browser clears the landing page but in-page PDF fetch still returns 403. Headed
  browser (non-headless Chromium on a real desktop) is the expected solution — NOT
  YET VERIFIED. Use when an agent needs chemistry preprint metadata; flag PDF as requiring
  human drop or headed-browser node until Mode D is confirmed.'
category: methodology
keywords: chemrxiv, chemistry, preprint, doi, metadata, literature, fetch, cloudflare,
  browser, library-card
version: 1.0.0
last-updated: '2026-06-11'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
  agent: UViiVe
created: '2026-06-11'
changelog:
- version: 1.0.0
  date: 2026-06-11
  agent: UViiVe
  event: authored
  summary: Initial skill. Metadata via Crossref verified. PDF requires headed browser
    or manual download due to Cloudflare protection.
description: 'Library card skill — programmatic access to ChemRxiv chemistry preprints.
  Metadata via Crossref API (fully open). PDFs are behind Cloudflare: headless browser
  clears the landing page but in-page PDF fetch still returns 403. Headed browser
  (non-headless Chromium on a real desktop) is the expected solution — NOT YET VERIFIED.
  Use when an agent needs chemistry preprint metadata; flag PDF as requiring human
  drop or headed-browser node until Mode D is confirmed.'
---

# ChemRxiv Paper Fetcher

Library card for ChemRxiv chemistry preprints (DOI prefix `10.26434`).

**Metadata:** ✅ Crossref API — fully open, always works.
**PDFs:** ⚠️ Cloudflare-walled — stricter than bioRxiv. Headless browser not sufficient.
Headed browser (Mode D) is the expected fix but **not yet verified on ChemRxiv**.

---

## When to Use This Skill

✅ **Use when:**
- You need metadata/abstract for a ChemRxiv preprint
- Searching for chemistry preprints (organic, medicinal, computational, materials)
- Literature review touching drug design, synthesis routes, or computational chemistry

❌ **Do not use when:**
- The paper is a published journal article → `ncbi-api`
- The paper is bioRxiv/medRxiv → `biorxiv-fetch`
- The paper is arXiv → `arxiv-fetch` (fully open, much simpler)

---

## ⚠️ CRITICAL: ASK THIS FIRST

**Do you need metadata only, or the actual PDF?**

| Need | Path |
|---|---|
| Metadata + abstract | **Mode A** (Crossref — works now) |
| Full PDF | Mode D headed browser **(unverified on ChemRxiv)** or human drop |

Extract DOI from URL: `re.search(r'10\.26434/[^\s?#/]+(?:/v\d+)?', url).group()`

---

## Wall Status (tested 2026-06-11)

| Endpoint | Direct requests | Headless browser | Headed browser |
|---|---|---|---|
| Landing page | ❌ 403 | ✅ renders (hasChallenge=False) | Expected ✅ |
| PDF (`/doi/pdf/...`) | ❌ 403 | ❌ 403 from in-page fetch | **Untested — needs field report** |

ChemRxiv is stricter than bioRxiv. Headless clears the JS challenge but Cloudflare still blocks
the PDF fetch — likely checking TLS fingerprint or IP reputation beyond `cf_clearance`.
**A headed browser on a physical desktop is the recommended approach for PDF access.**

---

## Mode A — Metadata via Crossref (always works)

```python
import requests, re

def fetch_chemrxiv_meta(doi: str) -> dict:
    """
    Fetch ChemRxiv paper metadata via Crossref API.
    doi: e.g. '10.26434/chemrxiv-2025-xjkxl' or full URL
    """
    doi = re.sub(r'^https?://(?:doi\.org|dx\.doi\.org|chemrxiv\.org/[^?]+\?doi=)/?', '', doi.strip())
    doi = doi.rstrip('/')

    r = requests.get(
        f"https://api.crossref.org/works/{doi}",
        # Use your own email here — this opts you into Crossref's "polite pool"
        # for better rate limits and priority access.
        headers={"User-Agent": "my-agent/1.0 (mailto:your-email@example.com)"},
        timeout=15
    )
    r.raise_for_status()
    work = r.json().get("message", {})

    title    = work.get("title", ["N/A"])[0]
    authors  = [f"{a.get('given','')} {a.get('family','')}".strip()
                for a in work.get("author", [])]
    abstract = re.sub(r'<[^>]+>', '', work.get("abstract", "N/A"))
    posted   = work.get("posted", {}).get("date-parts", [[]])[0]
    date_str = "-".join(str(p) for p in posted) if posted else "N/A"

    return {
        "title":        title,
        "authors":      authors,
        "abstract":     abstract,
        "doi":          doi,
        "date":         date_str,
        "subject":      work.get("subject", []),
        "pdf_url":      f"https://chemrxiv.org/doi/pdf/{doi}",
        "landing_url":  f"https://chemrxiv.org/engage/chemrxiv/view-preprint/{doi}",
        "source":       "crossref",
    }
```

---

## Mode B — Search via Crossref

```python
def search_chemrxiv(query: str, max_results: int = 10) -> list:
    """Search ChemRxiv preprints via Crossref. Filter: DOI prefix 10.26434."""
    r = requests.get(
        "https://api.crossref.org/works",
        params={
            "query": query,
            "filter": "prefix:10.26434",
            "rows": max_results * 3,
            "sort": "published", "order": "desc",
        },
        # Use your own email here — this opts you into Crossref's "polite pool"
        # for better rate limits and priority access.
        headers={"User-Agent": "my-agent/1.0 (mailto:your-email@example.com)"},
        timeout=20,
    )
    r.raise_for_status()
    items = r.json().get("message", {}).get("items", [])

    results = []
    for work in items:
        doi = work.get("DOI", "")
        if "chemrxiv" not in doi.lower():
            continue
        posted = work.get("posted", {}).get("date-parts", [[]])[0]
        results.append({
            "doi":     doi,
            "title":   work.get("title", ["N/A"])[0],
            "authors": [f"{a.get('given','')} {a.get('family','')}".strip()
                        for a in work.get("author", [])],
            "date":    "-".join(str(p) for p in posted) if posted else "N/A",
            "pdf_url": f"https://chemrxiv.org/doi/pdf/{doi}",
        })
        if len(results) >= max_results:
            break

    print(f"✓ chemrxiv_search: {len(results)} results for '{query}'")
    return results
```

---

## Mode D — PDF via Headed Browser (UNVERIFIED on ChemRxiv)

⚠️ **This has not been confirmed to work on ChemRxiv.** Headless browser fails at the PDF fetch step (403) even after clearing the landing page challenge. Headed Chromium is the next test.

If you have a headed browser node available, adapt `biorxiv-fetch` Mode D:

```
Landing page: https://chemrxiv.org/engage/chemrxiv/view-preprint/{doi}
PDF URL:      https://chemrxiv.org/doi/pdf/{doi}
```

Same steps as biorxiv-fetch Mode D:
1. `browser(action=start)` — **must be headed (non-headless)**
2. Navigate to landing page, wait 5–8s for Cloudflare to clear
3. Verify `hasChallenge: false` AND that the real article title is visible
4. In-page `fetch(pdf_url, {credentials: 'include'})` inside the browser session
5. Verify: `status=200`, `Content-Type=application/pdf`, `>100KB`, magic `%PDF`
6. Persist via base64 → `file_write`

**If this works, add a field report. This is the key open question for ChemRxiv.**

**Fallback (reliable):** Ask the human to drop the PDF into Discord. 100% works.

---

## Common Issues

| Issue | Cause | Fix |
|---|---|---|
| Crossref 404 | DOI not yet indexed (~24h after posting) | Wait and retry; try `web_fetch` on chemrxiv.org instead |
| Abstract has JATS XML tags | Crossref raw response | `re.sub(r'<[^>]+>', '', abstract)` |
| PDF 403 from direct requests | Cloudflare wall | Mode D headed browser (unverified) or human drop |
| PDF 403 from headless browser in-page fetch | Cloudflare deeper check | Use headed browser on physical desktop node |
| DOI prefix mismatch in search | Old-format ChemRxiv DOIs | Some older DOIs use `chemrxiv.org/` instead of `10.26434/` |
| Landing URL 404 | URL format changed | Inspect chemrxiv.org network tab for current pattern |
| Author list incomplete | Crossref truncation | Full list in `work['author']` array |
| Abstract missing | Not deposited by author | `web_fetch` the landing page for scraped abstract |
| Rate limit on Crossref | >50 req/s | Add `mailto:` in User-Agent (already included) for polite pool |
| DOI version suffix handling | `/v1` suffix in some DOIs | Include `/v1` in requests — Crossref handles it |

---

## Suggested Next Steps
- **`docling-parse`** — parse PDF into structured text (once PDF is obtained)
- **`paperclip`** — deposit into the lab library
- **`arxiv-fetch`** — for arXiv preprints (fully open, no walls)
- **`biorxiv-fetch`** — for bioRxiv/medRxiv (Mode D PDF verified)
- **`semantic-scholar-api`** — citation graph

---

## References
- ChemRxiv: https://chemrxiv.org
- Crossref API: https://api.crossref.org/swagger-ui/index.html
- ChemRxiv DOI prefix: `10.26434`
- biorxiv-fetch Mode D (PDF browser method): `skills/biorxiv-fetch/SKILL.md`

---

## 🔄 Learning Loop — Field Reports

**Open question: does Mode D work with a headed browser on ChemRxiv?**
A headed browser on a physical desktop is the recommended approach for PDF access.

If you run a field test, add a report below and contribute your findings back to the skills library.

### 2026-06-11 — Initial field test
**Result:** ⚠️ PARTIAL
**Verified:** Crossref metadata API ✅ (metadata/abstracts returned correctly for multiple ChemRxiv DOIs).
**PDF status:** ❌ Direct requests → 403 (Cloudflare HTML). Headless browser → landing page renders (`hasChallenge=False`, real title), but in-page PDF fetch still 403.
**Finding:** Cloudflare appears to fingerprint beyond `cf_clearance` — TLS fingerprint or IP reputation check. Headed browser on a real desktop likely needed.
**Next step:** Test Mode D with headed Chromium on a physical desktop. Use DOI `10.26434/chemrxiv.15004614/v1` as a test case.

<!-- Field report template:
### [DATE] — [ENVIRONMENT]
**Result:** ✅ SUCCESS | ⚠️ PARTIAL | ❌ FAILED
**Time:** <how long>  **New blockers:** <new>  **Resolved by:** <fix>
**What to skip:** <unnecessary steps>  **Next note:** <one line>
-->
