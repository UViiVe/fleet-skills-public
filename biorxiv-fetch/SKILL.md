---
id: "skill_f269db8bc6a24ab780a04f60d459e5e8"
name: "biorxiv-fetch"
short-description: "Library card skill - programmatic access to bioRxiv/medRxiv preprint METADATA (title, authors, abstract, DOI, date) via the official REST API, AND a VERIFIED browser-automation method to download t..."
display-name: "bioRxiv / medRxiv Paper Fetcher"
description: "Library card skill - programmatic access to bioRxiv/medRxiv preprint METADATA (title, authors, abstract, DOI, date) via the official REST API, AND a VERIFIED browser-automation method to download the full PDF past Cloudflare. The metadata API and PDF host are separate servers; only the PDF host is Cloudflare-walled. Direct curl/requests of .full.pdf ALWAYS returns ~5KB Cloudflare HTML (silent-failure trap). v2 adds the proven PDF path: headed browser clears the JS challenge on the landing page, then an IN-PAGE fetch(credentials:include) pulls real PDF bytes. Use when an agent needs to discover, identify, summarize, OR download a preprint."
category: "methodology"
keywords: "biorxiv, medrxiv, preprint, doi, metadata, literature, fetch, pdf, cloudflare, browser, library-card"
version: "2.0.1"
created: "2026-06-11"
last-updated: "2026-07-07"
provenance:
  wraps: "bioRxiv / medRxiv (Cold Spring Harbor Laboratory)"
  api_docs: "https://api.biorxiv.org"
  citation: "Abdill RJ, Blekhman R. \"Tracking the popularity and outcomes of all bioRxiv preprints.\" eLife. 2019;8:e45133. doi:10.7554/eLife.45133"
changelog:
  - version: 2.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section."
  - version: 2.0.0
    date: 2026-06-11
    event: major
    summary: >-
      Added a VERIFIED, first-class PDF download method (Mode D) using browser
      automation. Documented a browser-session PDF download path for Cloudflare-protected PDF pages and the validation checks needed to avoid saving HTML error pages as PDFs.
  - version: 1.1.0
    date: 2026-06-11
    event: correction
    summary: >-
      Corrected scope: this skill fetches METADATA only, not PDFs (at v1).
      Direct PDF requests can return Cloudflare HTML; metadata remains available through the official API.
  - version: 1.0.0
    date: 2026-06-11
    event: authored
    summary: >-
      Initial skill. Wraps the official bioRxiv REST API for metadata/discovery.
      Initial public metadata workflow.
author: "Fleet Skills contributors"
status: active
promoted: true
license: MIT
---

# bioRxiv / medRxiv Paper Fetcher (v2)

Fetch preprints from bioRxiv or medRxiv. **Metadata via the official REST API** (always works),
and **full PDF via verified browser automation** (clears Cloudflare).

> **v2 headline:** The PDF download problem is SOLVED. See **Mode D**. Do NOT use the
> deprecated curl/requests `.full.pdf` approach - it silently returns ~5KB of Cloudflare
> HTML that looks like a download but isn't a PDF.

---

## When to Use This Skill

 **Use when:**
- You have a DOI / biorxiv URL and need metadata, abstract, OR the full PDF
- You need to search for recent preprints by keyword or date range
- A user pastes a biorxiv.org URL and wants a summary, analysis, or the file
- Running a literature scan for a project
- Depositing a preprint into a literature library

 **Do not use when:**
- The paper is in PubMed / a published journal - use `ncbi-api` or `semantic-scholar-api`
- You need bulk corpus access (100k+ papers) - use the S3 TDM bucket `s3://biorxiv-src-monthly`

---

##  ASK THIS FIRST

**What do you have, and do you need the PDF or just the text?**

| You have… | You need… | Mode |
|---|---|---|
| DOI or biorxiv URL | metadata + abstract | **A** (API) |
| DOI or biorxiv URL | full text (no figures) | **A** then `web_fetch` the `.full` HTML |
| DOI or biorxiv URL | the actual PDF file | **D** (browser - verified) |
| search terms / topic | discovery | **B** (search) |
| date range | batch metadata | **C** (range) |

Extract DOI from a URL: `re.search(r'10\.1101/[^\s/v]+', url).group()` → strip version with `re.sub(r'v\d+$','',doi)`.

---

## Access Notes

bioRxiv exposes papers through **two completely separate servers**:

1. **`api.biorxiv.org`** - metadata/discovery REST API. NOT Cloudflare-walled. Returns title, authors, abstract, DOI, date, category, version as JSON. **Always works from curl/requests.**
2. **`www.biorxiv.org/content/...full.pdf`** - the PDF host. **Cloudflare bot-walled.** Any direct curl/requests/web_fetch returns ~5KB of `Just a moment...` HTML with a 200 or 403 - never the PDF. The `cf_clearance` cookie that authorizes access is set by solving a JS challenge in a real browser and is **httpOnly** (so you cannot scrape it out and reuse it in curl).

➡️ Metadata = API. PDF = browser. Do not cross the streams.

---

## Mode A - Metadata + Abstract by DOI (always works)

```python
import requests, re

def fetch_biorxiv_meta(doi: str) -> dict:
    doi_clean = re.sub(r'v\d+$', '', doi.strip())
    for server in ["biorxiv", "medrxiv"]:
        r = requests.get(f"https://api.biorxiv.org/details/{server}/{doi_clean}/na/json", timeout=15)
        collection = r.json().get("collection", [])
        if collection:
            p = collection[-1]
            v = p.get("version", "1")
            return {"title": p.get("title"), "authors": p.get("authors"),
                    "abstract": p.get("abstract"), "doi": doi_clean, "date": p.get("date"),
                    "category": p.get("category"), "version": v, "server": server,
                    "html_url": f"https://www.biorxiv.org/content/{doi_clean}v{v}.full",
                    "pdf_url":  f"https://www.biorxiv.org/content/{doi_clean}v{v}.full.pdf"}
    raise ValueError(f"No paper found for DOI: {doi_clean}")
```

For **full text without figures**, take `html_url` and use the `web_fetch` tool
(`extractMode="markdown"`) - the HTML full-text page is generally NOT walled even though the PDF is.

---

## Mode B - Search by Query  •  Mode C - Batch by Date Range

(unchanged from v1 - both use `api.biorxiv.org/details/{server}/{start}/{end}/{cursor}/json`,
client-side keyword filter for search. See `scripts/biorxiv_api.py`.)

---

## Mode D - Download the PDF via Browser (VERIFIED 2026-06-11)

**This is the canonical agent path for getting the actual PDF.** Verified end-to-end on
Overath 2025 (`10.1101/2025.08.14.670059`): status 200, `application/pdf`, 19,173,446 bytes, `%PDF-1.5`.

### Mechanism (why it works)
1. A **headed** Chromium session loads the **landing page** (`.../{DOI}v{N}` - NOT the PDF URL). This triggers and *solves* the Cloudflare JS challenge, setting the httpOnly `cf_clearance` cookie on the session.
2. Once the landing page renders the real article (`hasChallenge:false`), an **in-page `fetch(pdfUrl, {credentials:'include'})`** runs *inside* that authorized session and returns the real PDF bytes.
3. Persist the bytes (base64 → `file_write`, or the browser tool's `pdf`/download action).

>  The curl-with-cookies shortcut FAILS: `cf_clearance` is httpOnly, so `document.cookie` does not expose it and curl cannot replay it. The fetch must happen **in the browser**.

### Step-by-step with the `browser` tool

```
# 1. Start headed browser (headed clears Cloudflare more reliably than headless)
browser(action=start)

# 2. Navigate to the LANDING PAGE (not the .pdf) to clear the challenge
browser(action=navigate, url="https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1",
        loadState="domcontentloaded")
browser(action=act, request={"kind":"wait","timeMs":4000})   # let challenge resolve

# 3. Confirm the challenge cleared + grab the canonical PDF link
browser(action=act, targetId=<id>, request={"kind":"evaluate",
  "fn":"() => ({ title: document.title, hasChallenge: document.body.innerText.includes('Just a moment') || document.body.innerText.includes('Verifying'), pdfLink: document.querySelector('a[href*=\".full.pdf\"]')?.href || null })"})
#   -> expect hasChallenge:false and a real title before proceeding

# 4. In-page fetch the PDF and hand bytes back as base64 (chunk to avoid call-stack overflow)
browser(action=act, targetId=<id>, request={"kind":"evaluate",
  "fn":"async () => { const r = await fetch('https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1.full.pdf',{credentials:'include'}); const b=await r.arrayBuffer(); const u=new Uint8Array(b); let s=''; const c=8192; for(let i=0;i<u.length;i+=c){ s+=String.fromCharCode.apply(null,u.subarray(i,i+c)); } return { status:r.status, type:r.headers.get('content-type'), bytes:b.byteLength, b64: btoa(s) }; }"})
#   -> verify status:200, type:'application/pdf', magic %PDF (first bytes of b64 decode), THEN:

# 5. Persist: pass result.b64 to file_write(contentBase64=..., path="/tmp/papers/<doi>.pdf")
#    OR use browser(action=pdf, ...) if you only need a print-render.

# 6. Stop the browser
browser(action=stop)
```

### Verification gate (MANDATORY - First Law)
A real PDF is **`application/pdf`, magic bytes `%PDF`, and typically >100 KB.**
If you get `text/html` or <20 KB, you got the Cloudflare page - the download FAILED.
Never report a paper as fetched without checking magic bytes + size. (This is the exact
silent-failure that bit the v1 approach: a 5,666-byte "PDF" that was an error page.)

### Fallbacks (in order)
1. **Mode D browser** (above) - primary agent path.
2. **Human-in-the-loop drop** - ask the human to drop the PDF into chat. 100% reliable; best for one-offs or when no browser node is available.
3. **Europe PMC / PMC mirror** - works but only AFTER indexing (days-to-weeks lag). Useless for fresh preprints (e.g. a paper posted this week).

---

## Common Issues

| Issue | Cause | Fix |
|---|---|---|
| `collection: []` | DOI not found on that server | Try the other server; strip `vN` suffix |
| "PDF" is ~5KB HTML / `Just a moment` | Hit the Cloudflare PDF host directly | Use **Mode D** (browser in-page fetch). NEVER curl the `.full.pdf`. |
| Mode D still shows challenge | JS not finished | Increase wait to 6–8s; ensure **headed** (not headless); re-check `hasChallenge` before fetching |
| `btoa` / call-stack error on large PDF | String too long | Chunk the Uint8Array (8192) before `btoa` (shown in Mode D) |
| curl with `document.cookie` fails | `cf_clearance` is httpOnly | Cannot scrape it - fetch must run in-browser |
| HTML full-text empty | Paper just posted | Wait 24–48h after posting |
| Rate limit (429) on API | Too many metadata calls | `time.sleep(1)` between calls; S3 bucket for bulk |

---

## Suggested Next Steps
- **`docling-parse`** - parse the fetched PDF into structured text
- **`paperclip`** - deposit into the lab library with metadata
- **`arxiv-search`** / **`ncbi-api`** / **`semantic-scholar-api`** - other sources

---

## References

- **bioRxiv**: Abdill RJ, Blekhman R. "Tracking the popularity and outcomes of all bioRxiv preprints." *eLife*. 2019;8:e45133. doi:10.7554/eLife.45133
- **API documentation**: https://api.biorxiv.org
- **bioRxiv TDM S3 bucket**: https://www.biorxiv.org/tdm
- **PDF URL pattern**: `https://www.biorxiv.org/content/{DOI}v{N}.full.pdf` (browser-only access)
- **Browser automation tool**: navigate → evaluate(in-page fetch) → file_write(base64)

---

## 🔄 Learning Loop - Field Reports

Every agent that runs this skill in production appends a field report, bumps `version:`,
adds a `changelog:` entry when the note is reusable.

<!-- Field report template:
### [AGENT] [DATE] - [ENVIRONMENT]
**Result:**  SUCCESS |  PARTIAL |  FAILED
**Time:** <how long>  **New blockers:** <new>  **Resolved by:** <fix>
**What to skip:** <unnecessary steps>  **Next agent note:** <one line>
-->
