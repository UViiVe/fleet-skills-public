---
id: "skill_8f1937e238824ba4a18e1d1476427a83"
name: "arxiv-fetch"
short-description: "Library card skill - programmatic access to arXiv preprints. Both metadata (Atom API) and PDFs are fully open: no Cloudflare, no auth, no browser required. Direct HTTP always works. Use when an age..."
display-name: "arXiv Paper Fetcher"
description: "Library card skill - programmatic access to arXiv preprints. Both metadata (Atom API) and PDFs are fully open: no Cloudflare, no auth, no browser required. Direct HTTP always works. Use when an agent needs any arXiv paper by ID, wants to search recent preprints, or needs to download a PDF."
category: "methodology"
keywords: "arxiv, preprint, pdf, metadata, literature, fetch, library-card, cs, physics, math, quantbio, open-access"
version: "1.0.1"
created: "2026-06-11"
last-updated: "2026-07-07"
provenance:
  wraps: "arXiv (Cornell University)"
  api_docs: "https://info.arxiv.org/help/api/index.html"
  citation: "arXiv is maintained by Cornell University. No primary paper - see https://arxiv.org"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section."
  - version: 1.0.0
    date: 2026-06-11
    event: authored
    summary: >-
      Initial skill. Verified direct metadata and PDF access using public arXiv endpoints.
author: "Fleet Skills contributors"
status: active
promoted: true
license: MIT
---

# arXiv Paper Fetcher

Library card for arXiv. **Fully open - no walls anywhere.** Direct HTTP works.
No browser automation needed (though Mode D also works if you have it running).

---

## Learning Loop - Field Reports

Append field reports only when they add reusable public information, such as a changed API behavior, rate-limit change, or parser update. Do not include personal machine names, private paths, or internal project context.
