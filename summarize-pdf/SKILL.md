---
id: skill_3e6fd55322f833b10f2ca195072268db
name: summarize-pdf
display-name: Summarize PDF
short-description: Summarize PDF reusable workflow. Use when this workflow is needed.
category: methodology
keywords: PDF, text extraction, summarize, citation verification, literature review
version: 1.0.0
last-updated: '2026-04-26'
status: active
license: MIT
author:
  human: Stephen C. Ekker, PhD / UViiVe
provenance:
  tool: steipete/summarize
  repo: https://github.com/steipete/summarize
  version-tested: 0.13.1 (brew) / 0.14.0 (adds local PDF --extract without LLM)
  steve-note: PDF extraction was the lead for me — 2026-04-26
  validated-on: arm64, macOS Sequoia
changelog:
- version: 1.0.0
  date: '2026-07-07'
  event: schema-compliance
  summary: 'Backfilled public schema. De-agentified: removed fleet deployment table
    and internal agent references. Vienna public subset.'
---

# summarize-pdf — PDF Text Extraction Skill

**Version:** 1.0.0 | **Date:** 2026-04-26

## What This Does

Extracts clean markdown or plain text from PDF files (local or URL) using `summarize` CLI — without requiring an LLM call. Useful as a cheap pre-processing step before any agent reads a paper, verifies a citation, or feeds text into a pipeline.

**Primary use cases:**
- Citation verification gate (extract → verify title/author/abstract)
- Periodic literature scouting pipeline paper intake (extract → summarize with LLM)
- Literature review pre-processing (extract → feed to analysis)
- Any agent reading papers — saves tokens vs. sending raw PDF to LLM

---

## Installation

```bash
# macOS (Homebrew) — installs v0.13.x
brew install summarize

# For v0.14.0+ (adds local PDF --extract without LLM):
npm install -g @steipete/summarize
# or download binary from:
# https://github.com/steipete/summarize/releases/tag/v0.14.0
```

**Verify:**
```bash
summarize --version
# → 0.13.1 or higher
```

**No API key needed for extraction-only use.** API keys only required for LLM summarization features.

---

## Core Usage

### Extract text from a PDF URL (no LLM)
```bash
summarize "https://arxiv.org/pdf/PAPER_ID" --extract --format md
# Returns: clean markdown text of the paper
```

### Extract text from a local PDF file (v0.14.0+)
```bash
summarize --extract /path/to/paper.pdf
# Returns: clean text without any LLM call
```

### Summarize a PDF with LLM (uses ANTHROPIC_API_KEY or similar)
```bash
summarize "https://arxiv.org/pdf/PAPER_ID" --format md --model anthropic/claude-sonnet-4-5
```

### Pipe output to file
```bash
summarize "https://arxiv.org/pdf/PAPER_ID" --extract --format md > /tmp/paper_text.md
```

---

## Agent Pattern — Paper Intake

Use this pattern in any agent reading a paper:

```bash
# Step 1: Extract text (cheap, no LLM)
summarize "https://arxiv.org/pdf/PMCID_OR_ARXIV_ID" --extract --format md > /tmp/paper.md

# Step 2: Verify citation metadata
head -50 /tmp/paper.md  # Check title, authors, abstract

# Step 3: Feed to LLM analysis (now you control what goes in)
# Pass /tmp/paper.md content to your analysis prompt
```

**Why this pattern beats sending raw PDF to LLM:**
- Extraction is free (no tokens)
- You see the text before sending — catch bad PDFs early
- Smaller, cleaner input → better LLM output
- Works offline for local files (v0.14.0+)

---

## Integration Points

### Citation Verification Gate
```bash
# Before confirming any citation:
summarize "$PDF_URL" --extract --format md | head -100
# → verify: title matches, authors present, abstract relevant
# → if mismatch: EXCLUDE, never silently drop
```

### Periodic Literature Scouting Pipeline Paper Intake
```bash
# In the literature scouting pipeline, replace direct PDF-to-LLM with:
TEXT=$(summarize "$PDF_URL" --extract --format md)
# Then pass $TEXT to summarization prompt
```

### Fetching papers by PMID / arXiv ID
```bash
# After fetching a PMID or arXiv ID, extract full text:
summarize "https://arxiv.org/pdf/${ARXIV_ID}" --extract --format md
# or for PMC papers:
summarize "https://ftp.ncbi.nlm.nih.gov/pub/pmc/${PMC_PATH}" --extract --format md
```

---

## Version Notes

| Version | Key difference |
|---------|---------------|
| 0.13.x (brew) | URL PDFs via `--extract --format md` ✅ · Local PDF `--extract` ❌ (media files only) |
| 0.14.0+ (npm/binary) | Local PDF `--extract` ✅ · Reddit thread JSON extraction · GPT-5.5 Fast presets |

**For most use cases (URL-based paper extraction): brew v0.13.x is sufficient.**
Install v0.14.0 only if you need local PDF extraction without LLM.

---

## Known Limitations

- **Scanned PDFs** (image-only): extraction returns garbage — needs OCR. Use `--slides-ocr` flag or fall back to the `pdf` tool.
- **Very large PDFs** (>20MB): may timeout on slow connections. Download first, then extract.
- **PMC HTML stubs**: some PMC PDF URLs return HTML not PDF (~1,817 bytes). Validate with `wc -c` > 10KB before extracting.
- **v0.13.x local files**: `--extract` only works on media files (MP3/MP4/WAV), NOT PDFs. Upgrade to v0.14.0 for local PDF support.

---

## Why This Exists

Steve's note (2026-04-26): *"PDF extraction was the lead for me... basically any agent reading papers would find it useful."*

Before this skill, agents were either:
1. Sending raw PDF bytes to LLM (expensive, lossy)
2. Using the `pdf` tool (good but synchronous, token-heavy for extraction-only tasks)
3. Fetching HTML stubs that weren't actual PDFs

`summarize --extract` gives a clean, free pre-processing step that any paper-reading agent can use as a first pass before deciding whether to burn LLM tokens on analysis.
