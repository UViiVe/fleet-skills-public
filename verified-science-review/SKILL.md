---
name: verified-science-review
id: "skill_15cd3d478c724870a69fa60e5fc06415"
display-name: "Verified Science Review"
short-description: "High-confidence systematic literature review pipeline for biomedical research. Multi-pass architecture: retrieve → verify citations → synthesize → audit. Every PMID verified by direct PubMed fetch...."
category: methodology
keywords: [science review, verification, citations, methodology, gate]
last-updated: "2026-07-07"
version: "2.0.0"
status: active
license: MIT
description: "High-confidence systematic literature review pipeline for biomedical research. Multi-pass architecture: retrieve → verify citations → synthesize → audit. Every PMID verified by direct PubMed fetch. Use when building disease reviews, grant background sections, drug landscape analyses, or any scientific literature synthesis where citation accuracy is critical."
author: "Fleet Skills contributors"
provenance:
  wraps: "PubMed E-utilities, Europe PMC, Semantic Scholar"
  upstream_source: "original work"
  citation: "N/A — original pipeline. Cite underlying databases: Sayers et al., NAR 2022 (doi:10.1093/nar/gkab1112)"
changelog:
  - version: "2.0.0"
    date: "2026-07-07"
    event: public-release
    summary: >-
      Generalized for public release. Core verification pipeline preserved
      intact. New general-purpose templates added for biomedical review use
      cases. MIT license.
  - version: "1.0.0"
    date: "2026-06-01"
    event: authored
    summary: "Initial systematic review pipeline."
promoted: true
---

# Verified Science Review Pipeline

A multi-pass literature review pipeline built on one rule: **every citation must be verified before it appears in the output.** No hallucinated PMIDs. No fabricated drug approvals. No plausible-sounding fiction.

Works with any biomedical topic — disease reviews, grant backgrounds, drug landscapes, gene function summaries, collaborator briefings.

---

## When to Use This Skill

✅ **Use when:**
- Building a comprehensive literature review that will be read by scientists or clinicians
- Writing grant background sections where citation accuracy matters
- Summarizing the therapeutic landscape for a disease or target
- Any synthesis task where a hallucinated PMID would be a real problem

❌ **Do not use when:**
- You need a quick 2-paragraph summary (use `summarize-pdf` or a direct LLM query instead)
- Speed matters more than accuracy
- The topic has fewer than ~10 papers (just fetch them directly)

---

## ⚠️ CRITICAL: ASK THIS FIRST

**What is the scope?** A single gene, a disease, a drug class, or a broader topic?

Scope determines which template to use and how many subagent passes to run. Broad scope + single pass = citation errors. Narrow scope + single pass = fine.

---

## Core Principle — The Verification Rule

**Every PMID must be verified by direct HTTP fetch before inclusion.**

```
Fetch: https://pubmed.ncbi.nlm.nih.gov/[PMID]/
Confirm: title and first author match the claim being cited
```

If the paper doesn't match — reject it. Flag with ⚠️ if unverified. Never include a PMID you haven't fetched.

This rule eliminates the most common and most damaging AI literature failure mode: confident fabrication of plausible-sounding citations.

---

## Pipeline Overview

Four passes. Each catches what the previous missed:

```
1. Briefing       → structured context injected before the review starts
2. Review pass    → full literature synthesis with per-citation verification
3. Audit pass     → independent check: are all cited PMIDs real and on-topic?
4. Corrections    → log any errors caught; improve future briefings
```

For high-stakes reviews, run passes 2 and 3 as separate subagents (parallel or sequential).

---

## Quick Start

### Single-topic review
```
Task for review subagent:
"Write a comprehensive literature review of [TOPIC].
Verify every PMID by fetching https://pubmed.ncbi.nlm.nih.gov/[PMID]/.
Confirm title and authors match before including.
Flag any unverified citations with ⚠️.
Do NOT fabricate PMIDs. Do NOT claim drug approvals without verification.
Use the structure: Overview → Genetics/Mechanism → Models → Therapeutics → Open Questions."
```

### Batch reviews (multiple genes/diseases)
Spawn up to 5 review subagents in parallel. Each gets the same verification instructions. Collect and audit outputs.

---

## Templates

See `references/` for ready-to-use briefing templates:

| Template | Use for |
|---|---|
| `general-review-template.md` | Any gene, disease, or biological topic |
| `grant-section-template.md` | Grant background / significance sections |
| `drug-landscape-template.md` | Therapeutic landscape for a disease or target |
| `quality-protocol.md` | Full PMID verification procedure + common hallucination patterns |

---

## Standard Workflow

### Step 1 — Fill the briefing
Choose a template from `references/`. Fill in:
- Gene/disease name, key identifiers (OMIM, gene ID, etc.)
- Known context: inheritance, prevalence, key papers, relevant models
- Scope: how deep, how many citations expected
- Always include the verification warning block (copy from template)

### Step 2 — Run the review pass
Spawn a review subagent with the filled briefing as its task. For long reviews, set `runTimeoutSeconds: 300`.

### Step 3 — Audit pass
After the review lands, run a second pass:
```
"Audit the following literature review. For each cited PMID,
fetch https://pubmed.ncbi.nlm.nih.gov/[PMID]/ and confirm
the title matches the claim. List any mismatches or unverified citations."
```

### Step 4 — Log corrections
Any errors caught → note them. Future briefings on the same topic can inject known corrections upfront, reducing errors progressively.

---

## PMID Verification Protocol

See `references/quality-protocol.md` for full procedure.

**Quick version:**

| Status | Meaning |
|---|---|
| ✓ Verified | Title/authors confirmed by direct PubMed fetch |
| ⚠️ Secondary | Cited from a review paper, not directly fetched |
| ❌ Rejected | Fetched paper is wrong topic — state what it actually is |

**Never include a PMID with ❌ status. Always state what the paper actually covers.**

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| Subagent fabricates PMIDs | No verification instruction | Always include the PMID warning block from template |
| PubMed fetch returns 404 | PMID doesn't exist | ❌ Reject immediately — hallucinated citation |
| Title mismatch | Real PMID, wrong paper | ❌ Reject — model associated wrong paper with claim |
| Drug approval claim unverified | Model confabulating | Verify at FDA.gov or EMA.europa.eu before including |
| Review too broad, low quality | Scope not defined | Narrow scope in briefing; one gene/disease per pass |
| Subagent times out on large review | Too many citations to verify | Split into sections; verify in batches of 20 PMIDs |
| Duplicate citations across sections | No deduplication pass | Add dedup step: "List all PMIDs cited more than once" |
| Secondary sources dominating | Insufficient primary lit | Add to briefing: "Prefer primary research over review articles" |
| Missing recent papers | Training cutoff | Explicitly instruct: "Search PubMed for papers from 2023–2025" |
| Inconsistent citation format | Mixed styles | Specify format in briefing: "Use [PMID: XXXXXXXX] inline" |

---

## Suggested Next Steps

- **`abstract-lookup`** — waterfall abstract fetch before running the review
- **`ncbi-api`** — PubMed search to seed the PMID list before verification
- **`europepmc-search`** — full-text search to catch papers PubMed misses
- **`semantic-scholar-api`** — citation graph to find highly-cited related work
- **`pdf-report-generation`** — format the completed review as a PDF

---

## References

**Bundled:**
- `references/general-review-template.md` — general-purpose briefing template
- `references/grant-section-template.md` — grant background section template
- `references/drug-landscape-template.md` — therapeutic landscape template
- `references/quality-protocol.md` — PMID verification procedure + hallucination patterns

**Primary citations:**
- **PubMed / NCBI E-utilities**: Sayers EW, et al. "Database resources of the National Center for Biotechnology Information." *Nucleic Acids Research*. 2022;50:D20–D26. doi:10.1093/nar/gkab1112
- **Europe PMC**: Europe PMC Consortium. "Europe PMC: a full-text literature database for the life sciences." *Nucleic Acids Research*. 2015;43:D1042–D1048. doi:10.1093/nar/gku1061

---

## 🔄 Learning Loop — Field Reports

<!--
### [Lab/Agent] [Date] — [Topic reviewed]
**Result:** ✅ SUCCESS | ⚠️ PARTIAL | ❌ FAILED
**PMIDs verified:** <N> / <total cited>
**Errors caught in audit pass:** <N>
**Most common error type:** <fabricated PMID | title mismatch | wrong drug approval | other>
**Next reviewer note:** <one line>
-->
