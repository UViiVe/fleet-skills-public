# Skill Template — The 14-Section Checklist

Every agentic-science skill must contain a `SKILL.md` with at least eight of the fourteen sections below, including the seven marked mandatory.

---

## Required Frontmatter

```yaml
---
name: "<kebab-case>"                # must match parent directory name exactly
display-name: "<Human Readable Title>"
description: "<1-2 sentences ending with 'Use when...'>"
category: "<one of: transcriptomics | epigenomics | genetics | proteomics | multi_omics | functional_genomics | clinical | wet_lab | reporting | methodology>"
keywords: "comma, separated, search terms"
version: "1.0.0"                    # semantic versioning — quoted string
status: active                      # active | deprecated | experimental
license: "MIT"                      # MIT | Apache-2.0 | other SPDX identifier
author: "Fleet Skills contributors" # public-safe author or organization
created: "YYYY-MM-DD"
last-updated: "YYYY-MM-DD"
provenance:
  wraps: "<underlying tool or platform>"        # e.g. DESeq2, Scanpy, arXiv API
  upstream_source: "<GitHub URL or original work>"  # upstream repo if wrapping
  citation: "<Author et al., Journal Year, doi:XX>" # paper to cite in methods
changelog:
  - version: "1.0.0"
    date: "YYYY-MM-DD"
    event: authored
    summary: "Initial skill authored for <purpose>."
---
```

Notes:
- `name` must match the parent directory name exactly.
- `version` must be a quoted string: `"1.0.0"` not `1.0.0`.
- `provenance.wraps` — if this skill wraps an existing tool (DESeq2, an API, etc.), name it here. Use `"original work"` for skills with no upstream.
- `provenance.citation` — the paper a user should cite in their methods section. Required for any skill wrapping a published tool.
- `last-updated` is the ISO date of the most recent substantive edit.

---

## Section checklist

The MANDATORY label marks sections that must be present.

### 1. When to Use This Skill — MANDATORY ✅

Includes both ✅ "use when..." cases and ❌ "do not use when..." anti-patterns. Concrete; reads like an instruction to a teammate.

### 2. Quick Start

A 10–30 line worked example showing the typical invocation, end-to-end. Names runtime estimate.

### 3. Installation

Table of packages: `Package | Version | License | Commercial use | Install command`. Use a one-liner install for the common case.

### 4. Inputs — MANDATORY

What the skill expects. Distinguish required vs. optional. Specify formats (CSV, RDS, h5ad, etc.).

### 5. Outputs — MANDATORY

What the skill produces. File names, formats, and which downstream skills consume each artifact.

### 6. Clarification Questions — MANDATORY

A short list (4–7 items) the agent should ask before starting. Q1 must be labeled `⚠️ CRITICAL: ASK THIS FIRST` — this exact phrase enables agent recognition.

### 7. Standard Workflow

Numbered steps. Each step shows the command(s) to run with the exact script invocation pattern. Use code blocks; don't paraphrase.

### 8. When Scripts Fail

The four-level escalation hierarchy:

1. Fix and retry (90%).
2. Modify the script (5%) — with documentation.
3. Use the script as reference (4%) — adapt and cite.
4. Write from scratch (1%) — explain why.

### 9. Decision Points

Table of the major decisions, their defaults, and override conditions. Three columns: Decision, Default, Override when.

### 10. Common Issues — MANDATORY

Table with at least 10 rows: `Issue | Cause | Solution`. Pull from real-world failure modes, not hypotheticals.

### 11. Best Practices

Bullet list of recommendations. Each item is concise and actionable.

### 12. Suggested Next Steps — MANDATORY

Named downstream skills with a one-line rationale per skill. This is what the next agent reads first when figuring out what to do after this skill.

### 13. Related Skills

Upstream, downstream, and companion skills with brief notes.

### 14. References — MANDATORY

Three subsections:
- **Bundled** — list every file under `references/` with a one-line description.
- **Scripts** — list every file under `scripts/` with a one-line description.
- **Primary citations** — academic papers / external docs the skill builds on.

---

## Minimum acceptable skill

Eight of the fourteen sections, including all seven mandatory ones. Drop the optional sections only when they truly don't apply.

A "stub" skill that the rewrite team is filling in should aim for 250–400 lines of `SKILL.md`. Less than 200 lines usually signals missing depth.

---

## Style notes

- Tables outperform prose for decision content.
- Code blocks should be runnable; comments above explain the why.
- Avoid first-person pronouns; the skill is impersonal.
- Avoid emojis except the ✅ / ❌ / ⚠️ tokens used by convention.
- Match the tone of existing skills — opinionated, terse, practical.

---

## Template example (skeleton only)

```markdown
---
name: "<name>"
display-name: "<Title>"
description: "<1-2 sentences>. Use when..."
category: "<category>"
keywords: "<terms>"
version: "1.0.0"
status: active
license: "MIT"
author: "Fleet Skills contributors"
created: "YYYY-MM-DD"
last-updated: "YYYY-MM-DD"
provenance:
  wraps: "<tool or 'original work'>"
  upstream_source: "<URL or N/A>"
  citation: "<Author et al., Journal Year, doi:XX or N/A>"
changelog:
  - version: "1.0.0"
    date: "YYYY-MM-DD"
    event: authored
    summary: "Initial skill."
---

# <Display Title>

<One-paragraph overview.>

## When to Use This Skill
- ✅ ...
- ❌ ...

## Quick Start
<runnable example>

## Installation
<table>

## Inputs
**Required:** ...
**Optional:** ...

## Outputs
<table>

## Clarification Questions
### 1. ⚠️ CRITICAL: ASK THIS FIRST — ...
### 2. ...

## Standard Workflow
### Step 1 — ...
### Step 2 — ...

## When Scripts Fail
1. Fix and retry (90%).
2. ...

## Decision Points
<table>

## Common Issues
<table with ≥ 10 rows>

## Best Practices
- ...

## Suggested Next Steps
1. **<skill-name>** — <rationale>.
2. ...

## Related Skills
- ...

## References
**Bundled:**
- [references/<file>.md](references/<file>.md) — ...

**Scripts:**
- `scripts/<file>.py` — ...

**Primary citations:**
- ...
```
