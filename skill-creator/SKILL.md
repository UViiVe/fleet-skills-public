---
name: skill-creator
id: "skill_2df7aeaf7e8f4ef292e5219c1ca2de25"
display-name: "Skill Creator"
short-description: "Draft, validate, and package a reusable agentic-science skill with public-safe metadata, references, scripts, and validation checks. Use when designing a new workflow."
category: methodology
keywords: [skill, create, template, workflow, methodology]
version: "1.4.1"
author: "Fleet Skills contributors"
created: 2026-05-20
last-updated: 2026-07-07
status: active
license: MIT
changelog:
  - version: "1.4.1"
    date: 2026-07-07
    event: template-update
    summary: >-
      Template updated to require provenance fields in all new skills: wraps,
      upstream_source, and citation. Added status, license, and public-safe
      author fields.
  - version: "1.4.0"
    date: 2026-07-07
    event: public-schema-update
    summary: >-
      Frontmatter template updated for public release: added status, license,
      public author, and provenance fields.
  - version: "1.3.0"
    date: 2026-07-06
    event: architecture-update
    summary: >-
      Optional shared skill library location documented for teams that maintain
      an internal deployment target.
  - version: 1.2.0
    date: 2026-06-05
    event: architecture-update
    summary: >-
      Changelog update added as a mandatory step in the skill creation workflow.
  - version: 1.1.0
    date: 2026-06-05
    event: architecture-update
    summary: >-
      Added Learning Loop system as a mandatory skill component. Every skill
      now ships with a versioned changelog in frontmatter and a Learning Loop
      section for field reports.
  - version: 1.0.0
    date: 2026-05-20
    event: authored
    summary: Initial import.
description: "Draft, validate, and package a new agentic-science skill with consistent frontmatter, mandatory sections, reference docs, scripts, and validation checks. Use when designing a reusable workflow that other skills will depend on."
source: https://github.com/UViiVe/fleet-skills-public
canonical-skills-library:
  host: <your-skills-server>
  path: <your-skills-library>/
  install-cmd: "scp -r <skill-name>/ <your-skills-library>/<skill-name>/"
  verify-cmd: "ssh <your-skills-server> 'ls <your-skills-library>/<skill-name>/'"
  github-fallback: https://github.com/UViiVe/fleet-skills-public
promoted: true
---

# Skill Creator

> **Shared skills library:** `<your-skills-library>/`
> Teams that maintain a shared skill library can deploy completed skills there.
> Public source: https://github.com/UViiVe/fleet-skills-public

A meta-skill for designing and packaging new agentic-science skills. The output is a folder under `<skills-root>/<skill-name>/` containing a `SKILL.md`, a `references/` directory, and a `scripts/` directory - built to the project quality bar and passing validation checks.

Use this when the work you're doing has the shape of a reusable workflow (same inputs, same steps, same outputs across cases) and a future agent would benefit from a saved playbook.

---

## When to Use This Skill

✅ **Use when:**
- A workflow you've now done twice looks like it'll happen again.
- A teammate has asked for a documented version of an ad-hoc procedure.
- You want to package institutional knowledge — best practices, gotchas, references — into a single discoverable artifact.
- An upstream skill needs a downstream skill that doesn't exist yet.
- You're standardizing how a category of analyses is done across the team.

❌ **Do not use when:**
- The work is one-off and unlikely to recur.
- A skill already exists that mostly covers the need (extend the existing one).
- The user is exploring; a `discuss-phase` conversation is more appropriate than a new skill.
- The workflow is in flux — write notes first, formalize later.

**Core principle:** the skill is the contract. If two different agents pick it up six months apart, they should produce comparable artifacts.

---

## Quick Start

Scaffold a new skill in three commands:

```bash
# 1. Create the directory structure
python scripts/scaffold_skill.py --name my-new-skill --category methodology

# 2. Edit SKILL.md, references/, scripts/ — fill in the template

# 3. Validate before declaring done
python scripts/validate_skill.py --skill-dir skills/my-new-skill --public-release
```

The scaffold writes a `SKILL.md` template with all mandatory sections, an empty `references/` and `scripts/`, and a checklist of what to fill in.

---

## Installation

| Package | Purpose |
|---|---|
| python ≥ 3.11 | Required for f-string and walrus features |
| pyyaml | Parsing SKILL.md frontmatter |
| pathlib (stdlib) | Path handling |

```bash
pip install pyyaml
```

No external dependencies beyond the validator's needs.

---

## Inputs

**Required for scaffolding:**
- `name`: kebab-case directory name (must match the `name:` field in frontmatter).
- `category`: one of `transcriptomics | epigenomics | genetics | proteomics | multi_omics | functional_genomics | clinical | wet_lab | reporting | methodology`.

**Required for the skill itself (to be authored):**
- A clear `short-description` ending with "Use when..."
- The full When-to-Use list.
- Inputs and Outputs with formats.
- A `CRITICAL: ASK THIS FIRST` clarification question.
- A Common Issues table with ≥ 10 rows.
- Named Suggested Next Steps.
- A References section with primary citations.

---

## Outputs

A fully-formed skill directory:

```
skills/<skill-name>/
├── SKILL.md                  # The skill itself (frontmatter + 8+ sections)
├── references/               # Long-form reference docs
│   └── <topic>.md
└── scripts/                  # Executable code
    ├── <thing>.py
    └── <thing>.R
```

The `validate_skill.py` script enforces:

1. No banned tokens or public-release leaks (see `scripts/validate_skill.py` for the authoritative checks and `references/banned_terms.md` for category guidance).
2. Frontmatter parses and has all required fields.
3. Every `.py` script compiles cleanly.
4. Every `.sh` script parses cleanly.

---

## Clarification Questions

### 1. CRITICAL: ASK THIS FIRST - Is this actually a skill?

Before writing anything, verify:

- Is the workflow procedural (same steps every time)?
- Is it likely to be invoked again?
- Are inputs, outputs, and failure modes stable enough to document?
- Is there a clear "done" state?

If two answers are "no", this is probably a *project task* (use a plan) or a *one-off* (use the conversation). Don't bake noise into the skill library.

### 2. Category and scope

| Category | Use for |
|---|---|
| `methodology` | Cross-cutting analysis discipline (this skill is one) |
| `transcriptomics` | Bulk or single-cell RNA workflows |
| `epigenomics` | ATAC, ChIP, methylation |
| `genetics` | GWAS, MR, variant annotation |
| `proteomics` | Mass-spec proteomics |
| `multi_omics` | Integrative workflows |
| `functional_genomics` | CRISPR screens, gene essentiality |
| `clinical` | Survival, trials, clinical-data analyses |
| `wet_lab` | Bench-facing protocols |
| `reporting` | PDF / PPTX / DOCX generation |

Pick one. Cross-cutting skills usually belong in `methodology`.

### 3. Naming

- Kebab-case (`my-new-skill`).
- Specific enough that the name encodes purpose (`crispr-screen-mageck` not `screens`).
- No vendor names unless central (allowed: `pooled-crispr-screens`, `gwas-to-function-twas`).
- No banned tokens (see Banned Terms below).

### 4. What references and scripts will it need?

Sketch the references/ files (long-form domain explanations) and the scripts/ files (executable code) before writing. The skeleton helps focus the SKILL.md.

---

## Standard Workflow

### Step 1 — Decide whether the skill should exist

Apply the "is this a skill" test. If unsure, draft a half-page outline and gut-check.

### Step 2 — Scaffold

```bash
python scripts/scaffold_skill.py --name <name> --category <category>
```

Generates the directory tree and a `SKILL.md` template with all mandatory sections marked.

### Step 3 — Fill the frontmatter

```yaml
---
name: "<kebab-case>"
display-name: "<Human Readable Title>"
description: "<1-2 sentences ending with 'Use when...'>"
category: "<category>"
keywords: "comma, separated, search terms"
version: "1.0.0"
status: active
license: "MIT"                    # MIT for public release; change if proprietary
author: "Fleet Skills contributors"
created: "YYYY-MM-DD"
last-updated: "YYYY-MM-DD"
provenance:
  wraps: "<underlying tool or platform, if applicable>"
  upstream_source: "<GitHub URL or 'original work'>"  # upstream repo if wrapping existing tool
  citation: "<Author et al., Journal Year, doi:XX>"    # paper to cite in methods sections
changelog:
  - version: "1.0.0"
    date: "YYYY-MM-DD"
    event: authored
    summary: "Initial skill authored for <purpose>."
---
```

Generate the `id` with:

```bash
python3 -c "import uuid; print('skill_' + uuid.uuid4().hex)"
```

**Version format:** `MAJOR.MINOR.PATCH`
- PATCH: typo fixes, clarifications, new field report in Learning Loop
- MINOR: new blockers documented, new steps added, behavioral changes
- MAJOR: skill redesigned or deprecated path removed

**Changelog rule:** Every version bump requires a `changelog:` entry. This keeps the skill history durable and reviewable.

### Step 4 — Draft SKILL.md

Use the 14-section checklist in [references/skill_template.md](references/skill_template.md). Mandatory sections (must be present):

1. When to Use This Skill (with ✅ and ❌).
2. Inputs (required + optional, with formats).
3. Outputs (file names, formats, downstream consumers).
4. Clarification Questions (Q1 labeled `⚠️ CRITICAL: ASK THIS FIRST`).
5. Common Issues (≥ 10 rows in a table).
6. Suggested Next Steps (named downstream skills with rationale).
7. References (bundled docs + primary citations).
8. **Learning Loop** — see below. This is now mandatory for any skill that involves installation, deployment, or execution steps.

Recommended additions: Quick Start, Installation, Standard Workflow, Decision Points, Best Practices, When Scripts Fail.

### Step 4b — Write the Learning Loop section

Every skill involving installation, deployment, or execution steps **must** include a Learning Loop section. This is the mechanism by which fleet agents teach each other through real experience — not just document what *should* work, but record what *did* work and what broke.

**Scaffold for the Learning Loop section:**

```markdown
## Learning Loop - Field Reports

This section is the skill's living memory. Every agent that runs this skill
in production can append a field report. The frontmatter `changelog:` tracks
versions formally; this section provides the human-readable narrative.

**Protocol for the executing agent:**
1. Complete the task and verify outputs
2. Write a field report below (template follows)
3. Bump `version:` in frontmatter (patch for minor, minor for new blockers)
4. Add a `changelog:` entry with `date`, `event`, and `summary`
5. If your team uses a shared skills library, sync the completed skill there.

**Field report template:**
<!-- Copy this block, fill in, remove the comment tags -->
<!--
### [AGENT] [DATE] — [ENVIRONMENT DESCRIPTION]
**Result:** ✅ SUCCESS | ⚠️ PARTIAL | ❌ FAILED
**Time:** <how long>
**New blockers:** <anything not in prior reports>
**Resolved by:** <fix>
**What to skip:** <steps that proved unnecessary>
**Next agent note:** <one-line advice>
-->

*No field reports yet - first agent to run this skill: please add yours.*
```

**Why this matters:** The gap between "documented" and "works" is where institutional knowledge lives. A skill that only records the happy path fails the second agent who hits a different environment. Field reports close that gap. The versioning system makes it durable - each agent's experience is preserved, not overwritten.

**Scope:** Required for skills involving any of: installation, environment setup, external API calls, SSH/network operations, model downloads, service deployment. Optional but encouraged for pure analysis skills.

### Step 5 — Write references/

Each reference file should be a self-contained explanation of a concept the skill assumes. Examples from existing skills:

- `references/missing_data.md` — concept-level reference on MCAR/MAR/MNAR.
- `references/cv_tissue_database.csv` — data file consumed by scripts.
- `references/multiple_testing_guide.md` — decision-tree-style guide.

Avoid: bullet-list dumps with no narrative; restating SKILL.md; placeholder text.

### Step 6 — Write scripts/

Every script must satisfy all seven rules in [references/script_style_guide.md](references/script_style_guide.md):

1. Module docstring stating purpose, inputs, outputs, example invocation.
2. Type hints on every public function.
3. `argparse` with `--help` for CLI scripts; `def main()` entry point.
4. Success line: `print("✓ <thing> completed: <n> rows / <path>")`.
5. Graceful ImportError: `try: import X; except ImportError: sys.exit("Install with: pip install X")`.
6. `if __name__ == "__main__": main()` guard.
7. No hardcoded absolute paths.

### Step 7 — Validate

```bash
python scripts/validate_skill.py --skill-dir skills/<name>
```

Runs the four required checks. Zero failures is the bar.

### Step 8 - Update the changelog

A skill that ships without a changelog entry is harder to review and maintain.

```bash
# Append to your project changelog
echo "[$(date +%Y-%m-%d)] [skill] [<name> v1.0.0] <what it does and why it was created>" >> CHANGELOG.md
```

**Entry format:**
```
[YYYY-MM-DD] [skill] [skill-name vX.Y.Z] - <one paragraph: what the skill does,
why it was created, what problem it solves, any key design decisions>
```

For **revisions** (not just new skills): same format, but lead with what changed and why.

### Step 9 — Cross-link

Update related skills' `Suggested Next Steps` to point at the new skill if appropriate. Discoverability is a feature, not a bonus.

---

## When Scripts Fail

1. **Fix and retry (90%).** Missing dependency or typo in frontmatter. Re-run after fix.
2. **Modify the validator (5%).** If a check is overly strict for a legitimate case, update the validator with a justification.
3. **Use as reference (4%).** For unusual skill structures (multi-language, multi-environment), read the validator and adapt the check manually.
4. **Write from scratch (1%).** Reserved for the rare case where the existing scaffold genuinely doesn't fit.

---

## Decision Points

| Decision | Default | Override when |
|---|---|---|
| Category | most specific match | Cross-cutting → `methodology` |
| References count | 2–6 files | Stub skills (< 100 lines source) can start with 2; full skills need 4+ |
| Scripts count | 3–7 files | Pure-reference skills (no code) can have 0; complex pipelines may have 10+ |
| Common Issues rows | ≥ 10 | Smaller skills can have 8–10; pad with realistic cases |
| Q1 phrasing | "CRITICAL: ASK THIS FIRST" | Don't override - this token enables recognition |
| Frontmatter format | YAML | Don't override |

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| Banned token check fails on a dashed license abbreviation | Word boundary at the dash matches a banned token | Rephrase as "Creative Commons Attribution 4.0" |
| Frontmatter doesn't parse | Stray colon or quote in description | Wrap string values in double quotes |
| Duplicate `id:` across skills | Re-used template | Always regenerate UUID for new skills |
| `name:` doesn't match directory | Renamed directory but not frontmatter | Keep them synchronized — script enforces |
| Missing `⚠️ CRITICAL: ASK THIS FIRST` | Skipped the rule | Add to Q1; agent dispatching relies on this token |
| Common Issues table < 10 rows | Author ran out of cases | Brainstorm with the workflow's failure modes; ask "what could go wrong here?" |
| Scripts hardcode paths | Pasted from a Jupyter notebook | Refactor to take input/output via `argparse` |
| Reference files don't exist but are linked | Broken plan during drafting | Validator catches via link existence checks; create or remove the link |
| Quick Start uses imports that aren't installed | Forgot the Installation section | Make the install list match every import in Quick Start |
| Categories don't match the global enum | Made up a new category | Use the enum in clarification Q2; if new category needed, propose first |
| Stub source vs expanded skill | Source has 50 lines, new skill needs 250 | Expand the stub fully; pull from related skills + domain knowledge |
| Cross-references to non-existent skills | Linked to a skill not yet written | Use `[[name]]` notation in references or fix link |
| `last-updated:` not refreshed | Copied template without update | Always set to today's ISO date |
| Changelog not updated | Treated as optional | Add a concise entry in Step 8 |
| `changelog:` missing from frontmatter | Used old template | Every skill needs `changelog:` from v1.0.0 onward |
| Learning Loop section missing | Skipped as optional | Required for any skill with install/deploy/execution steps |
| Field report not written after first run | Agent forgot the protocol | The protocol is in the Learning Loop scaffold — read it before closing the task |
| Banned-token substring inside a longer scientific word | False positive on substring | The check uses word boundaries (`\b`), so substrings inside longer scientific words are fine |
| Scripts pass but Python version too old | Used 3.11+ syntax in a 3.10 project | Pin Python ≥ 3.11 in requirements or downgrade syntax |

---

## Best Practices

- **Start from a stub, fail fast.** Drop a placeholder SKILL.md with the 8 mandatory sections empty, then fill in. Iterating on structure is cheaper than iterating on prose.
- **Write the Inputs and Outputs first.** They constrain everything else. If you can't write them, the skill isn't ready.
- **Question 1 is for the agent.** It should be the single most-important thing to know before doing anything. Make it specific.
- **Tables beat bullet lists.** Most decisions benefit from a 3-column "case / default / override" table.
- **Reference docs are for narrative depth.** SKILL.md is the contract; references are the supporting evidence.
- **Scripts should run unmodified.** If a junior teammate clones the repo and runs `python scripts/<thing>.py --help`, it should work.
- **The Learning Loop is not optional for deploy skills.** If your skill has installation or execution steps and you've just run it, you owe the fleet a field report. Two paragraphs. What broke, what fixed it, what to skip. Future agents will thank you.
- **Version numbers communicate intent.** A PATCH says "same behavior, cleaner docs." A MINOR says "new things work now or new things broke." Be precise.
- **Read the validator.** It encodes project conventions you might not have absorbed yet.

---

## Suggested Next Steps

After creating a new skill:

1. **Validate** with `validate_skill.py`. Fix any failures.
2. **Cross-link** with related skills' Suggested Next Steps.
3. **Publish or sync** if your project maintains a shared skill library:
   ```bash
   ssh <your-skills-server> "mkdir -p <your-skills-library>/<skill-name>"
   scp -r <skill-name>/ <your-skills-library>/<skill-name>/
   ```
   Confirm with: `ssh <your-skills-server> "ls <your-skills-library>/<skill-name>/"`
   Then update your project changelog.
4. **Smoke-test** by following your own Quick Start as if you were the agent — does it produce something usable?
5. **Index** if your project maintains a skill catalog or knowledge graph (e.g., `gsd-graphify` if available).
6. **Iterate** based on first real use. The second invocation surfaces gaps the first hides.

---

## Related Skills

- All skills in `skills/` are examples — read 2–3 of them before starting.
- **Reference:** [data-analysis-best-practices](../data-analysis-best-practices/SKILL.md) is a good example of a methodology-category skill.
- **Reference:** [experimental-design-statistics](../experimental-design-statistics/SKILL.md) is a heavier reference for what a multi-script skill looks like.

---

## References

**Bundled:**
- [references/skill_template.md](references/skill_template.md) — the 14-section checklist with explanations.
- [references/script_style_guide.md](references/script_style_guide.md) — the seven script rules with rationale.
- [references/validation_checks.md](references/validation_checks.md) — what `validate_skill.py` checks and why each rule exists.
- [references/banned_terms.md](references/banned_terms.md) — the project-wide do-not-use list and how to rephrase common offenders.

**Scripts:**
- `scripts/scaffold_skill.py` — create the directory tree and `SKILL.md` template.
- `scripts/validate_skill.py` — run all four validation checks.
- `scripts/check_frontmatter.py` — verify YAML frontmatter shape only.

**External references:**
- The existing skill catalog under `skills/` — read these for live examples.
