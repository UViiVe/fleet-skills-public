---
name: peer-review-gate
id: "skill_10a9e9c4904a4f8e9a3413be95ee4910"
display-name: "Peer Review Gate"
short-description: "Require independent review before scaled, expensive, irreversible, or high-risk scientific work. Use when a second-review checkpoint should happen before execution."
keywords: [peer review, validation, gate, science, methodology]
last-updated: "2026-07-09"
version: "1.0.2"
author: "Fleet Skills contributors"
created: 2026-06-14
description: "Require independent review before scaled, expensive, irreversible, or high-risk scientific work. Use when a second-review checkpoint should happen before execution."
category: methodology
changelog:
  - version: "1.0.2"
    date: 2026-07-09
    event: public-cleanup
    notes: "Replaced internal incident details, private channel names, and private logging examples with generic public review guidance."
  - version: "1.0.1"
    date: 2026-06-14
    event: clarification
    notes: "Clarified scope for broad scientific review."
  - version: "1.0.0"
    date: 2026-06-14
    event: authored
    notes: "Initial release."
status: active
promoted: true
license: MIT
---

# Peer Review Gate

Use this skill to require a second reviewer before work that is expensive,
irreversible, externally visible, or scientifically high risk. The goal is not
bureaucracy; it is to catch mistakes before they scale.

---

## When This Applies

Always require an independent review before:

- A compute or design run expected to consume substantial time or cost.
- Any irreversible file write, API submission, message send, or external action at volume.
- A new structure, reference genome, annotation set, or model is used as the basis for downstream work.
- A result will support a manuscript, grant, regulatory, clinical, or commercial claim.
- A workflow depends on manual mapping between numbering systems, coordinate systems, identifiers, or units.

Use lighter review for:

- Small exploratory runs.
- Local-only analysis drafts.
- Routine reruns with unchanged inputs and reviewed configs.

---

## Reviewer Rules

- The reviewer must be independent from the producer.
- Self-review does not count.
- The reviewer should be qualified for the failure mode: domain reviewer for biology/science, technical reviewer for code/config, and project owner for scope or claims.
- A blocked review is not a failure. It is the gate doing its job.

---

## Minimum Review Checklist

### Inputs

- [ ] Input files exist and match the intended dataset, organism, build, sample set, or structure.
- [ ] Identifiers are mapped correctly across systems.
- [ ] Units and coordinate systems are explicit.
- [ ] Any manual curation step is reproducible from notes or scripts.

### Configuration

- [ ] Config paths exist on the machine that will run the job.
- [ ] Output paths have enough disk space and will not overwrite canonical data.
- [ ] Parameters match the stated objective.
- [ ] The command can be rerun from a clean shell.

### Scientific Claim

- [ ] The result supports the claim being made.
- [ ] Alternative explanations or confounders are named.
- [ ] Required controls or comparison groups are present.
- [ ] Any regulatory, clinical, IP, or commercial language is routed to the right human reviewer before use.

---

## Sign-Off Format

Post or save this before the run starts:

```text
PEER REVIEW - [Reviewer] reviewed [artifact/config name] for [Producer]
Date: YYYY-MM-DD
Artifact: [filename, commit, or description]
Scope: [T1 scaled run | T2 scientific claim | T3 project assessment]
Verified:
  - Input identity:
  - Identifier/coordinate mapping:
  - Config and output paths:
  - Claim scope:
Result: CLEARED / BLOCKED
If blocked: [specific issue and required fix]
```

If cleared, the producer may proceed within the reviewed scope. If blocked, the
producer fixes the issue and resubmits the corrected artifact for review.

---

## Common Failure Modes

| Failure | Detection | Fix |
|---|---|---|
| Identifier mismatch | IDs do not match the source table or annotation build | Re-map from the source and document the mapping step |
| Coordinate mismatch | Results shift by build, chain, or indexing convention | State the coordinate system and verify against source records |
| Wrong input file | File name looks right but content is stale or from another run | Check hashes, row counts, sample IDs, and timestamps |
| Overbroad claim | Result supports association but text claims mechanism | Downgrade language or add experimental validation |
| Missing control | No baseline, negative control, or comparator | Add control or mark result exploratory |
| Unsafe output path | Run overwrites canonical or shared data | Write to a versioned output directory |
| Scale mismatch | Small pilot config reused for full production | Recalculate resources and review production parameters |

---

## Records

Keep review records somewhere durable for the project: issue tracker, pull
request, lab notebook, run log, or project database. The public skill does not
prescribe a private schema.

---

## Related Skills

- `data-analysis-best-practices` - data validation and reproducibility checks.
- `verified-science-review` - citation-backed review for scientific claims.
- `experimental-design-statistics` - study design and power checks.
