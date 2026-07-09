"""Scaffold a new Fleet Skills skill: directory tree + SKILL.md template.

Inputs:
    --name STR              Kebab-case skill name (must match directory name).
    --category STR          One of the allowed categories.
    --display-name STR      Human title (default derived from --name).
    --skills-root DIR       Root directory under which the new skill goes
                            (default: 'skills' in the current directory).

Outputs:
    <skills-root>/<name>/SKILL.md
    <skills-root>/<name>/references/
    <skills-root>/<name>/scripts/

Example:
    python scaffold_skill.py --name pcr-primer-design --category wet_lab
"""

from __future__ import annotations

import argparse
import sys
import uuid
from datetime import date
from pathlib import Path


CATEGORIES = {
    "transcriptomics",
    "epigenomics",
    "genetics",
    "proteomics",
    "multi_omics",
    "functional_genomics",
    "clinical",
    "wet_lab",
    "reporting",
    "methodology",
}


SKILL_TEMPLATE = """---
id: "{skill_id}"
name: "{name}"
display-name: "{display_name}"
short-description: "<1-2 sentences ending with 'Use when...'>"
category: "{category}"
keywords: "<comma, separated, terms>"
version: "1.0.0"
status: active
license: "MIT"
author: "Fleet Skills contributors"
created: "{today}"
last-updated: "{today}"
provenance:
  wraps: "original work"
  upstream_source: "original work"
  citation: "N/A"
changelog:
  - version: "1.0.0"
    date: "{today}"
    event: authored
    summary: "Initial skill."
---

# {display_name}

<One-paragraph overview describing what the skill does and why someone would use it.>

---

## When to Use This Skill

✅ **Use when:**
- <case 1>
- <case 2>

❌ **Do not use when:**
- <case 1>
- <case 2>

**Core principle:** <one sentence>.

---

## Quick Start

```bash
# runnable example
```

Runtime: <estimate>.

---

## Installation

| Package | Version | License | Commercial use | Install |
|---|---|---|---|---|
| <pkg> | <ver> | <license> | yes | `pip install <pkg>` |

---

## Inputs

**Required:**
- <input>

**Optional:**
- <input>

---

## Outputs

| File | Description |
|---|---|
| <name> | <description> |

---

## Clarification Questions

### 1. ⚠️ CRITICAL: ASK THIS FIRST — <the most important question>

<Explanation of why this is asked first.>

### 2. <next question>

<context>

---

## Standard Workflow

### Step 1 — <step name>

```bash
# command
```

### Step 2 — <step name>

```bash
# command
```

---

## When Scripts Fail

1. **Fix and retry (90%).** <how>
2. **Modify the script (5%).** <when>
3. **Use as reference (4%).** <when>
4. **Write from scratch (1%).** <when>

---

## Decision Points

| Decision | Default | Override when |
|---|---|---|
| <decision> | <default> | <override condition> |

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| <issue 1> | <cause> | <solution> |
| <issue 2> | <cause> | <solution> |
| <issue 3> | <cause> | <solution> |
| <issue 4> | <cause> | <solution> |
| <issue 5> | <cause> | <solution> |
| <issue 6> | <cause> | <solution> |
| <issue 7> | <cause> | <solution> |
| <issue 8> | <cause> | <solution> |
| <issue 9> | <cause> | <solution> |
| <issue 10> | <cause> | <solution> |

---

## Best Practices

- <practice>
- <practice>

---

## Suggested Next Steps

1. **<skill-name>** — <one-line rationale>.

---

## Related Skills

- <other skill>

---

## References

**Bundled:**
- [references/<file>.md](references/<file>.md) — <one-line description>.

**Scripts:**
- `scripts/<file>.py` — <one-line description>.

**Primary citations:**
- <author year, journal>.
"""


def scaffold(
    name: str,
    category: str,
    display_name: str | None = None,
    skills_root: str | Path = "skills",
) -> Path:
    """Create the skill directory tree and write a SKILL.md template."""
    if category not in CATEGORIES:
        raise ValueError(
            f"category must be one of {sorted(CATEGORIES)}, got {category!r}"
        )

    if display_name is None:
        display_name = " ".join(part.capitalize() for part in name.split("-"))

    skill_id = "skill_" + uuid.uuid4().hex
    root = Path(skills_root) / name
    if root.exists():
        sys.exit(f"Directory already exists: {root}")

    (root / "references").mkdir(parents=True)
    (root / "scripts").mkdir(parents=True)

    skill_md = SKILL_TEMPLATE.format(
        skill_id=skill_id,
        name=name,
        display_name=display_name,
        category=category,
        today=date.today().isoformat(),
    )
    (root / "SKILL.md").write_text(skill_md)

    print(f"✓ Scaffolded skill at: {root}")
    print(f"   id: {skill_id}")
    print(f"   category: {category}")
    print(f"   Next: edit SKILL.md and fill in references/ + scripts/")
    return root


def main() -> None:
    parser = argparse.ArgumentParser(description="Scaffold a new claude-science skill.")
    parser.add_argument("--name", required=True, help="Kebab-case directory name.")
    parser.add_argument(
        "--category",
        required=True,
        choices=sorted(CATEGORIES),
        help="Skill category.",
    )
    parser.add_argument("--display-name", default=None,
                        help="Human-readable title (default derived from --name).")
    parser.add_argument("--skills-root", default="skills",
                        help="Root directory under which the new skill goes.")
    args = parser.parse_args()

    scaffold(
        name=args.name,
        category=args.category,
        display_name=args.display_name,
        skills_root=args.skills_root,
    )


if __name__ == "__main__":
    main()
