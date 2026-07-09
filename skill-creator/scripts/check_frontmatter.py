"""Quick frontmatter-only validation for a SKILL.md.

A focused check that complements validate_skill.py — useful for iterating on
frontmatter without paying the cost of full validation.

Inputs:
    --skill-md PATH     Path to a SKILL.md file.

Outputs:
    Prints the parsed frontmatter and any missing required keys.

Example:
    python check_frontmatter.py --skill-md skills/my-skill/SKILL.md
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("Install pyyaml with: pip install pyyaml")


REQUIRED_KEYS = (
    "id", "name", "display-name", "short-description",
    "category", "keywords", "version", "last-updated",
)


def parse_frontmatter(skill_md_path: str | Path) -> dict[str, object]:
    """Parse YAML frontmatter from a SKILL.md file."""
    text = Path(skill_md_path).read_text()
    parts = text.split("---")
    if len(parts) < 3:
        raise ValueError("SKILL.md does not contain frontmatter delimited by '---'")
    data = yaml.safe_load(parts[1])
    if not isinstance(data, dict):
        raise ValueError("Frontmatter is not a YAML mapping")
    return data


def check(skill_md_path: str | Path) -> tuple[bool, list[str]]:
    """Return (ok, list of error messages)."""
    p = Path(skill_md_path)
    if not p.exists():
        return False, [f"file not found: {p}"]

    try:
        data = parse_frontmatter(p)
    except (ValueError, yaml.YAMLError) as e:
        return False, [str(e)]

    errors: list[str] = []
    missing = [k for k in REQUIRED_KEYS if k not in data]
    if missing:
        errors.append(f"missing keys: {missing}")

    expected_name = p.parent.name
    if "name" in data and data["name"] != expected_name:
        errors.append(
            f"name '{data['name']}' does not match directory '{expected_name}'"
        )

    if "id" in data and not isinstance(data["id"], str):
        errors.append("id must be a string")
    if "id" in data and isinstance(data["id"], str) and not data["id"].startswith("skill_"):
        errors.append("id should start with 'skill_'")

    return (len(errors) == 0), errors


def main() -> None:
    parser = argparse.ArgumentParser(description="Quick frontmatter validation.")
    parser.add_argument("--skill-md", required=True,
                        help="Path to a SKILL.md file.")
    args = parser.parse_args()

    ok, errors = check(args.skill_md)
    if ok:
        data = parse_frontmatter(args.skill_md)
        print("✓ Frontmatter validated: all required keys present")
        for k in REQUIRED_KEYS:
            print(f"   {k}: {data.get(k, '<missing>')}")
    else:
        print("✗ Frontmatter validation failed:")
        for e in errors:
            print(f"   - {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
