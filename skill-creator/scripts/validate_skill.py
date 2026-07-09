"""Run the four validation checks on a skill directory.

Inputs:
    --skill-dir DIR     Path to the skill directory (e.g., skills/my-skill).

Outputs:
    Prints pass/fail per check; exits non-zero on any failure.

Checks performed:
    1. No banned tokens.
    2. SKILL.md frontmatter parses and has all required keys.
    3. Every .py script compiles cleanly.
    4. Every .sh script parses cleanly.
    5. Optional public-release scan for internal metadata and likely secrets.

Example:
    python validate_skill.py --skill-dir skills/data-analysis-best-practices
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    sys.exit("Install pyyaml with: pip install pyyaml")


# Build the banned-token regex from string fragments so this source file does
# not itself contain any of the rejected word-boundary tokens.
_BANNED_TOKENS = (
    "biomn" + "i",
    "Biomn" + "i",
    "phy" + "lo",
    "Phy" + "lo",
    "blatan" + "t",
    "Blatan" + "t",
    "B" + "Y",
    "Blatan" + "t-Why",
)
BANNED_PATTERN = re.compile(r"\b(" + "|".join(_BANNED_TOKENS) + r")\b")

REQUIRED_FRONTMATTER_KEYS = (
    "id", "name", "display-name", "short-description",
    "category", "keywords", "version", "last-updated",
)

SCANNABLE_SUFFIXES = (".md", ".py", ".sh", ".yaml", ".yml", ".json", ".R", ".csv")

PUBLIC_RELEASE_TOKENS = (
    "Ches" + "ter",
    "Ze" + "vo",
    "At" + "las",
    "Gal" + "en",
    "Open" + "Claw",
    "open" + "claw",
    "St" + "eve" + " " + "Ek" + "ker",
    "001" + "TMF",
    "Vienna" + " public subset",
)
PUBLIC_RELEASE_PATTERNS = (
    re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}"),
    re.compile(r"\b(?:AIza|AKIA|ASIA)[A-Za-z0-9_-]{16,}\b"),
    re.compile(r"\bgh[pousr]_[A-Za-z0-9_]{20,}\b"),
    re.compile(r"\bhf_[A-Za-z0-9]{20,}\b"),
    re.compile(r"\bplane_api_[A-Za-z0-9_]{20,}\b"),
    re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
    re.compile(r"\b1[A-Za-z0-9_-]{25,}\b"),
)

PUBLIC_RELEASE_ALLOWED_EMAILS = {
    "you@example.com",
}
PUBLIC_RELEASE_ALLOWED_EMAIL_SUFFIXES = (
    "@meta.data",
    "@data.signaling",
)
PUBLIC_RELEASE_ALLOWED_ATLAS_CONTEXTS = (
    "ChIP-Atlas",
    "Human Cell Atlas",
    "Human Primary Cell Atlas",
    "humancellatlas.org",
)


def _frontmatter_text(skill_md: Path) -> str:
    text = skill_md.read_text()
    parts = text.split("---")
    if len(parts) < 3:
        return ""
    return parts[1]


def check_banned_tokens(skill_dir: Path) -> list[str]:
    """Return a list of 'file:line:text' hits for banned tokens."""
    hits: list[str] = []
    for path in skill_dir.rglob("*"):
        if not path.is_file() or path.suffix not in SCANNABLE_SUFFIXES:
            continue
        try:
            text = path.read_text(errors="ignore")
        except OSError:
            continue
        for lineno, line in enumerate(text.splitlines(), 1):
            if BANNED_PATTERN.search(line):
                hits.append(f"{path}:{lineno}:{line.strip()}")
    return hits


def check_duplicate_frontmatter_keys(skill_dir: Path) -> list[str]:
    """Return duplicate top-level YAML keys from SKILL.md frontmatter."""
    skill_md = skill_dir / "SKILL.md"
    if not skill_md.exists():
        return []
    seen: set[str] = set()
    duplicates: list[str] = []
    for line in _frontmatter_text(skill_md).splitlines():
        match = re.match(r"^([A-Za-z0-9_-]+):", line)
        if not match:
            continue
        key = match.group(1)
        if key in seen and key not in duplicates:
            duplicates.append(key)
        seen.add(key)
    return duplicates


def check_public_release(skill_dir: Path) -> list[str]:
    """Return public-release hygiene hits for internal tokens and likely secrets."""
    hits: list[str] = []
    token_pattern = re.compile(
        r"\b(" + "|".join(re.escape(t) for t in PUBLIC_RELEASE_TOKENS) + r")\b"
    )

    def allowed_public_release_match(match_text: str, line: str) -> bool:
        if match_text in PUBLIC_RELEASE_ALLOWED_EMAILS:
            return True
        if any(match_text.endswith(suffix) for suffix in PUBLIC_RELEASE_ALLOWED_EMAIL_SUFFIXES):
            return True
        if match_text == "At" + "las" and any(
            context in line for context in PUBLIC_RELEASE_ALLOWED_ATLAS_CONTEXTS
        ):
            return True
        return False

    for path in skill_dir.rglob("*"):
        if not path.is_file() or path.suffix not in SCANNABLE_SUFFIXES:
            continue
        try:
            lines = path.read_text(errors="ignore").splitlines()
        except OSError:
            continue
        for lineno, line in enumerate(lines, 1):
            token_hit = any(
                not allowed_public_release_match(match.group(0), line)
                for match in token_pattern.finditer(line)
            )
            pattern_hit = any(
                not allowed_public_release_match(match.group(0), line)
                for pattern in PUBLIC_RELEASE_PATTERNS
                for match in pattern.finditer(line)
            )
            if token_hit or pattern_hit:
                hits.append(f"{path}:{lineno}:{line.strip()}")
    return hits


def check_frontmatter(skill_dir: Path) -> list[str]:
    """Return a list of validation errors; empty list means OK."""
    skill_md = skill_dir / "SKILL.md"
    if not skill_md.exists():
        return [f"missing SKILL.md at {skill_md}"]

    text = skill_md.read_text()
    parts = text.split("---")
    if len(parts) < 3:
        return ["SKILL.md does not contain frontmatter delimited by '---'"]

    try:
        data = yaml.safe_load(parts[1])
    except yaml.YAMLError as e:
        return [f"frontmatter YAML parse error: {e}"]
    if not isinstance(data, dict):
        return ["frontmatter is not a mapping"]

    missing = [k for k in REQUIRED_FRONTMATTER_KEYS if k not in data]
    errors: list[str] = []
    if missing:
        errors.append(f"missing keys: {missing}")
    duplicates = check_duplicate_frontmatter_keys(skill_dir)
    if duplicates:
        errors.append(f"duplicate top-level frontmatter keys: {duplicates}")
    if "name" in data and data["name"] != skill_dir.name:
        errors.append(
            f"name '{data['name']}' does not match directory '{skill_dir.name}'"
        )
    return errors


def check_python_scripts(skill_dir: Path) -> list[str]:
    """Compile each .py with py_compile; return list of failing files."""
    failures: list[str] = []
    for p in (skill_dir / "scripts").glob("*.py"):
        result = subprocess.run(
            [sys.executable, "-m", "py_compile", str(p)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            failures.append(f"{p}: {result.stderr.strip()}")
    return failures


def check_shell_scripts(skill_dir: Path) -> list[str]:
    """Parse each .sh with bash -n; return list of failing files."""
    failures: list[str] = []
    for p in (skill_dir / "scripts").glob("*.sh"):
        result = subprocess.run(
            ["bash", "-n", str(p)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            failures.append(f"{p}: {result.stderr.strip()}")
    return failures


def validate(skill_dir: str | Path, public_release: bool = False) -> tuple[int, list[str]]:
    """Run all four checks; return (failure_count, messages)."""
    skill_dir = Path(skill_dir)
    if not skill_dir.is_dir():
        return 1, [f"not a directory: {skill_dir}"]

    failures = 0
    messages: list[str] = []

    banned_hits = check_banned_tokens(skill_dir)
    if banned_hits:
        failures += 1
        messages.append("CHECK 1 (banned tokens): FAIL")
        messages.extend(f"  - {h}" for h in banned_hits[:10])
    else:
        messages.append("CHECK 1 (banned tokens): pass")

    fm_errors = check_frontmatter(skill_dir)
    if fm_errors:
        failures += 1
        messages.append("CHECK 2 (frontmatter): FAIL")
        messages.extend(f"  - {e}" for e in fm_errors)
    else:
        messages.append("CHECK 2 (frontmatter): pass")

    py_failures = check_python_scripts(skill_dir)
    if py_failures:
        failures += 1
        messages.append("CHECK 3 (python compile): FAIL")
        messages.extend(f"  - {f}" for f in py_failures)
    else:
        messages.append("CHECK 3 (python compile): pass")

    sh_failures = check_shell_scripts(skill_dir)
    if sh_failures:
        failures += 1
        messages.append("CHECK 4 (shell parse): FAIL")
        messages.extend(f"  - {f}" for f in sh_failures)
    else:
        messages.append("CHECK 4 (shell parse): pass")

    if public_release:
        public_hits = check_public_release(skill_dir)
        if public_hits:
            failures += 1
            messages.append("CHECK 5 (public release): FAIL")
            messages.extend(f"  - {h}" for h in public_hits[:20])
        else:
            messages.append("CHECK 5 (public release): pass")

    return failures, messages


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a skill directory.")
    parser.add_argument("--skill-dir", required=True,
                        help="Path to the skill directory.")
    parser.add_argument("--public-release", action="store_true",
                        help="Also scan for internal metadata and likely secrets.")
    args = parser.parse_args()

    failures, messages = validate(args.skill_dir, public_release=args.public_release)
    for m in messages:
        print(m)

    if failures > 0:
        print(f"\n✗ Validation completed: {failures} check(s) failed.")
        sys.exit(1)
    else:
        n_checks = 5 if args.public_release else 4
        print(f"\n✓ Validation completed: all {n_checks} checks passed.")


if __name__ == "__main__":
    main()
