# Validation Checks — What and Why

`validate_skill.py` runs four checks. Zero failures is the bar for declaring a skill complete.

---

## Check 1 — No banned tokens

Searches every file in the skill directory for tokens that conflict with project identity or refer to other projects.

The authoritative list of banned tokens lives in `scripts/validate_skill.py` as the `BANNED_PATTERN` compiled regex. See [banned_terms.md](banned_terms.md) for the categories of strings that are disallowed (without spelling them out — that would make this file itself fail the check).

Why: a skill is a public artifact. Stale or wrong project names confuse future readers and create discoverability issues.

How the check works:

```python
import re
BANNED_PATTERN = re.compile(r"\b(...)\b")   # see validate_skill.py
hits = [
    (path, lineno, line)
    for path in skill_dir.rglob("*")
    for lineno, line in enumerate(path.read_text().splitlines(), 1)
    if BANNED_PATTERN.search(line)
]
```

Word boundaries (`\b`) protect substrings; the regex is case-sensitive.

The most common false-positive is when a string contains an internal two-letter abbreviation after a hyphen — the hyphen creates a word boundary. The fix is to write the affected term in its full long form (see [banned_terms.md](banned_terms.md)).

---

## Check 2 — Frontmatter parses and is complete

Loads the YAML between the first two `---` delimiters and verifies these keys exist:

- `id`
- `name`
- `display-name`
- `short-description`
- `category`
- `keywords`
- `version`
- `last-updated`

```python
import yaml
content = open(f'skills/{name}/SKILL.md').read()
fm = content.split('---')[1]
data = yaml.safe_load(fm)
required = ['id','name','display-name','short-description','category',
            'keywords','version','last-updated']
missing = [k for k in required if k not in data]
assert not missing
```

Why: frontmatter is the contract. Missing fields break agent dispatch and search.

Common failure modes:
- Stray colon inside a description not wrapped in quotes.
- Multiple `---` delimiters (e.g., one in a description) confusing the split.
- Indented YAML lists where YAML expects flat strings.

Fix: quote all string values; keep the frontmatter top-of-file with no blank lines between the opening `---` and the first key.

---

## Check 3 — Every Python script compiles

```bash
for f in skills/<name>/scripts/*.py; do python3 -m py_compile "$f" || echo "FAIL $f"; done
```

Why: a broken script is worse than no script — it suggests something works when it doesn't.

`py_compile` catches:
- Syntax errors.
- Indentation errors.
- Encoding errors.

It doesn't catch:
- Import errors (different machine may have the package).
- Runtime errors.
- Type errors.

For a deeper check, the maintainer can run `mypy` or `pyright`.

---

## Check 4 — Every shell script parses

```bash
for f in skills/<name>/scripts/*.sh; do bash -n "$f" || echo "FAIL $f"; done
```

`bash -n` parses without executing. Same rationale as the Python check.

If a skill has no shell scripts, the check is a no-op and trivially passes.

---

## Optional checks (not enforced but recommended)

### Link integrity

Every `[text](references/file.md)` or `[text](../other-skill/SKILL.md)` should point to a file that exists.

```python
import re
links = re.findall(r"\]\(([^)]+\.md)\)", skill_md_text)
for link in links:
    assert (skill_dir / link).exists(), link
```

### Common Issues row count

```bash
grep -c "^|" skills/<name>/SKILL.md
```

If the Common Issues table has fewer than 10 rows, consider augmenting.

### Reference file size

A useful skill has substantive references. A check for `wc -l` on each reference file flags placeholder content (< 30 lines is usually a placeholder).

### Q1 marker presence

```bash
grep "⚠️ CRITICAL: ASK THIS FIRST" skills/<name>/SKILL.md
```

Should match exactly once. Used by the agent for dispatch.

---

## When checks fail

The validator prints a diagnostic and exits non-zero. Common remediations:

| Failure | Fix |
|---|---|
| Banned token | Rephrase. See [banned_terms.md](banned_terms.md) for category guidance. |
| Frontmatter missing key | Add the key; copy from a similar skill. |
| YAML parse error | Quote string values; remove stray colons or `---`. |
| Python compile error | Run `python -m py_compile <file>` and fix the syntax. |
| Shell parse error | Run `bash -n <file>` and fix. |

Re-run the validator after each fix.

---

## Adding new checks

To propose a new check, edit `validate_skill.py` and add tests. Each check should:

- Be deterministic (same input → same result).
- Print a clear pass/fail line per file/skill.
- Exit non-zero on any failure.
- Be fast (no network, no heavy parsing).

Slow or flaky checks belong in CI, not in the standard validator.
