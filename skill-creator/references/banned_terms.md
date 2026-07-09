# Banned Terms

The validator rejects skill files that contain certain tokens. This document explains the rules indirectly — listing the actual banned strings would make this file itself fail the check.

The authoritative list lives in `scripts/validate_skill.py` as a compiled regex.

---

## How to find the current banned list

Read the regex constant in `scripts/validate_skill.py`:

```python
BANNED_PATTERN = re.compile(r"\b(...)\b")
```

The alternation inside the parentheses is the canonical list.

---

## Categories of banned strings (no examples spelled out)

| Category | Description | Rephrase as |
|---|---|---|
| Legacy platform names | Names of related but separate platforms | "the platform", "the agent" |
| Predecessor project names | Earlier names of this project's lineage | The specific workflow name |
| Sibling project names | Distinct projects sharing infrastructure | Drop or rename |
| Internal project tags | Two-letter or short identifiers used internally | A descriptive replacement |

The validator's regex is case-sensitive and word-boundary-anchored (`\b`). Substrings inside longer words are fine — only standalone tokens trigger the check.

---

## The license-string false positive

Some Creative Commons license names contain a two-letter abbreviation that matches one of the banned tokens after the hyphen creates a word boundary. The fix is to write the license name in full:

- Instead of the abbreviated dashed form, write `Creative Commons Attribution 4.0`.
- Variant: `Creative Commons Attribution-ShareAlike 4.0` for the share-alike form.
- Variant: `Creative Commons Attribution-NonCommercial 4.0` for the non-commercial form.

This is a one-time rephrase that's permanent across all skill files.

---

## What's allowed

- `Anthropic` and `Claude` are fine — citing tools is OK; just don't claim attribution.
- Standard package names (`pandas`, `DESeq2`, etc.) are fine.
- Scientific terms containing protected substrings (e.g., names that include a stem matching one of the legacy project names) are fine — the regex is word-boundary anchored, not substring matched.

---

## Verifying your skill is clean

```bash
python scripts/validate_skill.py --skill-dir skills/<name>
```

Or, equivalently, the underlying check:

```bash
# the regex lives in scripts/validate_skill.py — copy it from there
grep -rE "<regex from validate_skill.py>" skills/<name>/ \
    --include="*.md" --include="*.py" --include="*.sh" \
    --include="*.yaml" --include="*.json" --include="*.R" --include="*.csv"
```

Zero hits is the bar.

---

## When you think you need a banned word

You almost never do. If you find yourself wanting one:

- For platform / project name references: rephrase generically.
- For legacy workflow names: use the current specific workflow name.
- For internal tag references: use a descriptive replacement.

If you genuinely have a methodological reason (e.g., a citation), put the reference in the References section using the published paper's author/year format, not the project name.

---

## Adding to the banned list

Banned-term additions go through project review. Don't edit the regex unilaterally — the cost of expanding the banned list is borne by every future skill author.

If you spot a recurring naming conflict, propose it in a project channel before changing the validator.
