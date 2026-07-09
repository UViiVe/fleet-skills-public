# Sign Conventions and Identifier Quirks

The single largest source of bugs in DepMap analyses is sign confusion.
DepMap gene-effect scores point the opposite direction from what people
intuitively expect when they hear "essentiality." This reference is
defensive: if you read only one document in this skill, read this one.

---

## The Sign Convention

In a DepMap gene-effect matrix:

| Value | Meaning |
|---|---|
| **< 0** | Knockout *reduces* fitness — gene is essential |
| **0** | Knockout has no effect |
| **> 0** | Knockout *increases* fitness — rare; usually a growth suppressor |

The reference scale:

```
strongly essential  ←———————————————————————————→  pro-growth (suppressor KO)
       -2.0          -1.0          0           +0.5
                       ↑
            median of pan-essentials
```

`-1.0` is calibrated to the median of pan-essential genes (RPS / RPL,
POLR2A, EIF3A, …). It is the "yes this is real essentiality" sanity-check
floor. Anything more negative is highly essential.

---

## Why This Is Confusing

In a CRISPR fitness screen, the raw measurement is "how many cells with this
guide are left at the end relative to the start." A guide knocking out an
essential gene leaves *fewer* cells, so log fold-change is negative. CERES
and Chronos preserve that sign — they just normalize and de-confound.

But when people talk about "essentiality scores," they often mean
"high score = highly essential." DepMap's convention is the *physics* sign:
negative effect on fitness. That is internally consistent and survives every
downstream calculation correctly *as long as you remember it.*

---

## Always Invert Before Correlating

To compute a Spearman correlation with another quantity such that "more X →
more essential" comes out as a positive correlation, **invert the gene-effect
score** before correlating:

```python
essentiality = -gene_effect             # invert once
corr, pval = spearmanr(expression, essentiality)
```

Interpretation after inversion:
- **positive correlation** → higher expression ↔ higher essentiality
- **negative correlation** → higher expression ↔ lower essentiality

If you do not invert, your interpretation must be flipped:
- **positive correlation** → higher expression ↔ *lower* essentiality
- **negative correlation** → higher expression ↔ *higher* essentiality

Both are correct math; only one matches the intuition behind the natural-
language description of the result. Pick the convention that matches your
writing style and apply it consistently. The recommendation here is to
invert at load time and never think about it again.

---

## Common Mistakes

### Looking for the "most negative" raw correlation when you want strongest dependency

The intuition is: "find the gene whose expression most strongly predicts
essentiality." If you correlate against the **raw** (un-inverted) score, the
gene you want is the one with the *most negative* correlation. If you
correlate against the **inverted** score, it is the *most positive*. People
mix these up.

```python
# WRONG: looking for "most positive correlation" with raw scores
# returns genes whose expression goes UP when essentiality DROPS — opposite of intent
corr = expression.corrwith(gene_effect_raw, method="spearman")
top = corr.nlargest(20)

# RIGHT: invert first
corr = expression.corrwith(-gene_effect_raw, method="spearman")
top = corr.nlargest(20)
```

### Pan-essential noise dominating the ranking

If you do not pre-filter pan-essentials, the top "selective" hits are
overwhelmingly RPS / RPL / proteasome subunits, because their essentiality
correlates with anything that correlates with general cell health (cell
size, doubling time, etc.). Always filter against the pan-essential set.

### Sign-flipping across publications

Some papers report `-gene_effect` and call it "essentiality score." Others
report the raw score. Some use a `gene_effect_normalized` quantity that is
re-scaled. When reproducing prior work, read the methods carefully — the
plot direction is usually a giveaway (highly essential genes at the top vs
bottom of a waterfall).

---

## Sanity Checks at Load Time

Run these on every fresh load:

```python
known_pan_essentials = ["RPS14", "RPL5", "POLR2A", "EIF3A"]
for sym in known_pan_essentials:
    matching = gene_effect.index[gene_effect.index.str.startswith(f"{sym} (")]
    if len(matching):
        med = gene_effect.loc[matching[0]].median()
        assert med < -0.8, f"sign convention may be off for {sym}: median={med:.2f}"
```

If any of these assertions fails, the data is inverted, mislabelled, or you
have hit a CERES-vs-Chronos mismatch. Stop and figure it out.

---

## Gene Identifiers

Current DepMap releases use a composite row index: `SYMBOL (entrez_id)`,
e.g. `MLL1 (4297)`, `TP53 (7157)`.

Why the composite? Symbol-only is not unique across releases — when HGNC
renames a gene, two rows for the same Entrez ID can show up under different
symbols across releases. The Entrez ID is the stable key.

To match by symbol only:

```python
def find_gene(matrix, symbol):
    candidates = matrix.index[matrix.index.str.startswith(f"{symbol} (")]
    if len(candidates) == 0:
        raise KeyError(f"{symbol!r} not in matrix")
    if len(candidates) > 1:
        print(f"warning: multiple matches for {symbol!r}: {list(candidates)}")
    return matrix.loc[candidates[0]]
```

To match by Entrez ID exactly:

```python
def find_by_entrez(matrix, entrez):
    candidates = matrix.index[matrix.index.str.endswith(f" ({entrez})")]
    if len(candidates) == 1:
        return matrix.loc[candidates[0]]
    raise KeyError(f"no/multiple matches for Entrez {entrez}")
```

---

## Cell-Line Identifiers

- **`ModelID`** (e.g. `ACH-000001`) — the stable, recommended key
- **`CellLineName`** (e.g. `A549_LUNG`) — human-readable; **not unique** across
  databases
- **`StrippedCellLineName`** (e.g. `A549`) — alphanumeric; also not unique

Use `ModelID` for joins. Use the names for display only.

In releases before 22Q1, the file `sample_info.csv` exposed `DepMap_ID` —
that is the same thing as `ModelID`.

---

## A One-Page Decision Card

When you receive a DepMap-derived analysis or are about to write one:

1. ☐ Did the data come from a tagged DepMap release? Which one?
2. ☐ Is the score CERES (≤ 20Q4) or Chronos (≥ 21Q1)?
3. ☐ Have known pan-essentials been verified to be < -0.8 (sign check)?
4. ☐ When correlating, was the score inverted before correlation?
5. ☐ Are pan-essentials excluded before claiming selectivity?
6. ☐ Is the cohort size adequate for the test used (≥10 for asymptotic,
      otherwise permute)?
7. ☐ Is the figure caption stating the score direction in plain English?

If you cannot answer all seven, the result is not yet reliable.
