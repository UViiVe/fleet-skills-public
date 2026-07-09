# PMID Verification Protocol

## Mandatory check for every citation

```
Fetch: https://pubmed.ncbi.nlm.nih.gov/[PMID]/
Parse: title, first author, journal, year
Confirm: matches the claim being cited
```

**Status codes:**
- ✓ Verified — title/authors confirmed by direct fetch
- ⚠️ Secondary — cited from a review paper, not directly verified
- ❌ REJECTED — fetched paper is wrong topic (state what it actually covers)

---

## Common Hallucination Patterns — Check for These Actively

| Pattern | How to catch |
|---|---|
| Fabricated PMID | PubMed 404 or redirects to unrelated paper |
| Real PMID, wrong claim | Title/abstract don't support the statement being cited |
| Wrong drug approval status | Says "FDA approved" for a drug in Phase 2 — verify at FDA.gov |
| Wrong approval date | Off by years — verify at FDA label database |
| Wrong approval agency | Confuses EMA and FDA approvals |
| Invented trial (NCT number) | Verify at clinicaltrials.gov before including any NCT |
| Overstated effect size | "Dramatically reduces" when paper shows modest effect — fetch abstract |
| Wrong gene/protein name | Especially common with paralogs — confirm HGNC symbol |
| Wrong species | Claims human data from mouse study |
| Missing negative results | Only cites positive studies — ask explicitly for null results |

---

## Verification workflow for a 20-citation review

```bash
# For each PMID in the review:
curl -s "https://pubmed.ncbi.nlm.nih.gov/[PMID]/" | grep -o '<title>[^<]*'
# Confirm the title matches your citation claim
# Flag ⚠️ if not fetched; ❌ if mismatched
```

Or via ncbi-api skill:
```
ncbi-api fetch_abstract --pmid [PMID]
# Returns: title, authors, abstract, journal, year
```

---

## Audit pass prompt template

```
Audit the following literature review.
For each cited PMID, fetch https://pubmed.ncbi.nlm.nih.gov/[PMID]/
and confirm the title matches the claim being cited.

Return a table:
| PMID | Claimed topic | Actual title | Status (✓/⚠️/❌) |

Flag any drug approval claims and state the source used to verify them.
```
