# General Review Briefing Template

Fill in [BRACKETED] fields before spawning the review subagent.

---

## Briefing

**Topic:** [Gene name / disease name / biological process]
**Key identifiers:** [OMIM ID / gene ID / MeSH term / ICD-10 — whatever applies]
**Scope:** [narrow: single gene | medium: disease + therapeutics | broad: full landscape]

**Context to inject:**
[2-4 sentences: what is known, what the key papers are, any specific angles to cover,
any known errors in prior AI-generated reviews of this topic to avoid]

**⚠️ VERIFICATION WARNING — include this verbatim in every briefing:**
Verify every PMID by direct fetch from https://pubmed.ncbi.nlm.nih.gov/[PMID]/ —
confirm title and authors match before including. Do NOT fabricate PMIDs.
Flag all unverified citations with ⚠️. Do NOT claim drug approvals without
independent verification at FDA.gov or EMA.europa.eu.

---

## Task: [TOPIC] Comprehensive Literature Review

Write a comprehensive literature review of [TOPIC] using the following structure:

### 1. Overview
- Definition, prevalence, inheritance (if genetic), clinical features
- Why this topic matters

### 2. Molecular Mechanism / Biology
- What is known about the gene/protein/pathway
- Key functional studies (verified PMIDs)

### 3. Disease Models
- Animal models, cell lines, organoids
- What they have and haven't shown

### 4. Therapeutic Landscape
- Approved therapies (verify approval status independently)
- Clinical pipeline (Phase 1/2/3 — verify at clinicaltrials.gov)
- Preclinical approaches

### 5. Open Questions
- Key gaps in the field
- Controversies or conflicting data

### Bibliography
For each citation: Authors. "Title." *Journal* Year;Vol:Pages. PMID: XXXXXXXX [✓/⚠️/❌]
