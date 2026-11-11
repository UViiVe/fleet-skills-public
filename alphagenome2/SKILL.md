---
id: skill_6fb29df52df543219eb75348e507b5ac
name: alphagenome2
display-name: "AlphaGenome 2"
description: "Discover supported conditions and run a bounded public-reference regulatory-genomics prediction with AlphaGenome 2. Use when evaluating noncommercial research access or planning condition-aware variant analyses."
short-description: "Condition-aware regulatory genomics with a bring-your-own API key. Use when testing AlphaGenome 2 research access."
category: genetics
keywords: [alphagenome2, genomics, regulatory, variant, condition, expression]
version: "1.0.0"
status: active
license: MIT
author: "Stephen C. Ekker / Ekker Lab; scientific software assistance"
created: "2026-11-11"
last-updated: "2026-11-11"
provenance:
  wraps: "Google DeepMind AlphaGenome 2 API"
  upstream_source: "https://github.com/google-deepmind/alphagenome"
  citation: "AlphaGenome 2 official documentation, accessed 2026-11-11; see References."
---

# AlphaGenome 2

A portable, bring-your-own-key research workflow. The supplied helper lists reference
conditions or executes **one** 16,384-bp human reference-interval RNA-seq prediction.
It deliberately accepts neither private sequences nor single-cell uploads and does not
compute variant pathogenicity. The detailed workflow below explains how to extend an
analysis through the official v2 tutorials rather than silently using v1 interfaces.

## When to Use This Skill

- ✅ Check authenticated AlphaGenome 2 metadata access before planning research.
- ✅ Discover tissue/cell conditions and distinguish observed from imputed outputs.
- ✅ Run a small technical smoke prediction with a documented public human locus.
- ✅ Plan condition-aware noncoding-variant research with separate scientific validation.
- ❌ Clinical diagnosis, treatment selection, or claims of established pathogenicity.
- ❌ Commercial activity under a noncommercial research entitlement.
- ❌ Zebrafish or other species extrapolation presented as supported inference.
- ❌ Unapproved patient data or collaborator single-cell uploads.

## Clarification Questions

### 1. ⚠️ CRITICAL: ASK THIS FIRST
Is the intended use eligible under the current noncommercial API terms, and may every
input be sent to the external service? Stop if either answer is unclear.

2. Which species, genome assembly, tissue/cell condition and output are relevant?
3. Is this a technical access check or a biological hypothesis test?
4. Are observed tracks necessary, or are imputed/exploratory tracks acceptable?
5. What independent evidence will test any proposed biological interpretation?

## Installation

Install this directory in a generic `skills/` directory, or `~/.claude/skills/` for
Claude Code. Run the following from the installed skill directory.

| Component | Tested version | Licensing / eligibility |
|---|---|---|
| Python | 3.13 | Python 3.11+ recommended; exact tested interpreter recorded in receipt |
| alphagenome | 2.0.0 | Upstream SDK license is separate from API use terms |
| This integration | 1.0.0 | MIT; does not license the model, service or returned data |

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
```

No local GPU or model-weight download is needed. SDK dependencies need disk/RAM;
allow for API network latency and provider availability. Keep old project environments
unchanged. The SDK pin is tested; transitive dependencies are not a universal lockfile.

## Quick Start

Obtain **your own** key at https://alphagenome.google/api and review the current terms.
Do not put a key in a script, command argument, notebook output, commit or support log.
In a Bash shell, read it silently into a process environment:

```bash
read -r -s -p 'AlphaGenome API key: ' ALPHAGENOME_API_KEY; printf '\n'
export ALPHAGENOME_API_KEY
python scripts/alphagenome2.py conditions --output results/conditions --confirm-noncommercial
# Optional: exactly one public human reference prediction, not a variant-effect test.
python scripts/alphagenome2.py smoke --output results/smoke --confirm-noncommercial
unset ALPHAGENOME_API_KEY
```

Metadata is a separate RPC from prediction. A successful conditions call proves only
that this credential can access v2 metadata. Inspect `receipt.json` for what actually
ran. Latency/quota are demand-dependent; there is no promised completion time.

## Inputs

**Required:** `ALPHAGENOME_API_KEY` process environment variable, `conditions` or
`smoke` mode, a new/empty `--output` directory, and `--confirm-noncommercial`.

**Optional:** `--organism human|mouse` for metadata (default human),
`--timeout` seconds (default 60; allowed >0 to 300). Smoke is human only.

The fixed smoke uses the public tutorial GLB1L locus (chr2:219243701, GRCh38),
resized to 16,384 bp and blood CD14-positive monocytes
`UBERON:0000178_CL:0001054`. It predicts a reference interval only; the tutorial
C→T variant is used solely to center that interval. No alternate allele is submitted.
No input file, batch execution or custom expression condition is supported by the CLI.

## Outputs

| Mode | Files | Consumer / interpretation |
|---|---|---|
| conditions | `conditions.tsv`, `receipt.json` | Condition selection; includes output-specific imputation labels |
| smoke | `rna_seq.npy`, `metadata.tsv`, `receipt.json` | Numerical tracks and matching channel annotations for technical inspection |

The receipt records SDK/Python versions, API namespace, organism and execution scope.
The SDK server-default model is used; the helper does **not** claim an immutable model
version. Preserve conditions, assemblies and provenance in any downstream analysis.
Output arrays are saved without pickle. Partial files from filesystem failures must
not be treated as completion; require successful exit plus receipt and expected files.

## Standard Workflow

1. Confirm entitlement and data-sharing permission; no consent is inferred from a key.
2. Install the pinned SDK in isolation. `from alphagenome.models.v2 import dna_client`
   is required; the old v1 `ModelVersion` enumeration is not an AG2 selector.
3. Run `conditions`; search `biosample_name` in the TSV and inspect output-specific
   `observed`, `imputed_supported`, and `imputed_exploratory` labels.
4. Optionally run the single public smoke. Require nonempty finite values and inspect
   shape, channel annotations and interval in the receipt. This is not a validation of
   accuracy, clinical utility, cross-species transfer or effect direction.
5. For a real variant project, use official v2 `predict_variant` / `score_variant`
   tutorials with verified assembly/alleles, selected condition and relevant scorers.
   The helper intentionally does not expose these higher-scope operations.
6. For expression conditioning, consult the dedicated tutorial: up to 256 **human**
   cells with raw integer UMI counts and Ensembl gene IDs. Review privacy separately;
   do not normalize counts or upload collaborator data without authorization.
7. Keep prediction, effect scoring, gene mapping and experimental confirmation separate.
   Record exact parameters, model selection and observed/imputed status in a methods log.

## Rate Limits and Error Handling

Calls are serial and user-triggered; there is no automatic loop or script-level retry.
The SDK may have internal transport behavior; `--timeout` is an SDK timeout, not a
hard process-wide deadline. On HTTP 429/gRPC RESOURCE_EXHAUSTED, stop. Honor any
provider Retry-After/retry guidance, then retry manually only within the allowed quota.
If no wait is specified, wait at least 60 seconds and check current limits before one
manual retry. Never add parallel batch traffic as a workaround. Authentication failures
must not trigger retries, new keys or credential-scope changes automatically.

## Decision Points

| Decision | Default | Change only when |
|---|---|---|
| API namespace | v2 | Reproducing an explicitly v1 study separately |
| Condition | Exact supported CURIE | A separately approved expression-conditioning study exists |
| Organism | Human | Metadata for mouse is relevant; no fish support is asserted |
| Biological conclusion | None from smoke | Independent research analysis and validation exist |
| Retries | None automatic | Provider wait guidance and quota permit manual retry |

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| v2 import missing | Old SDK | Use isolated pinned 2.0.0 install |
| Missing environment key | Key not exported | Set own key privately; never commit it |
| UNAUTHENTICATED | Invalid/revoked key | Check provider account; do not echo key |
| PERMISSION_DENIED | Entitlement or policy | Resolve with provider; metadata success is not prediction proof |
| RESOURCE_EXHAUSTED / 429 | Demand/quota | Stop and honor retry guidance |
| DEADLINE_EXCEEDED | Latency | Check status; consider one longer bounded manual attempt |
| Unsupported condition | Wrong CURIE/species | Use live supported_conditions for that species |
| Observed versus imputed confusion | Assay coverage differs | Inspect labels per output, not just condition name |
| Wrong coordinates | Assembly or coordinate convention | Verify GRCh38 / GRCm39 and official interval/variant conventions |
| Mouse smoke rejected | Fixed example is human | Use metadata only; design mouse analysis separately |
| Output directory rejected | Existing nonempty directory | Choose a fresh directory; preserve previous receipts |
| No RNA-seq or nonfinite values | Service/output mismatch | Treat as failed, not a scientific null |
| Cell conditioning fails | Wrong counts/IDs/cell count | Follow official human scRNA tutorial; not accepted by helper |
| Zebrafish requested | No established supported species | Do not extrapolate human/mouse predictions as validated fish biology |

## When Scripts Fail

Use the reported error category/code without printing raw credential-bearing request
objects. Fix configuration locally; retry only when policy permits. If SDK interfaces
change, review official docs, update tests and version the integration. Do not quietly
fall back to AG1 or replace an absent result with plausible predictions.

## Best Practices

Keep secrets and output directories out of source control. Use only public data for
smoke checks. Track which scientific claims have not been validated. MIT applies to
this authored workflow only; upstream API terms govern service use, including any
commercial agreement. A noncommercial key is not blanket permission for company work.

## Suggested Next Steps

- **genetic-variant-annotation** — independently verify variant identity and assembly.
- **experimental-design-statistics** — plan controls and validation for candidate effects.
- **data-analysis-best-practices** — preserve provenance and avoid interpreting a smoke test as evidence.

These are optional companion skills; this workflow has no dependency on them.

## References

**Bundled:** [API contract](references/api-contract.md),
[Safety and interpretation](references/safety.md),
[Reproducibility](references/reproducibility.md),
[Examples](references/examples.md).

**Scripts:** `scripts/alphagenome2.py` — bounded metadata/smoke CLI;
`scripts/test_alphagenome2.py` — offline unit tests.

**Primary sources (accessed 2026-11-11):**
- https://www.alphagenomedocs.com/getting_started.html — current AG2 workflow.
- https://www.alphagenomedocs.com/colabs/v2/choosing_a_condition.html — metadata contract.
- https://www.alphagenomedocs.com/api/models_v2.html — v2 API reference.
- https://www.alphagenomedocs.com/colabs/v2/conditioning_on_single_cell_data.html — human expression conditioning.
- https://deepmind.google/blog/alphagenome-2-deciphering-how-cells-read-dna-to-better-understand-disease/ — AG2 announcement, species and research scope.
- https://github.com/google-deepmind/alphagenome — upstream SDK and license.

## Learning Loop — Field Reports

After use, record interpreter/SDK, actual access versus prediction outcomes and blockers;
never include credentials, private sequences or identifying cell metadata. Bump the
patch version for new field reports and add a body Changelog entry.

### 2026-11-11 — Initial isolated SDK qualification
SDK 2.0.0 provides the v2 namespace and documented signatures. Live qualification
outcomes are recorded in `references/reproducibility.md`; they describe one environment,
not universal entitlement or reproducibility of an unversioned server model.

## Changelog

### v1.0.0 — 2026-11-11 — Ekker Lab — authored
Initial public BYOK workflow with isolated SDK setup, bounded metadata and reference
prediction helpers, offline tests, explicit noncommercial/privacy and interpretation gates.
