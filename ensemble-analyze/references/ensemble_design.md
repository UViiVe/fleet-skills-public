# Ensemble Design — Architecture Reference

## Why an ensemble?

A single frontier model at 95–98% accuracy on USMLE still gets ~10–25 questions wrong per 500.
More importantly, for *open-ended scientific interpretation*, there is no ground truth — only
well-reasoned positions. The question isn't "who is right" but "where do smart models disagree,
and what does that disagreement tell us about scientific uncertainty?"

## The dissent map is the product

The consensus answer is useful. The dissent map is the insight.

When Gemini, GPT, and Claude all agree that "ND3 mutations cause secondary CIV loss via ROS,"
that's high confidence. When Qwen disagrees and says "CIV reduction is via supercomplex
destabilization, not ROS," that's a research question worth investigating.

In scientific analysis, model disagreement = uncertainty = where new experiments are needed.

## Weighting rationale

Equal weights (1.0) by default because model-specific strengths vary by domain
and should be validated on the user's own task distribution. Published and vendor
benchmarks can inform the default order of precedence in tie-breaking, but not
hard weights.

Domain-specific adjustment guidelines:
- Quantitative/structured data → upweight Gemini
- Literature synthesis / mechanistic argument → upweight Claude  
- Clinical reasoning → upweight GPT
- Local/cost-sensitive sanity check → Qwen

## Synthesis model choice

Claude Opus 4.8 is the default synthesis model for one reason: it is better than
the other models at holding conflicting positions simultaneously and producing
coherent prose that honestly represents the disagreement rather than papering
over it. This is an empirical observation from extended use, not a formal benchmark.

## Context window management

The synthesis pass receives all individual model responses — potentially 4× 2000 tokens
of input. At 8000 tokens of panel output + 500 tokens of dissent map + 500 tokens of
question/data, the synthesis input is ~9000 tokens. All models in the panel have
≥ 200K context windows so this is not a constraint.

For data context: summarize to < 8000 chars before passing. The models don't need
raw FASTQ or full VCF files — they need structured summaries (top DE genes, enriched
pathways, key statistics).

## Latency profile (empirical, 2026-06-19)

| Model | Avg time/call | 500Q MedQA accuracy |
|---|---|---|
| Gemini 3.1 Pro Preview | 5.2s | 98.2% |
| GPT-5.5 | 5.5s | 97.2% |
| Claude Opus 4.8 | 1.7s | 95.6% |
| Qwen3.6-35B (local) | 2.8s | 84.8% |
| Qwen3.6-35B (OpenRouter) | ~4.5s | ~85% (estimated) |

For a 4-model ensemble with parallel execution, total panel time ≈ max(individual times)
≈ 5–6 seconds. Synthesis adds ~3–5s. Total: ~10–12 seconds per ensemble call.
