---
id: "skill_9faaf5708b1346cf9286e01b47d29936"
name: "ensemble-analyze"
short-description: "Run a scientific question or dataset analysis through a configurable panel of frontier LLMs (Gemini, GPT-5.5, Claude Opus, Qwen local) in parallel, synthesize a consensus answer, and surface model ..."
display-name: "Ensemble Scientific Analysis"
description: "Run a scientific question or dataset analysis through a configurable panel of frontier LLMs (Gemini, GPT-5.5, Claude Opus, Qwen local) in parallel, synthesize a consensus answer, and surface model disagreements as signal. Use when you want higher-confidence scientific interpretation than any single model provides, or when model disagreement itself is informative."
category: "methodology"
keywords: "ensemble, multi-model, consensus, scientific analysis, LLM panel, disagreement detection, omics, interpretation"
version: "1.0.0"
created: "2026-06-19"
last-updated: "2026-06-19"
changelog:
  - version: 1.0.0
    date: 2026-06-19
    event: authored
    summary: >-
      Initial skill authored following a multi-model benchmarking review.
      Inspired by OpenRouter Fusion model concepts applied to domain-specific scientific data.
      Reusable for scientific data interpretation and methods review.
author: "Fleet Skills contributors"
status: active
promoted: true
license: MIT
---

# Ensemble Scientific Analysis

Run any scientific question or data interpretation task through a panel of frontier LLMs simultaneously. Each model reasons independently, then a synthesis pass produces a consensus answer plus a **dissent map** - where models disagreed is often where the science is most uncertain or most interesting.

Latency of 2–5 minutes per ensemble call is acceptable for deep analysis tasks. Nobody cares if a pathway interpretation takes 5 minutes instead of 2.

---

## When to Use This Skill

✅ **Use when:**
- You want higher-confidence interpretation than a single model provides
- You're analyzing a dataset where model disagreement would be scientifically informative
- The question is open-ended (pathway analysis, mechanistic interpretation, differential expression narrative)
- You're running a flagship analysis and want a "second opinion" audit trail
- Latency is acceptable (analysis tasks, not real-time queries)
- You want to know *where* uncertainty lives in an interpretation, not just the answer

❌ **Do not use when:**
- The task is simple/unambiguous (single-model is fine)
- You need a response in < 30 seconds (use single frontier model)
- The question requires live tool calls (database lookups, PubMed search) — run those first, feed results into ensemble
- Cost is tightly constrained (4× API calls vs 1×)

---

## Inputs

**Required:**
- `question` (str): The scientific question or analysis task. Can be long — include all relevant context.
- `data_context` (str | path): The dataset, results table, or evidence to reason over. Paste directly or provide file path.

**Optional:**
- `models` (list): Which models to include. Default: `["gemini", "gpt", "claude", "qwen"]`
- `synthesis_model` (str): Which model runs the final synthesis pass. Default: `"claude"` (strong at synthesis)
- `output_path` (str): Where to save the ensemble report. Default: `outputs/ensemble_<timestamp>.json`
- `weights` (dict): Per-model weights for consensus scoring. Default: equal weights.
- `dissent_threshold` (float): Fraction of models that must disagree to flag a point as contested. Default: `0.5`

---

## Outputs

```
outputs/
  ensemble_<timestamp>.json       # Full machine-readable output
  ensemble_<timestamp>_report.md  # Human-readable synthesis report
```

**JSON structure:**
```json
{
  "question": "...",
  "timestamp": "2026-06-19T...",
  "models_used": ["gemini-3.1-pro-preview", "gpt-5.5", "claude-opus-4.8", "qwen-local"],
  "individual_responses": {
    "gemini": {"response": "...", "key_claims": [...], "confidence": 0.92, "time_s": 4.2},
    "gpt":    {"response": "...", "key_claims": [...], "confidence": 0.88, "time_s": 5.1},
    "claude": {"response": "...", "key_claims": [...], "confidence": 0.91, "time_s": 1.8},
    "qwen":   {"response": "...", "key_claims": [...], "confidence": 0.79, "time_s": 12.4}
  },
  "consensus": {
    "answer": "...",
    "agreement_score": 0.87,
    "high_confidence_claims": [...],
    "contested_claims": [...]
  },
  "dissent_map": [
    {
      "claim": "...",
      "majority_position": "...",
      "dissenting_models": ["qwen"],
      "dissent_reason": "..."
    }
  ],
  "synthesis_model": "claude-opus-4.8"
}
```

---

## Quick Start

```bash
cd ensemble-analyze

# Run on a summarized differential-expression dataset
python3 scripts/ensemble_analyze.py \
  --question "What are the most novel mechanistic insights in this dataset?" \
  --data outputs/example_de_results.csv \
  --models gemini gpt claude qwen \
  --output outputs/ensemble_$(date +%Y%m%d).json

# Run on a text question only (no data file)
python3 scripts/ensemble_analyze.py \
  --question "Interpret the evidence for pathway-level regulation in this result table" \
  --models gemini gpt claude \
  --output outputs/pathway_ensemble.json
```

---

## Standard Workflow

### Step 1 — Prepare the question + data context
Frame the question specifically. Include:
- What the data represents
- What you already know / the prior interpretation
- What you specifically want the models to address
- Any constraints ("focus on novel findings", "avoid restating known background")

### Step 2 — Run parallel model calls
The script fires all model calls simultaneously via `ThreadPoolExecutor`. Each model:
- Receives identical prompt: `[SYSTEM ROLE] + [DATA CONTEXT] + [QUESTION]`
- Is blind to other models' responses
- Returns: full response + extracted key claims + self-reported confidence

Model routing:
```
gemini  → OpenRouter: google/gemini-3.1-pro-preview
gpt     → OpenRouter: openai/gpt-5.5
claude  → OpenRouter: anthropic/claude-opus-4.8
qwen    → Local model server: $LOCAL_LLM_URL OR OpenRouter: qwen/qwen3.6-35b-a3b
```

### Step 3 — Extract key claims
A lightweight pass extracts structured claims from each model's freetext response. Claims are normalized to enable agreement/disagreement detection.

### Step 4 — Build dissent map
For each claim, compute what fraction of models agree. Claims below `dissent_threshold` are flagged as **contested**. These are the scientifically interesting disagreements.

### Step 5 — Synthesis pass
The synthesis model (default: Claude) receives:
- All individual responses
- The dissent map
- Instruction: "Synthesize a consensus answer. Explicitly address contested claims. Flag where uncertainty remains."

### Step 6 — Output report
JSON + markdown report. The markdown report is designed to be pasted directly into a manuscript discussion or lab notebook entry.

---

## Model Roles (Default Weighting Rationale)

| Model | Strength in scientific context | Default weight |
|---|---|---|
| Gemini 3.1 Pro Preview | Structured data, tables, quantitative reasoning | 1.0 |
| GPT-5.5 | Clinical reasoning, broad biomedical knowledge | 1.0 |
| Claude Opus 4.8 | Literature synthesis, nuanced mechanistic argument | 1.0 (synthesis default) |
| Qwen3.6-35B | Free local sanity check, independent prior | 0.8 |

Weights are adjustable. For omics analysis, consider upweighting Gemini. For mechanistic interpretation, upweight Claude.

---

## Clarification Questions

### ⚠️ CRITICAL: ASK THIS FIRST — Is the data context ready?

The ensemble is only as good as the context it reasons over. Before running:
- Is the data summarized to a form a language model can reason over? (CSV, markdown table, bullet list of findings — not raw FASTQ)
- Have upstream analyses (DE, pathway enrichment, etc.) already been run?
- Is the question specific enough that model disagreement is meaningful?

If the data isn't ready, run the appropriate upstream skill first (e.g., `bulk-rnaseq-differential-expression`, `pathway-enrichment`) and feed those outputs in.

### Q2 — How many models?
3–4 is the sweet spot. Below 3, "consensus" is just one model agreeing with another. Above 4, synthesis becomes noisy and cost climbs. Default 4 (Gemini + GPT + Claude + Qwen) is recommended.

### Q3 — What synthesis model?
Claude Opus 4.8 is the default synthesis model because it excels at integrating conflicting positions into coherent prose. Override with Gemini if the analysis is heavily quantitative.

### Q4 — What do you do with the dissent map?
Contested claims are the output you didn't expect. Options:
- Flag them as "requires further experimental validation" in the manuscript
- Run targeted literature search on the contested point
- Design an experiment to resolve the disagreement
- Report them as an honest uncertainty statement

---

## Framework Integration

This skill is a **framework-neutral analysis helper** — any agent can call it on any scientific question.

```python
# From any agent session:
# 1. Prepare context
question = "What novel biology does this dataset reveal?"
data_context = open("results/de_genes.csv").read()

# 2. Run ensemble
# python3 ensemble-analyze/scripts/ensemble_analyze.py \
#   --question "..." --data results/de_genes.csv --output results/ensemble.json
```

**Model availability:**
- Qwen local: set `LOCAL_LLM_URL` to use a local OpenAI-compatible endpoint; otherwise the script falls back to OpenRouter `qwen/qwen3.6-35b-a3b`
- All other models: OpenRouter (always available)
- OpenRouter key: set `OPENROUTER_API_KEY` in the environment

---

## Common Issues

| Issue | Cause | Solution |
|---|---|---|
| Qwen local not responding | Fleet node offline / MLX server hung | Fall back to `qwen/qwen3.6-35b-a3b` via OpenRouter automatically |
| Models give near-identical responses | Question too simple or too constrained | Reframe to invite genuine interpretation; add "what is surprising or unexpected?" |
| Synthesis too long | Claude verbose on complex datasets | Add `--max_synthesis_tokens 800` to constrain output |
| High dissent on everything | Data context too ambiguous | Add more structured context (tables > prose summaries) |
| Rate limit 429 on OpenRouter | Running too many parallel threads | Reduce `--threads` to 2; add `--delay_s 2` between calls |
| Qwen local gives lower quality | Model quantization artifacts | Expected; downweight Qwen or use it as tie-breaker only |
| Claims extraction misses key points | Freetext parsing imperfect | Review raw individual responses; claims extraction is a heuristic |
| Synthesis model contradicts consensus | Synthesis model's own prior dominates | Switch synthesis model; or pass stricter instruction to stay faithful to majority |
| Cost too high | 4× frontier model calls on large context | Summarize data context to < 2000 tokens before passing; use Qwen local to reduce cost |
| Output JSON too large | Large individual responses stored in full | Use `--slim` flag to store only key_claims, not full responses |
| Dissent threshold wrong | Default 0.5 misses subtle disagreements | Lower to 0.25 for sensitive analyses; raise to 0.75 to report only strong disagreements |

---

## Suggested Next Steps

After ensemble analysis:
- **`pathway-enrichment`** — if ensemble flags novel pathway involvement, run formal enrichment
- **`literature-preclinical`** — search primary literature to resolve contested claims
- **`pdf-report-generation`** — convert ensemble report to a shareable PDF
- **`science-project`** — if ensemble reveals a fundable novel finding, open a science project

---

## References

**Bundled:**
- [references/ensemble_design.md](references/ensemble_design.md) — rationale for ensemble architecture, weighting, and dissent detection
- [references/prompt_templates.md](references/prompt_templates.md) — system prompts and question framing templates for scientific analysis

**External:**
- Lapuente-Brun et al. (2013) *Science* PMID 23812712 — inspiration: dedicated pools, not shared averages
- Published and vendor model evaluations can inform panel selection, but local validation should drive final model choice.

---

## Learning Loop - Field Reports

This section is the skill's living memory. Every agent that runs this skill in production appends a field report.

**Protocol:**
1. Complete the task and verify outputs
2. Write a field report below
3. Bump `version:` in frontmatter
4. Add `changelog:` entry
5. If your team uses a shared skills library, sync the completed skill there.

*No field reports yet - first agent to run this skill: please add yours.*
