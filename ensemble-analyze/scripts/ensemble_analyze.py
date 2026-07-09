#!/usr/bin/env python3
"""
ensemble_analyze.py - Multi-model scientific analysis ensemble

Fires a scientific question at a panel of frontier LLMs in parallel,
extracts key claims from each response, builds a dissent map, then runs
a synthesis pass to produce a consensus answer.

Inputs:
    --question      Scientific question or analysis task (str)
    --data          Path to data file (CSV, TSV, MD, TXT) OR inline text
    --models        Space-separated model keys: gemini gpt claude qwen (default: all four)
    --synthesis     Model to run synthesis pass (default: claude)
    --output        Output JSON path (default: outputs/ensemble_<timestamp>.json)
    --threads       Parallel threads (default: 4)
    --dissent       Dissent threshold 0.0–1.0 (default: 0.5)
    --slim          Store only key_claims, not full responses (flag)
    --max_synthesis_tokens  Max tokens for synthesis response (default: 1200)

Example:
    python3 ensemble_analyze.py \\
        --question "What novel biology does this differential-expression dataset reveal?" \\
        --data outputs/de_genes.csv \\
        --models gemini gpt claude qwen \\
        --output outputs/ensemble_example.json
"""

import argparse
import json
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

try:
    import requests
except ImportError:
    sys.exit("Install with: pip install requests")

# Model Configuration

OPENROUTER_API_KEY = os.environ.get(
    "OPENROUTER_API_KEY",
    ""
)
OPENROUTER_URL = "https://openrouter.ai/api/v1/chat/completions"

MODEL_MAP = {
    "gemini": {
        "slug": "google/gemini-3.1-pro-preview",
        "backend": "openrouter",
        "display": "Gemini 3.1 Pro Preview",
    },
    "gpt": {
        "slug": "openai/gpt-5.5",
        "backend": "openrouter",
        "display": "GPT-5.5",
    },
    "claude": {
        "slug": "anthropic/claude-opus-4.8",
        "backend": "openrouter",
        "display": "Claude Opus 4.8",
    },
    "qwen": {
        "slug": "qwen/qwen3.6-35b-a3b",
        "backend": "local",
        "local_url": os.environ.get("LOCAL_LLM_URL", ""),
        "local_model": os.environ.get("LOCAL_LLM_MODEL", "qwen/qwen3.6-35b-a3b"),
        "fallback_slug": "qwen/qwen3.6-35b-a3b",
        "display": "Qwen3.6-35B (local / OpenRouter fallback)",
    },
}

SYSTEM_PROMPT = """You are an expert scientist performing rigorous data analysis.
You will be given a scientific question and supporting data or context.
Your job:
1. Answer the question thoroughly and precisely based on the provided data
2. Identify 3-7 specific KEY CLAIMS your answer makes (label them clearly as "KEY CLAIM 1:", etc.)
3. For each key claim, state your confidence: HIGH / MEDIUM / LOW
4. Flag anything surprising, novel, or that challenges conventional understanding
5. Be honest about uncertainty - do not overstate what the data supports

Do NOT simply summarize what is already known. Focus on what this specific data reveals."""

SYNTHESIS_PROMPT = """You are synthesizing the outputs of a multi-model scientific analysis panel.

You have received independent responses from {n_models} models on the same question and dataset.
Your job:
1. Identify claims that ALL or MOST models agree on - these are HIGH CONFIDENCE findings
2. Identify claims where models DISAGREE - these are scientifically contested and must be flagged
3. Write a synthesis paragraph that represents the consensus view
4. For contested claims, explicitly state "Models disagree on: [X]. [Model A] argues [Y], [Model B] argues [Z]."
5. End with: "AREAS REQUIRING FURTHER EVIDENCE: [list]"

Be faithful to what the models said. Do not introduce new claims.
Prioritize surfacing disagreement - that is where the science is most uncertain and most interesting."""


# ─── Core Functions ──────────────────────────────────────────────────────────

def load_data_context(data_arg: str | None) -> str:
    """Load data from file path or return inline string."""
    if data_arg is None:
        return ""
    p = Path(data_arg)
    if p.exists():
        content = p.read_text(encoding="utf-8")
        # Truncate very large files to ~8000 chars
        if len(content) > 8000:
            content = content[:8000] + f"\n\n[TRUNCATED — {len(p.read_text())} total chars]"
        return f"\n\nDATA CONTEXT:\n```\n{content}\n```"
    else:
        # Treat as inline text
        return f"\n\nDATA CONTEXT:\n{data_arg}"


def call_openrouter(model_slug: str, messages: list[dict], max_tokens: int = 2000) -> tuple[str, float]:
    """Call OpenRouter API. Returns (response_text, elapsed_seconds)."""
    t0 = time.time()
    headers = {
        "Authorization": f"Bearer {OPENROUTER_API_KEY}",
        "Content-Type": "application/json",
    }
    if os.environ.get("OPENROUTER_SITE_URL"):
        headers["HTTP-Referer"] = os.environ["OPENROUTER_SITE_URL"]
    if os.environ.get("OPENROUTER_APP_TITLE"):
        headers["X-Title"] = os.environ["OPENROUTER_APP_TITLE"]
    payload = {
        "model": model_slug,
        "messages": messages,
        "max_tokens": max_tokens,
        "temperature": 0.3,
    }
    r = requests.post(OPENROUTER_URL, headers=headers, json=payload, timeout=120)
    r.raise_for_status()
    elapsed = time.time() - t0
    return r.json()["choices"][0]["message"]["content"], elapsed


def call_local(local_url: str, fallback_slug: str, messages: list[dict], max_tokens: int = 2000) -> tuple[str, float, str]:
    """Call local inference node; fall back to OpenRouter if unavailable."""
    if not local_url:
        text, elapsed = call_openrouter(fallback_slug, messages, max_tokens)
        return text, elapsed, "openrouter_fallback"
    t0 = time.time()
    payload = {
        "model": MODEL_MAP["qwen"].get("local_model", fallback_slug),
        "messages": messages,
        "temperature": 0.0,
        "max_tokens": max_tokens,
    }
    try:
        r = requests.post(local_url, json=payload, timeout=60)
        r.raise_for_status()
        elapsed = time.time() - t0
        return r.json()["choices"][0]["message"]["content"], elapsed, "local"
    except Exception as e:
        print(f"  [qwen] Local inference failed ({e}), falling back to OpenRouter...")
        text, elapsed = call_openrouter(fallback_slug, messages, max_tokens)
        return text, elapsed, "openrouter_fallback"


def extract_key_claims(response_text: str) -> list[dict]:
    """Extract KEY CLAIM entries from model response."""
    claims = []
    lines = response_text.split("\n")
    for line in lines:
        line = line.strip()
        if "KEY CLAIM" in line.upper():
            # Extract claim text after the label
            parts = line.split(":", 1)
            if len(parts) > 1:
                claim_text = parts[1].strip()
                confidence = "MEDIUM"
                for conf in ["HIGH", "MEDIUM", "LOW"]:
                    if conf in line.upper():
                        confidence = conf
                        break
                if claim_text:
                    claims.append({"claim": claim_text, "confidence": confidence})
    return claims


def run_model(model_key: str, messages: list[dict], max_tokens: int) -> dict:
    """Run a single model and return structured result."""
    cfg = MODEL_MAP[model_key]
    print(f"  [{cfg['display']}] Starting...")
    try:
        if cfg["backend"] == "local":
            text, elapsed, route = call_local(
                cfg["local_url"], cfg["fallback_slug"], messages, max_tokens
            )
        else:
            text, elapsed = call_openrouter(cfg["slug"], messages, max_tokens)
            route = "openrouter"

        claims = extract_key_claims(text)
        print(f"  [{cfg['display']}] ✓ Done in {elapsed:.1f}s — {len(claims)} key claims")
        return {
            "model_key": model_key,
            "model_display": cfg["display"],
            "model_slug": cfg.get("slug") or cfg.get("fallback_slug"),
            "route": route,
            "response": text,
            "key_claims": claims,
            "response_time_s": round(elapsed, 2),
            "error": None,
        }
    except Exception as e:
        print(f"  [{cfg['display']}] ✗ Failed: {e}")
        return {
            "model_key": model_key,
            "model_display": cfg["display"],
            "response": None,
            "key_claims": [],
            "response_time_s": None,
            "error": str(e),
        }


def build_dissent_map(individual_responses: dict, dissent_threshold: float) -> list[dict]:
    """
    Compare key claims across models. Flag claims where models diverge.
    Simple heuristic: if a claim appears in < (1-dissent_threshold) of models, it's contested.
    """
    # Collect all claims with their source model
    all_claims = []
    for model_key, result in individual_responses.items():
        for claim in result.get("key_claims", []):
            all_claims.append({"model": model_key, "claim": claim["claim"], "confidence": claim["confidence"]})

    # Models that contributed
    contributing_models = [k for k, v in individual_responses.items() if not v.get("error")]
    n_models = len(contributing_models)

    if n_models < 2:
        return []

    # Group similar claims (naive: just return all with model attribution)
    # A more sophisticated version would use embedding similarity
    dissent_entries = []
    seen_models_per_claim: dict[str, list[str]] = {}

    for entry in all_claims:
        # Use first 60 chars as approximate deduplication key
        key = entry["claim"][:60].lower().strip()
        if key not in seen_models_per_claim:
            seen_models_per_claim[key] = []
        seen_models_per_claim[key].append(entry["model"])

    # Claims appearing in only one model are potentially unique/contested
    for claim_key, models_with_claim in seen_models_per_claim.items():
        coverage = len(set(models_with_claim)) / n_models
        if coverage < (1.0 - dissent_threshold):
            # Find the full claim text
            full_claim = next(
                (e["claim"] for e in all_claims if e["claim"][:60].lower().strip() == claim_key),
                claim_key
            )
            dissent_entries.append({
                "claim": full_claim,
                "models_that_made_claim": list(set(models_with_claim)),
                "models_silent_on_claim": [m for m in contributing_models if m not in models_with_claim],
                "coverage_fraction": round(coverage, 2),
                "flag": "UNIQUE_TO_SUBSET" if coverage < 0.5 else "MAJORITY_ONLY",
            })

    return sorted(dissent_entries, key=lambda x: x["coverage_fraction"])


def run_synthesis(
    question: str,
    data_context: str,
    individual_responses: dict,
    dissent_map: list[dict],
    synthesis_model_key: str,
    max_tokens: int,
) -> str:
    """Run synthesis pass over all model responses."""
    cfg = MODEL_MAP[synthesis_model_key]

    # Build synthesis input
    responses_text = ""
    for model_key, result in individual_responses.items():
        if result.get("response"):
            responses_text += f"\n\n=== {result['model_display']} ===\n{result['response']}"

    dissent_text = ""
    if dissent_map:
        dissent_text = "\n\nDISSENT MAP (claims not shared across all models):\n"
        for d in dissent_map[:10]:  # top 10 most contested
            dissent_text += f"- Claim: '{d['claim'][:100]}...'\n"
            dissent_text += f"  Made by: {d['models_that_made_claim']} | Silent: {d['models_silent_on_claim']}\n"

    synthesis_input = (
        f"ORIGINAL QUESTION: {question}\n"
        f"{data_context}\n"
        f"INDIVIDUAL MODEL RESPONSES:{responses_text}\n"
        f"{dissent_text}"
    )

    n_models = len([v for v in individual_responses.values() if not v.get("error")])
    system = SYNTHESIS_PROMPT.format(n_models=n_models)

    messages = [
        {"role": "system", "content": system},
        {"role": "user", "content": synthesis_input},
    ]

    print(f"\n  [Synthesis — {cfg['display']}] Running...")
    if cfg["backend"] == "local":
        text, elapsed, _ = call_local(cfg["local_url"], cfg["fallback_slug"], messages, max_tokens)
    else:
        text, elapsed = call_openrouter(cfg["slug"], messages, max_tokens)
    print(f"  [Synthesis] ✓ Done in {elapsed:.1f}s")
    return text


# Main

def main():
    parser = argparse.ArgumentParser(
        description="Multi-model scientific analysis ensemble"
    )
    parser.add_argument("--question", required=True, help="Scientific question or analysis task")
    parser.add_argument("--data", default=None, help="Data file path or inline data context")
    parser.add_argument("--models", nargs="+", default=["gemini", "gpt", "claude", "qwen"],
                        choices=list(MODEL_MAP.keys()), help="Models to include in panel")
    parser.add_argument("--synthesis", default="claude", choices=list(MODEL_MAP.keys()),
                        help="Model to run synthesis pass")
    parser.add_argument("--output", default=None, help="Output JSON path")
    parser.add_argument("--threads", type=int, default=4, help="Parallel threads")
    parser.add_argument("--dissent", type=float, default=0.5, help="Dissent threshold 0.0-1.0")
    parser.add_argument("--slim", action="store_true", help="Store only key_claims, not full responses")
    parser.add_argument("--max_tokens", type=int, default=2000, help="Max tokens per model response")
    parser.add_argument("--max_synthesis_tokens", type=int, default=1200, help="Max tokens for synthesis")
    args = parser.parse_args()

    # Output path
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = Path(args.output) if args.output else Path(f"outputs/ensemble_{ts}.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    report_path = output_path.with_suffix("").with_name(output_path.stem + "_report.md")

    print("=" * 60)
    print("ENSEMBLE ANALYZE")
    print(f"  Question : {args.question[:80]}...")
    print(f"  Models   : {args.models}")
    print(f"  Synthesis: {args.synthesis}")
    print(f"  Output   : {output_path}")
    print("=" * 60)

    # Load data context
    data_context = load_data_context(args.data)

    # Build messages for panel models
    messages = [
        {"role": "system", "content": SYSTEM_PROMPT},
        {"role": "user", "content": f"QUESTION: {args.question}{data_context}"},
    ]

    # Run panel in parallel
    print(f"\nRunning {len(args.models)} models in parallel...")
    t_start = time.time()
    individual_responses = {}

    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(run_model, model_key, messages, args.max_tokens): model_key
            for model_key in args.models
        }
        for future in as_completed(futures):
            result = future.result()
            individual_responses[result["model_key"]] = result

    panel_elapsed = time.time() - t_start
    n_successful = sum(1 for v in individual_responses.values() if not v.get("error"))
    print(f"\nPanel complete: {n_successful}/{len(args.models)} models succeeded in {panel_elapsed:.1f}s")

    # Build dissent map
    dissent_map = build_dissent_map(individual_responses, args.dissent)
    print(f"Dissent map: {len(dissent_map)} contested claims identified")

    # Synthesis pass
    synthesis_text = run_synthesis(
        args.question, data_context, individual_responses,
        dissent_map, args.synthesis, args.max_synthesis_tokens
    )

    # Build output
    if args.slim:
        for v in individual_responses.values():
            v.pop("response", None)

    output = {
        "question": args.question,
        "data_source": str(args.data) if args.data else None,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "models_used": args.models,
        "synthesis_model": args.synthesis,
        "panel_elapsed_s": round(panel_elapsed, 2),
        "individual_responses": individual_responses,
        "dissent_map": dissent_map,
        "consensus": {
            "synthesis": synthesis_text,
            "n_contested_claims": len(dissent_map),
            "models_succeeded": n_successful,
        },
    }

    # Save JSON
    output_path.write_text(json.dumps(output, indent=2, ensure_ascii=False))
    print(f"\n✓ Ensemble analysis complete")
    print(f"  JSON   : {output_path}")

    # Save markdown report
    md = f"""# Ensemble Analysis Report
**Question:** {args.question}
**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M')}
**Models:** {', '.join(args.models)}
**Synthesis model:** {args.synthesis}

---

## Consensus

{synthesis_text}

---

## Dissent Map ({len(dissent_map)} contested claims)

"""
    for d in dissent_map[:10]:
        md += f"**Claim:** {d['claim'][:150]}\n"
        md += f"- Made by: {', '.join(d['models_that_made_claim'])}\n"
        md += f"- Silent: {', '.join(d['models_silent_on_claim'])}\n"
        md += f"- Coverage: {d['coverage_fraction']:.0%} of models\n\n"

    md += "---\n\n## Individual Model Responses\n\n"
    for model_key, result in individual_responses.items():
        md += f"### {result['model_display']}\n"
        if result.get("error"):
            md += f"*Error: {result['error']}*\n\n"
        else:
            md += f"*{result['response_time_s']:.1f}s | {len(result['key_claims'])} key claims*\n\n"
            if result.get("response"):
                md += result["response"] + "\n\n"

    md += f"\n---\n*Generated by ensemble-analyze v1.0.0 | {datetime.now().strftime('%Y-%m-%d')}*\n"
    report_path.write_text(md)
    print(f"  Report : {report_path}")


if __name__ == "__main__":
    main()
