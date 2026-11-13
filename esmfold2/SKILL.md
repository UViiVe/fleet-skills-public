---
id: skill_35509fd9be05412296d338ac2ea5cb59
name: esmfold2
display-name: ESMFold2 structure prediction
description: "Predict a short protein structure using Biohub ESMFold2. Use when an isolated, explicitly approved local folding trial is needed."
short-description: "Source-reviewed local folding workflow; inference qualification pending. Use when testing Biohub ESMFold2."
category: methodology
keywords: [esmfold2, protein, structure, mps]
version: "3.0.0-rc.1"
created: "2026-11-13"
last-updated: "2026-11-13"
status: candidate
license: MIT
author: UViiVe
source: https://github.com/Biohub/esm
provenance: {wraps: Biohub ESMFold2, upstream_source: "https://github.com/Biohub/esm", citation: "Biohub ESM repository and model card; see references/SOURCES.md"}
---

# ESMFold2 structure prediction

This is Biohub ESMFold2, not Meta ESMFold1. This complete-folder candidate contains a monomer smoke runner, unit tests and source notes. It has **not** passed model inference. Do not equate copying this folder with installing dependencies, downloading weights, or validating a device.

## When to Use This Skill

- Use for an explicitly approved, short monomer feasibility trial.
- Do not use as a validated production pipeline or binder-affinity assay.
- Do not extrapolate a short smoke result to long proteins or complexes.

## ⚠️ CRITICAL: ASK THIS FIRST

Has the operator approved package installation, ~26GB checkpoint download, compute use and disk/memory budget? Is sequence confidentiality compatible with the chosen deployment? This runner performs local inference but initial package/model retrieval uses external services. It sends no protein sequence to a hosted folding API.

## Installation

Use a separate Python>=3.12 environment. Review dependency resolution for the target platform before approving installation:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install 'esm==3.4.3'
```

This pins the main package, not its transitive dependency lock. GPU-specific installation differs by platform. Do not use CUDA/Triton kernels on Apple Silicon. No API key is needed for the public ungated checkpoint. Do not copy credentials into this folder.

## Inputs

The bundled runner uses public ubiquitin76aa, the pinned Biohub checkpoint revision, FP32, one loop, ten sampling steps and one sample. This is a technical smoke, not a quality benchmark. Select `--device mps`, `cuda`, or `cpu`; silent CPU fallback is disabled.

## Outputs

A newly created output directory with `structure.cif` and `metrics.json`. Existing output directories are refused. Acceptance requires finite nonempty coordinates, finite pLDDT/pTM and a CIF. Monomer iPTM is not an acceptance criterion. Report wall time separately from download/load time; memory and latency remain target-dependent.

## Workflow

1. Review sources and approve isolated installation/download separately.
2. Check device availability and free memory. Never raise MPS allocation limits to force a fit.
3. Run tests without downloading weights: `python -m unittest discover -s scripts -p 'test_*.py'`.
4. After approval, run `python scripts/smoke_fold.py --device mps --output smoke-output`.
5. Monitor memory pressure; stop on sustained swap growth or an agreed timeout. No automatic retries.
6. Inspect coordinates, confidence, elapsed time and process memory; record actual outcomes.

## Common Issues

| Issue | Response |
|---|---|
| Python too old | Create Python>=3.12 environment. |
| Wrong ESMFold generation | Use Biohub/esm, not the Meta model. |
| MPS unavailable | Stop; choose another explicitly approved device. |
| CUDA unavailable | Stop instead of silently using CPU. |
| Out of memory | Stop; shorten trial or select a larger device. |
| Fast-model assumption | Fast weights remain about26GB; not a small-memory model. |
| FP8 accelerator error | This candidate explicitly selects FP32. |
| Unsupported MPS operation | Report traceback; do not claim inference pass. |
| SVD failure | Review actual version; CPU SVD cannot accept CUDA-only gesvd driver. |
| Missing current API | Check pinned esm version; old wrapper APIs may differ. |
| Nonfinite output | Fail and preserve logs for diagnosis. |
| Long runtime | Enforce operator-agreed time limit externally. |
| Existing output | Select a fresh directory; do not overwrite results. |

## Suggested Next Steps

After a real monomer pass, evaluate intended sequence length and quality settings on the actual device. Complex workflows require separate validation. Confidence scores do not demonstrate binding or biological function.

## References

See [sources and limitations](references/SOURCES.md).

## Changelog

- 3.0.0-rc.1, 2026-11-13: sanitized candidate reconstructed from current upstream interfaces; no private infrastructure/project details or legacy package duplicates. Source and dependency-free tests only; model inference pending.
