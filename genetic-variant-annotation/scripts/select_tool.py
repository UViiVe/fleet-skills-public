#!/usr/bin/env python3
"""Decide between VEP and SNPEff based on organism, use case, and resources.

Public entry points:
  - ``select_annotation_tool(organism, use_case, resources, annotation_priorities)``
  - ``get_vep_recommended_config(use_case, annotation_priorities)``
  - ``get_snpeff_recommended_config(use_case)``
  - ``compare_tools(organism, use_case)``

Example:
    python scripts/select_tool.py --organism human --use-case clinical \\
        --resources standard --priorities clinical_significance pathogenicity_predictions
"""
from __future__ import annotations

import argparse
from typing import Any

_MODEL_ORGANISMS = {
    "human", "mouse", "rat", "zebrafish", "drosophila", "celegans",
}


def select_annotation_tool(
    organism: str,
    use_case: str,
    resources: str,
    annotation_priorities: list[str],
) -> str:
    """Return ``'vep'`` or ``'snpeff'``.

    Rules (first matching rule wins):
      1. Non-model organism → SNPEff (38k+ genomes).
      2. Limited compute (< 16 GB RAM) → SNPEff.
      3. Clinical use case with ClinVar priority → VEP.
      4. Regulatory-impact priority → VEP (ENCODE plugins).
      5. Two or more pathogenicity-style priorities → VEP.
      6. Default: VEP for human, SNPEff otherwise.
    """
    org = organism.lower().strip()
    if org not in _MODEL_ORGANISMS:
        return "snpeff"
    if resources.lower() == "limited":
        return "snpeff"
    if use_case == "clinical" and "clinical_significance" in annotation_priorities:
        return "vep"
    if "regulatory_impacts" in annotation_priorities:
        return "vep"
    advanced = {"pathogenicity_predictions", "protein_domains",
                "conservation_scores"}
    if len(advanced.intersection(annotation_priorities)) >= 2:
        return "vep"
    return "vep" if org == "human" else "snpeff"


def get_vep_recommended_config(
    use_case: str, annotation_priorities: list[str],
) -> dict[str, Any]:
    """Return a VEP CLI-style config dict for the chosen use case."""
    cfg: dict[str, Any] = {
        "everything":      True,
        "vcf":             True,
        "force_overwrite": True,
        "fork":            4,
        "buffer_size":     5000,
    }
    if use_case == "clinical":
        cfg["plugins"] = ["CADD", "dbNSFP,ALL", "REVEL"]
        cfg["clinical_annotations"] = True
        cfg["check_existing"] = True
    if "conservation_scores" in annotation_priorities:
        cfg.setdefault("plugins", []).extend(["Conservation", "Blosum62"])
    if "regulatory_impacts" in annotation_priorities:
        cfg["regulatory"] = True
    if "population_frequencies" in annotation_priorities:
        cfg["max_af"] = "gnomAD"
    return cfg


def get_snpeff_recommended_config(use_case: str) -> dict[str, Any]:
    cfg: dict[str, Any] = {
        "stats":          "snpeff_summary.html",
        "csv_stats":      "snpeff_stats.csv",
        "format_eff":     False,   # ANN format (modern)
        "canon":          True,
        "hgvs":           True,
        "threads":        4,
    }
    if use_case in ("clinical", "research"):
        cfg["lof"]            = True
        cfg["no_downstream"]  = True
        cfg["no_upstream"]    = True
        cfg["no_intergenic"]  = True
    return cfg


def compare_tools(organism: str = "human",
                  use_case: str = "clinical") -> dict[str, Any]:
    """Return a side-by-side pros/cons comparison for the chosen context."""
    out: dict[str, Any] = {
        "vep": {
            "pros": [
                "Comprehensive clinical databases (ClinVar, COSMIC, HGMD links)",
                "Multiple pathogenicity predictors (SIFT, PolyPhen, CADD, REVEL)",
                "Extensive regulatory annotation (ENCODE)",
                "Quarterly Ensembl release cadence",
                "VCF output compatible with most downstream tools",
            ],
            "cons": [
                "Large cache (~15–20 GB for human)",
                "Slower than SNPEff",
                "More complex initial setup",
            ],
        },
        "snpeff": {
            "pros": [
                "Fast annotation, lightweight setup (~2–5 GB per genome)",
                "Simple installation via conda",
                "38,000+ supported genome databases",
                "Excellent integration with GATK pipelines",
            ],
            "cons": [
                "Fewer pathogenicity predictors out of the box",
                "Less comprehensive clinical databases",
                "Limited regulatory annotation",
            ],
        },
    }
    if organism.lower() not in {"human", "mouse", "rat"}:
        out["snpeff"]["pros"].insert(
            0, f"Strong support for {organism}")
        out["vep"]["cons"].append(f"Limited support for {organism}")
    if use_case == "clinical":
        out["vep"]["recommended"] = True
    elif use_case == "research":
        out["snpeff"]["recommended"] = True
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--organism", required=True)
    parser.add_argument("--use-case", required=True,
                        choices=["clinical", "population", "cancer", "research"])
    parser.add_argument("--resources", required=True,
                        choices=["high-performance", "standard", "limited"])
    parser.add_argument("--priorities", nargs="+", default=["basic"])
    args = parser.parse_args()

    tool = select_annotation_tool(
        args.organism, args.use_case, args.resources, args.priorities)
    print(f"  organism : {args.organism}")
    print(f"  use case : {args.use_case}")
    print(f"  resources: {args.resources}")
    print(f"  priorities: {', '.join(args.priorities)}")
    print(f"✓ Tool selection completed: {tool}")


if __name__ == "__main__":
    main()
