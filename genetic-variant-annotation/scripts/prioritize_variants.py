#!/usr/bin/env python3
"""Assign an ACMG/AMP-style classification label to each variant.

The label is a deterministic function of:
  - ClinVar significance (if present),
  - Impact (HIGH / MODERATE / ...),
  - Pathogenicity scores (CADD, REVEL),
  - Allele frequency (gnomAD).

Five buckets: ``Pathogenic``, ``Likely_Pathogenic``, ``VUS``, ``Likely_Benign``,
``Benign``. The classifier is purposely transparent — it implements a small,
inspectable rule set, not the full ACMG/AMP scoring matrix. For clinical-grade
reporting, layer manual review on top.

Public function: ``prioritize_variants(df)``.

Example:
    python scripts/prioritize_variants.py --input all_variants.csv \\
        --output prioritized.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    import pandas as pd
except ImportError:
    sys.exit("Install required packages with: pip install pandas")


_RANK = {"Pathogenic": 5, "Likely_Pathogenic": 4,
         "VUS": 3, "Likely_Benign": 2, "Benign": 1}


def _clinvar_terms(sig: str) -> set[str]:
    """Normalize ClinVar significance strings into controlled terms."""
    normalized = sig.lower().replace(" ", "_").replace("-", "_")
    return {part.strip("_") for part in re.split(r"[,;/|]+", normalized) if part.strip("_")}


def _classify_row(row: pd.Series) -> tuple[str, list[str]]:
    """Return (label, evidence_list)."""
    evidence: list[str] = []
    sig = str(row.get("clinvar_sig", "")).lower()
    impact = str(row.get("impact", "")).upper()
    cadd   = pd.to_numeric(row.get("cadd_phred"), errors="coerce")
    revel  = pd.to_numeric(row.get("revel"),       errors="coerce")
    af     = pd.to_numeric(row.get("gnomad_af",
                                   row.get("af")),  errors="coerce")

    # ClinVar dominates only when the controlled term is unambiguous.
    terms = _clinvar_terms(sig)
    uncertain_terms = {
        "conflicting_interpretations_of_pathogenicity",
        "uncertain_significance",
        "not_provided",
        "other",
    }
    if terms & uncertain_terms:
        return "VUS", [f"ClinVar_{next(iter(terms & uncertain_terms))}"]

    has_pathogenic = "pathogenic" in terms
    has_likely_pathogenic = "likely_pathogenic" in terms
    has_benign = "benign" in terms
    has_likely_benign = "likely_benign" in terms

    if (has_pathogenic or has_likely_pathogenic) and (has_benign or has_likely_benign):
        return "VUS", ["ClinVar_mixed_assertions"]
    if has_pathogenic and not has_likely_pathogenic:
        return "Pathogenic", ["ClinVar_pathogenic"]
    if has_likely_pathogenic:
        return "Likely_Pathogenic", ["ClinVar_likely_pathogenic"]
    if has_likely_benign:
        return "Likely_Benign", ["ClinVar_likely_benign"]
    if has_benign:
        return "Benign", ["ClinVar_benign"]

    # No ClinVar — synthesise from impact + scores + frequency.
    if impact == "HIGH":
        evidence.append("HIGH_impact (loss-of-function-like)")
        if pd.notna(af) and af < 1e-4:
            evidence.append("AF<1e-4 (ultra-rare)")
            return "Likely_Pathogenic", evidence
        return "VUS", evidence

    if impact == "MODERATE":
        if pd.notna(cadd) and cadd >= 25:
            evidence.append(f"CADD>=25 ({cadd:.1f})")
        if pd.notna(revel) and revel >= 0.7:
            evidence.append(f"REVEL>=0.7 ({revel:.2f})")
        if pd.notna(af) and af < 1e-4:
            evidence.append("AF<1e-4 (ultra-rare)")
        if len(evidence) >= 2:
            return "Likely_Pathogenic", evidence
        if len(evidence) == 1:
            return "VUS", evidence
        if pd.notna(af) and af > 0.05:
            return "Likely_Benign", ["AF>5%"]
        return "VUS", ["MODERATE_no_supporting_evidence"]

    if impact in {"LOW", "MODIFIER"}:
        if pd.notna(af) and af > 0.05:
            return "Benign", ["LOW/MODIFIER_with_AF>5%"]
        return "Likely_Benign", ["LOW/MODIFIER_impact"]

    return "VUS", ["unclassified"]


def prioritize_variants(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``classification`` and ``evidence`` columns; sort by rank."""
    if df.empty:
        return df
    labels, evidences = [], []
    for _, row in df.iterrows():
        label, ev = _classify_row(row)
        labels.append(label); evidences.append(";".join(ev))
    out = df.copy()
    out["classification"]      = labels
    out["classification_rank"] = [_RANK[l] for l in labels]
    out["evidence"]            = evidences
    out = out.sort_values(["classification_rank", "impact_rank"],
                          ascending=[False, False],
                          kind="mergesort", na_position="last")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    print("=== Variant prioritisation ===")
    df = pd.read_csv(args.input)
    out = prioritize_variants(df)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    counts = out["classification"].value_counts().to_dict()
    for label in ("Pathogenic", "Likely_Pathogenic", "VUS",
                  "Likely_Benign", "Benign"):
        print(f"  {label}: {counts.get(label, 0):,}")
    print(f"✓ Variant prioritisation completed: {len(out):,} variants / "
          f"{args.output}")


if __name__ == "__main__":
    main()
