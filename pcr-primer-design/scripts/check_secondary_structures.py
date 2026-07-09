"""Detect hairpins, self-complementarity, GC clamp, and poly-runs in primers.

Hairpin ΔG below −2 kcal/mol or self-dimer ΔG below −5 kcal/mol blocks
primer extension. The 3' end deserves extra attention because that's where
the polymerase initiates synthesis.

Inputs: a single primer sequence (CLI: positional argument).
Outputs: dict with hairpin/self-dimer ΔG, GC-clamp status, poly-run length.

Example invocation
------------------
    python -m scripts.check_secondary_structures ATGCGTACGATCGATCGAA --temp 60
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict

try:
    import primer3
except ImportError:
    sys.exit("Install with: pip install primer3-py")

try:
    from Bio.Seq import Seq
except ImportError:
    sys.exit("Install with: pip install biopython")


HAIRPIN_DG_THRESHOLD = -2.0
SELF_DIMER_DG_THRESHOLD = -5.0
THREE_PRIME_COMP_LIMIT = 8
TOTAL_COMP_LIMIT = 12
MAX_POLY_RUN = 4


def _validate(primer: str) -> str:
    primer = primer.upper().strip()
    if not set(primer).issubset({"A", "T", "C", "G", "N"}):
        raise ValueError(f"Invalid bases in primer: {primer!r}")
    return primer


def analyze_secondary_structures(
    primer: str,
    temperature: float = 60.0,
    mv_conc: float = 50.0,
    dv_conc: float = 0.0,
    dntp_conc: float = 0.8,
    dna_conc: float = 50.0,
) -> Dict[str, object]:
    """Run hairpin and self-dimer thermodynamic checks via primer3."""
    primer = _validate(primer)

    hp = primer3.calc_hairpin(
        primer,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        temp_c=temperature,
        output_structure=True,
    )
    hp_dg = hp.dg / 1000.0

    sd = primer3.calc_homodimer(
        primer,
        mv_conc=mv_conc,
        dv_conc=dv_conc,
        dntp_conc=dntp_conc,
        dna_conc=dna_conc,
        temp_c=temperature,
        output_structure=True,
    )
    sd_dg = sd.dg / 1000.0

    three_prime = check_3prime_self_complementarity(primer)
    has_issues = (
        hp_dg < HAIRPIN_DG_THRESHOLD
        or sd_dg < SELF_DIMER_DG_THRESHOLD
        or three_prime["problematic"]
    )

    return {
        "primer": primer,
        "length": len(primer),
        "hairpin": {
            "dg_kcal_mol": round(hp_dg, 2),
            "tm_C": round(hp.tm, 2) if hp.tm is not None else None,
            "structure": getattr(hp, "ascii_structure", None),
            "problematic": hp_dg < HAIRPIN_DG_THRESHOLD,
            "threshold": HAIRPIN_DG_THRESHOLD,
        },
        "self_dimer": {
            "dg_kcal_mol": round(sd_dg, 2),
            "tm_C": round(sd.tm, 2) if sd.tm is not None else None,
            "problematic": sd_dg < SELF_DIMER_DG_THRESHOLD,
            "threshold": SELF_DIMER_DG_THRESHOLD,
        },
        "three_prime_complementarity": three_prime,
        "has_issues": has_issues,
        "parameters": {
            "temp_C": temperature,
            "mv_mM": mv_conc,
            "dv_mM": dv_conc,
            "dntp_mM": dntp_conc,
            "dna_nM": dna_conc,
        },
    }


def check_3prime_self_complementarity(primer: str, window: int = 5) -> Dict[str, object]:
    """Count complementary 3' bases the primer makes with itself."""
    primer = _validate(primer)
    tail = primer[-window:]
    tail_rc = str(Seq(tail).reverse_complement())

    matches = 0
    for i in range(min(len(tail), len(tail_rc))):
        if tail[-(i + 1)] == tail_rc[i]:
            matches += 1
        else:
            break

    total = 0
    for i, ba in enumerate(primer):
        for j, bb in enumerate(primer):
            if i < j and ba == str(Seq(bb).complement()):
                total += 1

    return {
        "three_prime_window": tail,
        "three_prime_matches": matches,
        "total_self_complementary": total,
        "problematic": matches > THREE_PRIME_COMP_LIMIT or total > TOTAL_COMP_LIMIT,
        "threshold_3prime": THREE_PRIME_COMP_LIMIT,
        "threshold_total": TOTAL_COMP_LIMIT,
    }


def check_gc_clamp(primer: str) -> Dict[str, object]:
    """Verify the 3' GC-clamp is in the 1-2 sweet spot."""
    primer = _validate(primer)
    tail5 = primer[-5:]
    tail3 = primer[-3:]
    gc_in_last5 = sum(1 for b in tail5 if b in "GC")
    all_gc_last3 = all(b in "GC" for b in tail3)

    if gc_in_last5 == 0:
        verdict = "missing — add 1–2 G/C to the 3' end"
    elif all_gc_last3:
        verdict = "excessive — GGG/GCG/CCG at 3' end risks non-specific priming"
    elif 1 <= gc_in_last5 <= 2 and not all_gc_last3:
        verdict = "optimal"
    else:
        verdict = "borderline"

    return {
        "tail_last5": tail5,
        "gc_in_last5": gc_in_last5,
        "all_gc_in_last3": all_gc_last3,
        "verdict": verdict,
        "optimal": verdict == "optimal",
    }


def check_poly_runs(primer: str, max_run: int = MAX_POLY_RUN) -> Dict[str, object]:
    """Detect any A/T/G/C run longer than `max_run` bases."""
    primer = _validate(primer)
    runs = {b: 0 for b in "ATCG"}
    for base in "ATCG":
        current = best = 0
        for c in primer:
            if c == base:
                current += 1
                if current > best:
                    best = current
            else:
                current = 0
        runs[base] = best

    longest_base = max(runs, key=runs.get)
    longest = runs[longest_base]
    too_long = longest > max_run

    return {
        "per_base_runs": runs,
        "longest_base": longest_base,
        "longest_run": longest,
        "max_acceptable": max_run,
        "problematic": too_long,
        "message": (
            f"Run of {longest} {longest_base}s exceeds the {max_run}-base limit"
            if too_long
            else f"Longest run is {longest} {longest_base}s — within tolerance"
        ),
    }


def comprehensive_qc_check(primer: str, temperature: float = 60.0) -> Dict[str, object]:
    """Roll up structural QC into a single pass/fail with a list of issues."""
    sec = analyze_secondary_structures(primer, temperature=temperature)
    clamp = check_gc_clamp(primer)
    runs = check_poly_runs(primer)

    issues = []
    if sec["hairpin"]["problematic"]:
        issues.append(f"hairpin ΔG = {sec['hairpin']['dg_kcal_mol']} kcal/mol")
    if sec["self_dimer"]["problematic"]:
        issues.append(f"self-dimer ΔG = {sec['self_dimer']['dg_kcal_mol']} kcal/mol")
    if sec["three_prime_complementarity"]["problematic"]:
        issues.append("3' self-complementarity exceeds threshold")
    if not clamp["optimal"]:
        issues.append(f"GC clamp: {clamp['verdict']}")
    if runs["problematic"]:
        issues.append(runs["message"])

    return {
        "primer": primer,
        "passes": not issues,
        "issues": issues,
        "secondary_structures": sec,
        "gc_clamp": clamp,
        "poly_runs": runs,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("primer", help="Primer sequence (5' to 3')")
    p.add_argument("--temp", type=float, default=60.0)
    p.add_argument("--json", action="store_true")
    args = p.parse_args()

    res = comprehensive_qc_check(args.primer, temperature=args.temp)
    if args.json:
        print(json.dumps(res, indent=2, default=str))
    else:
        verdict = "PASS" if res["passes"] else "FAIL"
        print(f"[{verdict}] {args.primer}")
        for issue in res["issues"]:
            print(f"  - {issue}")
    print(f"\n✓ Secondary-structure QC completed for {args.primer}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
