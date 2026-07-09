"""Screen primers for homodimer, heterodimer, and 3'-end pairing.

For each primer pair (and optionally a multiplex panel), compute the
thermodynamic ΔG of homodimer and heterodimer formation at the annealing
temperature. ΔG below −5 kcal/mol marks a problematic dimer.

Inputs: a list of primer sequences (CLI: --primers ATGC... CGTA...).
Outputs: structured dict with ΔG per interaction and a `problematic` flag.

Example invocation
------------------
    python -m scripts.check_dimers --primers ATGCGTACGATCGATCGAA CGATCGATCGAAGGCCATC --temp 60
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List

try:
    import primer3
except ImportError:
    sys.exit("Install with: pip install primer3-py")

try:
    from Bio.Seq import Seq
except ImportError:
    sys.exit("Install with: pip install biopython")


def analyze_dimers(
    primers: List[str],
    temperature: float = 60.0,
    dg_threshold: float = -5.0,
    mv_conc: float = 50.0,
    dv_conc: float = 0.0,
    dntp_conc: float = 0.8,
    dna_conc: float = 50.0,
) -> Dict[str, object]:
    """Compute homo- and heterodimer ΔG for every pair of primers.

    A primer set with no problematic interaction has every ΔG above
    `dg_threshold` (default −5 kcal/mol).
    """
    if len(primers) < 2:
        raise ValueError("Need at least two primers for dimer analysis.")
    primers = [p.upper().strip() for p in primers]
    for p in primers:
        if not set(p).issubset({"A", "T", "C", "G", "N"}):
            raise ValueError(f"Invalid primer: {p!r}")

    interactions = []
    for i, a in enumerate(primers):
        for j, b in enumerate(primers):
            if i > j:
                continue
            same = i == j
            calc = primer3.calc_homodimer if same else primer3.calc_heterodimer
            args = (a,) if same else (a, b)
            result = calc(
                *args,
                mv_conc=mv_conc,
                dv_conc=dv_conc,
                dntp_conc=dntp_conc,
                dna_conc=dna_conc,
                temp_c=temperature,
                output_structure=True,
            )
            dg_kcal = result.dg / 1000.0
            interactions.append(
                {
                    "type": "self-dimer" if same else "cross-dimer",
                    "primer_a": f"P{i + 1}",
                    "primer_b": f"P{j + 1}",
                    "seq_a": a,
                    "seq_b": b,
                    "dg_kcal_mol": round(dg_kcal, 2),
                    "tm_C": round(result.tm, 2) if result.tm is not None else None,
                    "structure": getattr(result, "ascii_structure", None),
                    "problematic": dg_kcal < dg_threshold,
                }
            )

    problematic = [x for x in interactions if x["problematic"]]
    return {
        "n_primers": len(primers),
        "n_interactions": len(interactions),
        "n_problematic": len(problematic),
        "has_issues": bool(problematic),
        "interactions": interactions,
        "problematic_dimers": problematic,
        "parameters": {
            "temp_C": temperature,
            "dg_threshold": dg_threshold,
            "mv_mM": mv_conc,
            "dv_mM": dv_conc,
            "dntp_mM": dntp_conc,
            "dna_nM": dna_conc,
        },
    }


def check_3prime_complementarity(
    primer_a: str,
    primer_b: str,
    window: int = 5,
) -> Dict[str, object]:
    """Count complementary bases at the 3' end of two primers.

    Three or more bases pairing at the 3' end will extend in solution.
    """
    a = primer_a.upper()[-window:]
    b = primer_b.upper()[-window:]
    b_rc = str(Seq(b).reverse_complement())
    n = 0
    for i in range(min(len(a), len(b_rc))):
        if a[-(i + 1)] == b_rc[i]:
            n += 1
        else:
            break
    return {
        "primer_a_3prime": a,
        "primer_b_3prime": b,
        "complementary_bases": n,
        "is_problematic": n >= 3,
    }


def analyze_multiplex_compatibility(
    pairs: List[Dict[str, object]],
    temperature: float = 60.0,
    dg_threshold: float = -5.0,
) -> Dict[str, object]:
    """Screen a panel of primer pairs for cross-dimers and Tm spread.

    `pairs` is a list of dicts with `forward_seq`, `reverse_seq`, and
    optionally `forward_tm`, `reverse_tm`.
    """
    flat: List[str] = []
    for pair in pairs:
        flat.append(pair["forward_seq"])
        flat.append(pair["reverse_seq"])

    dim = analyze_dimers(flat, temperature=temperature, dg_threshold=dg_threshold)

    tms = []
    for pair in pairs:
        if "forward_tm" in pair and "reverse_tm" in pair:
            tms.extend([pair["forward_tm"], pair["reverse_tm"]])

    tm_range = round(max(tms) - min(tms), 2) if tms else None
    tm_ok = tm_range is None or tm_range <= 5.0
    compatible = not dim["has_issues"] and tm_ok

    recs = []
    if dim["has_issues"]:
        recs.append(
            f"Redesign primers in {dim['n_problematic']} problematic interactions "
            f"(ΔG < {dg_threshold} kcal/mol)."
        )
    if not tm_ok:
        recs.append(f"Tm spread {tm_range}°C exceeds the 5°C multiplex window.")
    if compatible:
        recs.append("Panel is multiplex-compatible.")

    return {
        "is_compatible": compatible,
        "n_pairs": len(pairs),
        "n_primers": len(flat),
        "tm_range": tm_range,
        "tm_compatible": tm_ok,
        "dimer_analysis": dim,
        "recommendations": recs,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--primers", nargs="+", required=True,
                   help="Two or more primer sequences")
    p.add_argument("--temp", type=float, default=60.0,
                   help="Annealing temperature in °C")
    p.add_argument("--dg", type=float, default=-5.0,
                   help="ΔG threshold in kcal/mol")
    p.add_argument("--json", action="store_true",
                   help="Emit full JSON instead of a summary")
    args = p.parse_args()

    result = analyze_dimers(args.primers, temperature=args.temp, dg_threshold=args.dg)
    if args.json:
        print(json.dumps(result, indent=2, default=str))
    else:
        for it in result["interactions"]:
            flag = "FAIL" if it["problematic"] else "ok"
            print(f"[{flag}] {it['type']:>12s}  "
                  f"{it['primer_a']} x {it['primer_b']}  "
                  f"ΔG = {it['dg_kcal_mol']:>6.2f} kcal/mol")
    print(f"\n✓ Dimer screening completed: "
          f"{result['n_interactions']} interactions, "
          f"{result['n_problematic']} problematic")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
