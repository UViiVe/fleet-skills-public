"""Design MIQE-compliant qPCR primers.

Uses primer3 with the MIQE 2.0 constraint set: 70-140 bp amplicon, primer
length 18-25 nt, Tm 58-62°C matched within 2°C, GC 40-60%, GC clamp = 1,
and strict 3' end matching.

Inputs:
    --sequence TARGET_DNA  (or --fasta path)
Outputs:
    JSON of primer pairs with a `miqe_compliant` flag per pair.

Example invocation
------------------
    python -m scripts.design_qpcr_primers --fasta target.fasta --n 5
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, Optional, Tuple

try:
    import primer3
except ImportError:
    sys.exit("Install with: pip install primer3-py")

try:
    from Bio import SeqIO
except ImportError:
    sys.exit("Install with: pip install biopython")


MIQE_AMPLICON = (70, 140)
MIQE_TM = (58.0, 62.0)
MIQE_DTM = 2.0
MIQE_GC = (40.0, 60.0)
MIQE_LENGTH = (18, 25)


def design_qpcr_primers(
    sequence: str,
    target_region: Optional[Tuple[int, int]] = None,
    amplicon_size_range: Tuple[int, int] = MIQE_AMPLICON,
    tm_match_threshold: float = MIQE_DTM,
    primer_tm_range: Tuple[float, float] = MIQE_TM,
    avoid_3prime_mismatch: bool = True,
    num_return: int = 5,
) -> Dict[str, object]:
    """Design qPCR primers and tag MIQE compliance per pair."""
    sequence = sequence.upper().strip()
    if not set(sequence).issubset({"A", "T", "C", "G", "N"}):
        raise ValueError("Sequence must contain only A/T/C/G/N.")
    if len(sequence) < amplicon_size_range[0]:
        raise ValueError(
            f"Sequence length {len(sequence)} below minimum amplicon "
            f"{amplicon_size_range[0]} bp."
        )

    opt_tm = (primer_tm_range[0] + primer_tm_range[1]) / 2.0
    params = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": primer_tm_range[0],
        "PRIMER_MAX_TM": primer_tm_range[1],
        "PRIMER_PAIR_MAX_DIFF_TM": tm_match_threshold,
        "PRIMER_MIN_GC": MIQE_GC[0],
        "PRIMER_MAX_GC": MIQE_GC[1],
        "PRIMER_PRODUCT_SIZE_RANGE": [list(amplicon_size_range)],
        "PRIMER_NUM_RETURN": max(num_return * 2, 10),
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_GC_CLAMP": 1,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY": 8,
        "PRIMER_MAX_SELF_END": 3 if avoid_3prime_mismatch else 5,
        "PRIMER_PAIR_MAX_COMPL_ANY": 8,
        "PRIMER_PAIR_MAX_COMPL_END": 3 if avoid_3prime_mismatch else 5,
        "PRIMER_MAX_HAIRPIN_TH": 47.0,
    }

    seq_args = {"SEQUENCE_TEMPLATE": sequence}
    if target_region is not None:
        start, end = target_region
        if start < 0 or end > len(sequence) or start >= end:
            raise ValueError(f"Invalid target_region {target_region}.")
        seq_args["SEQUENCE_TARGET"] = [start, end - start]

    raw = primer3.design_primers(seq_args, params)
    n_returned = raw.get("PRIMER_PAIR_NUM_RETURNED", 0)

    pairs = []
    for i in range(n_returned):
        fwd_pos, fwd_len = raw[f"PRIMER_LEFT_{i}"]
        rev_pos, rev_len = raw[f"PRIMER_RIGHT_{i}"]
        fwd_tm = raw[f"PRIMER_LEFT_{i}_TM"]
        rev_tm = raw[f"PRIMER_RIGHT_{i}_TM"]
        fwd_gc = raw[f"PRIMER_LEFT_{i}_GC_PERCENT"]
        rev_gc = raw[f"PRIMER_RIGHT_{i}_GC_PERCENT"]
        amp = raw[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
        dt = abs(fwd_tm - rev_tm)
        miqe = (
            amplicon_size_range[0] <= amp <= amplicon_size_range[1]
            and dt <= tm_match_threshold
            and MIQE_GC[0] <= fwd_gc <= MIQE_GC[1]
            and MIQE_GC[0] <= rev_gc <= MIQE_GC[1]
            and MIQE_LENGTH[0] <= fwd_len <= MIQE_LENGTH[1]
            and MIQE_LENGTH[0] <= rev_len <= MIQE_LENGTH[1]
        )
        pairs.append(
            {
                "pair_id": i,
                "forward_seq": raw[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "reverse_seq": raw[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "forward_tm": round(fwd_tm, 2),
                "reverse_tm": round(rev_tm, 2),
                "tm_diff": round(dt, 2),
                "forward_gc": round(fwd_gc, 2),
                "reverse_gc": round(rev_gc, 2),
                "forward_length": fwd_len,
                "reverse_length": rev_len,
                "forward_pos": fwd_pos,
                "reverse_pos": rev_pos,
                "amplicon_size": amp,
                "penalty": round(raw[f"PRIMER_PAIR_{i}_PENALTY"], 3),
                "miqe_compliant": miqe,
            }
        )

    # Compliant pairs first, then by primer3 penalty
    pairs.sort(key=lambda p: (not p["miqe_compliant"], p["penalty"]))
    pairs = pairs[:num_return]

    return {
        "primers": pairs,
        "miqe_compliant_present": any(p["miqe_compliant"] for p in pairs),
        "parameters": {
            "amplicon_size_range": amplicon_size_range,
            "tm_match_threshold": tm_match_threshold,
            "primer_tm_range": primer_tm_range,
            "avoid_3prime_mismatch": avoid_3prime_mismatch,
        },
        "sequence_length": len(sequence),
        "num_primers_found": len(pairs),
        "num_miqe_compliant": sum(1 for p in pairs if p["miqe_compliant"]),
    }


def validate_miqe_compliance(pair: Dict[str, object]) -> Dict[str, object]:
    """Per-pair MIQE audit."""
    checks = {
        "amplicon_70_140": MIQE_AMPLICON[0] <= pair["amplicon_size"] <= MIQE_AMPLICON[1],
        "tm_diff_le_2": pair["tm_diff"] <= MIQE_DTM,
        "tm_in_58_62": (
            MIQE_TM[0] <= pair["forward_tm"] <= MIQE_TM[1]
            and MIQE_TM[0] <= pair["reverse_tm"] <= MIQE_TM[1]
        ),
        "gc_40_60": (
            MIQE_GC[0] <= pair["forward_gc"] <= MIQE_GC[1]
            and MIQE_GC[0] <= pair["reverse_gc"] <= MIQE_GC[1]
        ),
        "length_18_25": (
            MIQE_LENGTH[0] <= pair["forward_length"] <= MIQE_LENGTH[1]
            and MIQE_LENGTH[0] <= pair["reverse_length"] <= MIQE_LENGTH[1]
        ),
    }
    return {
        "miqe_compliant": all(checks.values()),
        "checks": checks,
        "failed_checks": [k for k, v in checks.items() if not v],
    }


def _load_sequence(args) -> str:
    if args.fasta:
        return str(SeqIO.read(args.fasta, "fasta").seq)
    if args.sequence:
        return args.sequence
    raise SystemExit("Provide --sequence or --fasta")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--sequence")
    g.add_argument("--fasta")
    p.add_argument("--amplicon-min", type=int, default=MIQE_AMPLICON[0])
    p.add_argument("--amplicon-max", type=int, default=MIQE_AMPLICON[1])
    p.add_argument("--tm-min", type=float, default=MIQE_TM[0])
    p.add_argument("--tm-max", type=float, default=MIQE_TM[1])
    p.add_argument("--n", type=int, default=5, dest="num_return")
    p.add_argument("--out")
    args = p.parse_args()

    seq = _load_sequence(args)
    result = design_qpcr_primers(
        sequence=seq,
        amplicon_size_range=(args.amplicon_min, args.amplicon_max),
        primer_tm_range=(args.tm_min, args.tm_max),
        num_return=args.num_return,
    )

    if args.out:
        with open(args.out, "w") as f:
            json.dump(result, f, indent=2)
        print(f"✓ qPCR primers written: {args.out} "
              f"({result['num_miqe_compliant']}/{result['num_primers_found']} MIQE)")
    else:
        print(json.dumps(result, indent=2))
        print(f"\n✓ qPCR design completed: "
              f"{result['num_miqe_compliant']}/{result['num_primers_found']} MIQE-compliant")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
