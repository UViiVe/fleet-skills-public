"""Design standard PCR primers via primer3.

Standard PCR allows wider Tm and amplicon ranges than qPCR — use this when
designing for cloning, screening, or genotyping where MIQE compliance is
not required.

Inputs:
    --sequence TARGET_DNA  (or --fasta path)
Outputs:
    JSON of primer pairs with Tm, GC%, amplicon size, and primer3 penalty.

Example invocation
------------------
    python -m scripts.design_standard_primers \\
        --sequence ATGCGT...AAAAA \\
        --amplicon-min 200 --amplicon-max 500 \\
        --tm-min 55 --tm-max 65 --n 5
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


def design_pcr_primers(
    sequence: str,
    target_region: Optional[Tuple[int, int]] = None,
    primer_size_range: Tuple[int, int] = (18, 25),
    tm_range: Tuple[float, float] = (55.0, 65.0),
    gc_range: Tuple[float, float] = (40.0, 60.0),
    amplicon_size_range: Tuple[int, int] = (100, 1000),
    num_return: int = 5,
) -> Dict[str, object]:
    """Run primer3 with standard-PCR defaults and parse the result."""
    sequence = sequence.upper().strip()
    if not set(sequence).issubset({"A", "T", "C", "G", "N"}):
        raise ValueError("Sequence must contain only A/T/C/G/N.")
    if len(sequence) < amplicon_size_range[0]:
        raise ValueError(
            f"Sequence length {len(sequence)} is shorter than the minimum "
            f"amplicon size {amplicon_size_range[0]}."
        )

    opt_tm = (tm_range[0] + tm_range[1]) / 2.0
    params = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": primer_size_range[0],
        "PRIMER_MAX_SIZE": primer_size_range[1],
        "PRIMER_OPT_TM": opt_tm,
        "PRIMER_MIN_TM": tm_range[0],
        "PRIMER_MAX_TM": tm_range[1],
        "PRIMER_MIN_GC": gc_range[0],
        "PRIMER_MAX_GC": gc_range[1],
        "PRIMER_PRODUCT_SIZE_RANGE": [list(amplicon_size_range)],
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_GC_CLAMP": 1,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY": 8,
        "PRIMER_MAX_SELF_END": 3,
        "PRIMER_PAIR_MAX_COMPL_ANY": 8,
        "PRIMER_PAIR_MAX_COMPL_END": 3,
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
        pair = {
            "pair_id": i,
            "forward_seq": raw[f"PRIMER_LEFT_{i}_SEQUENCE"],
            "reverse_seq": raw[f"PRIMER_RIGHT_{i}_SEQUENCE"],
            "forward_tm": round(fwd_tm, 2),
            "reverse_tm": round(rev_tm, 2),
            "tm_diff": round(abs(fwd_tm - rev_tm), 2),
            "forward_gc": round(raw[f"PRIMER_LEFT_{i}_GC_PERCENT"], 2),
            "reverse_gc": round(raw[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
            "forward_length": fwd_len,
            "reverse_length": rev_len,
            "forward_pos": fwd_pos,
            "reverse_pos": rev_pos,
            "amplicon_size": raw[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
            "penalty": round(raw[f"PRIMER_PAIR_{i}_PENALTY"], 3),
        }
        pairs.append(pair)

    return {
        "primers": pairs,
        "parameters": {
            "primer_size_range": primer_size_range,
            "tm_range": tm_range,
            "gc_range": gc_range,
            "amplicon_size_range": amplicon_size_range,
            "target_region": target_region,
        },
        "sequence_length": len(sequence),
        "num_primers_found": n_returned,
    }


def _load_sequence(args) -> str:
    if args.fasta:
        record = SeqIO.read(args.fasta, "fasta")
        return str(record.seq)
    if args.sequence:
        return args.sequence
    raise SystemExit("Provide --sequence or --fasta")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--sequence", help="Raw DNA target sequence")
    g.add_argument("--fasta", help="Path to FASTA file with a single record")
    p.add_argument("--amplicon-min", type=int, default=100)
    p.add_argument("--amplicon-max", type=int, default=1000)
    p.add_argument("--tm-min", type=float, default=55.0)
    p.add_argument("--tm-max", type=float, default=65.0)
    p.add_argument("--gc-min", type=float, default=40.0)
    p.add_argument("--gc-max", type=float, default=60.0)
    p.add_argument("--n", type=int, default=5, dest="num_return",
                   help="Number of primer pairs to return")
    p.add_argument("--out", help="Optional JSON output path")
    args = p.parse_args()

    sequence = _load_sequence(args)
    result = design_pcr_primers(
        sequence=sequence,
        primer_size_range=(18, 25),
        tm_range=(args.tm_min, args.tm_max),
        gc_range=(args.gc_min, args.gc_max),
        amplicon_size_range=(args.amplicon_min, args.amplicon_max),
        num_return=args.num_return,
    )

    if args.out:
        with open(args.out, "w") as f:
            json.dump(result, f, indent=2)
        print(f"✓ Standard PCR primers written: {args.out} "
              f"({result['num_primers_found']} pairs)")
    else:
        print(json.dumps(result, indent=2))
        print(f"\n✓ Standard PCR design completed: "
              f"{result['num_primers_found']} pairs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
