"""Design a TaqMan assay (primers + internal probe) via primer3.

Probe Tm is set 5-10°C above primer Tm by default. Probes starting with G
are flagged because the 5' G quenches the fluorophore.

Inputs:
    --sequence TARGET_DNA  (or --fasta path)
Outputs:
    JSON of assays (primer pair + probe) with validation flags.

Example invocation
------------------
    python -m scripts.design_taqman_probes --fasta target.fasta --probe-offset 8.0
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


def design_taqman_assay(
    sequence: str,
    target_region: Optional[Tuple[int, int]] = None,
    probe_tm_offset: float = 8.0,
    primer_tm_range: Tuple[float, float] = (58.0, 62.0),
    amplicon_size_range: Tuple[int, int] = (70, 140),
    num_return: int = 5,
) -> Dict[str, object]:
    """Run primer3 with internal-oligo (probe) picking enabled."""
    sequence = sequence.upper().strip()
    if not set(sequence).issubset({"A", "T", "C", "G", "N"}):
        raise ValueError("Sequence must contain only A/T/C/G/N.")
    if len(sequence) < amplicon_size_range[0]:
        raise ValueError("Target sequence shorter than minimum amplicon.")

    opt_primer = (primer_tm_range[0] + primer_tm_range[1]) / 2.0
    opt_probe = opt_primer + probe_tm_offset

    params = {
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": opt_primer,
        "PRIMER_MIN_TM": primer_tm_range[0],
        "PRIMER_MAX_TM": primer_tm_range[1],
        "PRIMER_PAIR_MAX_DIFF_TM": 2.0,
        "PRIMER_MIN_GC": 40.0,
        "PRIMER_MAX_GC": 60.0,
        "PRIMER_PRODUCT_SIZE_RANGE": [list(amplicon_size_range)],
        "PRIMER_NUM_RETURN": num_return,
        "PRIMER_MAX_POLY_X": 4,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_PICK_INTERNAL_OLIGO": 1,
        "PRIMER_INTERNAL_OPT_SIZE": 24,
        "PRIMER_INTERNAL_MIN_SIZE": 18,
        "PRIMER_INTERNAL_MAX_SIZE": 30,
        "PRIMER_INTERNAL_OPT_TM": opt_probe,
        "PRIMER_INTERNAL_MIN_TM": opt_probe - 3,
        "PRIMER_INTERNAL_MAX_TM": opt_probe + 3,
        "PRIMER_INTERNAL_MIN_GC": 40.0,
        "PRIMER_INTERNAL_MAX_GC": 60.0,
        "PRIMER_INTERNAL_MAX_POLY_X": 4,
        "PRIMER_INTERNAL_MAX_SELF_ANY": 8,
        "PRIMER_INTERNAL_MAX_SELF_END": 3,
    }

    seq_args = {"SEQUENCE_TEMPLATE": sequence}
    if target_region is not None:
        start, end = target_region
        if start < 0 or end > len(sequence) or start >= end:
            raise ValueError(f"Invalid target_region {target_region}.")
        seq_args["SEQUENCE_TARGET"] = [start, end - start]

    raw = primer3.design_primers(seq_args, params)
    n_returned = raw.get("PRIMER_PAIR_NUM_RETURNED", 0)

    assays = []
    for i in range(n_returned):
        probe_key = f"PRIMER_INTERNAL_{i}_SEQUENCE"
        if probe_key not in raw:
            continue
        fwd_pos, fwd_len = raw[f"PRIMER_LEFT_{i}"]
        rev_pos, rev_len = raw[f"PRIMER_RIGHT_{i}"]
        probe_pos, probe_len = raw[f"PRIMER_INTERNAL_{i}"]
        fwd_tm = raw[f"PRIMER_LEFT_{i}_TM"]
        rev_tm = raw[f"PRIMER_RIGHT_{i}_TM"]
        probe_tm = raw[f"PRIMER_INTERNAL_{i}_TM"]
        probe_seq = raw[probe_key]
        probe_starts_g = probe_seq.startswith("G")

        assays.append(
            {
                "assay_id": i,
                "forward_seq": raw[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "reverse_seq": raw[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "probe_seq": probe_seq,
                "forward_tm": round(fwd_tm, 2),
                "reverse_tm": round(rev_tm, 2),
                "probe_tm": round(probe_tm, 2),
                "primer_tm_diff": round(abs(fwd_tm - rev_tm), 2),
                "probe_primer_tm_diff": round(probe_tm - (fwd_tm + rev_tm) / 2.0, 2),
                "forward_gc": round(raw[f"PRIMER_LEFT_{i}_GC_PERCENT"], 2),
                "reverse_gc": round(raw[f"PRIMER_RIGHT_{i}_GC_PERCENT"], 2),
                "probe_gc": round(raw[f"PRIMER_INTERNAL_{i}_GC_PERCENT"], 2),
                "forward_length": fwd_len,
                "reverse_length": rev_len,
                "probe_length": probe_len,
                "forward_pos": fwd_pos,
                "reverse_pos": rev_pos,
                "probe_pos": probe_pos,
                "amplicon_size": raw[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                "penalty": round(raw[f"PRIMER_PAIR_{i}_PENALTY"], 3),
                "probe_starts_with_g": probe_starts_g,
                "warning": "Probe starts with G — fluorophore quenching risk"
                if probe_starts_g
                else None,
            }
        )

    return {
        "assays": assays,
        "parameters": {
            "probe_tm_offset": probe_tm_offset,
            "primer_tm_range": primer_tm_range,
            "amplicon_size_range": amplicon_size_range,
        },
        "sequence_length": len(sequence),
        "num_assays_found": len(assays),
    }


def validate_taqman_assay(assay: Dict[str, object]) -> Dict[str, object]:
    """Per-assay validation against published TaqMan recommendations."""
    checks = {
        "amplicon_70_140": 70 <= assay["amplicon_size"] <= 140,
        "primer_tm_diff_le_2": assay["primer_tm_diff"] <= 2.0,
        "probe_offset_5_10": 5 <= assay["probe_primer_tm_diff"] <= 10,
        "probe_length_18_30": 18 <= assay["probe_length"] <= 30,
        "probe_gc_40_60": 40 <= assay["probe_gc"] <= 60,
        "probe_does_not_start_g": not assay["probe_starts_with_g"],
    }
    warnings = []
    if assay["warning"]:
        warnings.append(assay["warning"])
    return {
        "valid": all(checks.values()),
        "checks": checks,
        "failed_checks": [k for k, v in checks.items() if not v],
        "warnings": warnings,
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
    p.add_argument("--probe-offset", type=float, default=8.0)
    p.add_argument("--tm-min", type=float, default=58.0)
    p.add_argument("--tm-max", type=float, default=62.0)
    p.add_argument("--amplicon-min", type=int, default=70)
    p.add_argument("--amplicon-max", type=int, default=140)
    p.add_argument("--n", type=int, default=5, dest="num_return")
    p.add_argument("--out")
    args = p.parse_args()

    seq = _load_sequence(args)
    result = design_taqman_assay(
        sequence=seq,
        probe_tm_offset=args.probe_offset,
        primer_tm_range=(args.tm_min, args.tm_max),
        amplicon_size_range=(args.amplicon_min, args.amplicon_max),
        num_return=args.num_return,
    )

    if args.out:
        with open(args.out, "w") as f:
            json.dump(result, f, indent=2)
        print(f"✓ TaqMan assays written: {args.out} "
              f"({result['num_assays_found']} assays)")
    else:
        print(json.dumps(result, indent=2))
        print(f"\n✓ TaqMan design completed: {result['num_assays_found']} assays")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
