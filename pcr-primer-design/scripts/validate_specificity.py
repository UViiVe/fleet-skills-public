"""Local in-silico PCR specificity checks for primer pairs.

Searches a sequence (or an entire FASTA) for forward and reverse-complement
primer matches and reports predicted PCR products. Approximate matching
allows up to `max_mismatches` per primer.

Inputs:
    --forward / --reverse primers + --sequence or --fasta
Outputs:
    JSON list of predicted products with positions and product size.

Example invocation
------------------
    python -m scripts.validate_specificity \\
        --forward ATGCGTACGATCGATCGAA --reverse TAGCTAGGCATAGCTACGT \\
        --fasta target.fasta --max-mismatches 1
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError:
    sys.exit("Install with: pip install biopython")


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def _approximate_find_all(sequence: str, query: str, max_mismatches: int) -> List[int]:
    """Return positions where `query` matches `sequence` within mismatch budget."""
    if max_mismatches == 0:
        starts, i = [], 0
        while True:
            j = sequence.find(query, i)
            if j == -1:
                return starts
            starts.append(j)
            i = j + 1
    n, q = len(sequence), len(query)
    return [
        i
        for i in range(n - q + 1)
        if _hamming(sequence[i : i + q], query) <= max_mismatches
    ]


def in_silico_pcr(
    forward_primer: str,
    reverse_primer: str,
    sequence: str,
    max_product_size: int = 5000,
    max_mismatches: int = 0,
) -> List[Dict[str, object]]:
    """Predict PCR products formed by `forward` and `reverse` on `sequence`."""
    fwd = forward_primer.upper()
    rev = reverse_primer.upper()
    sequence = sequence.upper()
    rev_rc = str(Seq(rev).reverse_complement())

    fwd_starts = _approximate_find_all(sequence, fwd, max_mismatches)
    rev_starts = _approximate_find_all(sequence, rev_rc, max_mismatches)
    rev_ends = [s + len(rev_rc) for s in rev_starts]

    products = []
    for fs in fwd_starts:
        for re in rev_ends:
            if fs < re:
                size = re - fs
                if size <= max_product_size:
                    products.append(
                        {
                            "forward_start": fs,
                            "reverse_end": re,
                            "product_size": size,
                            "product_sequence": sequence[fs:re],
                        }
                    )
    return products


def check_primer_specificity(
    forward_primer: str,
    reverse_primer: str,
    fasta_path: str,
    max_product_size: int = 5000,
    max_mismatches: int = 1,
) -> Dict[str, object]:
    """Scan every FASTA record for in-silico PCR hits.

    Useful as an offline alternative to NCBI Primer-BLAST when a local
    transcriptome or genome is available.
    """
    all_products = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        hits = in_silico_pcr(
            forward_primer,
            reverse_primer,
            str(record.seq),
            max_product_size=max_product_size,
            max_mismatches=max_mismatches,
        )
        for h in hits:
            h["template"] = record.id
            all_products.append(h)

    return {
        "n_products": len(all_products),
        "is_specific": len(all_products) == 1,
        "products": all_products,
        "parameters": {
            "max_product_size": max_product_size,
            "max_mismatches": max_mismatches,
        },
    }


def _load_sequence(args) -> str:
    if args.fasta:
        return str(SeqIO.read(args.fasta, "fasta").seq)
    if args.sequence:
        return args.sequence
    raise SystemExit("Provide --sequence or --fasta")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--forward", required=True)
    p.add_argument("--reverse", required=True)
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--sequence")
    g.add_argument("--fasta")
    p.add_argument("--max-product-size", type=int, default=5000)
    p.add_argument("--max-mismatches", type=int, default=0)
    p.add_argument("--scan-all-records", action="store_true",
                   help="Scan every FASTA record (requires --fasta)")
    p.add_argument("--out")
    args = p.parse_args()

    if args.scan_all_records and not args.fasta:
        raise SystemExit("--scan-all-records requires --fasta")

    if args.scan_all_records:
        result = check_primer_specificity(
            args.forward, args.reverse, args.fasta,
            max_product_size=args.max_product_size,
            max_mismatches=args.max_mismatches,
        )
    else:
        sequence = _load_sequence(args)
        products = in_silico_pcr(
            args.forward, args.reverse, sequence,
            max_product_size=args.max_product_size,
            max_mismatches=args.max_mismatches,
        )
        result = {
            "n_products": len(products),
            "is_specific": len(products) == 1,
            "products": products,
            "parameters": {
                "max_product_size": args.max_product_size,
                "max_mismatches": args.max_mismatches,
            },
        }

    if args.out:
        with open(args.out, "w") as f:
            json.dump(result, f, indent=2)
        print(f"✓ Specificity check written: {args.out} "
              f"({result['n_products']} products)")
    else:
        print(json.dumps(result, indent=2))
        print(f"\n✓ Specificity check completed: "
              f"{result['n_products']} predicted product(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
