"""Compute primer melting temperatures using several thermodynamic models.

The default model is the SantaLucia (1998) nearest-neighbor implementation
exposed by primer3-py. Salt-adjusted Wallace and pure-%GC formulas are
provided for sanity checks and short oligos.

Inputs: a primer string (A/T/C/G/N).
Outputs: a dict with `tm_nn`, `tm_salt_adjusted`, `tm_gc`, and the
parameters used.

Example invocation
------------------
    python -m scripts.calculate_tm ATGCGTACGATCGATCGAA --method all --mg 3.0
"""

from __future__ import annotations

import argparse
import math
import sys
from typing import Dict, List

try:
    import primer3
except ImportError:
    sys.exit("Install with: pip install primer3-py")

try:
    from Bio.Seq import Seq
    from Bio.SeqUtils import MeltingTemp as mt
except ImportError:
    sys.exit("Install with: pip install biopython")


VALID_BASES = set("ATCGN")


def _validate(primer: str) -> str:
    primer = primer.upper().strip()
    if len(primer) < 10:
        raise ValueError(
            f"Primer length {len(primer)} is below the 10-nt minimum for "
            "meaningful Tm calculation."
        )
    bad = set(primer) - VALID_BASES
    if bad:
        raise ValueError(f"Primer contains non-DNA characters: {sorted(bad)}")
    return primer


def calculate_tm_comprehensive(
    primer: str,
    method: str = "nearest_neighbor",
    salt_conc: float = 50.0,
    dna_conc: float = 200.0,
    mg_conc: float = 0.0,
    dntp_conc: float = 0.0,
) -> Dict[str, object]:
    """Return Tm computed by the requested method(s).

    Parameters
    ----------
    primer : str
        Primer sequence (5' to 3').
    method : str
        One of `"nearest_neighbor"`, `"salt_adjusted"`, `"gc_content"`, `"all"`.
    salt_conc : float
        Monovalent cation concentration in mM (Na+ + K+).
    dna_conc : float
        DNA/primer concentration in nM.
    mg_conc, dntp_conc : float
        Divalent cation and dNTP concentrations (mM) for the Owczarzy correction.

    Returns
    -------
    dict
        Tm values, the method used, and the input parameters.
    """
    primer = _validate(primer)
    out: Dict[str, object] = {"primer": primer, "length": len(primer)}

    if method in ("nearest_neighbor", "all"):
        tm_p3 = primer3.calc_tm(
            primer,
            mv_conc=salt_conc,
            dv_conc=mg_conc,
            dntp_conc=dntp_conc,
            dna_conc=dna_conc,
        )
        out["tm_nn"] = round(tm_p3, 2)
        try:
            tm_bio = mt.Tm_NN(
                Seq(primer),
                Na=salt_conc,
                Mg=mg_conc,
                dNTPs=dntp_conc,
                saltcorr=7,
            )
            out["tm_nn_biopython"] = round(tm_bio, 2)
        except Exception:
            pass

    gc = sum(1 for b in primer if b in "GC")
    at = sum(1 for b in primer if b in "AT")
    gc_pct = 100.0 * gc / len(primer)

    if method in ("salt_adjusted", "all"):
        tm_wallace = 2 * at + 4 * gc
        salt_term = 16.6 * math.log10(salt_conc / 1000.0)
        gc_term = 0.41 * gc_pct
        out["tm_salt_adjusted"] = round(tm_wallace - salt_term + gc_term, 2)

    if method in ("gc_content", "all"):
        out["tm_gc"] = round(2 * at + 4 * gc, 2)

    if method == "all" and len(primer) > 13:
        out["tm_percent"] = round(64.9 + 41 * (gc_pct - 16.4) / 100, 2)

    out["parameters"] = {
        "salt_mM": salt_conc,
        "mg_mM": mg_conc,
        "dntp_mM": dntp_conc,
        "dna_nM": dna_conc,
    }
    return out


def calculate_tm_range(primers: List[str], **kwargs) -> Dict[str, object]:
    """Return the min/max/mean Tm across a set of primers."""
    if len(primers) < 2:
        raise ValueError("Need at least two primers to compute a range.")
    tms = []
    rows = []
    for p in primers:
        r = calculate_tm_comprehensive(p, method="nearest_neighbor", **kwargs)
        tms.append(r["tm_nn"])
        rows.append({"sequence": p, "tm": r["tm_nn"]})
    return {
        "n_primers": len(primers),
        "primer_tms": rows,
        "tm_min": round(min(tms), 2),
        "tm_max": round(max(tms), 2),
        "tm_mean": round(sum(tms) / len(tms), 2),
        "tm_range": round(max(tms) - min(tms), 2),
        "compatible_within_5C": (max(tms) - min(tms)) <= 5.0,
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("primer", help="Primer sequence (5' to 3')")
    p.add_argument(
        "--method",
        default="nearest_neighbor",
        choices=["nearest_neighbor", "salt_adjusted", "gc_content", "all"],
    )
    p.add_argument("--salt", type=float, default=50.0, help="Na+ in mM")
    p.add_argument("--mg", type=float, default=0.0, help="Mg2+ in mM")
    p.add_argument("--dntp", type=float, default=0.0, help="dNTP in mM")
    p.add_argument("--dna", type=float, default=200.0, help="primer in nM")
    args = p.parse_args()

    res = calculate_tm_comprehensive(
        args.primer,
        method=args.method,
        salt_conc=args.salt,
        mg_conc=args.mg,
        dntp_conc=args.dntp,
        dna_conc=args.dna,
    )
    for k, v in res.items():
        print(f"{k}: {v}")
    print(f"\n✓ Tm calculation completed for {args.primer}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
