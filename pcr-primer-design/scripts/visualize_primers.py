"""Render primer-design figures with matplotlib.

Three plots: (1) primer binding-site map on the target sequence,
(2) Tm distribution across primer pairs, (3) amplicon-size bars vs
target range. Outputs default to 300 DPI PNG, sized for inclusion in a
US-Letter document at ~80% page width (8.4 × 4 in).

Inputs:
    Dict from a design_* function (positions and Tm required).
Outputs:
    PNG/SVG file at `output_file`.

Example invocation
------------------
    python -m scripts.visualize_primers --input primers.json --kind alignment \\
        --sequence ATGCGT...AAAAA --out alignment.png
"""

from __future__ import annotations

import argparse
import json
import sys
from typing import Dict, List, Optional, Tuple

try:
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("Install with: pip install matplotlib")


FIGSIZE_LETTER = (8.4, 4.0)
DPI = 300


def plot_primer_alignment(
    sequence: str,
    pair: Dict[str, object],
    output_file: str,
    figsize: Tuple[float, float] = FIGSIZE_LETTER,
) -> str:
    """Draw forward, reverse, and amplicon spans on a coordinate axis."""
    fig, ax = plt.subplots(figsize=figsize)
    seq_length = len(sequence)

    fwd_start = pair["forward_pos"]
    fwd_end = fwd_start + pair["forward_length"]
    rev_end = pair["reverse_pos"]
    rev_start = rev_end - pair["reverse_length"]

    ax.hlines(0.5, fwd_start, rev_end, color="#cbd5e1", lw=8, label="Amplicon")
    ax.hlines(1.0, fwd_start, fwd_end, color="#dc2626", lw=10,
              label=f"Forward Tm {pair['forward_tm']:.1f}°C")
    ax.hlines(1.0, rev_start, rev_end, color="#2563eb", lw=10,
              label=f"Reverse Tm {pair['reverse_tm']:.1f}°C")
    ax.text((fwd_start + rev_end) / 2, 0.35,
            f"{pair['amplicon_size']} bp", ha="center", fontsize=10)

    ax.set_xlim(0, seq_length)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_xlabel("Position (bp)")
    ax.set_title("Primer binding sites")
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    plt.tight_layout()
    fig.savefig(output_file, dpi=DPI)
    plt.close(fig)
    return output_file


def plot_tm_distribution(
    pairs: List[Dict[str, object]],
    output_file: str,
    figsize: Tuple[float, float] = FIGSIZE_LETTER,
) -> str:
    """Bar chart of forward and reverse Tm per pair, with a 60°C reference line."""
    if not pairs:
        raise ValueError("No primer pairs to plot.")

    fig, ax = plt.subplots(figsize=figsize)
    x = list(range(1, len(pairs) + 1))
    width = 0.4
    fwd = [p["forward_tm"] for p in pairs]
    rev = [p["reverse_tm"] for p in pairs]
    ax.bar([i - width / 2 for i in x], fwd, width=width, color="#dc2626", label="Forward")
    ax.bar([i + width / 2 for i in x], rev, width=width, color="#2563eb", label="Reverse")
    ax.axhline(60.0, ls="--", color="gray", alpha=0.5, label="60°C reference")
    ax.set_xticks(x)
    ax.set_xticklabels([f"Pair {i}" for i in x])
    ax.set_ylabel("Tm (°C)")
    ax.set_title("Tm distribution across primer pairs")
    ax.legend(loc="upper right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_file, dpi=DPI)
    plt.close(fig)
    return output_file


def plot_amplicon_sizes(
    pairs: List[Dict[str, object]],
    output_file: str,
    target_range: Optional[Tuple[int, int]] = None,
    figsize: Tuple[float, float] = FIGSIZE_LETTER,
) -> str:
    """Bar chart of amplicon size with optional target band shading."""
    if not pairs:
        raise ValueError("No primer pairs to plot.")
    fig, ax = plt.subplots(figsize=figsize)
    x = list(range(1, len(pairs) + 1))
    sizes = [p["amplicon_size"] for p in pairs]
    if target_range:
        in_range = [target_range[0] <= s <= target_range[1] for s in sizes]
        colors = ["#16a34a" if k else "#dc2626" for k in in_range]
        ax.axhspan(target_range[0], target_range[1], color="#a7f3d0", alpha=0.4,
                   label=f"Target {target_range[0]}–{target_range[1]} bp")
    else:
        colors = ["#3b82f6"] * len(sizes)

    ax.bar(x, sizes, color=colors)
    ax.set_xticks(x)
    ax.set_xticklabels([f"Pair {i}" for i in x])
    ax.set_ylabel("Amplicon size (bp)")
    ax.set_title("Amplicon size by pair")
    if target_range:
        ax.legend(loc="upper right", fontsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()
    fig.savefig(output_file, dpi=DPI)
    plt.close(fig)
    return output_file


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input", required=True, help="JSON from a design_* run")
    p.add_argument("--kind", required=True,
                   choices=["alignment", "tm", "amplicon"])
    p.add_argument("--out", required=True)
    p.add_argument("--sequence", help="Required for kind=alignment")
    p.add_argument("--target-min", type=int)
    p.add_argument("--target-max", type=int)
    args = p.parse_args()

    with open(args.input) as f:
        result = json.load(f)
    pairs = result.get("primers", [])

    if args.kind == "alignment":
        if not args.sequence:
            raise SystemExit("--sequence is required for kind=alignment")
        if not pairs:
            raise SystemExit("No primer pairs in input file")
        plot_primer_alignment(args.sequence, pairs[0], args.out)
    elif args.kind == "tm":
        plot_tm_distribution(pairs, args.out)
    elif args.kind == "amplicon":
        target = None
        if args.target_min is not None and args.target_max is not None:
            target = (args.target_min, args.target_max)
        plot_amplicon_sizes(pairs, args.out, target_range=target)

    print(f"✓ Figure saved: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
