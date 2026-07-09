"""Tier-2 on-target validation: did each perturbation move its own gene?

Given DE results from the screening tier, check that the target gene's own
log2FC sits in the expected direction at the expected magnitude. CRISPRi
expects down (negative log2FC), CRISPRa expects up. Anything else is treated
as an inactive guide.

Example:
    python validate_perturbations.py --de-dir results/initial_screening/ \\
        --direction down --out results/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Literal

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas numpy")


def validate_target_effect(
    de_results: Dict[str, "pd.DataFrame"],
    expected_direction: Literal["down", "up"] = "down",
    min_log2fc: float = -0.5,
    max_pval: float = 0.05,
    verbose: bool = True,
) -> "pd.DataFrame":
    """Per-perturbation validation by target log2FC direction + significance.

    For CRISPRi pass ``min_log2fc = -0.5`` (target should drop ≥30%);
    for CRISPRa pass ``min_log2fc = 0.5`` (target should rise ≥40%).
    """
    rows = []
    for gene, df in de_results.items():
        if gene not in df["gene"].values:
            rows.append({
                "perturbation": gene, "target_affected": False,
                "target_log2fc": np.nan, "target_pval": np.nan,
                "reason": "target_not_measured",
            })
            continue
        target = df[df["gene"] == gene].iloc[0]
        lfc = float(target["log2fc"])
        pval = float(target["pval"])
        if expected_direction == "down":
            ok = (lfc < min_log2fc) and (pval < max_pval)
            if ok:
                reason = "valid"
            elif lfc >= 0:
                reason = "wrong_direction_upregulated"
            elif pval >= max_pval:
                reason = "not_significant"
            else:
                reason = "insufficient_effect"
        elif expected_direction == "up":
            ok = (lfc > min_log2fc) and (pval < max_pval)
            if ok:
                reason = "valid"
            elif lfc <= 0:
                reason = "wrong_direction_downregulated"
            elif pval >= max_pval:
                reason = "not_significant"
            else:
                reason = "insufficient_effect"
        else:
            raise ValueError(f"expected_direction must be 'down' or 'up'")
        rows.append({
            "perturbation": gene, "target_affected": ok,
            "target_log2fc": lfc, "target_pval": pval, "reason": reason,
        })

    val = pd.DataFrame(rows)

    if verbose:
        n = len(val)
        ok = int(val["target_affected"].sum())
        not_meas = int((val["reason"] == "target_not_measured").sum())
        insuf = int((val["reason"] == "insufficient_effect").sum())
        not_sig = int((val["reason"] == "not_significant").sum())
        wrong = int(val["reason"].str.startswith("wrong_direction").sum())
        print(f"Validation — direction={expected_direction}")
        print(f"  total: {n}")
        print(f"  valid:                  {ok}  ({ok/n:.1%})")
        print(f"  target not measured:    {not_meas}")
        print(f"  insufficient magnitude: {insuf}")
        print(f"  not significant:        {not_sig}")
        print(f"  wrong direction:        {wrong}")
        if n:
            rate = ok / n
            grade = "good" if rate >= 0.6 else "acceptable" if rate >= 0.5 else "poor"
            print(f"  → validation rate {rate:.1%} ({grade})")
    return val


def filter_to_validated_perturbations(
    de_results: Dict[str, "pd.DataFrame"],
    validation_df: "pd.DataFrame",
    output_dir: str = "results/validated_perturbations/",
) -> Dict[str, "pd.DataFrame"]:
    """Drop DE results that failed validation; write validated CSVs."""
    valid = validation_df[validation_df["target_affected"]]["perturbation"].tolist()
    filtered = {g: de_results[g] for g in valid if g in de_results}
    print(f"Filtered to {len(filtered)}/{len(de_results)} validated perturbations")
    if output_dir:
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        for g, df in filtered.items():
            df.to_csv(out / f"{g}_validated.csv", index=False)
    return filtered


def generate_validation_report(
    validation_df: "pd.DataFrame",
    output_file: str = "results/validation_report.txt",
) -> None:
    """Plain-text per-perturbation validation status report."""
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "CRISPR Screen — Target-Gene Validation Report",
        "=" * 50,
        f"Total perturbations: {len(validation_df)}",
        "",
        "Status counts:",
    ]
    for reason, n in validation_df["reason"].value_counts().items():
        lines.append(f"  {reason}: {n} ({n/len(validation_df):.1%})")
    lines.extend(["", "Valid perturbations (sorted by target log2FC):", "-" * 50])
    valid = validation_df[validation_df["target_affected"]].sort_values("target_log2fc")
    for _, r in valid.iterrows():
        lines.append(
            f"  {r['perturbation']:20s}  log2FC={r['target_log2fc']:6.2f}  p={r['target_pval']:.2e}"
        )
    lines.extend(["", "Failed perturbations by reason:", "-" * 50])
    failed = validation_df[~validation_df["target_affected"]]
    for reason in failed["reason"].unique():
        names = failed[failed["reason"] == reason]["perturbation"].tolist()
        lines.append(f"  {reason} ({len(names)}):")
        for nm in names[:20]:
            lines.append(f"    - {nm}")
        if len(names) > 20:
            lines.append(f"    … and {len(names) - 20} more")
    Path(output_file).write_text("\n".join(lines))
    print(f"✓ Wrote {output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--de-dir", required=True,
                        help="Directory of <gene>_ttest.csv from the screening tier")
    parser.add_argument("--direction", choices=["down", "up"], required=True)
    parser.add_argument("--min-log2fc", type=float, default=None)
    parser.add_argument("--max-pval", type=float, default=0.05)
    parser.add_argument("--out", default="results/")
    args = parser.parse_args()

    if args.min_log2fc is None:
        args.min_log2fc = -0.5 if args.direction == "down" else 0.5

    de_results: Dict[str, pd.DataFrame] = {}
    for csv in Path(args.de_dir).glob("*_ttest.csv"):
        gene = csv.stem.replace("_ttest", "")
        de_results[gene] = pd.read_csv(csv)
    print(f"Loaded DE results for {len(de_results)} perturbations")

    val = validate_target_effect(
        de_results,
        expected_direction=args.direction,
        min_log2fc=args.min_log2fc,
        max_pval=args.max_pval,
    )
    out_dir = Path(args.out); out_dir.mkdir(parents=True, exist_ok=True)
    val.to_csv(out_dir / "validation_results.csv", index=False)
    generate_validation_report(val, output_file=str(out_dir / "validation_report.txt"))
    print(f"✓ Wrote validation outputs to {out_dir}")


if __name__ == "__main__":
    main()
