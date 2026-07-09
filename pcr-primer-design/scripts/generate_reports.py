"""Render primer-design results into Markdown, plain text, or HTML reports.

Wraps the design output, validation findings, and (optionally) a MIQE
checklist into a single human-readable document.

Inputs:
    Dict from a design_* function + optional validation_results.
Outputs:
    Report string in the chosen format (also writes to file when --out is set).

Example invocation
------------------
    python -m scripts.generate_reports --input primers.json --fmt markdown --out report.md
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime
from typing import Dict, List, Optional


def generate_primer_report(
    primers: Dict[str, object],
    validation_results: Optional[Dict[str, object]] = None,
    output_format: str = "markdown",
    include_miqe_checklist: bool = False,
) -> str:
    """Render the report in the requested format."""
    renderers = {
        "markdown": _render_markdown,
        "text": _render_text,
        "html": _render_html,
    }
    if output_format not in renderers:
        raise ValueError(f"Unknown format {output_format!r}.")
    return renderers[output_format](primers, validation_results, include_miqe_checklist)


def _render_markdown(primers, validation, include_miqe):
    L: List[str] = []
    L.append("# PCR Primer Design Report")
    L.append("")
    L.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"Sequence length: {primers.get('sequence_length', 'NA')} bp")
    L.append(f"Pairs found: {primers.get('num_primers_found', 0)}")
    L.append("")

    L.append("## Design Parameters")
    for k, v in (primers.get("parameters") or {}).items():
        if v is not None:
            L.append(f"- {k.replace('_', ' ')}: {v}")
    L.append("")

    L.append("## Primer Pairs")
    for i, pair in enumerate(primers.get("primers", [])[:5], 1):
        L.append(f"### Pair {i}")
        if "miqe_compliant" in pair:
            L.append(f"MIQE compliant: {pair['miqe_compliant']}")
        L.append("")
        L.append(f"- Forward: `{pair['forward_seq']}`")
        L.append(f"- Reverse: `{pair['reverse_seq']}`")
        if "probe_seq" in pair:
            L.append(f"- Probe:   `{pair['probe_seq']}`")
        L.append("")
        L.append(f"- Amplicon: {pair['amplicon_size']} bp")
        L.append(f"- Forward Tm: {pair['forward_tm']}°C (GC {pair['forward_gc']}%, "
                 f"{pair['forward_length']} nt)")
        L.append(f"- Reverse Tm: {pair['reverse_tm']}°C (GC {pair['reverse_gc']}%, "
                 f"{pair['reverse_length']} nt)")
        L.append(f"- ΔTm: {pair['tm_diff']}°C")
        if "probe_tm" in pair:
            L.append(f"- Probe Tm: {pair['probe_tm']}°C (GC {pair['probe_gc']}%)")
        L.append(f"- primer3 penalty: {pair['penalty']:.3f}")
        L.append("")

    if validation:
        L.append("## Validation")
        spec = validation.get("specificity")
        if spec is not None:
            L.append("### Specificity")
            L.append(f"- Single product expected: {spec.get('is_specific')}")
            L.append(f"- On-target hits: {spec.get('on_target_hits', 'NA')}")
            L.append(f"- Off-target hits: {spec.get('off_target_hits', 'NA')}")
            L.append("")

        dim = validation.get("dimers")
        if dim is not None:
            L.append("### Primer dimers")
            if dim.get("has_issues") or dim.get("n_problematic", 0) > 0:
                L.append(f"- Problematic interactions: "
                         f"{dim.get('n_problematic', len(dim.get('problematic_dimers', [])))}")
                for it in dim.get("problematic_dimers", []):
                    L.append(f"  - {it['type']}: ΔG = {it['dg_kcal_mol']} kcal/mol")
            else:
                L.append("- No problematic dimers")
            L.append("")

        sec = validation.get("secondary_structures")
        if sec is not None:
            L.append("### Secondary structures")
            if sec.get("has_issues"):
                if sec.get("hairpin", {}).get("problematic"):
                    L.append(f"- Hairpin ΔG = {sec['hairpin']['dg_kcal_mol']} kcal/mol")
                if sec.get("self_dimer", {}).get("problematic"):
                    L.append(f"- Self-dimer ΔG = {sec['self_dimer']['dg_kcal_mol']} kcal/mol")
            else:
                L.append("- No problematic structures")
            L.append("")

    if include_miqe:
        L.append("## MIQE 2.0 Checklist (skeleton)")
        L.extend(_miqe_checklist_lines())
        L.append("")

    L.append("## Recommendations")
    L.extend(_recommendations(primers, validation))
    return "\n".join(L)


def _render_text(primers, validation, include_miqe):
    L = []
    bar = "=" * 60
    L.append(bar)
    L.append("PCR PRIMER DESIGN REPORT")
    L.append(bar)
    L.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"Sequence length: {primers.get('sequence_length', 'NA')} bp")
    L.append(f"Pairs found: {primers.get('num_primers_found', 0)}")
    L.append("")

    for i, pair in enumerate(primers.get("primers", [])[:3], 1):
        L.append(f"Pair {i}:")
        L.append(f"  Forward: {pair['forward_seq']}")
        L.append(f"  Reverse: {pair['reverse_seq']}")
        if "probe_seq" in pair:
            L.append(f"  Probe:   {pair['probe_seq']}")
        L.append(f"  Amplicon: {pair['amplicon_size']} bp")
        L.append(f"  Tm: F {pair['forward_tm']:.1f}, R {pair['reverse_tm']:.1f}, "
                 f"ΔTm {pair['tm_diff']:.1f}")
        L.append("")
    return "\n".join(L)


def _render_html(primers, validation, include_miqe):
    rows = []
    for i, pair in enumerate(primers.get("primers", [])[:5], 1):
        rows.append(
            f"<tr><td>{i}</td>"
            f"<td><code>{pair['forward_seq']}</code></td>"
            f"<td><code>{pair['reverse_seq']}</code></td>"
            f"<td>{pair['amplicon_size']}</td>"
            f"<td>{pair['forward_tm']:.1f}</td>"
            f"<td>{pair['reverse_tm']:.1f}</td>"
            f"<td>{pair['tm_diff']:.1f}</td></tr>"
        )

    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        "<title>PCR Primer Design Report</title>"
        "<style>"
        "body{font-family:system-ui,sans-serif;margin:24px}"
        "table{border-collapse:collapse;width:100%}"
        "th,td{border:1px solid #ddd;padding:6px}"
        "th{background:#3a5fcd;color:#fff;text-align:left}"
        "code{font-family:ui-monospace,Menlo,monospace;background:#f1f3f5;padding:2px 4px}"
        "</style></head><body>"
        "<h1>PCR Primer Design Report</h1>"
        f"<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>"
        f"<p>Sequence length: {primers.get('sequence_length', 'NA')} bp · "
        f"Pairs found: {primers.get('num_primers_found', 0)}</p>"
        "<table><thead><tr>"
        "<th>Pair</th><th>Forward</th><th>Reverse</th><th>Amplicon (bp)</th>"
        "<th>Fwd Tm (°C)</th><th>Rev Tm (°C)</th><th>ΔTm (°C)</th>"
        "</tr></thead><tbody>"
        + "".join(rows)
        + "</tbody></table></body></html>"
    )


def _miqe_checklist_lines() -> List[str]:
    return [
        "### Experimental Design",
        "- [ ] qPCR purpose stated",
        "- [ ] ≥ 3 biological replicates",
        "- [ ] ≥ 2 technical replicates per sample",
        "",
        "### Sample",
        "- [ ] RNA/DNA quality verified (RIN, 260/280)",
        "- [ ] Storage conditions documented",
        "- [ ] RT method documented (RT-qPCR)",
        "",
        "### Assay Validation",
        "- [x] Primer sequences provided",
        "- [ ] Standard curve (R² ≥ 0.98)",
        "- [ ] Efficiency in 0.90–1.10",
        "- [ ] Linear range ≥ 5 logs",
        "- [ ] Single melt-curve peak",
        "",
        "### Protocol",
        "- [ ] Reaction volume specified",
        "- [ ] Primer concentration specified",
        "- [ ] Cycling parameters specified",
        "- [ ] Master mix vendor / lot",
    ]


def _recommendations(primers, validation) -> List[str]:
    recs: List[str] = []
    best = (primers.get("primers") or [{}])[0]
    if best.get("tm_diff", 0) > 2.0:
        recs.append("- Top pair has ΔTm > 2°C; trim a base to balance Tm or pick a deeper pair.")
    if best.get("amplicon_size", 0) > 200:
        recs.append("- For qPCR, redesign with amplicon ≤ 200 bp to keep efficiency in range.")
    if validation:
        if validation.get("dimers", {}).get("has_issues"):
            recs.append("- Redesign to remove problematic dimers (ΔG < −5 kcal/mol).")
        if validation.get("specificity", {}).get("is_specific") is False:
            recs.append("- In-silico PCR returned > 1 product; check off-targets.")
    if not recs:
        recs.append("- All checks passed; proceed to wet-lab validation.")
    return recs


def generate_summary_table(designs: List[Dict[str, object]]) -> str:
    """Render a markdown table comparing multiple design runs."""
    L = ["| Set | Forward | Reverse | Amplicon | F Tm | R Tm | ΔTm |",
         "|-----|---------|---------|----------|------|------|-----|"]
    for i, d in enumerate(designs, 1):
        if not d.get("primers"):
            continue
        b = d["primers"][0]
        L.append(
            f"| {i} | `{b['forward_seq'][:15]}…` | `{b['reverse_seq'][:15]}…` | "
            f"{b['amplicon_size']} bp | {b['forward_tm']:.1f} | "
            f"{b['reverse_tm']:.1f} | {b['tm_diff']:.1f} |"
        )
    return "\n".join(L)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input", required=True)
    p.add_argument("--validation")
    p.add_argument("--fmt", default="markdown", choices=["markdown", "text", "html"])
    p.add_argument("--miqe", action="store_true")
    p.add_argument("--out")
    args = p.parse_args()

    with open(args.input) as f:
        result = json.load(f)
    val = None
    if args.validation:
        with open(args.validation) as f:
            val = json.load(f)

    report = generate_primer_report(
        result, val, output_format=args.fmt, include_miqe_checklist=args.miqe,
    )
    if args.out:
        with open(args.out, "w") as f:
            f.write(report)
        print(f"✓ Report written: {args.out}")
    else:
        sys.stdout.write(report)
        print(f"\n✓ Report rendered in {args.fmt}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
