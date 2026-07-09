"""Export primer-design results to CSV, Excel, JSON, IDT order sheet, or MIQE.

Single entry point `export_primers(result, fmt=..., output_file=...)` dispatches
to format-specific writers. Excel and MIQE checklist require pandas + openpyxl.

Inputs:
    A dict produced by design_* functions plus optional validation_results.
Outputs:
    File on disk in the requested format.

Example invocation
------------------
    python -m scripts.export_results --input primers.json --fmt csv --out primers.csv
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import sys
from datetime import datetime
from typing import Dict, List, Optional

VALID_FORMATS = {"csv", "excel", "json", "idt_order", "miqe_checklist", "benchling"}


def export_primers(
    result: Dict[str, object],
    fmt: str = "csv",
    output_file: Optional[str] = None,
    include_validation: bool = True,
    validation_results: Optional[Dict[str, object]] = None,
) -> str:
    """Dispatch to the requested writer.

    Returns the output path when `output_file` is given, otherwise the
    serialized payload as a string.
    """
    fmt = fmt.lower()
    if fmt not in VALID_FORMATS:
        raise ValueError(f"Unknown format {fmt!r}. Choose from {sorted(VALID_FORMATS)}.")
    writers = {
        "csv": _write_csv,
        "excel": _write_excel,
        "json": _write_json,
        "idt_order": _write_idt,
        "miqe_checklist": _write_miqe,
        "benchling": _write_benchling,
    }
    return writers[fmt](result, output_file, include_validation, validation_results)


def _flat_rows(result: Dict[str, object]) -> List[Dict[str, object]]:
    rows = []
    for pair in result.get("primers", []):
        rows.append(
            {
                "Pair_ID": pair.get("pair_id"),
                "Forward_Sequence": pair["forward_seq"],
                "Reverse_Sequence": pair["reverse_seq"],
                "Amplicon_bp": pair["amplicon_size"],
                "Forward_Tm_C": pair["forward_tm"],
                "Reverse_Tm_C": pair["reverse_tm"],
                "Tm_Diff_C": pair["tm_diff"],
                "Forward_GC_pct": pair["forward_gc"],
                "Reverse_GC_pct": pair["reverse_gc"],
                "Forward_Length": pair["forward_length"],
                "Reverse_Length": pair["reverse_length"],
                "Forward_Position": pair["forward_pos"],
                "Reverse_Position": pair["reverse_pos"],
                "Primer3_Penalty": pair["penalty"],
                "MIQE_Compliant": pair.get("miqe_compliant"),
            }
        )
    return rows


def _write_csv(result, output_file, include_validation, validation_results):
    rows = _flat_rows(result)
    if include_validation and validation_results:
        for r in rows:
            r["Specific"] = validation_results.get("specificity", {}).get("is_specific")
            r["Has_Dimers"] = validation_results.get("dimers", {}).get("has_issues")
            r["Has_Structures"] = validation_results.get("secondary_structures", {}).get(
                "has_issues"
            )
    cols = list(rows[0].keys()) if rows else []
    if output_file:
        with open(output_file, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols)
            w.writeheader()
            w.writerows(rows)
        return output_file
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue()


def _write_excel(result, output_file, include_validation, validation_results):
    try:
        import pandas as pd
    except ImportError:
        sys.exit("Install with: pip install pandas openpyxl")
    if not output_file:
        raise ValueError("Excel export requires output_file.")

    primers_df = pd.DataFrame(_flat_rows(result))
    params_df = pd.DataFrame([result.get("parameters", {})])

    with pd.ExcelWriter(output_file, engine="openpyxl") as xw:
        primers_df.to_excel(xw, sheet_name="Primers", index=False)
        params_df.to_excel(xw, sheet_name="Parameters", index=False)
        if include_validation and validation_results:
            val_rows = []
            for category, payload in validation_results.items():
                if isinstance(payload, dict):
                    for key, value in payload.items():
                        val_rows.append(
                            {"Category": category, "Check": key, "Result": str(value)}
                        )
            if val_rows:
                pd.DataFrame(val_rows).to_excel(xw, sheet_name="Validation", index=False)
    return output_file


def _write_json(result, output_file, include_validation, validation_results):
    payload = {
        "metadata": {
            "generated": datetime.now().isoformat(),
            "sequence_length": result.get("sequence_length"),
            "num_primers_found": result.get("num_primers_found"),
        },
        "parameters": result.get("parameters", {}),
        "primers": result.get("primers", []),
    }
    if include_validation and validation_results:
        payload["validation"] = validation_results
    text = json.dumps(payload, indent=2, default=str)
    if output_file:
        with open(output_file, "w") as f:
            f.write(text)
        return output_file
    return text


def _write_idt(result, output_file, include_validation, validation_results):
    lines = ["Name\tSequence\tScale\tPurification"]
    for i, pair in enumerate(result.get("primers", []), 1):
        lines.append(f"Pair{i}_Forward\t{pair['forward_seq']}\t25nm\tSTD")
        lines.append(f"Pair{i}_Reverse\t{pair['reverse_seq']}\t25nm\tSTD")
        if "probe_seq" in pair:
            lines.append(f"Pair{i}_Probe\t{pair['probe_seq']}\t25nm\tHPLC")
    text = "\n".join(lines) + "\n"
    if output_file:
        with open(output_file, "w") as f:
            f.write(text)
        return output_file
    return text


def _write_miqe(result, output_file, include_validation, validation_results):
    try:
        import pandas as pd
    except ImportError:
        sys.exit("Install with: pip install pandas openpyxl")
    if not output_file:
        raise ValueError("MIQE checklist export requires output_file.")

    first = result.get("primers", [{}])[0] if result.get("primers") else {}
    items = [
        ("Experimental Design", "qPCR purpose", "", "Essential"),
        ("Experimental Design", "Experimental design", "", "Essential"),
        ("Experimental Design", "Sample description", "", "Essential"),
        ("Sample", "RNA/DNA quality (RIN, 260/280)", "", "Essential"),
        ("Sample", "Storage conditions", "", "Essential"),
        ("Sample", "RT enzyme and priming (RT-qPCR)", "", "Essential"),
        ("Primers", "Forward sequence", first.get("forward_seq", ""), "Essential"),
        ("Primers", "Reverse sequence", first.get("reverse_seq", ""), "Essential"),
        (
            "Primers",
            "Target positions",
            f"{first.get('forward_pos', 'NA')}-{first.get('reverse_pos', 'NA')}",
            "Essential",
        ),
        ("Primers", "Amplicon length", f"{first.get('amplicon_size', 'NA')} bp", "Essential"),
        ("Primers", "In-silico specificity check", "Performed", "Essential"),
        ("Primers", "Empirical specificity check", "", "Essential"),
        ("Protocol", "Reaction volume", "", "Essential"),
        ("Protocol", "Primer concentration", "", "Essential"),
        ("Protocol", "Cycling parameters", "", "Essential"),
        ("Validation", "Specificity (melt curve)", "", "Essential"),
        ("Validation", "Efficiency E", "", "Essential"),
        ("Validation", "R-squared", "", "Essential"),
        ("Validation", "Linear dynamic range", "", "Essential"),
        ("Validation", "LOD", "", "Desirable"),
        ("Data Analysis", "Cq method", "", "Essential"),
        ("Data Analysis", "Outlier handling", "", "Essential"),
        ("Data Analysis", "NTC results", "", "Essential"),
        ("Data Analysis", "Reference gene validation", "", "Essential"),
        ("Data Analysis", "Statistical methods", "", "Essential"),
    ]
    df_check = pd.DataFrame(items, columns=["Category", "Item", "Value", "MIQE Level"])

    primer_rows = []
    for i, pair in enumerate(result.get("primers", []), 1):
        primer_rows.append(
            {
                "Pair": i,
                "Forward": pair["forward_seq"],
                "Reverse": pair["reverse_seq"],
                "Amplicon_bp": pair["amplicon_size"],
                "Tm_F": pair["forward_tm"],
                "Tm_R": pair["reverse_tm"],
            }
        )
    df_primers = pd.DataFrame(primer_rows)

    with pd.ExcelWriter(output_file, engine="openpyxl") as xw:
        df_check.to_excel(xw, sheet_name="MIQE_Checklist", index=False)
        df_primers.to_excel(xw, sheet_name="Primers", index=False)
    return output_file


def _write_benchling(result, output_file, include_validation, validation_results):
    rows = [("Name", "Bases", "Description")]
    for i, pair in enumerate(result.get("primers", []), 1):
        rows.append(
            (
                f"Pair{i}_Forward",
                pair["forward_seq"],
                f"Forward Tm={pair['forward_tm']}C GC={pair['forward_gc']}%",
            )
        )
        rows.append(
            (
                f"Pair{i}_Reverse",
                pair["reverse_seq"],
                f"Reverse Tm={pair['reverse_tm']}C GC={pair['reverse_gc']}%",
            )
        )
    if output_file:
        with open(output_file, "w", newline="") as f:
            csv.writer(f).writerows(rows)
        return output_file
    buf = io.StringIO()
    csv.writer(buf).writerows(rows)
    return buf.getvalue()


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--input", required=True, help="JSON file from a design_* run")
    p.add_argument("--fmt", required=True, choices=sorted(VALID_FORMATS))
    p.add_argument("--out", required=True)
    p.add_argument("--validation", help="Optional JSON file with validation_results")
    args = p.parse_args()

    with open(args.input) as f:
        result = json.load(f)
    val = None
    if args.validation:
        with open(args.validation) as f:
            val = json.load(f)

    path = export_primers(result, fmt=args.fmt, output_file=args.out,
                          include_validation=val is not None,
                          validation_results=val)
    print(f"✓ Export completed: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
