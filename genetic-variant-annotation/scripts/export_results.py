#!/usr/bin/env python3
"""Persist all variant-annotation outputs (CSVs, filtered VCF, Excel, pickle).

Public function: ``export_all(df, gene_df, original_vcf, output_dir, tool_name)``.

The script writes:
  - all_variants.csv, high_impact_variants.csv, moderate_impact_variants.csv,
    rare_variants.csv
  - filtered_high_impact.vcf.gz
  - gene_summary.csv
  - summary_report.xlsx (multi-sheet)
  - analysis_object.pkl

Example:
    python scripts/export_results.py \\
        --variants all_variants.csv --genes gene_summary.csv \\
        --original-vcf input.vcf.gz --output-dir results/ --tool snpeff
"""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path
from typing import Any

try:
    import pandas as pd
except ImportError:
    sys.exit("Install required packages with: pip install pandas openpyxl pysam")


def _filter_subset(df: pd.DataFrame, predicate) -> pd.DataFrame:
    if df.empty:
        return df
    return df.loc[predicate].copy()


def _write_filtered_vcf(
    original_vcf: Path, kept_variants: pd.DataFrame, output_vcf: Path,
) -> None:
    """Subset the original VCF to exact variant records present in ``kept_variants``."""
    import cyvcf2
    output_vcf.parent.mkdir(parents=True, exist_ok=True)
    key_cols = {"chrom", "pos", "ref", "alt"}
    if key_cols.issubset(kept_variants.columns):
        keep_keys = set(
            zip(
                kept_variants["chrom"].astype(str),
                kept_variants["pos"].astype(int),
                kept_variants["ref"].astype(str),
                kept_variants["alt"].astype(str),
            )
        )
        use_full_key = True
    else:
        keep_keys = set(
            zip(kept_variants["chrom"].astype(str), kept_variants["pos"].astype(int))
        )
        use_full_key = False

    vcf = cyvcf2.VCF(str(original_vcf))
    writer = cyvcf2.Writer(str(output_vcf), vcf)
    for record in vcf:
        if use_full_key:
            alts = record.ALT or []
            record_keys = {(record.CHROM, record.POS, record.REF, alt) for alt in alts}
            should_keep = bool(record_keys & keep_keys)
        else:
            should_keep = (record.CHROM, record.POS) in keep_keys
        if should_keep:
            writer.write_record(record)
    writer.close(); vcf.close()


def _write_excel_workbook(
    output: Path, *,
    all_variants: pd.DataFrame, gene_summary: pd.DataFrame,
    high_impact: pd.DataFrame, rare: pd.DataFrame, tool_name: str,
) -> None:
    """Multi-sheet workbook with summary + slices."""
    output.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(output, engine="openpyxl") as xw:
        # Summary sheet
        summary_rows = [
            ("annotation_tool",   tool_name),
            ("total variants",    int(len(all_variants))),
            ("HIGH impact",       int((all_variants["impact"]
                                       == "HIGH").sum()) if "impact" in all_variants.columns else 0),
            ("MODERATE impact",   int((all_variants["impact"]
                                       == "MODERATE").sum()) if "impact" in all_variants.columns else 0),
            ("rare (AF < 0.01)",  int(len(rare))),
            ("unique genes",      int(gene_summary["gene_symbol"].nunique())
                                     if not gene_summary.empty else 0),
        ]
        pd.DataFrame(summary_rows, columns=["metric", "value"]).to_excel(
            xw, sheet_name="Summary", index=False)

        all_variants.head(500_000).to_excel(xw, sheet_name="Variants",
                                            index=False)
        if not high_impact.empty:
            high_impact.to_excel(xw, sheet_name="High_Impact", index=False)
        if not rare.empty:
            rare.head(200_000).to_excel(xw, sheet_name="Rare_Variants",
                                        index=False)
        if not gene_summary.empty:
            gene_summary.to_excel(xw, sheet_name="Gene_Summary", index=False)


def export_all(
    df: pd.DataFrame,
    gene_df: pd.DataFrame,
    original_vcf: str | Path,
    output_dir: str | Path,
    tool_name: str = "snpeff",
) -> dict[str, Path]:
    """Persist every artefact. Returns a mapping of output role -> path."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    paths: dict[str, Path] = {}

    # All variants
    p = output_dir / "all_variants.csv"
    df.to_csv(p, index=False); paths["all_variants"] = p
    print(f"  -> {p}")

    if "impact" in df.columns:
        impact_upper = df["impact"].astype(str).str.upper()
        high = _filter_subset(df, impact_upper == "HIGH")
        moderate = _filter_subset(df, impact_upper.isin({"HIGH", "MODERATE"}))
    else:
        high = moderate = df.iloc[0:0]

    p = output_dir / "high_impact_variants.csv"
    high.to_csv(p, index=False); paths["high_impact"] = p
    print(f"  -> {p}")

    p = output_dir / "moderate_impact_variants.csv"
    moderate.to_csv(p, index=False); paths["moderate_impact"] = p
    print(f"  -> {p}")

    if "gnomad_af" in df.columns or "af" in df.columns:
        col = "gnomad_af" if "gnomad_af" in df.columns else "af"
        af = pd.to_numeric(df[col], errors="coerce")
        rare = df.loc[af.notna() & (af < 0.01)].copy()
        unknown_af = df.loc[af.isna()].copy()
    else:
        rare = df.iloc[0:0]
        unknown_af = df.iloc[0:0]
    p = output_dir / "rare_variants.csv"
    rare.to_csv(p, index=False); paths["rare"] = p
    print(f"  -> {p}")

    if not unknown_af.empty:
        p = output_dir / "unknown_frequency_variants.csv"
        unknown_af.to_csv(p, index=False); paths["unknown_frequency"] = p
        print(f"  -> {p}")

    # Gene summary
    if gene_df is not None and not gene_df.empty:
        p = output_dir / "gene_summary.csv"
        gene_df.to_csv(p, index=False); paths["gene_summary"] = p
        print(f"  -> {p}")

    # Filtered VCF (HIGH impact only)
    if not high.empty:
        p = output_dir / "filtered_high_impact.vcf.gz"
        try:
            _write_filtered_vcf(Path(original_vcf), high, p)
            paths["filtered_vcf"] = p
            print(f"  -> {p}")
        except Exception as exc:
            print(f"  ⚠ filtered VCF skipped ({exc})")

    # Excel workbook
    p = output_dir / "summary_report.xlsx"
    _write_excel_workbook(p, all_variants=df, gene_summary=gene_df,
                          high_impact=high, rare=rare, tool_name=tool_name)
    paths["excel"] = p
    print(f"  -> {p}")

    # Pickled analysis object
    p = output_dir / "analysis_object.pkl"
    payload: dict[str, Any] = {
        "tool":     tool_name,
        "variants": df, "gene_summary": gene_df,
        "high_impact": high, "rare": rare,
        "n_variants": len(df), "n_genes": (
            int(gene_df["gene_symbol"].nunique()) if gene_df is not None
            and not gene_df.empty else 0),
    }
    with open(p, "wb") as fh:
        pickle.dump(payload, fh)
    paths["pickle"] = p
    print(f"  -> {p}")

    print("\n=== Export Complete ===")
    return paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--variants",       type=Path, required=True)
    parser.add_argument("--genes",          type=Path, default=None)
    parser.add_argument("--original-vcf",   type=Path, required=True)
    parser.add_argument("--output-dir",     type=Path, required=True)
    parser.add_argument("--tool", default="snpeff",
                        choices=["vep", "snpeff"])
    args = parser.parse_args()

    print("=== Exporting variant-annotation results ===")
    df = pd.read_csv(args.variants)
    gene_df = pd.read_csv(args.genes) if args.genes and args.genes.exists() \
              else pd.DataFrame()
    export_all(df, gene_df, args.original_vcf, args.output_dir,
               tool_name=args.tool)
    print(f"✓ Export completed: {args.output_dir}")


if __name__ == "__main__":
    main()
