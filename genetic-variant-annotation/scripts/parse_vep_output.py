#!/usr/bin/env python3
"""Parse a VEP-annotated VCF into a tidy pandas DataFrame.

VEP writes annotations inside the INFO/CSQ field, with a pipe-separated
schema declared in the header (``##INFO=<ID=CSQ,...,Description="...|Allele|...">``).
This module reads the declared schema, splits per-transcript entries, and
returns one row per (variant × transcript) pair with canonical column names.

Public entry point: ``parse_vep_vcf(vcf_path)``.

Example:
    python scripts/parse_vep_output.py --input annotated.vep.vcf.gz \\
        --output all_variants.csv
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

try:
    import cyvcf2
    import pandas as pd
except ImportError:
    sys.exit("Install required packages with: pip install cyvcf2 pandas")


_FIELDS_RE = re.compile(r"Format:\s*(.+?)\s*\"?>?$")
_IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}


def _read_csq_schema(vcf: cyvcf2.VCF) -> list[str]:
    """Extract the CSQ field list from the VCF header."""
    header_str = vcf.raw_header
    for line in header_str.splitlines():
        if line.startswith("##INFO=<ID=CSQ"):
            match = _FIELDS_RE.search(line)
            if match:
                return [s.strip() for s in match.group(1).split("|")]
    raise ValueError("CSQ INFO field not declared in VCF header")


def _rank_impact(impact: str | None) -> int:
    return _IMPACT_RANK.get((impact or "").upper(), 0)


def parse_vep_vcf(vcf_path: str | Path) -> pd.DataFrame:
    """Return a tidy DataFrame of VEP annotations (one row per transcript)."""
    vcf_path = Path(vcf_path)
    vcf = cyvcf2.VCF(str(vcf_path))
    schema = _read_csq_schema(vcf)

    rows = []
    for record in vcf:
        csq_raw = record.INFO.get("CSQ")
        if not csq_raw:
            continue
        entries = csq_raw.split(",")
        for entry in entries:
            parts = entry.split("|")
            if len(parts) < len(schema):
                parts += [""] * (len(schema) - len(parts))
            row = dict(zip(schema, parts[:len(schema)]))
            row["chrom"] = record.CHROM
            row["pos"]   = record.POS
            row["ref"]   = record.REF
            row["alt"]   = record.ALT[0] if record.ALT else ""
            row["id"]    = record.ID
            row["qual"]  = record.QUAL
            row["filter"] = record.FILTER
            rows.append(row)
    vcf.close()

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    # Canonical column aliases used downstream.
    rename = {"Allele": "allele", "Consequence": "consequence",
              "IMPACT": "impact", "SYMBOL": "gene_symbol",
              "Gene": "gene_id", "Feature": "transcript_id",
              "BIOTYPE": "biotype", "Existing_variation": "existing_variation",
              "MAX_AF": "max_af", "gnomAD_AF": "gnomad_af",
              "AF": "af", "ClinVar_CLNSIG": "clinvar_sig",
              "CLIN_SIG": "clinvar_sig",
              "CADD_PHRED": "cadd_phred", "REVEL": "revel",
              "SIFT": "sift", "PolyPhen": "polyphen",
              "HGVSc": "hgvsc", "HGVSp": "hgvsp",
              "EXON": "exon", "INTRON": "intron"}
    df = df.rename(columns={k: v for k, v in rename.items() if k in df.columns})
    if "impact" in df.columns:
        df["impact_rank"] = df["impact"].map(_rank_impact)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    print(f"=== Parsing VEP output: {args.input} ===")
    df = parse_vep_vcf(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"✓ VEP parsing completed: {len(df):,} transcript annotations / "
          f"{args.output}")


if __name__ == "__main__":
    main()
