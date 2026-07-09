#!/usr/bin/env python3
"""Parse a SNPEff-annotated VCF into a tidy pandas DataFrame.

SNPEff writes its modern annotation schema in INFO/ANN with this field order:
  Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID |
  Feature_Type | Feature_ID | Transcript_BioType | Rank/Total | HGVS.c |
  HGVS.p | cDNA_pos/cDNA_len | CDS_pos/CDS_len | AA_pos/AA_len | Distance |
  ERRORS / WARNINGS / INFO

This module reads each ANN entry and returns one row per (variant × transcript).

Public entry point: ``parse_snpeff_vcf(vcf_path)``.

Example:
    python scripts/parse_snpeff_output.py --input annotated.snpeff.vcf.gz \\
        --output all_variants.csv
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import cyvcf2
    import pandas as pd
except ImportError:
    sys.exit("Install required packages with: pip install cyvcf2 pandas")


_ANN_SCHEMA = [
    "allele", "consequence", "impact", "gene_symbol", "gene_id",
    "feature_type", "transcript_id", "biotype", "exon_rank",
    "hgvsc", "hgvsp", "cdna_position", "cds_position",
    "aa_position", "distance", "errors",
]

_IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}


def _rank_impact(impact: str | None) -> int:
    return _IMPACT_RANK.get((impact or "").upper(), 0)


def parse_snpeff_vcf(vcf_path: str | Path) -> pd.DataFrame:
    """Return a tidy DataFrame from a SNPEff-annotated VCF."""
    vcf_path = Path(vcf_path)
    vcf = cyvcf2.VCF(str(vcf_path))
    rows = []
    for record in vcf:
        ann_raw = record.INFO.get("ANN")
        if not ann_raw:
            continue
        entries = ann_raw.split(",")
        for entry in entries:
            parts = entry.split("|")
            if len(parts) < len(_ANN_SCHEMA):
                parts += [""] * (len(_ANN_SCHEMA) - len(parts))
            row = dict(zip(_ANN_SCHEMA, parts[:len(_ANN_SCHEMA)]))
            row["chrom"] = record.CHROM
            row["pos"]   = record.POS
            row["ref"]   = record.REF
            row["alt"]   = record.ALT[0] if record.ALT else ""
            row["id"]    = record.ID
            row["qual"]  = record.QUAL
            row["filter"] = record.FILTER

            # Optional auxiliary INFO fields commonly carried by GATK / bcftools.
            for info_field, col in (("AF", "af"), ("AC", "ac"), ("AN", "an"),
                                    ("DP", "dp"),
                                    ("CLNSIG", "clinvar_sig"),
                                    ("CADD_PHRED", "cadd_phred"),
                                    ("REVEL", "revel"),
                                    ("LOF", "lof_flag"),
                                    ("NMD", "nmd_flag")):
                if info_field in record.INFO:
                    row[col] = record.INFO[info_field]
            rows.append(row)
    vcf.close()

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
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

    print(f"=== Parsing SNPEff output: {args.input} ===")
    df = parse_snpeff_vcf(args.input)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(f"✓ SNPEff parsing completed: {len(df):,} annotations / {args.output}")


if __name__ == "__main__":
    main()
