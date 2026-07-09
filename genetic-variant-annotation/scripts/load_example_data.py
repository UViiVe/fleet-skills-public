#!/usr/bin/env python3
"""Generate a tiny example VCF with known ClinVar-pathogenic variants in BRCA1/BRCA2.

Used by the quick-start path and for smoke-testing installations. The VCF
contains ~10 SNVs spread across BRCA1 and BRCA2 with realistic ClinVar SIG
classifications. The genome coordinates use GRCh38.

Example:
    python scripts/load_example_data.py --output-dir example_data
"""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Any


_VCF_HEADER = """##fileformat=VCFv4.2
##source=load_example_data.py
##reference=GRCh38
##contig=<ID=chr13>
##contig=<ID=chr17>
##INFO=<ID=GENE,Number=1,Type=String,Description="Gene symbol">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="ClinVar significance">
##INFO=<ID=AF,Number=1,Type=Float,Description="Population allele frequency (gnomAD)">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


# A curated mini-set of ClinVar pathogenic / likely-pathogenic variants in
# BRCA1 (chr17) and BRCA2 (chr13), GRCh38.
_VARIANTS: list[dict[str, Any]] = [
    # BRCA1 (chr17, NM_007294.4)
    {"chrom": "chr17", "pos": 43093919, "id": "rs80357906",
     "ref": "G", "alt": "A", "gene": "BRCA1",
     "clnsig": "Pathogenic", "af": 0.00002},
    {"chrom": "chr17", "pos": 43094426, "id": "rs80357713",
     "ref": "G", "alt": "T", "gene": "BRCA1",
     "clnsig": "Pathogenic", "af": 0.00001},
    {"chrom": "chr17", "pos": 43095919, "id": "rs80357711",
     "ref": "C", "alt": "T", "gene": "BRCA1",
     "clnsig": "Pathogenic", "af": 0.00002},
    {"chrom": "chr17", "pos": 43097243, "id": "rs28897696",
     "ref": "A", "alt": "G", "gene": "BRCA1",
     "clnsig": "Likely_pathogenic", "af": 0.00001},
    {"chrom": "chr17", "pos": 43124027, "id": "rs80357709",
     "ref": "C", "alt": "G", "gene": "BRCA1",
     "clnsig": "Pathogenic", "af": 0.00001},
    # BRCA2 (chr13, NM_000059.4)
    {"chrom": "chr13", "pos": 32316461, "id": "rs80359550",
     "ref": "G", "alt": "T", "gene": "BRCA2",
     "clnsig": "Pathogenic", "af": 0.00002},
    {"chrom": "chr13", "pos": 32319076, "id": "rs80358971",
     "ref": "A", "alt": "G", "gene": "BRCA2",
     "clnsig": "Likely_pathogenic", "af": 0.00001},
    {"chrom": "chr13", "pos": 32319124, "id": "rs80359077",
     "ref": "C", "alt": "T", "gene": "BRCA2",
     "clnsig": "Pathogenic", "af": 0.00001},
    {"chrom": "chr13", "pos": 32340301, "id": "rs80359501",
     "ref": "T", "alt": "A", "gene": "BRCA2",
     "clnsig": "Pathogenic", "af": 0.00002},
    {"chrom": "chr13", "pos": 32340793, "id": "rs80359152",
     "ref": "G", "alt": "C", "gene": "BRCA2",
     "clnsig": "Pathogenic", "af": 0.00001},
]


def _format_record(v: dict[str, Any]) -> str:
    info = f"GENE={v['gene']};CLNSIG={v['clnsig']};AF={v['af']:.5g}"
    return (f"{v['chrom']}\t{v['pos']}\t{v['id']}\t{v['ref']}\t{v['alt']}"
            f"\t99.0\tPASS\t{info}")


def load_clinvar_pathogenic_sample(
    output_dir: str | Path = "example_data",
) -> dict[str, Any]:
    """Write the demo VCF and return a metadata dict.

    Returns
    -------
    dict with keys: ``vcf_path``, ``n_variants``, ``genes``, ``description``.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    vcf_path = out_dir / "clinvar_brca_demo.vcf.gz"

    body = "\n".join(_format_record(v) for v in _VARIANTS) + "\n"
    content = _VCF_HEADER + body

    with gzip.open(vcf_path, "wt") as fh:
        fh.write(content)
    print(f"  wrote {len(_VARIANTS)} variants to {vcf_path}")

    return {
        "vcf_path": str(vcf_path),
        "n_variants": len(_VARIANTS),
        "genes": sorted({v["gene"] for v in _VARIANTS}),
        "description": "BRCA1/BRCA2 ClinVar pathogenic demo (GRCh38)",
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--output-dir", type=Path, default="example_data")
    args = parser.parse_args()

    info = load_clinvar_pathogenic_sample(args.output_dir)
    print(f"✓ Demo VCF written: {info['vcf_path']} "
          f"({info['n_variants']} variants in {', '.join(info['genes'])})")


if __name__ == "__main__":
    main()
