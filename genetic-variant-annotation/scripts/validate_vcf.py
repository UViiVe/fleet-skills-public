#!/usr/bin/env python3
"""Validate a VCF file and report summary statistics.

Two public entry points:
  - ``validate_vcf(vcf_path, reference_genome=None)`` returns a dict with
    ``is_valid``, ``errors``, and ``warnings``.
  - ``vcf_summary_stats(vcf_path)`` returns counts of total variants, SNVs,
    indels, transitions/transversions, multi-allelic sites, and per-chromosome
    breakdown.

Example:
    python scripts/validate_vcf.py --vcf input.vcf.gz
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    import cyvcf2
except ImportError:
    sys.exit("Install required packages with: pip install cyvcf2 pysam pandas numpy")


_TRANSITIONS = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}


def _is_transition(ref: str, alt: str) -> bool:
    return (ref.upper(), alt.upper()) in _TRANSITIONS


def _validate_header(vcf: cyvcf2.VCF) -> dict[str, Any]:
    out: dict[str, Any] = {"is_valid": True, "errors": [], "warnings": []}
    if not vcf.seqnames:
        out["warnings"].append("no contig definitions in header (##contig lines)")
    if vcf.samples and "GT" not in vcf:
        out["warnings"].append("GT (genotype) FORMAT field not declared in header")
    for info_field in ("AC", "AF", "AN", "DP"):
        if info_field not in vcf:
            out["warnings"].append(f"INFO/{info_field} not declared in header")
    return out


def _validate_variants(
    vcf: cyvcf2.VCF, *, max_variants: int, reference_genome: Path | None,
) -> dict[str, Any]:
    out: dict[str, Any] = {"is_valid": True, "errors": [], "warnings": []}
    ref_fasta = None
    if reference_genome is not None:
        try:
            import pysam
            ref_fasta = pysam.FastaFile(str(reference_genome))
        except Exception as exc:
            out["warnings"].append(f"could not open reference FASTA: {exc}")

    seen_chroms: set[str] = set()
    inspected = 0
    for record in vcf:
        if inspected >= max_variants:
            break
        inspected += 1
        seen_chroms.add(record.CHROM)

        if not record.REF or set(record.REF.upper()) - set("ACGTN"):
            out["warnings"].append(
                f"non-ACGTN REF at {record.CHROM}:{record.POS}: {record.REF}")
        for alt in record.ALT:
            if alt and set(alt.upper()) - set("ACGTN,<>"):
                out["warnings"].append(
                    f"non-ACGTN ALT at {record.CHROM}:{record.POS}: {alt}")

        if ref_fasta is not None:
            try:
                ref_seq = ref_fasta.fetch(
                    record.CHROM, record.POS - 1,
                    record.POS - 1 + len(record.REF)).upper()
                if ref_seq != record.REF.upper():
                    out["errors"].append(
                        f"REF mismatch at {record.CHROM}:{record.POS}: "
                        f"VCF={record.REF}, FASTA={ref_seq}")
                    out["is_valid"] = False
            except Exception:
                pass  # contig missing — caller handles separately

    out["chromosomes_seen"] = sorted(seen_chroms)
    out["inspected_variants"] = inspected
    return out


def validate_vcf(
    vcf_path: str | Path, reference_genome: str | Path | None = None,
    *, max_variants: int = 10_000,
) -> dict[str, Any]:
    """Validate header + variants of a VCF file."""
    out: dict[str, Any] = {"is_valid": True, "errors": [], "warnings": []}
    vcf_path = Path(vcf_path)
    if not vcf_path.exists():
        out["is_valid"] = False
        out["errors"].append(f"file not found: {vcf_path}")
        return out
    if vcf_path.stat().st_size == 0:
        out["is_valid"] = False
        out["errors"].append("file is empty")
        return out

    try:
        vcf = cyvcf2.VCF(str(vcf_path))
    except Exception as exc:
        out["is_valid"] = False
        out["errors"].append(f"cannot parse VCF: {exc}")
        return out

    header = _validate_header(vcf)
    out["warnings"].extend(header["warnings"])
    if not header["is_valid"]:
        out["is_valid"] = False
        out["errors"].extend(header["errors"])

    body = _validate_variants(vcf, max_variants=max_variants,
                              reference_genome=Path(reference_genome)
                              if reference_genome else None)
    out["warnings"].extend(body["warnings"])
    out["errors"].extend(body["errors"])
    out["inspected_variants"] = body["inspected_variants"]
    out["chromosomes_seen"]   = body["chromosomes_seen"]
    if body["errors"]:
        out["is_valid"] = False

    vcf.close()
    return out


def vcf_summary_stats(vcf_path: str | Path) -> dict[str, Any]:
    """Counts of variants, SNVs, indels, transitions/transversions, etc."""
    vcf_path = Path(vcf_path)
    vcf = cyvcf2.VCF(str(vcf_path))
    counts: dict[str, Any] = {
        "total_variants": 0, "snvs": 0, "indels": 0,
        "multi_allelic": 0, "transitions": 0, "transversions": 0,
        "per_chrom": defaultdict(int),
    }
    try:
        for record in vcf:
            counts["total_variants"] += 1
            counts["per_chrom"][record.CHROM] += 1
            if len(record.ALT) > 1:
                counts["multi_allelic"] += 1
            for alt in record.ALT:
                if alt is None:
                    continue
                if len(record.REF) == 1 and len(alt) == 1:
                    counts["snvs"] += 1
                    if _is_transition(record.REF, alt):
                        counts["transitions"] += 1
                    else:
                        counts["transversions"] += 1
                else:
                    counts["indels"] += 1
    finally:
        vcf.close()

    ti, tv = counts["transitions"], counts["transversions"]
    counts["ti_tv_ratio"] = round(ti / tv, 3) if tv else float("nan")
    counts["per_chrom"] = dict(counts["per_chrom"])
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--vcf", type=Path, required=True)
    parser.add_argument("--reference", type=Path, default=None)
    parser.add_argument("--max-variants", type=int, default=10_000)
    args = parser.parse_args()

    print(f"=== Validating VCF: {args.vcf} ===")
    result = validate_vcf(args.vcf, args.reference,
                          max_variants=args.max_variants)
    for w in result["warnings"]:
        print(f"  WARN: {w}")
    for e in result["errors"]:
        print(f"  ERROR: {e}")
    print(f"  is_valid: {result['is_valid']}")

    stats = vcf_summary_stats(args.vcf)
    print(f"  total variants: {stats['total_variants']:,}")
    print(f"  SNVs / indels : {stats['snvs']:,} / {stats['indels']:,}")
    print(f"  Ti/Tv ratio   : {stats['ti_tv_ratio']}")
    print(f"  multi-allelic : {stats['multi_allelic']:,}")
    print(f"✓ VCF summary completed: {stats['total_variants']:,} variants")


if __name__ == "__main__":
    main()
