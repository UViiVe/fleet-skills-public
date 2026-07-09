#!/usr/bin/env python3
"""Run SNPEff against a VCF and return the annotated VCF path.

Public entry points:
  - ``check_snpeff_installation()`` returns ``(is_installed, snpeff_path, version)``.
  - ``download_snpeff_database(genome)``
  - ``run_snpeff(input_vcf, output_vcf, genome, ...)``

Example:
    python scripts/run_snpeff.py --input input.vcf.gz \\
        --output annotated.snpeff.vcf.gz --genome GRCh38.105
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def check_snpeff_installation() -> tuple[bool, str | None, str | None]:
    """Return (is_installed, path, version_string)."""
    path = shutil.which("snpEff")
    if not path:
        return False, None, None
    try:
        out = subprocess.run(["snpEff", "-version"], capture_output=True,
                             text=True, timeout=15)
        version = (out.stdout or out.stderr).strip().splitlines()[0]
        return True, path, version
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
        return False, path, None


def download_snpeff_database(genome: str) -> bool:
    """Download a SNPEff genome database."""
    cmd = ["snpEff", "download", "-v", genome]
    print(f"  $ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"✓ SNPEff database downloaded: {genome}")
    return True


def run_snpeff(
    input_vcf: str | Path,
    output_vcf: str | Path,
    *,
    genome: str = "GRCh38.105",
    canon: bool = True,
    hgvs: bool = True,
    no_intergenic: bool = False,
    no_downstream: bool = False,
    no_upstream: bool = False,
    lof: bool = False,
    threads: int = 4,
    stats_html: str | Path | None = None,
    csv_stats: str | Path | None = None,
    extra_args: list[str] | None = None,
) -> Path:
    """Run SNPEff and write a gzipped annotated VCF."""
    input_vcf  = Path(input_vcf)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["snpEff", "-t", str(threads)]
    if canon:          cmd.append("-canon")
    if hgvs:           cmd.append("-hgvs")
    if no_intergenic:  cmd.append("-noIntergenic")
    if no_downstream:  cmd.append("-no-downstream")
    if no_upstream:    cmd.append("-no-upstream")
    if lof:            cmd.append("-lof")
    if stats_html:     cmd += ["-stats", str(stats_html)]
    if csv_stats:      cmd += ["-csvStats", str(csv_stats)]
    if extra_args:     cmd += extra_args

    cmd += [genome, str(input_vcf)]
    print(f"  $ {' '.join(cmd)}  | gzip > {output_vcf}")

    try:
        with gzip.open(output_vcf, "wb") as out_fh:
            proc = subprocess.run(cmd, check=True, stdout=out_fh)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "snpEff not found. Install with: "
            "python scripts/install_tools.py --tool snpeff") from exc
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(f"SNPEff failed: {exc}\n")
        raise

    print(f"✓ SNPEff annotation completed: {output_vcf}")
    return output_vcf


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genome", default="GRCh38.105")
    parser.add_argument("--no-canon", action="store_true",
                        help="Disable -canon (use all transcripts)")
    parser.add_argument("--no-hgvs", action="store_true")
    parser.add_argument("--no-intergenic", action="store_true")
    parser.add_argument("--no-downstream", action="store_true")
    parser.add_argument("--no-upstream", action="store_true")
    parser.add_argument("--lof", action="store_true")
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--stats-html", type=Path, default=None)
    parser.add_argument("--csv-stats", type=Path, default=None)
    args = parser.parse_args()

    installed, path, version = check_snpeff_installation()
    if not installed:
        sys.exit("SNPEff not installed. Run: "
                 "python scripts/install_tools.py --tool snpeff")
    print(f"  SNPEff: {path} ({version})")

    run_snpeff(
        args.input, args.output,
        genome=args.genome,
        canon=not args.no_canon,
        hgvs=not args.no_hgvs,
        no_intergenic=args.no_intergenic,
        no_downstream=args.no_downstream,
        no_upstream=args.no_upstream,
        lof=args.lof,
        threads=args.threads,
        stats_html=args.stats_html,
        csv_stats=args.csv_stats,
    )


if __name__ == "__main__":
    main()
