#!/usr/bin/env python3
"""Run Ensembl VEP against a VCF and return the annotated VCF path.

Public entry points:
  - ``check_vep_installation()`` returns ``(is_installed, vep_path, version)``.
  - ``install_vep_cache(species, assembly, cache_dir)``
  - ``run_vep(input_vcf, output_vcf, genome, ...)``

Example:
    python scripts/run_vep.py --input input.vcf.gz --output annotated.vep.vcf.gz \\
        --genome GRCh38 --fork 4 --everything
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def check_vep_installation() -> tuple[bool, str | None, str | None]:
    """Return (is_installed, path, version_string)."""
    path = shutil.which("vep")
    if not path:
        return False, None, None
    try:
        out = subprocess.run(["vep", "--help"], capture_output=True,
                             text=True, timeout=15)
        version = None
        for line in out.stdout.splitlines():
            if "ensembl-vep" in line.lower():
                version = line.strip()
                break
        return True, path, version
    except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
        return False, path, None


def install_vep_cache(
    species: str = "homo_sapiens",
    assembly: str = "GRCh38",
    cache_dir: str | Path | None = None,
) -> str:
    """Download the VEP cache via ``vep_install``."""
    cache_dir = Path(cache_dir) if cache_dir else Path.home() / ".vep"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["vep_install", "--AUTO", "c",
           "--SPECIES", species, "--ASSEMBLY", assembly,
           "--CACHEDIR", str(cache_dir), "--NO_UPDATE"]
    print(f"  $ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    print(f"✓ VEP cache installed to {cache_dir}")
    return str(cache_dir)


def run_vep(
    input_vcf: str | Path,
    output_vcf: str | Path,
    *,
    genome: str = "GRCh38",
    species: str = "homo_sapiens",
    cache_dir: str | Path | None = None,
    everything: bool = True,
    fork: int = 4,
    buffer_size: int = 5_000,
    extra_args: list[str] | None = None,
) -> Path:
    """Run VEP on ``input_vcf`` and write ``output_vcf`` (gzipped)."""
    input_vcf  = Path(input_vcf)
    output_vcf = Path(output_vcf)
    output_vcf.parent.mkdir(parents=True, exist_ok=True)

    cmd = [
        "vep",
        "--input_file",  str(input_vcf),
        "--output_file", str(output_vcf),
        "--vcf", "--compress_output", "gzip",
        "--force_overwrite",
        "--species", species,
        "--assembly", genome,
        "--cache",
        "--fork", str(fork),
        "--buffer_size", str(buffer_size),
    ]
    if everything:
        cmd.append("--everything")
    if cache_dir:
        cmd += ["--dir_cache", str(cache_dir)]
    if extra_args:
        cmd += extra_args

    print(f"  $ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "vep binary not found. Install with: "
            "python scripts/install_tools.py --tool vep") from exc
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(f"VEP failed: {exc}\n")
        raise

    print(f"✓ VEP annotation completed: {output_vcf}")
    return output_vcf


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--genome", default="GRCh38")
    parser.add_argument("--species", default="homo_sapiens")
    parser.add_argument("--cache-dir", type=Path, default=None)
    parser.add_argument("--no-everything", action="store_true",
                        help="Disable --everything (basic annotation only)")
    parser.add_argument("--fork", type=int, default=4)
    parser.add_argument("--buffer-size", type=int, default=5_000)
    parser.add_argument("extra", nargs="*",
                        help="Additional VEP flags forwarded verbatim")
    args = parser.parse_args()

    installed, path, version = check_vep_installation()
    if not installed:
        sys.exit("VEP not installed. Run: "
                 "python scripts/install_tools.py --tool vep")
    print(f"  VEP: {path} ({version})")

    run_vep(
        args.input, args.output,
        genome=args.genome,
        species=args.species,
        cache_dir=args.cache_dir,
        everything=not args.no_everything,
        fork=args.fork,
        buffer_size=args.buffer_size,
        extra_args=args.extra or None,
    )


if __name__ == "__main__":
    main()
