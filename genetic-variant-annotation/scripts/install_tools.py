#!/usr/bin/env python3
"""Install VEP or SNPEff via conda (or report a useful error if conda is not present).

Provides two convenience functions plus a CLI:
  - ``install_vep(channel='bioconda')``    — installs ``ensembl-vep`` via conda.
  - ``install_snpeff(channel='bioconda')`` — installs ``snpeff`` via conda.

The functions return ``True`` on success and raise ``RuntimeError`` otherwise.
Cache / database installation is a separate step — see ``installation_guide.md``.

Example:
    python scripts/install_tools.py --tool vep
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import sys


def _conda_available() -> str | None:
    """Return path to a usable conda-compatible binary, or None."""
    for candidate in ("mamba", "micromamba", "conda"):
        path = shutil.which(candidate)
        if path:
            return path
    return None


def _conda_install(package: str, channel: str) -> bool:
    binary = _conda_available()
    if not binary:
        raise RuntimeError(
            "no conda / mamba binary on PATH. "
            "Install Miniconda first, then rerun.")
    cmd = [binary, "install", "-y", "-c", channel, package]
    print(f"  $ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"conda install failed: {exc}") from exc
    return True


def install_vep(channel: str = "bioconda") -> bool:
    """Install Ensembl VEP via conda."""
    print("=== Installing Ensembl VEP ===")
    _conda_install("ensembl-vep", channel)
    if shutil.which("vep") is None:
        raise RuntimeError(
            "vep binary not on PATH after install; activate the conda env first.")
    print("✓ VEP install completed: " + shutil.which("vep"))
    print("  next: download cache with `vep_install -a c -s homo_sapiens -y GRCh38`")
    return True


def install_snpeff(channel: str = "bioconda") -> bool:
    """Install SNPEff via conda."""
    print("=== Installing SNPEff ===")
    _conda_install("snpeff", channel)
    if shutil.which("snpEff") is None:
        raise RuntimeError(
            "snpEff not on PATH after install; activate the conda env first.")
    print("✓ SNPEff install completed: " + shutil.which("snpEff"))
    print("  next: download database with `snpEff download GRCh38.105`")
    return True


def install_python_dependencies() -> bool:
    """Install the Python runtime dependencies via pip."""
    deps = ["pandas>=2.0", "numpy>=1.24", "pysam", "cyvcf2",
            "matplotlib>=3.7", "openpyxl>=3.1"]
    cmd = [sys.executable, "-m", "pip", "install", *deps]
    print(f"  $ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"pip install failed: {exc}") from exc
    print("✓ Python dependencies completed")
    return True


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--tool", choices=["vep", "snpeff", "deps", "all"],
                        required=True)
    parser.add_argument("--channel", default="bioconda")
    args = parser.parse_args()

    if args.tool in ("deps", "all"):
        install_python_dependencies()
    if args.tool in ("vep", "all"):
        install_vep(args.channel)
    if args.tool in ("snpeff", "all"):
        install_snpeff(args.channel)


if __name__ == "__main__":
    main()
