"""Reproducibility utilities — seed pinning and environment snapshots.

Inputs:
    --seed INT             Seed to write (default 42).
    --output-dir DIR       Where to write seeds.json and requirements.txt.

Outputs:
    <output_dir>/seeds.json
    <output_dir>/requirements.txt
    <output_dir>/python_version.txt

Example:
    python repro.py --seed 42 --output-dir manifest/
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import random
import subprocess
import sys
from pathlib import Path


def set_seeds(seed: int = 42) -> dict[str, int]:
    """Set Python, numpy, and (if available) torch / tensorflow seeds.

    Returns a dict naming which generators were set.
    """
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)

    set_generators: dict[str, int] = {"python_random": seed, "pythonhashseed": seed}

    try:
        import numpy as np
        np.random.seed(seed)
        set_generators["numpy"] = seed
    except ImportError:
        pass

    try:
        import torch
        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
        set_generators["torch"] = seed
    except ImportError:
        pass

    try:
        import tensorflow as tf
        tf.random.set_seed(seed)
        set_generators["tensorflow"] = seed
    except ImportError:
        pass

    return set_generators


def freeze_environment() -> str:
    """Return the output of `pip freeze` for the current interpreter."""
    try:
        out = subprocess.check_output(
            [sys.executable, "-m", "pip", "freeze"],
            stderr=subprocess.PIPE,
            text=True,
        )
    except subprocess.CalledProcessError as e:
        return f"# pip freeze failed: {e.stderr}\n"
    return out


def write_environment(output_path: str | Path) -> None:
    """Write pip freeze output to disk."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(freeze_environment())
    print(f"✓ Environment frozen: {p}")


def write_python_version(output_path: str | Path) -> None:
    """Write a brief Python + OS fingerprint."""
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"python_version: {sys.version}",
        f"executable: {sys.executable}",
        f"platform: {platform.platform()}",
        f"machine: {platform.machine()}",
        f"system: {platform.system()}",
    ]
    p.write_text("\n".join(lines) + "\n")
    print(f"✓ Python fingerprint: {p}")


def write_seeds(seeds: dict[str, int], output_path: str | Path) -> None:
    p = Path(output_path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w") as fh:
        json.dump(seeds, fh, indent=2)
    print(f"✓ Seeds written: {p}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Write reproducibility manifest.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", default="manifest")
    args = parser.parse_args()

    out = Path(args.output_dir)
    seeds = set_seeds(args.seed)
    write_seeds(seeds, out / "seeds.json")
    write_environment(out / "requirements.txt")
    write_python_version(out / "python_version.txt")
    print(f"✓ Reproducibility manifest written to: {out}")


if __name__ == "__main__":
    main()
