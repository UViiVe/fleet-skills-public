"""Explicit-device monomer smoke; imports heavy dependencies only after CLI parsing."""
import argparse
import json
import math
from pathlib import Path
import time

MODEL = "biohub/ESMFold2"
REVISION = "69869f737beffec5294845ede23db5fc0b4f509e"
SEQUENCE = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"


def require_finite(values):
    values = list(values)
    if not values or not all(math.isfinite(float(x)) for x in values):
        raise ValueError("Empty or nonfinite result")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", required=True, choices=["mps", "cuda", "cpu"])
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Output directory already exists")
    import numpy as np
    import torch
    from esm.models.esmfold2 import (
        EsmFold2Model, ESMFold2InputBuilder, ProteinInput, StructurePredictionInput,
    )
    if args.device == "mps" and not torch.backends.mps.is_available():
        raise RuntimeError("MPS unavailable; no silent CPU fallback")
    if args.device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA unavailable; no silent CPU fallback")
    model = EsmFold2Model.from_pretrained(
        MODEL, revision=REVISION, device=args.device,
        esmc_precision="fp32", dtype=torch.float32,
    ).float().eval()
    model.set_kernel_backend(None)
    spi = StructurePredictionInput(sequences=[ProteinInput(id="A", sequence=SEQUENCE)])
    if args.device == "mps":
        torch.mps.synchronize()
    elif args.device == "cuda":
        torch.cuda.synchronize()
    start = time.monotonic()
    with torch.inference_mode():
        result = ESMFold2InputBuilder().fold(
            model, spi, num_loops=1, num_sampling_steps=10,
            num_diffusion_samples=1, seed=42,
        )
    if isinstance(result, list):
        result = result[0]
    if args.device == "mps":
        torch.mps.synchronize()
    elif args.device == "cuda":
        torch.cuda.synchronize()
    elapsed = time.monotonic() - start
    coords = np.asarray(result.complex.atom_positions)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError("Invalid coordinate shape")
    require_finite(coords.reshape(-1))
    plddt = float(result.plddt.mean())
    ptm = float(result.ptm)
    require_finite([plddt, ptm])
    cif = result.complex.to_mmcif()
    if not cif.strip():
        raise ValueError("Empty CIF")
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "structure.cif").write_text(cif)
    metrics = dict(model=MODEL, revision=REVISION, device=args.device,
                   plddt_mean=plddt, ptm=ptm, fold_seconds=elapsed,
                   atom_count=int(coords.shape[0]), inference_pass=True)
    (args.output / "metrics.json").write_text(json.dumps(metrics, indent=2, allow_nan=False))
    print(json.dumps(metrics, allow_nan=False))


if __name__ == "__main__":
    main()
