"""Ambient RNA removal via CellBender (preferred) or SoupX (Python wrapper).

Inputs:  Raw 10X `.h5` (CellBender) or raw+filtered matrix dirs (SoupX).
Outputs: A corrected AnnData object; CellBender writes additional artifacts to disk.
Example:
    python remove_ambient_rna.py --method cellbender \\
        --raw raw_feature_bc_matrix.h5 \\
        --expected-cells 10000 --total-droplets 30000 \\
        --output-dir results/cellbender
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from anndata import AnnData


def run_cellbender(
    raw_h5: str,
    expected_cells: int,
    total_droplets: int,
    output_dir: str = "results/cellbender",
    epochs: int = 200,
    fpr: float = 0.01,
    use_cuda: bool = True,
) -> "AnnData":
    """Shell out to CellBender's CLI; return the corrected `.h5` as AnnData."""
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    if shutil.which("cellbender") is None:
        sys.exit("Install with: pip install cellbender")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    out_h5 = Path(output_dir) / "cellbender_output.h5"

    cmd = [
        "cellbender", "remove-background",
        "--input", str(raw_h5),
        "--output", str(out_h5),
        "--expected-cells", str(expected_cells),
        "--total-droplets-included", str(total_droplets),
        "--epochs", str(epochs),
        "--fpr", str(fpr),
    ]
    if use_cuda:
        cmd.append("--cuda")
    print(f"$ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    filtered = Path(output_dir) / "cellbender_output_filtered.h5"
    if not filtered.exists():
        filtered = out_h5
    adata = sc.read_10x_h5(str(filtered))
    adata.var_names_make_unique()
    print(f"✓ CellBender output loaded: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


def run_soupx_python(
    raw_matrix_dir: str,
    filtered_matrix_dir: str,
    contamination: Optional[float] = None,
) -> "AnnData":
    """Lightweight Python SoupX-style correction.

    Approximation: estimates the ambient profile from empty droplets and
    subtracts a contamination fraction × ambient profile from each cell.
    For production use, defer to the official R SoupX package.
    """
    try:
        import scanpy as sc
        import numpy as np
    except ImportError:
        sys.exit("Install with: pip install scanpy numpy")

    raw = sc.read_10x_mtx(raw_matrix_dir)
    raw.var_names_make_unique()
    filt = sc.read_10x_mtx(filtered_matrix_dir)
    filt.var_names_make_unique()

    # Estimate ambient profile from droplets with very few UMIs
    totals = np.asarray(raw.X.sum(axis=1)).ravel()
    threshold = np.quantile(totals, 0.50)
    empties = raw[totals < threshold, :]
    ambient = np.asarray(empties.X.sum(axis=0)).ravel()
    ambient = ambient / max(ambient.sum(), 1)

    if contamination is None:
        raise ValueError(
            "SoupX-style correction requires an explicit contamination fraction. "
            "Pass --contamination after reviewing CellBender/SoupX metrics."
        )
    print(f"  SoupX-py: contamination fraction = {contamination:.2f}")

    X = filt.X.toarray() if hasattr(filt.X, "toarray") else np.asarray(filt.X)
    cell_totals = X.sum(axis=1, keepdims=True)
    expected = contamination * cell_totals * ambient[None, :]
    corrected = np.maximum(X - expected, 0).astype(np.float32)
    filt.X = corrected
    print(f"✓ SoupX-py correction applied")
    return filt


def estimate_contamination_fraction(
    raw_adata: "AnnData",
    filtered_adata: "AnnData",
) -> float:
    """Heuristic contamination estimate: ratio of ambient counts to total."""
    try:
        import numpy as np
    except ImportError:
        sys.exit("Install with: pip install numpy")
    raw_totals = np.asarray(raw_adata.X.sum(axis=1)).ravel()
    filt_totals = np.asarray(filtered_adata.X.sum(axis=1)).ravel()
    estimate = 1 - filt_totals.mean() / max(raw_totals.mean(), 1)
    return float(max(0.0, min(0.5, estimate)))


def estimate_contamination(adata: "AnnData") -> float:
    """Estimate contamination from correction metrics stored on an AnnData object."""
    if "cellbender_metrics" in adata.uns:
        metrics = adata.uns["cellbender_metrics"]
        if "contamination_fraction" in metrics:
            return float(metrics["contamination_fraction"])
    raise ValueError(
        "No contamination metrics were found. Provide an explicit contamination "
        "fraction rather than relying on a guessed default."
    )


def compare_before_after(
    before: "AnnData",
    after: "AnnData",
    output_dir: str = "results/ambient",
    markers: Optional[list] = None,
) -> None:
    """Compare UMI distributions before and after correction; plot marker UMAPs."""
    try:
        import matplotlib.pyplot as plt
        import numpy as np
    except ImportError:
        sys.exit("Install with: pip install matplotlib numpy")
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(6, 4))
    b = np.asarray(before.X.sum(axis=1)).ravel()
    a = np.asarray(after.X.sum(axis=1)).ravel()
    ax.hist(b, bins=80, alpha=0.5, label="before", color="steelblue")
    ax.hist(a, bins=80, alpha=0.5, label="after", color="orange")
    ax.set_xscale("log")
    ax.set_xlabel("total UMIs / cell")
    ax.legend()
    plt.tight_layout()
    out = Path(output_dir) / "umi_before_after.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print(f"✓ Saved {out}")
    pct = (1 - a.mean() / max(b.mean(), 1)) * 100
    print(f"  mean UMI reduction: {pct:.1f}%")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--method", choices=["cellbender", "soupx"], default="cellbender")
    parser.add_argument("--raw", help="raw H5 (CellBender) or raw matrix dir (SoupX)")
    parser.add_argument("--filtered", help="filtered matrix dir (SoupX only)")
    parser.add_argument("--expected-cells", type=int, default=10000)
    parser.add_argument("--total-droplets", type=int, default=30000)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--output-dir", default="results/ambient")
    parser.add_argument("--no-cuda", action="store_true")
    parser.add_argument("--contamination", type=float, default=None)
    parser.add_argument("--output-h5ad", help="Optional: write corrected adata")
    args = parser.parse_args()

    if args.method == "cellbender":
        if not args.raw:
            sys.exit("--raw required")
        adata = run_cellbender(
            args.raw,
            expected_cells=args.expected_cells,
            total_droplets=args.total_droplets,
            output_dir=args.output_dir,
            epochs=args.epochs,
            use_cuda=not args.no_cuda,
        )
    else:
        if not (args.raw and args.filtered):
            sys.exit("--raw and --filtered both required for SoupX mode")
        adata = run_soupx_python(args.raw, args.filtered, contamination=args.contamination)

    if args.output_h5ad:
        adata.write_h5ad(args.output_h5ad)
        print(f"✓ Wrote {args.output_h5ad}")


if __name__ == "__main__":
    main()
