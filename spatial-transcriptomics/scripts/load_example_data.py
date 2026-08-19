"""Load 10x Visium spatial transcriptomics data.

Provides two entry points:

    load_visium_heart() -> AnnData
        Downloads (and caches) the V1_Human_Heart demonstration dataset from
        10x Genomics. Returns an AnnData with raw counts, spatial
        coordinates, and the H&E image attached.

    load_visium_data(path, library_id=None) -> AnnData
        Loads user-provided Visium data — either an .h5ad, a 10x .h5
        filtered matrix, or a Space Ranger output directory.

Inputs : optional path to user data; defaults to the demo dataset.
Outputs: anndata.AnnData with `.obsm['spatial']` populated.

Example:
    python load_example_data.py
    python load_example_data.py --path /data/visium/sample/
"""

from __future__ import annotations

import argparse
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Optional, Union


def load_visium_heart() -> "anndata.AnnData":
    """Load the V1_Human_Heart Visium demonstration dataset.

    Uses ``scanpy.datasets.visium_sge`` when available and falls back to a
    direct download from the 10x public CDN on environments where
    ``visium_sge`` is unusable (older Python tarfile lacks ``data_filter``).

    Returns
    -------
    anndata.AnnData
        Raw expression matrix with H&E image and spatial coordinates.
    """
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    print("Loading V1_Human_Heart Visium dataset")
    print("  Source: 10x Genomics Spatial Gene Expression")
    print("  Description: Human heart tissue section")
    print("  Platform: 10x Visium")
    print("  (First run downloads ~50 MB; subsequent runs use cache)")

    try:
        adata = sc.datasets.visium_sge(sample_id="V1_Human_Heart")
    except AttributeError:
        print("  scanpy.datasets.visium_sge unavailable; falling back to direct download")
        adata = _download_visium_manual("V1_Human_Heart")

    adata.var_names_make_unique()
    _print_summary(adata)
    return adata


def load_visium_data(
    path: Union[str, Path],
    library_id: Optional[str] = None,
) -> "anndata.AnnData":
    """Load a user-supplied Visium dataset.

    Parameters
    ----------
    path
        Path to a ``.h5ad`` file, a 10x ``.h5`` file, or a Space Ranger
        output directory.
    library_id
        Library identifier passed through to ``sc.read_visium`` when loading
        a Space Ranger directory with multiple captures.

    Returns
    -------
    anndata.AnnData
        Loaded spatial AnnData with ``.obsm['spatial']`` populated.
    """
    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Path not found: {path}")

    if path.suffix == ".h5ad":
        print(f"Reading AnnData from {path}")
        adata = sc.read_h5ad(path)
    elif path.suffix == ".h5":
        print(f"Reading 10x HDF5 from {path}")
        adata = sc.read_10x_h5(path)
    elif path.is_dir():
        print(f"Reading Space Ranger output from {path}")
        adata = sc.read_visium(path, library_id=library_id)
    else:
        raise ValueError(
            f"Unsupported input: {path}. Expected .h5ad, .h5, or a Space Ranger directory."
        )

    adata.var_names_make_unique()
    _print_summary(adata)
    return adata


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _print_summary(adata: "anndata.AnnData") -> None:
    has_image = "spatial" in adata.uns and len(adata.uns["spatial"]) > 0
    coords_shape = (
        adata.obsm["spatial"].shape if "spatial" in adata.obsm else "(missing)"
    )

    print("")
    print(f"✓ Data loaded successfully: {adata.n_obs} spots / {adata.n_vars} genes")
    print(f"  Spatial image attached: {'yes' if has_image else 'no'}")
    print(f"  Spatial coordinates   : {coords_shape}")


def _download_file(url: str, dest: Path) -> None:
    """Download a single file with a polite User-Agent header."""
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "Mozilla/5.0 (compatible; spatial-skill-downloader)"},
    )
    with urllib.request.urlopen(req) as response, open(dest, "wb") as out:
        while True:
            chunk = response.read(8192)
            if not chunk:
                break
            out.write(chunk)


def _download_visium_manual(sample_id: str) -> "anndata.AnnData":
    """Fallback Visium download for environments where ``visium_sge`` fails."""
    import scanpy as sc

    base_url = "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0"
    cache_root = Path(sc.settings.datasetdir) / sample_id
    cache_root.mkdir(parents=True, exist_ok=True)

    h5_path = cache_root / "filtered_feature_bc_matrix.h5"
    if not h5_path.exists():
        print("  Downloading filtered feature-barcode matrix...")
        _download_file(
            f"{base_url}/{sample_id}/{sample_id}_filtered_feature_bc_matrix.h5",
            h5_path,
        )

    spatial_dir = cache_root / "spatial"
    if not spatial_dir.exists():
        print("  Downloading spatial directory...")
        tar_path = cache_root / "spatial.tar.gz"
        _download_file(f"{base_url}/{sample_id}/{sample_id}_spatial.tar.gz", tar_path)
        with tarfile.open(tar_path, "r:gz") as tf:
            tf.extractall(cache_root)
        tar_path.unlink()

    return sc.read_visium(cache_root)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Load 10x Visium data (example dataset by default)."
    )
    parser.add_argument(
        "--path",
        default=None,
        help="Optional path to .h5ad, .h5, or Space Ranger directory.",
    )
    parser.add_argument(
        "--library-id",
        default=None,
        help="Library ID passed to sc.read_visium for multi-library directories.",
    )
    args = parser.parse_args()

    if args.path is None:
        load_visium_heart()
    else:
        load_visium_data(args.path, library_id=args.library_id)


if __name__ == "__main__":
    main()
