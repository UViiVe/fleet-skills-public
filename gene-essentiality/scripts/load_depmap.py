"""Download and cache DepMap release files.

Pulls a gene-effect matrix, dependency-probability matrix (if present), and
cell-line metadata from depmap.org. Files are cached under
``~/.cache/depmap/<release>/`` and converted to parquet on first access for
fast subsequent loads.

Example:
    python load_depmap.py --release 24Q2 --files gene_effect,dependency,model
"""

from __future__ import annotations

import argparse
import os
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional

try:
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install pandas pyarrow")


_BASE_URL = "https://depmap.org/portal/api/download/files"

_KNOWN_FILES = {
    "gene_effect":  "CRISPRGeneEffect.csv",
    "dependency":   "CRISPRGeneDependency.csv",
    "model":        "Model.csv",
    "expression":   "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
    "mutations":    "OmicsSomaticMutationsMatrixHotspot.csv",
    "copy_number":  "OmicsCNGene.csv",
    "proteomics":   "ProteinMatrixAveragedNormalizedAbundance.csv",
}


def _cache_dir(release: str) -> Path:
    home = Path(os.environ.get("DEPMAP_CACHE_DIR", Path.home() / ".cache" / "depmap"))
    out = home / release
    out.mkdir(parents=True, exist_ok=True)
    return out


def _download(url: str, dest: Path) -> Path:
    """Stream a file to disk; replace partials on failure."""
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    print(f"Downloading {url}")
    try:
        urllib.request.urlretrieve(url, str(tmp))
        tmp.rename(dest)
    except Exception:
        if tmp.exists():
            tmp.unlink()
        raise
    print(f"  → {dest} ({dest.stat().st_size / 1024 / 1024:.0f} MB)")
    return dest


def _ensure_file(release: str, key: str) -> Path:
    if key not in _KNOWN_FILES:
        raise ValueError(f"Unknown key {key!r}; choices: {sorted(_KNOWN_FILES)}")
    filename = _KNOWN_FILES[key]
    cache = _cache_dir(release)
    parquet = cache / (Path(filename).stem + ".parquet")
    csv = cache / filename

    if parquet.exists():
        return parquet
    if csv.exists():
        df = pd.read_csv(csv, index_col=0)
        df.to_parquet(parquet)
        return parquet

    url = f"{_BASE_URL}/{release}/{filename}"
    _download(url, csv)
    df = pd.read_csv(csv, index_col=0)
    df.to_parquet(parquet)
    return parquet


def load_gene_effect(release: str = "24Q2") -> "pd.DataFrame":
    """Genes × cell lines matrix of gene-effect scores."""
    path = _ensure_file(release, "gene_effect")
    df = pd.read_parquet(path)
    print(f"✓ Loaded gene effect: {df.shape[0]} genes × {df.shape[1]} cell lines (release {release})")
    return df


def load_gene_dependency(release: str = "24Q2") -> Optional["pd.DataFrame"]:
    """Posterior probability of dependency (Chronos releases only)."""
    try:
        path = _ensure_file(release, "dependency")
    except Exception as exc:
        print(f"⚠ Dependency-probability table not available for release {release}: {exc}")
        return None
    df = pd.read_parquet(path)
    print(f"✓ Loaded dependency prob: {df.shape[0]} genes × {df.shape[1]} cell lines")
    return df


def load_sample_info(release: str = "24Q2") -> "pd.DataFrame":
    """Cell-line metadata (lineage, subtype, mutations summary)."""
    path = _ensure_file(release, "model")
    df = pd.read_parquet(path)
    print(f"✓ Loaded sample info: {df.shape[0]} cell lines × {df.shape[1]} metadata columns")
    return df


def load_expression(release: str = "24Q2") -> "pd.DataFrame":
    """log2(TPM+1) protein-coding expression matrix."""
    path = _ensure_file(release, "expression")
    df = pd.read_parquet(path)
    print(f"✓ Loaded expression: {df.shape[0]} genes × {df.shape[1]} cell lines")
    return df


def load_mutations_hotspot(release: str = "24Q2") -> "pd.DataFrame":
    """Boolean hotspot-mutation matrix."""
    path = _ensure_file(release, "mutations")
    df = pd.read_parquet(path)
    print(f"✓ Loaded hotspot mutations: {df.shape[0]} genes × {df.shape[1]} cell lines")
    return df


def detect_scoring_algorithm(release_dir: str) -> str:
    """Heuristically identify CERES vs Chronos from cached files."""
    files = set(os.listdir(release_dir)) if os.path.exists(release_dir) else set()
    if "CRISPRGeneDependency.csv" in files or "CRISPRGeneDependency.parquet" in files:
        return "Chronos"
    if "Achilles_gene_effect.csv" in files and "Achilles_gene_dependency.csv" not in files:
        return "CERES"
    if "CRISPRGeneEffect.csv" in files or "CRISPRGeneEffect.parquet" in files:
        return "Chronos"
    return "unknown"


def verify_sign_convention(gene_effect: "pd.DataFrame", verbose: bool = True) -> bool:
    """Sanity check: known pan-essentials must have median < -0.8."""
    landmarks = ["RPS14", "RPL5", "POLR2A", "EIF3A"]
    ok = True
    for sym in landmarks:
        rows = [i for i in gene_effect.index if i.startswith(f"{sym} (")]
        if not rows:
            if verbose:
                print(f"⚠ {sym}: not in matrix")
            continue
        med = float(gene_effect.loc[rows[0]].median())
        if verbose:
            print(f"  {sym}: median gene effect = {med:.2f}")
        if med > -0.8:
            ok = False
            print(f"  ⚠ {sym} median > -0.8 — sign convention may be off or data is inverted")
    if ok and verbose:
        print("✓ Sign convention verified")
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--release", default="24Q2",
                        help="DepMap release tag, e.g. 24Q2, 23Q4, 21Q1")
    parser.add_argument("--files", default="gene_effect,model",
                        help="Comma-separated keys; choices: " + ", ".join(_KNOWN_FILES))
    parser.add_argument("--check-sign", action="store_true",
                        help="Verify sign convention via pan-essential medians")
    args = parser.parse_args()

    keys = [k.strip() for k in args.files.split(",") if k.strip()]
    paths: Dict[str, Path] = {}
    for key in keys:
        paths[key] = _ensure_file(args.release, key)

    if args.check_sign and "gene_effect" in paths:
        ge = pd.read_parquet(paths["gene_effect"])
        verify_sign_convention(ge)

    print("✓ Files cached:")
    for k, p in paths.items():
        print(f"  {k}: {p}")


if __name__ == "__main__":
    main()
