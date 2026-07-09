"""Fetch and reformat the Papalexi 2021 example Perturb-seq dataset.

Downloads a scPerturb-formatted h5ad from Zenodo (Papalexi & Satija 2021
ECCITE-seq CRISPRi screen in THP-1 cells, ~140 MB) and splits it into per-batch
AnnData objects with a matching sgRNA mapping file each — the form the rest
of the workflow expects. Also offers a small synthetic ``demo`` dataset for
offline testing.

Example:
    python load_example_data.py --dataset papalexi2021 --download-dir data/
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List, Tuple

try:
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install numpy pandas scanpy anndata")


_ZENODO_URL = (
    "https://zenodo.org/records/13350497/files/"
    "PapalexiSatija2021_eccite_RNA.h5ad?download=1"
)
_EXPECTED_SIZE_MB = 140


def load_example_data(dataset: str = "papalexi2021", download_dir: str = "data") -> Dict:
    """Return ``{adata_list, mapping_files, metadata}`` for testing."""
    if dataset == "papalexi2021":
        return _load_papalexi2021(download_dir)
    if dataset == "demo":
        return _generate_demo()
    raise ValueError(f"Unknown dataset {dataset!r}; choose 'papalexi2021' or 'demo'")


def _download_if_needed(download_dir: str) -> Path:
    dst = Path(download_dir)
    dst.mkdir(parents=True, exist_ok=True)
    target = dst / "papalexi2021_eccite_rna.h5ad"

    if target.exists():
        size_mb = target.stat().st_size / (1024 * 1024)
        if size_mb > 100:
            print(f"Cached: {target} ({size_mb:.0f} MB)")
            return target
        print(f"Cached file too small ({size_mb:.0f} MB) — refetching")

    print(f"Downloading Papalexi 2021 (~{_EXPECTED_SIZE_MB} MB) from Zenodo …")
    try:
        urllib.request.urlretrieve(_ZENODO_URL, str(target))
    except Exception as exc:
        if target.exists():
            target.unlink()
        raise RuntimeError(
            f"Download failed: {exc}\n"
            f"Manual fallback: https://zenodo.org/records/13350497\n"
            f"  file: PapalexiSatija2021_eccite_RNA.h5ad\n"
            f"  save to: {target}"
        ) from exc
    print(f"  → {target} ({target.stat().st_size / 1024 / 1024:.0f} MB)")
    return target


def _parse_perturbation(pert_id: str) -> Tuple[str, str]:
    """Strip the 'g<digit>' suffix used by scPerturb to encode the guide number."""
    if pert_id == "control":
        return "non-targeting", "0"
    m = re.match(r"(.+?)g(\d+)$", pert_id)
    if m:
        return m.group(1), m.group(2)
    return pert_id, "1"


def _reformat_sgrna(pert_id: str, guide_id: str) -> str:
    """Translate scPerturb naming to GENE_sgRNAn format used by the workflow."""
    if pert_id == "control":
        m = re.search(r"g(\d+)$", guide_id)
        n = m.group(1) if m else "0"
        return f"non-targeting_sgRNA{n}"
    gene, num = _parse_perturbation(pert_id)
    return f"{gene}_sgRNA{num}"


def _load_papalexi2021(download_dir: str) -> Dict:
    h5ad_path = _download_if_needed(download_dir)
    print("Loading and reformatting …")
    full = sc.read_h5ad(str(h5ad_path))

    for col in ("perturbation", "guide_id"):
        if col not in full.obs.columns:
            raise ValueError(
                f"Expected column '{col}' missing — available: {list(full.obs.columns)}"
            )

    batch_col = next(
        (c for c in ("hto", "batch", "replicate", "library") if c in full.obs.columns),
        None,
    )
    if batch_col is None:
        full.obs["_batch"] = "lib1"
        batch_col = "_batch"

    batches = sorted(full.obs[batch_col].unique())
    print(
        f"  cells={full.n_obs:,} genes={full.n_vars:,} "
        f"perturbations={full.obs['perturbation'].nunique()} "
        f"guides={full.obs['guide_id'].nunique()} batches={len(batches)}"
    )

    full.obs["_sgRNA_formatted"] = [
        _reformat_sgrna(p, g)
        for p, g in zip(full.obs["perturbation"], full.obs["guide_id"])
    ]
    full.obs["_gene_formatted"] = [
        _parse_perturbation(p)[0] for p in full.obs["perturbation"]
    ]
    n_ctrl = int((full.obs["_gene_formatted"] == "non-targeting").sum())
    print(f"  controls: {n_ctrl:,}  perturbed: {full.n_obs - n_ctrl:,}")

    min_batch = 50
    valid = [b for b in batches if int((full.obs[batch_col] == b).sum()) >= min_batch]
    skipped = set(batches) - set(valid)
    if skipped:
        print(f"  skipping small batches (<{min_batch} cells): {sorted(skipped)}")

    out_dir = Path(download_dir) / "papalexi2021_split"
    out_dir.mkdir(parents=True, exist_ok=True)

    adata_list: List[ad.AnnData] = []
    mapping_files: List[str] = []
    labels: List[str] = []

    for i, batch in enumerate(valid):
        mask = full.obs[batch_col] == batch
        a = full[mask].copy()

        # Strip metadata to mimic raw 10X output; the workflow will re-attach sgRNAs
        X = a.X.toarray() if hasattr(a.X, "toarray") else np.asarray(a.X)
        adata_clean = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=a.obs_names.copy()),
            var=pd.DataFrame(index=a.var_names.copy()),
        )
        adata_clean.obs["n_counts"] = X.sum(axis=1)
        adata_clean.obs["n_genes"] = (X > 0).sum(axis=1)
        mito = adata_clean.var_names.str.startswith("MT-")
        if mito.sum() > 0:
            adata_clean.obs["percent_mito"] = (
                X[:, mito].sum(axis=1) / (adata_clean.obs["n_counts"] + 1e-6)
            )
        else:
            adata_clean.obs["percent_mito"] = 0.0
        adata_list.append(adata_clean)

        label = f"lib{i + 1}"
        labels.append(label)
        mapping_path = out_dir / f"mapped_sgRNA_to_cell_{label}.txt"
        pd.DataFrame({
            "cell_barcode": a.obs_names,
            "sgRNA": a.obs["_sgRNA_formatted"].values,
        }).to_csv(mapping_path, sep="\t", index=False, header=False)
        mapping_files.append(str(mapping_path))
        print(
            f"  {label} ({batch}): {adata_clean.n_obs:,} cells, "
            f"{a.obs['_gene_formatted'].nunique()} perturbations"
        )

    target_genes = sorted(g for g in full.obs["_gene_formatted"].unique() if g != "non-targeting")
    meta = {
        "dataset": "papalexi2021",
        "description": "Papalexi & Satija 2021 ECCITE-seq CRISPRi screen",
        "cell_type": "THP-1 (human acute monocytic leukemia)",
        "screen_type": "CRISPRi",
        "organism": "human",
        "n_libraries": len(adata_list),
        "n_cells_total": full.n_obs,
        "n_genes_expression": full.n_vars,
        "n_target_genes": len(target_genes),
        "n_guides": int(full.obs["guide_id"].nunique()),
        "guides_per_gene": round(full.obs["guide_id"].nunique() / max(len(target_genes), 1), 1),
        "control_group": "non-targeting",
        "expected_direction": "down",
        "batch_labels": labels,
        "original_batches": list(batches),
        "target_genes": target_genes,
        "reference": "Papalexi et al. (2021) Nature Genetics 53:322-331",
        "source": "scPerturb (Peidli et al. 2024) via Zenodo",
    }
    print(f"✓ Example data loaded: {len(adata_list)} libraries, {full.n_obs:,} cells")
    return {"adata_list": adata_list, "mapping_files": mapping_files, "metadata": meta}


def _generate_demo() -> Dict:
    print("Generating synthetic demo dataset …")
    rng = np.random.default_rng(42)
    out_dir = Path("data/example_demo")
    out_dir.mkdir(parents=True, exist_ok=True)

    n_cells_per_lib = 500
    n_genes = 5000
    gene_names = [f"GENE{i:05d}" for i in range(n_genes)]

    strong = [f"TARGET{i:02d}" for i in range(5)]
    moderate = [f"TARGET{i:02d}" for i in range(5, 10)]
    weak = [f"TARGET{i:02d}" for i in range(10, 15)]
    all_perts = strong + moderate + weak + ["non-targeting"]

    sgrnas: Dict[str, List[str]] = {}
    for g in strong + moderate + weak:
        sgrnas[g] = [f"{g}_sgRNA{j+1}" for j in range(2)]
    sgrnas["non-targeting"] = [f"non-targeting_sgRNA{j+1}" for j in range(5)]

    adata_list: List[ad.AnnData] = []
    mapping_files: List[str] = []

    for lib in range(2):
        X = rng.negative_binomial(5, 0.5, size=(n_cells_per_lib, n_genes))
        barcodes = [f"LIB{lib+1}_{i:04d}" for i in range(n_cells_per_lib)]
        assignments = rng.choice(all_perts, size=n_cells_per_lib)

        for pert in all_perts:
            cells = np.where(assignments == pert)[0]
            if pert == "non-targeting" or len(cells) == 0:
                continue
            if pert in strong:
                n_de, effect = int(rng.integers(50, 100)), 2.0
            elif pert in moderate:
                n_de, effect = int(rng.integers(10, 30)), 1.0
            else:
                n_de, effect = int(rng.integers(0, 5)), 0.5
            if n_de > 0:
                de_idx = rng.choice(n_genes, size=n_de, replace=False)
                for gi in de_idx:
                    direction = -1 if rng.random() < 0.6 else 1
                    fc = 2.0 ** (direction * effect * rng.uniform(0.7, 1.3))
                    X[cells, gi] = (X[cells, gi] * fc).astype(int)

        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=barcodes),
            var=pd.DataFrame(index=gene_names),
        )
        adata.obs["n_counts"] = adata.X.sum(axis=1)
        adata.obs["n_genes"] = (adata.X > 0).sum(axis=1)
        adata_list.append(adata)

        path = out_dir / f"mapped_sgRNA_to_cell_lib{lib+1}.txt"
        rows = [[bc, rng.choice(sgrnas[pert])] for bc, pert in zip(barcodes, assignments)]
        pd.DataFrame(rows).to_csv(path, sep="\t", index=False, header=False)
        mapping_files.append(str(path))

    print(f"✓ Demo data: 2 libraries × {n_cells_per_lib} cells, {len(all_perts)} perturbations")
    return {
        "adata_list": adata_list,
        "mapping_files": mapping_files,
        "metadata": {
            "dataset": "demo",
            "n_libraries": 2,
            "n_cells_total": n_cells_per_lib * 2,
            "n_genes": n_genes,
            "n_perturbations": len(all_perts),
            "control_group": "non-targeting",
            "expected_direction": "down",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--dataset", choices=["papalexi2021", "demo"], default="papalexi2021")
    parser.add_argument("--download-dir", default="data")
    args = parser.parse_args()
    data = load_example_data(args.dataset, args.download_dir)
    print(f"✓ Loaded: {len(data['adata_list'])} libraries, "
          f"mapping_files={len(data['mapping_files'])}, "
          f"metadata keys={list(data['metadata'].keys())}")


if __name__ == "__main__":
    main()
