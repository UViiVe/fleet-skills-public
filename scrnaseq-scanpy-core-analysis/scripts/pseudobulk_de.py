"""Pseudobulk differential expression for multi-sample, multi-condition designs.

Aggregates raw counts per (sample, cell_type), then runs DESeq2 via pydeseq2.

Inputs:  AnnData with `layers['counts']`, plus `sample_key`, `celltype_key`, and
         a `condition` covariate in `obs`.
Outputs: per-celltype DE result CSVs and volcano/MA plots.
Example:
    python pseudobulk_de.py --input adata.h5ad --output-dir results/de \\
        --sample-key sample_id --celltype-key cell_type \\
        --formula "~ batch + condition" --contrast condition treated control
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def aggregate_to_pseudobulk(
    adata: "AnnData",
    sample_key: str = "sample_id",
    celltype_key: str = "cell_type",
    min_cells: int = 10,
    min_counts: int = 1,
) -> Dict[str, Dict[str, "pd.DataFrame"]]:
    """Return a dict of {cell_type: {'counts': DataFrame, 'meta': DataFrame}}."""
    if "counts" not in adata.layers:
        raise ValueError("adata.layers['counts'] required (raw counts).")
    if sample_key not in adata.obs.columns or celltype_key not in adata.obs.columns:
        raise KeyError(f"obs must contain `{sample_key}` and `{celltype_key}`")

    out: Dict[str, Dict[str, pd.DataFrame]] = {}
    samples = list(adata.obs[sample_key].astype(str).unique())
    celltypes = list(adata.obs[celltype_key].astype(str).unique())
    gene_names = list(adata.var_names)

    for ct in celltypes:
        counts_rows = []
        meta_rows = []
        for s in samples:
            mask = (adata.obs[sample_key].astype(str) == s) & (adata.obs[celltype_key].astype(str) == ct)
            n = int(mask.sum())
            if n < min_cells:
                continue
            X = adata[mask, :].layers["counts"]
            row = np.asarray(X.sum(axis=0)).ravel()
            if row.sum() < min_counts:
                continue
            counts_rows.append(row)
            meta = {"sample": s, "cell_type": ct, "n_cells": n}
            # Carry over per-sample condition / batch / donor metadata if uniquely defined
            for col in ("condition", "batch", "donor", "sex"):
                if col in adata.obs.columns:
                    vals = adata.obs.loc[mask, col].unique()
                    meta[col] = str(vals[0]) if len(vals) == 1 else ",".join(map(str, vals))
            meta_rows.append(meta)
        if not counts_rows:
            continue
        counts_df = pd.DataFrame(np.vstack(counts_rows), columns=gene_names,
                                  index=[m["sample"] for m in meta_rows])
        meta_df = pd.DataFrame(meta_rows).set_index("sample")
        out[ct] = {"counts": counts_df, "meta": meta_df}
        print(f"  {ct}: {counts_df.shape[0]} samples × {counts_df.shape[1]} genes")
    print(f"✓ Aggregated pseudobulk for {len(out)} cell types")
    return out


def validate_pseudobulk_design(
    meta: pd.DataFrame,
    condition: str,
    require_min: int = 2,
) -> bool:
    """Return True if every level of `condition` has ≥ `require_min` samples."""
    if condition not in meta.columns:
        print(f"  WARNING: condition column `{condition}` missing — cannot test")
        return False
    counts = meta[condition].value_counts()
    bad = counts[counts < require_min]
    if not bad.empty:
        print(f"  REFUSING to run DE — levels with N < {require_min}: {bad.to_dict()}")
        return False
    return True


def run_deseq2_analysis(
    pseudobulk: Dict[str, Dict[str, pd.DataFrame]],
    formula: str = "~ condition",
    contrast: Optional[List[str]] = None,
    celltype_key: str = "cell_type",
    alpha: float = 0.05,
    output_dir: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """Run DESeq2 via pydeseq2 per cell type. Skips cell types with N<2 per condition."""
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        sys.exit("Install with: pip install pydeseq2")

    if contrast is None:
        raise ValueError("contrast=[factor, level1, level0] required")
    factor, lvl_a, lvl_b = contrast
    design_factors = [t.strip() for t in formula.replace("~", "").split("+") if t.strip()]
    if factor not in design_factors:
        design_factors.append(factor)

    out: Dict[str, pd.DataFrame] = {}
    for ct, d in pseudobulk.items():
        counts = d["counts"]
        meta = d["meta"]
        if not validate_pseudobulk_design(meta, factor):
            continue
        print(f"-- {ct} --")
        dds = DeseqDataSet(counts=counts.astype(int), metadata=meta, design_factors=design_factors)
        dds.deseq2()
        stats = DeseqStats(dds, contrast=[factor, lvl_a, lvl_b], alpha=alpha)
        stats.summary()
        res = stats.results_df.copy()
        res["gene"] = res.index
        res["cell_type"] = ct
        out[ct] = res.reset_index(drop=True)
        print(f"  {ct}: {int((res['padj'] < alpha).sum())} sig at padj<{alpha}")

    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        for ct, res in out.items():
            safe = ct.replace("/", "_").replace(" ", "_")
            res.to_csv(Path(output_dir) / f"{safe}_de_results.csv", index=False)
        print(f"✓ Saved DE results to {output_dir}")
    return out


def export_de_results(results: Dict[str, pd.DataFrame], output_dir: str) -> None:
    """Write per-cell-type CSVs and a combined long-form CSV."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    combined = []
    for ct, df in results.items():
        safe = ct.replace("/", "_").replace(" ", "_")
        df.to_csv(out / f"{safe}_de_results.csv", index=False)
        combined.append(df)
    if combined:
        pd.concat(combined, ignore_index=True).to_csv(out / "all_de_results.csv", index=False)
    print(f"✓ Exported DE results to {out}")


def plot_volcano(
    results: pd.DataFrame,
    celltype: str = "",
    padj_threshold: float = 0.05,
    log2fc_threshold: float = 0.5,
    output_dir: str = "results/pseudobulk_de",
) -> None:
    """Save a volcano plot for one cell type."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    df = results.copy().dropna(subset=["padj", "log2FoldChange"])
    df["-log10_padj"] = -np.log10(np.maximum(df["padj"], 1e-300))
    sig = (df["padj"] < padj_threshold) & (df["log2FoldChange"].abs() > log2fc_threshold)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(df.loc[~sig, "log2FoldChange"], df.loc[~sig, "-log10_padj"],
               s=4, alpha=0.4, color="gray")
    ax.scatter(df.loc[sig, "log2FoldChange"], df.loc[sig, "-log10_padj"],
               s=6, alpha=0.7, color="crimson")
    ax.axhline(-np.log10(padj_threshold), color="k", lw=0.5, ls="--")
    ax.axvline(log2fc_threshold, color="k", lw=0.5, ls="--")
    ax.axvline(-log2fc_threshold, color="k", lw=0.5, ls="--")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    ax.set_title(f"Volcano — {celltype}")
    plt.tight_layout()
    safe = celltype.replace("/", "_").replace(" ", "_")
    out = Path(output_dir) / f"volcano_{safe}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved {out}")


def plot_ma(
    results: pd.DataFrame,
    celltype: str = "",
    output_dir: str = "results/pseudobulk_de",
) -> None:
    """Save an MA plot (baseMean vs log2FC) for one cell type."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        sys.exit("Install with: pip install matplotlib")
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    df = results.copy().dropna(subset=["baseMean", "log2FoldChange"])
    sig = df["padj"] < 0.05 if "padj" in df.columns else (np.abs(df["log2FoldChange"]) > 1)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(df.loc[~sig, "baseMean"], df.loc[~sig, "log2FoldChange"],
               s=4, alpha=0.4, color="gray")
    ax.scatter(df.loc[sig, "baseMean"], df.loc[sig, "log2FoldChange"],
               s=6, alpha=0.7, color="crimson")
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xscale("log")
    ax.set_xlabel("mean expression")
    ax.set_ylabel("log2 fold change")
    ax.set_title(f"MA — {celltype}")
    plt.tight_layout()
    safe = celltype.replace("/", "_").replace(" ", "_")
    out = Path(output_dir) / f"ma_{safe}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"✓ Saved {out}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--sample-key", default="sample_id")
    parser.add_argument("--celltype-key", default="cell_type")
    parser.add_argument("--formula", default="~ batch + condition")
    parser.add_argument("--contrast", nargs=3, required=True,
                        metavar=("FACTOR", "LEVEL_A", "LEVEL_B"))
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-cells", type=int, default=10)
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    pseudobulk = aggregate_to_pseudobulk(
        adata, sample_key=args.sample_key, celltype_key=args.celltype_key,
        min_cells=args.min_cells,
    )
    results = run_deseq2_analysis(
        pseudobulk, formula=args.formula, contrast=args.contrast,
        alpha=args.alpha, output_dir=args.output_dir,
    )
    for ct, res in results.items():
        plot_volcano(res, celltype=ct, padj_threshold=args.alpha, output_dir=args.output_dir)
        plot_ma(res, celltype=ct, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
