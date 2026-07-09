"""Tier-1 fast t-test screening for every perturbation vs control.

Loops over all perturbations, runs a Welch t-test against the non-targeting
control, FDR-corrects, attaches per-group means and cell counts, and writes
one CSV per perturbation.

Example:
    python screen_all_perturbations.py --in normalized.h5ad \\
        --out-dir results/initial_screening/
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional

try:
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import diffxpy.api as de
except ImportError:
    sys.exit("Install with: pip install scanpy diffxpy pandas numpy")


def screen_all_perturbations(
    adata,
    control_group: str = "non-targeting",
    gene_col: str = "gene",
    test_method: str = "t-test",
    min_cells_per_group: int = 20,
    output_dir: str = "results/initial_screening/",
    verbose: bool = True,
) -> Dict[str, "pd.DataFrame"]:
    """Run DE for every perturbation-vs-control comparison."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    perturbations = [g for g in adata.obs[gene_col].unique() if g != control_group]
    n_control = int((adata.obs[gene_col] == control_group).sum())
    if n_control < min_cells_per_group:
        raise ValueError(
            f"Only {n_control} control cells; need at least {min_cells_per_group}"
        )

    if verbose:
        print(f"Testing {len(perturbations)} perturbations vs '{control_group}' (n={n_control})")
        print(f"  method={test_method} min_cells_per_group={min_cells_per_group}")

    results: Dict[str, pd.DataFrame] = {}
    skipped: List[str] = []

    for i, gene in enumerate(perturbations, 1):
        if verbose and i % 50 == 0:
            print(f"  [{i}/{len(perturbations)}]")
        n_pert = int((adata.obs[gene_col] == gene).sum())
        if n_pert < min_cells_per_group:
            skipped.append(gene)
            continue

        sub = adata[adata.obs[gene_col].isin([gene, control_group])].copy()
        try:
            if test_method == "t-test":
                test = de.test.t_test(data=sub, grouping=gene_col, is_logged=True)
            elif test_method == "wilcoxon":
                test = de.test.wilcoxon(data=sub, grouping=gene_col)
            else:
                raise ValueError(f"Unknown test_method {test_method!r}")

            df = test.summary()
            ctrl = np.asarray(sub[sub.obs[gene_col] == control_group].X.mean(axis=0)).flatten()
            pert = np.asarray(sub[sub.obs[gene_col] == gene].X.mean(axis=0)).flatten()
            df["mean_control"] = ctrl
            df["mean_perturbed"] = pert
            df["n_cells_control"] = n_control
            df["n_cells_perturbed"] = n_pert
            df = df.sort_values("pval")
            df.to_csv(Path(output_dir) / f"{gene}_ttest.csv", index=False)
            results[gene] = df
        except Exception as exc:
            if verbose:
                print(f"  {gene}: error — {exc}")
            skipped.append(gene)

    if verbose:
        print(f"✓ Screening complete: {len(results)} perturbations tested, {len(skipped)} skipped")
    return results


def call_hits(
    de_results: Dict[str, "pd.DataFrame"],
    min_de_genes: int = 10,
    fdr_threshold: float = 0.05,
    min_abs_log2fc: float = 0.5,
    top_n: Optional[int] = None,
) -> "pd.DataFrame":
    """Filter perturbations by number of DE genes passing FDR + |log2FC|."""
    rows = []
    for gene, df in de_results.items():
        sig = df[(df["qval"] < fdr_threshold) & (np.abs(df["log2fc"]) > min_abs_log2fc)]
        top_gene = sig.iloc[0]["gene"] if len(sig) else None
        top_lfc = float(sig.iloc[0]["log2fc"]) if len(sig) else 0.0
        rows.append({
            "perturbation": gene,
            "n_de_genes": len(sig),
            "n_cells_perturbed": int(df["n_cells_perturbed"].iloc[0]) if len(df) else 0,
            "n_cells_control": int(df["n_cells_control"].iloc[0]) if len(df) else 0,
            "top_de_gene": top_gene,
            "top_de_log2fc": top_lfc,
        })
    out = pd.DataFrame(rows)
    out = out[out["n_de_genes"] >= min_de_genes].sort_values("n_de_genes", ascending=False)
    if top_n:
        out = out.head(top_n)
    return out


def summarize_screening_results(
    de_results: Dict[str, "pd.DataFrame"],
    output_file: str = "results/screening_summary.txt",
) -> None:
    """Distribution-of-DE-genes summary as plain text."""
    counts = np.array([
        int((df["qval"] < 0.05).sum()) for df in de_results.values()
    ])
    n = len(counts)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        f.write("CRISPR screen — Tier-1 screening summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Perturbations tested: {n}\n\n")
        if n:
            f.write("DE genes per perturbation:\n")
            f.write(f"  mean   {counts.mean():.1f}\n")
            f.write(f"  median {np.median(counts):.0f}\n")
            f.write(f"  range  [{counts.min()}, {counts.max()}]\n\n")
            f.write("Hit-rate thresholds:\n")
            for k in (10, 20, 50):
                hits = int((counts >= k).sum())
                f.write(f"  ≥{k:>3} DE genes: {hits}/{n} ({hits/n:.1%})\n")
    print(f"✓ Wrote {output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out-dir", default="results/initial_screening/")
    parser.add_argument("--control", default="non-targeting")
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--test", choices=["t-test", "wilcoxon"], default="t-test")
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--min-de-genes", type=int, default=10)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.in_path)
    de_results = screen_all_perturbations(
        adata, control_group=args.control, gene_col=args.gene_col,
        test_method=args.test, min_cells_per_group=args.min_cells,
        output_dir=args.out_dir,
    )
    hits = call_hits(de_results, min_de_genes=args.min_de_genes)
    hits.to_csv(Path(args.out_dir) / "hit_summary.csv", index=False)
    summarize_screening_results(de_results, output_file=str(Path(args.out_dir) / "screening_summary.txt"))
    print(f"✓ Tier-1 hits: {len(hits)}")


if __name__ == "__main__":
    main()
