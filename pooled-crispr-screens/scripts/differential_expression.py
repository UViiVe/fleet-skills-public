"""Per-perturbation differential expression on log-normalized data.

A thin wrapper around diffxpy that runs the chosen test for every perturbation
against the non-targeting control. Useful when the per-perturbation pickle is
not enough and you want raw DE CSVs in a directory.

Example:
    python differential_expression.py --in normalized.h5ad \\
        --out-dir DEG/ --test t-test
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, Optional

try:
    import pandas as pd
    import scanpy as sc
    import diffxpy.api as de
except ImportError:
    sys.exit("Install with: pip install scanpy diffxpy pandas")


def run_de_analysis(
    adata,
    control_group: str = "non-targeting",
    gene_col: str = "gene",
    test_method: str = "t-test",
    use_outliers_only: bool = False,
    outlier_cells: Optional[list] = None,
    output_dir: str = "DEG/",
    is_logged: bool = True,
) -> Dict[str, "pd.DataFrame"]:
    """Run DE per perturbation, writing one CSV per gene."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if use_outliers_only and outlier_cells:
        ctrl_cells = adata[adata.obs[gene_col] == control_group].obs_names.tolist()
        kept = set(outlier_cells) | set(ctrl_cells)
        adata = adata[adata.obs_names.isin(kept)].copy()
        print(f"Restricted to {len(outlier_cells)} outliers + {len(ctrl_cells)} controls")

    perturbations = [g for g in adata.obs[gene_col].unique() if g != control_group]
    print(f"Running {test_method} for {len(perturbations)} perturbations")

    out: Dict[str, pd.DataFrame] = {}
    for i, gene in enumerate(perturbations, 1):
        if i % 10 == 0:
            print(f"  [{i}/{len(perturbations)}] {gene}")
        sub = adata[adata.obs[gene_col].isin([control_group, gene])].copy()
        if test_method == "t-test":
            test = de.test.t_test(data=sub, grouping=gene_col, is_logged=is_logged)
        elif test_method == "wilcoxon":
            test = de.test.wilcoxon(data=sub, grouping=gene_col)
        elif test_method == "negative_binomial":
            test = de.test.wald(data=sub, formula_loc=f"~ 1 + {gene_col}",
                                factor_loc_totest=gene_col)
        else:
            raise ValueError(f"Unknown test_method {test_method!r}")
        df = test.summary().sort_values("pval")
        df.to_csv(os.path.join(output_dir, f"{gene}_{test_method}.csv"), index=False)
        out[gene] = df

    print(f"✓ DE complete: {len(out)} perturbations, results in {output_dir}")
    return out


def summarize_de_results(
    de_results: Dict[str, "pd.DataFrame"],
    pval_threshold: float = 0.05,
    log2fc_threshold: float = 0.5,
) -> "pd.DataFrame":
    """Compact summary table across all perturbations."""
    rows = []
    for gene, df in de_results.items():
        sig = df[df["pval"] < pval_threshold]
        n_up = n_down = 0
        if "log2fc" in df.columns:
            sigfc = sig[sig["log2fc"].abs() > log2fc_threshold]
            n_up = int((sigfc["log2fc"] > 0).sum())
            n_down = int((sigfc["log2fc"] < 0).sum())
        rows.append({
            "perturbation": gene,
            "n_de_genes": len(sig),
            "n_up": n_up,
            "n_down": n_down,
            "top_de_gene": sig.iloc[0]["gene"] if len(sig) else "none",
            "top_pval": float(sig.iloc[0]["pval"]) if len(sig) else 1.0,
        })
    return pd.DataFrame(rows).sort_values("n_de_genes", ascending=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out-dir", default="DEG/")
    parser.add_argument("--control", default="non-targeting")
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--test", choices=["t-test", "wilcoxon", "negative_binomial"], default="t-test")
    parser.add_argument("--outliers-csv", default=None,
                        help="Optional CSV with 'cell_barcode' column; restrict DE to those + controls")
    args = parser.parse_args()

    adata = sc.read_h5ad(args.in_path)
    outliers = None
    if args.outliers_csv:
        outliers = pd.read_csv(args.outliers_csv)["cell_barcode"].tolist()

    results = run_de_analysis(
        adata,
        control_group=args.control,
        gene_col=args.gene_col,
        test_method=args.test,
        use_outliers_only=bool(outliers),
        outlier_cells=outliers,
        output_dir=args.out_dir,
    )
    summary = summarize_de_results(results)
    summary.to_csv(os.path.join(args.out_dir, "de_summary.csv"), index=False)
    print(f"✓ Wrote summary table to {args.out_dir}/de_summary.csv")


if __name__ == "__main__":
    main()
