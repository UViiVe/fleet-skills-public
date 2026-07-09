"""Multiple-testing correction with method recommender.

Inputs:
    --input PATH         CSV with a p-value column.
    --pval-col COL       Name of the p-value column (default 'pvalue').
    --method {bh,bonferroni,by,holm,fdr_tsbh,recommend}
    --alpha FLOAT        Significance level (default 0.05).
    --output PATH        Output CSV with adjusted p-values and significance flag.

Outputs:
    CSV with two new columns:
        <pval_col>_adj    adjusted p-value (q-value for FDR methods)
        <pval_col>_sig    boolean significance flag

Example:
    python multiple_testing.py --input de_results.csv \
        --pval-col pvalue --method bh --alpha 0.05 --output adj.csv
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Literal

try:
    import pandas as pd
    import numpy as np
except ImportError:
    sys.exit("Install pandas and numpy with: pip install pandas numpy")

try:
    from statsmodels.stats.multitest import multipletests
except ImportError:
    sys.exit("Install statsmodels with: pip install statsmodels")


VALID_METHODS = {
    "bh": "fdr_bh",
    "bonferroni": "bonferroni",
    "by": "fdr_by",
    "holm": "holm",
    "holm-sidak": "holm-sidak",
    "fdr_tsbh": "fdr_tsbh",
}


def adjust_pvalues(
    pvalues: pd.Series | np.ndarray,
    method: str = "bh",
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    """Adjust p-values and return (reject_flags, adjusted_pvalues).

    Handles NA values by passing through (returns NaN for both).
    """
    if method not in VALID_METHODS:
        raise ValueError(
            f"method must be one of {list(VALID_METHODS)}, got {method!r}"
        )
    statsmodels_method = VALID_METHODS[method]

    pvals = pd.Series(pvalues)
    ok = pvals.notna()
    reject = np.full(len(pvals), False)
    adj = np.full(len(pvals), np.nan)

    if ok.sum() > 0:
        rej_ok, adj_ok, _, _ = multipletests(
            pvals[ok].values, alpha=alpha, method=statsmodels_method
        )
        reject[ok.values] = rej_ok
        adj[ok.values] = adj_ok

    return reject, adj


def recommend_method(
    n_tests: int,
    n_expected_true_positives: int | None = None,
    has_informative_covariate: bool = False,
    sample_size_per_group: int | None = None,
) -> dict[str, str]:
    """Recommend a multiple-testing method given context."""
    if n_tests < 100:
        return {
            "method": "bonferroni",
            "rationale": (
                f"With only {n_tests} tests, Bonferroni gives strong FWER control. "
                f"Threshold: p < {0.05 / n_tests:.2e}."
            ),
        }
    if n_tests > 500_000:
        return {
            "method": "bonferroni",
            "rationale": (
                "GWAS-scale: use the community-standard 5e-8 threshold."
            ),
        }
    if has_informative_covariate and n_tests > 5_000:
        return {
            "method": "IHW (use the R package — statsmodels does not implement it)",
            "rationale": (
                "Heterogeneous effects + informative covariate (e.g., mean expression) "
                "make IHW more powerful than BH-FDR."
            ),
        }
    if n_expected_true_positives and n_expected_true_positives > n_tests * 0.05:
        return {
            "method": "q-value (use the R 'qvalue' package)",
            "rationale": (
                f"With many expected true positives ({n_expected_true_positives}/{n_tests}), "
                "q-value is more powerful than BH."
            ),
        }
    if sample_size_per_group and sample_size_per_group < 6:
        return {
            "method": "permutation FDR",
            "rationale": (
                f"Sample size n={sample_size_per_group} is small — parametric "
                "p-values may be miscalibrated. Use a permutation null distribution."
            ),
        }
    return {
        "method": "bh",
        "rationale": (
            "BH-FDR is the genomics default — well understood by reviewers, "
            "broadly accepted, balances power and false-discovery control. "
            "Use q < 0.05."
        ),
    }


def _read_table(path: Path) -> pd.DataFrame:
    ext = path.suffix.lower()
    if ext in (".csv", ".txt"):
        return pd.read_csv(path)
    if ext == ".tsv":
        return pd.read_csv(path, sep="\t")
    if ext == ".parquet":
        return pd.read_parquet(path)
    raise ValueError(f"Unsupported extension: {ext}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Adjust p-values for multiple testing.")
    parser.add_argument("--input", required=True, help="CSV/TSV file with p-values")
    parser.add_argument("--pval-col", default="pvalue", help="Name of the p-value column")
    parser.add_argument(
        "--method",
        default="bh",
        choices=list(VALID_METHODS) + ["recommend"],
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--n-tests", type=int, default=None,
                        help="Override number of tests for recommend mode")
    parser.add_argument("--output", default=None, help="Output CSV path")
    args = parser.parse_args()

    df = _read_table(Path(args.input))
    if args.pval_col not in df.columns:
        sys.exit(f"Column '{args.pval_col}' not in input")

    if args.method == "recommend":
        n_tests = args.n_tests or len(df)
        rec = recommend_method(n_tests)
        print(f"recommended_method: {rec['method']}")
        print(f"rationale: {rec['rationale']}")
        return

    reject, adj = adjust_pvalues(df[args.pval_col], method=args.method, alpha=args.alpha)
    df[f"{args.pval_col}_adj"] = adj
    df[f"{args.pval_col}_sig"] = reject

    n_sig = int(np.nansum(reject))
    print(f"✓ Multiple testing completed: {n_sig}/{len(df)} significant "
          f"(method={args.method}, alpha={args.alpha})")

    if args.output:
        df.to_csv(args.output, index=False)
        print(f"✓ Adjusted output written: {args.output}")


if __name__ == "__main__":
    main()
