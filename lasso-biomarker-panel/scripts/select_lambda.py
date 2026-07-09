"""Select a LASSO regularization strength via the 1-SE rule.

Given a CSV of cross-validation deviance (or any loss metric) across a grid
of lambda values, pick lambda.min (best mean score) and lambda.1se (the
largest lambda whose score is within 1 standard error of the best). Also
produces a diagnostic plot of the CV curve with both lambda choices marked,
so the user can sanity-check the choice before locking the panel.

The 1-SE rule is the standard parsimony rule used by glmnet: it trades a
small amount of predictive accuracy for fewer non-zero coefficients, and
nearly always produces a more interpretable biomarker panel.

Inputs:
    --cv-curve PATH    CSV with columns: lambda, mean_score, sd_score, n.
                       Optional column `n_nonzero` is plotted on the
                       secondary axis if present.
    --metric STR       'deviance' (lower is better) or 'auc' (higher is
                       better). Default: 'auc'.
    --output-dir DIR   Directory to write outputs. Default: ./results.
    --plot-format STR  'png' or 'svg'. Default: 'png'.
    --no-plot          Skip the diagnostic plot.

Outputs:
    <output_dir>/lambda_selection.json    chosen lambda.min, lambda.1se
    <output_dir>/lambda_selection.<fmt>  diagnostic CV curve plot

Example:
    python select_lambda.py --cv-curve results/cv_curve.csv \\
        --metric auc --output-dir results/
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

try:
    import numpy as np
except ImportError:
    sys.exit("Install numpy with: pip install numpy")

try:
    import pandas as pd
except ImportError:
    sys.exit("Install pandas with: pip install pandas")


def _maybe_import_matplotlib() -> Any:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ImportError:
        sys.exit("Install matplotlib with: pip install matplotlib")


def _validate_curve(curve: pd.DataFrame) -> None:
    required = {"lambda", "mean_score", "sd_score"}
    missing = required - set(curve.columns)
    if missing:
        sys.exit(f"cv-curve CSV missing columns: {sorted(missing)}")
    if len(curve) < 2:
        sys.exit("cv-curve must contain at least 2 rows.")


def pick_lambda_one_se(
    curve: pd.DataFrame,
    metric: str,
) -> dict[str, Any]:
    """Return lambda.min, lambda.1se, and their CV positions."""
    metric = metric.lower()
    if metric not in {"deviance", "auc"}:
        sys.exit("metric must be 'deviance' or 'auc'.")

    higher_is_better = metric == "auc"

    df = curve.copy().sort_values("lambda").reset_index(drop=True)

    if higher_is_better:
        best_idx = int(df["mean_score"].idxmax())
    else:
        best_idx = int(df["mean_score"].idxmin())

    best_score = float(df.loc[best_idx, "mean_score"])
    best_se = float(df.loc[best_idx, "sd_score"])
    best_lambda = float(df.loc[best_idx, "lambda"])

    if higher_is_better:
        threshold = best_score - best_se
        valid = df[df["mean_score"] >= threshold]
    else:
        threshold = best_score + best_se
        valid = df[df["mean_score"] <= threshold]

    if valid.empty:
        one_se_idx = best_idx
    else:
        one_se_idx = int(valid["lambda"].idxmax())

    one_se_lambda = float(df.loc[one_se_idx, "lambda"])
    one_se_score = float(df.loc[one_se_idx, "mean_score"])

    return {
        "metric": metric,
        "lambda_min": best_lambda,
        "lambda_min_score": best_score,
        "lambda_min_idx": best_idx,
        "lambda_1se": one_se_lambda,
        "lambda_1se_score": one_se_score,
        "lambda_1se_idx": one_se_idx,
        "best_se": best_se,
        "threshold": float(threshold),
    }


def plot_cv_curve(
    curve: pd.DataFrame,
    selection: dict[str, Any],
    output_path: Path,
) -> None:
    plt = _maybe_import_matplotlib()
    df = curve.sort_values("lambda").reset_index(drop=True)
    log_lambda = np.log(df["lambda"].clip(lower=1e-12))

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(
        log_lambda,
        df["mean_score"],
        yerr=df["sd_score"],
        fmt="o-",
        color="#3498DB",
        ecolor="grey",
        elinewidth=1,
        capsize=2,
        markersize=4,
        label=f"mean +/- SE ({selection['metric']})",
    )
    ax.axvline(
        np.log(selection["lambda_min"]),
        color="#E74C3C",
        linestyle="--",
        label=f"lambda.min ({selection['lambda_min']:.4g})",
    )
    ax.axvline(
        np.log(selection["lambda_1se"]),
        color="#27AE60",
        linestyle="--",
        label=f"lambda.1se ({selection['lambda_1se']:.4g})",
    )
    ax.set_xlabel("log(lambda)")
    ax.set_ylabel(selection["metric"])
    ax.set_title("LASSO CV curve with 1-SE rule selection")
    ax.legend(loc="best", frameon=False)

    if "n_nonzero" in df.columns:
        ax2 = ax.twinx()
        ax2.plot(
            log_lambda,
            df["n_nonzero"],
            color="grey",
            alpha=0.5,
            linestyle=":",
            label="n nonzero",
        )
        ax2.set_ylabel("nonzero coefficients", color="grey")

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Pick LASSO lambda via 1-SE rule and plot the CV curve."
    )
    parser.add_argument("--cv-curve", required=True, type=Path)
    parser.add_argument("--metric", default="auc", choices=["auc", "deviance"])
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--plot-format", default="png", choices=["png", "svg"])
    parser.add_argument("--no-plot", action="store_true")
    args = parser.parse_args()

    if not args.cv_curve.exists():
        sys.exit(f"CV curve file not found: {args.cv_curve}")

    curve = pd.read_csv(args.cv_curve)
    _validate_curve(curve)

    selection = pick_lambda_one_se(curve, args.metric)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    json_path = args.output_dir / "lambda_selection.json"
    with open(json_path, "w") as fh:
        json.dump(selection, fh, indent=2)

    if not args.no_plot:
        plot_path = args.output_dir / f"lambda_selection.{args.plot_format}"
        plot_cv_curve(curve, selection, plot_path)
        print(f"  diagnostic plot: {plot_path}")

    delta = abs(selection["lambda_min_score"] - selection["lambda_1se_score"])
    print(
        "OK select_lambda: lambda.min=%.4g, lambda.1se=%.4g (delta %s = %.4f)"
        % (selection["lambda_min"], selection["lambda_1se"], args.metric, delta)
    )


if __name__ == "__main__":
    main()
