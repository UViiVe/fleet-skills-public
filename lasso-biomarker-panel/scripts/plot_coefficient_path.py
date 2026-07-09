"""Plot the LASSO / elastic-net coefficient regularization path.

Reads a feature matrix and binary outcome, fits a penalized logistic
regression along a grid of regularization strengths, and visualizes how
each feature's coefficient changes as lambda decreases. Optional vertical
lines mark lambda.min and lambda.1se if provided.

The regularization path is the standard diagnostic in glmnet's `plot.glmnet`
output: it shows which features "switch on" first as lambda shrinks (the
strongest signals enter early) and which are spurious (they oscillate
around zero before stabilizing). Useful for sanity-checking that the
features chosen by `fit_lasso_cv` are biologically meaningful and not
just numerically dominant at one specific lambda.

Inputs:
    --features PATH      CSV: samples x features + outcome column.
    --outcome-col STR    Outcome column name. Default: last column.
    --sample-col STR     Sample id column. Default: first column.
    --alpha FLOAT        Elastic-net mixing (1.0 = pure LASSO). Default 1.0.
    --n-lambdas INT      Lambda grid size (default 100).
    --top-n INT          Label the top-N features by max |coefficient|
                         (default 10; set 0 to label nothing).
    --highlight FILE     Optional text file with one feature name per line
                         to highlight on the plot (overrides --top-n).
    --lambda-min FLOAT   Optional vertical line at lambda.min.
    --lambda-1se FLOAT   Optional vertical line at lambda.1se.
    --output-dir DIR     Output directory. Default ./results.
    --plot-format STR    'png' or 'svg'. Default 'png'.

Outputs:
    <output_dir>/coefficient_path.<fmt>     coefficient trajectories plot
    <output_dir>/coefficient_path.csv        long-format path table

Example:
    python plot_coefficient_path.py --features data/discovery.csv \\
        --alpha 1.0 --top-n 12 --output-dir results/
"""

from __future__ import annotations

import argparse
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

try:
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
except ImportError:
    sys.exit("Install scikit-learn with: pip install scikit-learn")


def _maybe_import_matplotlib() -> Any:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        return plt
    except ImportError:
        sys.exit("Install matplotlib with: pip install matplotlib")


def _load_xy(
    path: Path,
    outcome_col: str | None,
    sample_col: str | None,
) -> tuple[pd.DataFrame, pd.Series]:
    df = pd.read_csv(path)
    sample_col = sample_col or str(df.columns[0])
    outcome_col = outcome_col or str(df.columns[-1])
    if outcome_col not in df.columns or sample_col not in df.columns:
        sys.exit("Sample or outcome column not present in features CSV.")
    y = df[outcome_col].astype(int)
    feature_cols = [c for c in df.columns if c not in (sample_col, outcome_col)]
    X = df[feature_cols].astype(float)
    return X, y


def _make_lambda_grid(X: np.ndarray, y: np.ndarray, n_lambdas: int) -> np.ndarray:
    """Make a logarithmic grid roughly matching glmnet's defaults.

    lambda_max is approximated as max |X^T (y - mean(y))| / n; the grid runs
    from lambda_max down to lambda_max * 1e-3 in log-spaced steps.
    """
    n = X.shape[0]
    y_centered = y - y.mean()
    lam_max = float(np.max(np.abs(X.T @ y_centered)) / n)
    lam_max = max(lam_max, 1e-3)
    return np.geomspace(lam_max, lam_max * 1e-3, n_lambdas)


def fit_path(
    X: pd.DataFrame,
    y: pd.Series,
    alpha: float,
    n_lambdas: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit the path. Returns (lambdas, coefs of shape (n_lambdas, n_features))."""
    X_arr = StandardScaler().fit_transform(X.to_numpy(dtype=float))
    y_arr = y.to_numpy(dtype=int)

    lambdas = _make_lambda_grid(X_arr, y_arr, n_lambdas)

    penalty = "l1" if alpha == 1.0 else "elasticnet"
    solver = "saga"
    l1_ratio: float | None = None if alpha == 1.0 else alpha

    coefs = np.zeros((len(lambdas), X_arr.shape[1]))
    for i, lam in enumerate(lambdas):
        kwargs: dict[str, Any] = {
            "penalty": penalty,
            "solver": solver,
            "C": 1.0 / max(lam * X_arr.shape[0], 1e-12),
            "max_iter": 5000,
            "fit_intercept": True,
        }
        if l1_ratio is not None:
            kwargs["l1_ratio"] = l1_ratio
        try:
            model = LogisticRegression(**kwargs)
            model.fit(X_arr, y_arr)
            coefs[i] = model.coef_[0]
        except ValueError:
            coefs[i] = 0.0

    return lambdas, coefs


def write_path_table(
    lambdas: np.ndarray,
    coefs: np.ndarray,
    feature_names: list[str],
    output_path: Path,
) -> None:
    rows: list[dict[str, Any]] = []
    for i, lam in enumerate(lambdas):
        for j, feat in enumerate(feature_names):
            rows.append(
                {
                    "lambda": float(lam),
                    "log_lambda": float(np.log(lam)),
                    "feature": feat,
                    "coefficient": float(coefs[i, j]),
                }
            )
    pd.DataFrame(rows).to_csv(output_path, index=False)


def plot_path(
    lambdas: np.ndarray,
    coefs: np.ndarray,
    feature_names: list[str],
    output_path: Path,
    top_n: int,
    highlight: list[str] | None,
    lambda_min: float | None,
    lambda_1se: float | None,
) -> None:
    plt = _maybe_import_matplotlib()
    log_lam = np.log(lambdas)

    if highlight:
        highlight_idx = [
            i for i, name in enumerate(feature_names) if name in set(highlight)
        ]
    elif top_n > 0:
        max_abs = np.max(np.abs(coefs), axis=0)
        highlight_idx = list(np.argsort(max_abs)[::-1][:top_n])
    else:
        highlight_idx = []

    fig, ax = plt.subplots(figsize=(9, 6))
    for j in range(coefs.shape[1]):
        if j in highlight_idx:
            continue
        ax.plot(log_lam, coefs[:, j], color="lightgrey", alpha=0.5, linewidth=0.8)

    palette = [
        "#E74C3C",
        "#3498DB",
        "#27AE60",
        "#F39C12",
        "#9B59B6",
        "#1ABC9C",
        "#E67E22",
        "#34495E",
        "#16A085",
        "#C0392B",
        "#2980B9",
        "#8E44AD",
    ]
    for k, j in enumerate(highlight_idx):
        color = palette[k % len(palette)]
        ax.plot(log_lam, coefs[:, j], color=color, linewidth=1.6,
                label=feature_names[j])
        idx_at_min = int(np.argmax(np.abs(coefs[:, j])))
        ax.text(
            log_lam[idx_at_min],
            coefs[idx_at_min, j],
            feature_names[j],
            fontsize=8,
            color=color,
        )

    if lambda_min is not None:
        ax.axvline(np.log(lambda_min), color="#E74C3C", linestyle="--",
                   label="lambda.min", alpha=0.6)
    if lambda_1se is not None:
        ax.axvline(np.log(lambda_1se), color="#27AE60", linestyle="--",
                   label="lambda.1se", alpha=0.6)

    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_xlabel("log(lambda)")
    ax.set_ylabel("coefficient")
    ax.set_title("LASSO coefficient regularization path")
    if highlight_idx:
        ax.legend(loc="best", fontsize=8, frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot the LASSO / elastic-net coefficient path."
    )
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--outcome-col", default=None)
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--n-lambdas", type=int, default=100)
    parser.add_argument("--top-n", type=int, default=10)
    parser.add_argument("--highlight", type=Path, default=None,
                        help="Optional text file: one feature name per line.")
    parser.add_argument("--lambda-min", type=float, default=None)
    parser.add_argument("--lambda-1se", type=float, default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--plot-format", default="png", choices=["png", "svg"])
    args = parser.parse_args()

    if not args.features.exists():
        sys.exit(f"Features file not found: {args.features}")
    if args.alpha <= 0 or args.alpha > 1:
        sys.exit("alpha must be in (0, 1].")

    X, y = _load_xy(args.features, args.outcome_col, args.sample_col)
    if set(int(v) for v in y.unique()) - {0, 1}:
        sys.exit("Outcome column must be binary 0/1.")

    highlight: list[str] | None = None
    if args.highlight is not None:
        if not args.highlight.exists():
            sys.exit(f"Highlight file not found: {args.highlight}")
        highlight = [
            line.strip()
            for line in args.highlight.read_text().splitlines()
            if line.strip()
        ]

    lambdas, coefs = fit_path(X, y, args.alpha, args.n_lambdas)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / "coefficient_path.csv"
    write_path_table(lambdas, coefs, list(X.columns), table_path)

    plot_path_file = args.output_dir / f"coefficient_path.{args.plot_format}"
    plot_path(
        lambdas=lambdas,
        coefs=coefs,
        feature_names=list(X.columns),
        output_path=plot_path_file,
        top_n=args.top_n,
        highlight=highlight,
        lambda_min=args.lambda_min,
        lambda_1se=args.lambda_1se,
    )

    nonzero_at_min = int(np.count_nonzero(coefs[-1]))
    print(
        "OK plot_coefficient_path: %d features over %d lambdas; %d nonzero at min lambda"
        % (len(X.columns), len(lambdas), nonzero_at_min)
    )
    print(f"  plot: {plot_path_file}")
    print(f"  table: {table_path}")


if __name__ == "__main__":
    main()
