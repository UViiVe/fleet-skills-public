"""Variable-stability analysis via bootstrap resampling of a LASSO fit.

Refits a penalized logistic regression on B bootstrap samples (sampled with
replacement) of the discovery cohort, and counts how often each feature is
selected (non-zero coefficient). Features with high selection probability
are robust to sampling fluctuations and form the "stable" biomarker panel.
This is a Meinshausen & Buhlmann (2010) -style stability check; the
companion R workflow already does repeated-CV stability, this Python
script provides bootstrap stability for cohorts where bootstrap is more
natural (e.g. case-control with very imbalanced groups).

Inputs:
    --features PATH        CSV: samples x features + outcome.
    --outcome-col STR      Outcome column name. Default: last column.
    --sample-col STR       Sample id column. Default: first column.
    --alpha FLOAT          Elastic-net mixing (1.0 = pure LASSO). Default 1.0.
    --C FLOAT              Inverse regularization strength (sklearn parameter).
                           If not given, picked via 5-fold CV on each bootstrap.
    --n-bootstraps INT     Number of bootstrap resamples (default 200).
    --stability-threshold  Min selection frequency to call a feature stable
                           (default 0.6 — typical Meinshausen-Buhlmann value).
    --stratify             Stratify each bootstrap by outcome (default on).
    --no-stratify          Sample without stratification.
    --seed INT             Random seed (default 42).
    --output-dir DIR       Output directory. Default ./results.
    --n-jobs INT           Parallel bootstraps (default 1).

Outputs:
    <output_dir>/bootstrap_stability.csv  feature x [selection_frequency,
                                          mean_coef, sd_coef, is_stable]
    <output_dir>/bootstrap_summary.json   parameters + stable-feature list

Example:
    python bootstrap_stability.py --features data/discovery.csv \\
        --n-bootstraps 200 --stability-threshold 0.6 --output-dir results/
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
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
    from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
    from sklearn.preprocessing import StandardScaler
except ImportError:
    sys.exit("Install scikit-learn with: pip install scikit-learn")


def _load_xy(
    path: Path,
    outcome_col: str | None,
    sample_col: str | None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    df = pd.read_csv(path)
    sample_col = sample_col or str(df.columns[0])
    outcome_col = outcome_col or str(df.columns[-1])
    if outcome_col not in df.columns or sample_col not in df.columns:
        sys.exit("Sample or outcome column not present in features CSV.")
    y = df[outcome_col].astype(int).to_numpy()
    feature_cols = [c for c in df.columns if c not in (sample_col, outcome_col)]
    X = df[feature_cols].astype(float).to_numpy()
    return X, y, feature_cols


def _bootstrap_indices(y: np.ndarray, stratify: bool, rng: np.random.Generator) -> np.ndarray:
    n = y.shape[0]
    if not stratify:
        return rng.integers(0, n, n)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    pos_sample = rng.choice(pos, size=len(pos), replace=True)
    neg_sample = rng.choice(neg, size=len(neg), replace=True)
    return np.concatenate([pos_sample, neg_sample])


def _fit_one(
    X: np.ndarray,
    y: np.ndarray,
    alpha: float,
    C: float | None,
    seed: int,
    stratify: bool,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    idx = _bootstrap_indices(y, stratify, rng)
    Xb = X[idx]
    yb = y[idx]

    if len(np.unique(yb)) < 2:
        return np.zeros(X.shape[1])

    Xb = StandardScaler().fit_transform(Xb)

    penalty = "l1" if alpha == 1.0 else "elasticnet"
    l1_ratio: float | None = None if alpha == 1.0 else alpha

    if C is None:
        kwargs: dict[str, Any] = {
            "Cs": 20,
            "cv": 5,
            "penalty": penalty,
            "solver": "saga",
            "scoring": "roc_auc",
            "max_iter": 5000,
            "random_state": seed,
            "n_jobs": 1,
            "refit": True,
        }
        if l1_ratio is not None:
            kwargs["l1_ratios"] = [l1_ratio]
        model = LogisticRegressionCV(**kwargs)
    else:
        kwargs = {
            "C": C,
            "penalty": penalty,
            "solver": "saga",
            "max_iter": 5000,
            "random_state": seed,
            "fit_intercept": True,
        }
        if l1_ratio is not None:
            kwargs["l1_ratio"] = l1_ratio
        model = LogisticRegression(**kwargs)

    try:
        model.fit(Xb, yb)
        return model.coef_[0]
    except (ValueError, RuntimeError):
        return np.zeros(X.shape[1])


def run_bootstrap(
    X: np.ndarray,
    y: np.ndarray,
    alpha: float,
    C: float | None,
    n_bootstraps: int,
    stratify: bool,
    seed: int,
    n_jobs: int,
) -> np.ndarray:
    """Returns coef matrix shape (n_bootstraps, n_features)."""
    coefs = np.zeros((n_bootstraps, X.shape[1]))

    if n_jobs <= 1:
        for b in range(n_bootstraps):
            coefs[b] = _fit_one(X, y, alpha, C, seed + b, stratify)
            if (b + 1) % max(1, n_bootstraps // 10) == 0:
                print(f"  bootstrap {b + 1}/{n_bootstraps}")
        return coefs

    with ProcessPoolExecutor(max_workers=n_jobs) as pool:
        futures = {
            pool.submit(_fit_one, X, y, alpha, C, seed + b, stratify): b
            for b in range(n_bootstraps)
        }
        done = 0
        for fut in as_completed(futures):
            b = futures[fut]
            coefs[b] = fut.result()
            done += 1
            if done % max(1, n_bootstraps // 10) == 0:
                print(f"  bootstrap {done}/{n_bootstraps}")
    return coefs


def summarize(
    coefs: np.ndarray,
    feature_names: list[str],
    stability_threshold: float,
) -> tuple[pd.DataFrame, list[str]]:
    selected = coefs != 0
    selection_freq = selected.mean(axis=0)
    mean_coef = np.where(selected.any(axis=0), coefs.sum(axis=0) / np.maximum(selected.sum(axis=0), 1), 0.0)
    sd_coef = np.std(coefs, axis=0, ddof=1) if coefs.shape[0] > 1 else np.zeros_like(selection_freq)
    is_stable = selection_freq >= stability_threshold

    df = pd.DataFrame(
        {
            "feature": feature_names,
            "selection_frequency": selection_freq,
            "mean_coefficient": mean_coef,
            "sd_coefficient": sd_coef,
            "is_stable": is_stable,
        }
    ).sort_values("selection_frequency", ascending=False)
    stable = df.loc[df["is_stable"], "feature"].tolist()
    return df, stable


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bootstrap stability selection for LASSO biomarker panels."
    )
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--outcome-col", default=None)
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--C", type=float, default=None)
    parser.add_argument("--n-bootstraps", type=int, default=200)
    parser.add_argument("--stability-threshold", type=float, default=0.6)
    strat_group = parser.add_mutually_exclusive_group()
    strat_group.add_argument(
        "--stratify", dest="stratify", action="store_true", default=True
    )
    strat_group.add_argument(
        "--no-stratify", dest="stratify", action="store_false"
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    if not args.features.exists():
        sys.exit(f"Features file not found: {args.features}")
    if args.alpha <= 0 or args.alpha > 1:
        sys.exit("alpha must be in (0, 1].")
    if not 0 < args.stability_threshold <= 1:
        sys.exit("stability-threshold must be in (0, 1].")

    X, y, feature_names = _load_xy(
        args.features, args.outcome_col, args.sample_col
    )

    print(
        "Running bootstrap stability: %d samples x %d features, B=%d, alpha=%.2f"
        % (X.shape[0], X.shape[1], args.n_bootstraps, args.alpha)
    )

    coefs = run_bootstrap(
        X=X,
        y=y,
        alpha=args.alpha,
        C=args.C,
        n_bootstraps=args.n_bootstraps,
        stratify=args.stratify,
        seed=args.seed,
        n_jobs=args.n_jobs,
    )

    df, stable = summarize(coefs, feature_names, args.stability_threshold)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_path = args.output_dir / "bootstrap_stability.csv"
    df.to_csv(table_path, index=False)

    summary = {
        "n_bootstraps": args.n_bootstraps,
        "alpha": args.alpha,
        "C": args.C,
        "stability_threshold": args.stability_threshold,
        "stratify": args.stratify,
        "seed": args.seed,
        "n_features": len(feature_names),
        "n_stable": len(stable),
        "stable_features": stable,
    }
    json_path = args.output_dir / "bootstrap_summary.json"
    with open(json_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(
        "OK bootstrap_stability: %d / %d features stable (>= %.2f selection frequency)"
        % (len(stable), len(feature_names), args.stability_threshold)
    )
    print(f"  table: {table_path}")
    print(f"  summary: {json_path}")


if __name__ == "__main__":
    main()
