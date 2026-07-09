"""Stratified k-fold cross-validated LASSO / elastic-net fit.

A standalone Python entry point for fitting a penalized logistic regression
on a samples-x-features matrix with a binary outcome. Uses scikit-learn's
LogisticRegressionCV with an L1 (LASSO) or elastic-net penalty, evaluated
under stratified k-fold cross-validation. Writes the fitted coefficients,
the per-fold AUCs, and a small JSON manifest with the chosen penalty
strength so downstream scripts can re-use the model.

Inputs:
    --features PATH    CSV with samples in rows, features in columns. The
                       first column is the sample id and the last column is
                       the binary outcome (0/1). Override with --outcome-col.
    --outcome-col STR  Name of the outcome column (default: last column).
    --sample-col STR   Name of the sample id column (default: first column).
    --alpha FLOAT      Elastic-net mixing (1.0 = pure LASSO, default 1.0,
                       0 < alpha <= 1; alpha = 0 = ridge is not allowed here
                       to keep selection behaviour).
    --n-folds INT      Number of stratified CV folds (default 10).
    --n-lambdas INT    Number of regularization strengths to search
                       (default 50).
    --scoring STR      Scoring metric (default roc_auc).
    --standardize      Standardize features before fitting (default on).
    --no-standardize   Disable standardization.
    --seed INT         Random seed (default 42).
    --output-dir DIR   Directory for outputs (default ./results).

Outputs:
    <output_dir>/lasso_cv_coefficients.csv   coefficient per feature
    <output_dir>/lasso_cv_fold_aucs.csv       per-fold AUC table
    <output_dir>/lasso_cv_manifest.json       parameters and chosen lambda

Example:
    python fit_lasso_cv.py --features data/discovery.csv \\
        --alpha 0.5 --n-folds 10 --output-dir results/
"""

from __future__ import annotations

import argparse
import json
import os
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
    from sklearn.linear_model import LogisticRegressionCV
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.preprocessing import StandardScaler
except ImportError:
    sys.exit("Install scikit-learn with: pip install scikit-learn")


def _load_xy(
    path: Path,
    outcome_col: str | None,
    sample_col: str | None,
) -> tuple[pd.DataFrame, pd.Series, list[str]]:
    """Return (X, y, sample_ids) from a flat CSV."""
    df = pd.read_csv(path)
    if df.shape[1] < 3:
        sys.exit("Feature CSV must have at least sample id + 1 feature + outcome.")

    sample_col = sample_col or str(df.columns[0])
    outcome_col = outcome_col or str(df.columns[-1])
    if sample_col not in df.columns:
        sys.exit(f"Sample column '{sample_col}' not found.")
    if outcome_col not in df.columns:
        sys.exit(f"Outcome column '{outcome_col}' not found.")

    sample_ids = df[sample_col].astype(str).tolist()
    y = df[outcome_col].astype(int)
    feature_cols = [c for c in df.columns if c not in (sample_col, outcome_col)]
    X = df[feature_cols].astype(float)
    return X, y, sample_ids


def fit_lasso(
    X: pd.DataFrame,
    y: pd.Series,
    alpha: float,
    n_folds: int,
    n_lambdas: int,
    scoring: str,
    standardize: bool,
    seed: int,
) -> dict[str, Any]:
    """Fit LogisticRegressionCV and capture per-fold AUCs.

    For alpha == 1.0 the penalty is L1; for 0 < alpha < 1 we use elasticnet
    with the saga solver. Returns coefficients, the chosen C, per-fold AUCs.
    """
    if alpha <= 0 or alpha > 1:
        sys.exit("alpha must be in (0, 1]; alpha=0 (ridge) is not supported.")

    X_arr = X.to_numpy(dtype=float)
    y_arr = y.to_numpy(dtype=int)

    if standardize:
        scaler = StandardScaler()
        X_arr = scaler.fit_transform(X_arr)
    else:
        scaler = None

    penalty = "l1" if alpha == 1.0 else "elasticnet"
    solver = "saga"
    l1_ratio = None if penalty == "l1" else alpha

    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)

    print(
        "Fitting LogisticRegressionCV: alpha=%.2f penalty=%s n_folds=%d"
        % (alpha, penalty, n_folds)
    )

    extra: dict[str, Any] = {}
    if l1_ratio is not None:
        extra["l1_ratios"] = [l1_ratio]

    model = LogisticRegressionCV(
        Cs=n_lambdas,
        cv=cv,
        penalty=penalty,
        solver=solver,
        scoring=scoring,
        max_iter=5000,
        random_state=seed,
        n_jobs=1,
        refit=True,
        **extra,
    )
    model.fit(X_arr, y_arr)

    chosen_C = float(model.C_[0])
    coefs = model.coef_[0]
    intercept = float(model.intercept_[0])

    fold_aucs: list[float] = []
    for fold_i, (train_idx, test_idx) in enumerate(cv.split(X_arr, y_arr), start=1):
        sub = LogisticRegressionCV(
            Cs=[chosen_C],
            cv=2,
            penalty=penalty,
            solver=solver,
            scoring=scoring,
            max_iter=5000,
            random_state=seed,
            refit=True,
            **extra,
        )
        sub.fit(X_arr[train_idx], y_arr[train_idx])
        probs = sub.predict_proba(X_arr[test_idx])[:, 1]
        auc = float(roc_auc_score(y_arr[test_idx], probs))
        fold_aucs.append(auc)
        print(f"  fold {fold_i}/{n_folds}  AUC={auc:.3f}")

    nonzero = int(np.count_nonzero(coefs))
    print(f"chosen C: {chosen_C:.4g}  nonzero coefficients: {nonzero}")

    return {
        "coefficients": coefs,
        "intercept": intercept,
        "chosen_C": chosen_C,
        "fold_aucs": fold_aucs,
        "scaler_mean": scaler.mean_.tolist() if scaler is not None else None,
        "scaler_scale": scaler.scale_.tolist() if scaler is not None else None,
        "penalty": penalty,
        "l1_ratio": l1_ratio,
    }


def write_outputs(
    fit: dict[str, Any],
    feature_names: list[str],
    output_dir: Path,
    alpha: float,
    n_folds: int,
    seed: int,
    standardize: bool,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    coef_df = pd.DataFrame(
        {
            "feature": feature_names,
            "coefficient": fit["coefficients"],
            "abs_coefficient": np.abs(fit["coefficients"]),
            "nonzero": fit["coefficients"] != 0,
        }
    ).sort_values("abs_coefficient", ascending=False)
    coef_df.to_csv(output_dir / "lasso_cv_coefficients.csv", index=False)

    fold_df = pd.DataFrame(
        {
            "fold": list(range(1, len(fit["fold_aucs"]) + 1)),
            "auc": fit["fold_aucs"],
        }
    )
    fold_df.to_csv(output_dir / "lasso_cv_fold_aucs.csv", index=False)

    manifest = {
        "alpha": alpha,
        "n_folds": n_folds,
        "seed": seed,
        "standardize": standardize,
        "penalty": fit["penalty"],
        "l1_ratio": fit["l1_ratio"],
        "chosen_C": fit["chosen_C"],
        "intercept": fit["intercept"],
        "mean_fold_auc": float(np.mean(fit["fold_aucs"])),
        "sd_fold_auc": float(np.std(fit["fold_aucs"], ddof=1))
        if len(fit["fold_aucs"]) > 1
        else 0.0,
        "n_nonzero": int(np.count_nonzero(fit["coefficients"])),
        "n_features": len(feature_names),
        "scaler_mean": fit["scaler_mean"],
        "scaler_scale": fit["scaler_scale"],
        "feature_names": feature_names,
    }
    with open(output_dir / "lasso_cv_manifest.json", "w") as fh:
        json.dump(manifest, fh, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Stratified k-fold cross-validated LASSO / elastic-net fit."
    )
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--outcome-col", default=None)
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--alpha", type=float, default=1.0)
    parser.add_argument("--n-folds", type=int, default=10)
    parser.add_argument("--n-lambdas", type=int, default=50)
    parser.add_argument("--scoring", default="roc_auc")
    standardize_group = parser.add_mutually_exclusive_group()
    standardize_group.add_argument(
        "--standardize", dest="standardize", action="store_true", default=True
    )
    standardize_group.add_argument(
        "--no-standardize", dest="standardize", action="store_false"
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    if not args.features.exists():
        sys.exit(f"Features file not found: {args.features}")

    X, y, _ = _load_xy(args.features, args.outcome_col, args.sample_col)

    if set(int(v) for v in y.unique()) - {0, 1}:
        sys.exit("Outcome column must be binary 0/1.")

    fit = fit_lasso(
        X=X,
        y=y,
        alpha=args.alpha,
        n_folds=args.n_folds,
        n_lambdas=args.n_lambdas,
        scoring=args.scoring,
        standardize=args.standardize,
        seed=args.seed,
    )

    write_outputs(
        fit=fit,
        feature_names=list(X.columns),
        output_dir=args.output_dir,
        alpha=args.alpha,
        n_folds=args.n_folds,
        seed=args.seed,
        standardize=args.standardize,
    )

    mean_auc = float(np.mean(fit["fold_aucs"]))
    print(
        "OK fit_lasso_cv: %d features, %d nonzero, mean CV AUC = %.3f"
        % (X.shape[1], int(np.count_nonzero(fit["coefficients"])), mean_auc)
    )
    print(f"  outputs in {os.fspath(args.output_dir)}")


if __name__ == "__main__":
    main()
