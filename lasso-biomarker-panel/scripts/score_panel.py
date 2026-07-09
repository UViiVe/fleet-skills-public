"""Apply a locked biomarker panel to new samples and emit classifications.

Loads a panel manifest (coefficients + intercept + standardization stats)
produced by `fit_lasso_cv.py`, projects it onto a new feature matrix, and
writes predicted probabilities plus a threshold-based class call. Designed
to be the deployment-time companion of the discovery workflow: once a
panel is locked, the only thing a clinical user needs is this script and
the manifest JSON.

The threshold can be supplied directly, or set to one of three named
strategies:

    'youden'  -- Youden's J (sens + spec - 1) maximum from a labelled
                 calibration set (requires --calibration-features).
    'prevalence' -- match the positive-class prevalence of a labelled
                 calibration set (requires --calibration-features).
    'fixed-specificity' -- pick the threshold that gives at least
                 --target-specificity specificity on calibration (requires
                 --calibration-features and --target-specificity).

Inputs:
    --manifest PATH         lasso_cv_manifest.json from fit_lasso_cv.py.
    --features PATH         CSV of samples x features to score. Sample id
                            in the first column.
    --sample-col STR        Sample id column name. Default: first column.
    --threshold FLOAT       Decision threshold on the probability scale.
                            Default 0.5. Mutually exclusive with --threshold-method.
    --threshold-method STR  'youden' | 'prevalence' | 'fixed-specificity'.
    --target-specificity F  Required for 'fixed-specificity' (default 0.9).
    --calibration-features PATH  CSV with the same columns + an outcome
                                 column, used to derive a non-fixed threshold.
    --calibration-outcome-col STR  Outcome column in calibration CSV.
    --output-dir DIR        Output directory. Default ./results.

Outputs:
    <output_dir>/panel_scores.csv  sample_id, predicted_prob, predicted_class
    <output_dir>/panel_score_summary.json  threshold, n calls, etc.

Example:
    python score_panel.py --manifest results/lasso_cv_manifest.json \\
        --features data/new_cohort.csv --threshold 0.4 --output-dir results/
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


def _load_manifest(path: Path) -> dict[str, Any]:
    if not path.exists():
        sys.exit(f"Manifest not found: {path}")
    with open(path) as fh:
        manifest = json.load(fh)
    required = {"feature_names", "intercept"}
    if missing := required - set(manifest):
        sys.exit(f"Manifest missing keys: {sorted(missing)}")
    if "coefficients" not in manifest:
        # fit_lasso_cv.py writes coefficients to a sibling CSV; lift them
        # from there if the manifest doesn't carry them.
        coef_csv = path.with_name("lasso_cv_coefficients.csv")
        if not coef_csv.exists():
            sys.exit(
                "Manifest has no 'coefficients' and no companion "
                f"{coef_csv.name} was found."
            )
        coef_df = pd.read_csv(coef_csv)
        coef_map = dict(zip(coef_df["feature"], coef_df["coefficient"]))
        manifest["coefficients"] = [
            float(coef_map.get(f, 0.0)) for f in manifest["feature_names"]
        ]
    return manifest


def _load_features(
    path: Path,
    sample_col: str | None,
    feature_names: list[str],
) -> tuple[np.ndarray, list[str], list[str]]:
    if not path.exists():
        sys.exit(f"Features file not found: {path}")
    df = pd.read_csv(path)
    sample_col = sample_col or str(df.columns[0])
    if sample_col not in df.columns:
        sys.exit(f"Sample column '{sample_col}' not found.")

    sample_ids = df[sample_col].astype(str).tolist()
    missing = [f for f in feature_names if f not in df.columns]
    available = [f for f in feature_names if f in df.columns]

    if not available:
        sys.exit("None of the panel features are present in the input matrix.")

    X = np.zeros((len(df), len(feature_names)))
    for j, name in enumerate(feature_names):
        if name in df.columns:
            X[:, j] = df[name].astype(float).to_numpy()
        # missing features stay at 0, which under standardization == mean

    return X, sample_ids, missing


def _standardize_with(
    X: np.ndarray,
    scaler_mean: list[float] | None,
    scaler_scale: list[float] | None,
) -> np.ndarray:
    if scaler_mean is None or scaler_scale is None:
        return X
    mean = np.asarray(scaler_mean)
    scale = np.asarray(scaler_scale)
    scale = np.where(scale == 0, 1.0, scale)
    return (X - mean) / scale


def _predict_proba(
    X: np.ndarray,
    coefs: np.ndarray,
    intercept: float,
) -> np.ndarray:
    z = X @ coefs + intercept
    # numerically stable logistic
    out = np.where(
        z >= 0,
        1.0 / (1.0 + np.exp(-z)),
        np.exp(z) / (1.0 + np.exp(z)),
    )
    return out


def _derive_threshold(
    method: str,
    target_spec: float,
    cal_features: Path,
    cal_outcome_col: str | None,
    sample_col: str | None,
    feature_names: list[str],
    scaler_mean: list[float] | None,
    scaler_scale: list[float] | None,
    coefs: np.ndarray,
    intercept: float,
) -> tuple[float, dict[str, Any]]:
    df = pd.read_csv(cal_features)
    sample_col = sample_col or str(df.columns[0])
    cal_outcome_col = cal_outcome_col or str(df.columns[-1])
    if cal_outcome_col not in df.columns:
        sys.exit(f"Calibration outcome column '{cal_outcome_col}' not present.")
    y = df[cal_outcome_col].astype(int).to_numpy()
    if set(np.unique(y).tolist()) - {0, 1}:
        sys.exit("Calibration outcome must be 0/1.")

    X = np.zeros((len(df), len(feature_names)))
    for j, name in enumerate(feature_names):
        if name in df.columns:
            X[:, j] = df[name].astype(float).to_numpy()
    Xs = _standardize_with(X, scaler_mean, scaler_scale)
    probs = _predict_proba(Xs, coefs, intercept)

    if method == "prevalence":
        prev = float(y.mean())
        thr = float(np.quantile(probs, 1.0 - prev))
        return thr, {"strategy": "prevalence", "prevalence": prev}

    sort_idx = np.argsort(probs)[::-1]
    sorted_probs = probs[sort_idx]
    sorted_y = y[sort_idx]

    tps = np.cumsum(sorted_y == 1)
    fps = np.cumsum(sorted_y == 0)
    P = max(int((y == 1).sum()), 1)
    N = max(int((y == 0).sum()), 1)
    sens = tps / P
    spec = 1 - (fps / N)

    if method == "youden":
        j_stat = sens + spec - 1
        best = int(np.argmax(j_stat))
        return float(sorted_probs[best]), {
            "strategy": "youden",
            "sensitivity": float(sens[best]),
            "specificity": float(spec[best]),
        }
    if method == "fixed-specificity":
        eligible = np.where(spec >= target_spec)[0]
        if eligible.size == 0:
            sys.exit(
                f"No threshold reaches target specificity {target_spec}; "
                "consider lowering --target-specificity."
            )
        idx = int(eligible[np.argmax(sens[eligible])])
        return float(sorted_probs[idx]), {
            "strategy": "fixed-specificity",
            "target_specificity": target_spec,
            "achieved_specificity": float(spec[idx]),
            "sensitivity": float(sens[idx]),
        }
    sys.exit(f"Unknown threshold method: {method}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply a locked LASSO biomarker panel to new samples."
    )
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--features", required=True, type=Path)
    parser.add_argument("--sample-col", default=None)
    parser.add_argument("--threshold", type=float, default=None)
    parser.add_argument(
        "--threshold-method",
        default=None,
        choices=["youden", "prevalence", "fixed-specificity"],
    )
    parser.add_argument("--target-specificity", type=float, default=0.9)
    parser.add_argument("--calibration-features", type=Path, default=None)
    parser.add_argument("--calibration-outcome-col", default=None)
    parser.add_argument("--output-dir", type=Path, default=Path("results"))
    args = parser.parse_args()

    manifest = _load_manifest(args.manifest)
    feature_names: list[str] = list(manifest["feature_names"])
    coefs = np.asarray(manifest["coefficients"], dtype=float)
    intercept = float(manifest["intercept"])
    scaler_mean = manifest.get("scaler_mean")
    scaler_scale = manifest.get("scaler_scale")

    X, sample_ids, missing = _load_features(
        args.features, args.sample_col, feature_names
    )
    if missing:
        print(
            f"WARNING: {len(missing)} / {len(feature_names)} panel features missing "
            f"from input; imputing as mean (post-standardization zero)."
        )

    Xs = _standardize_with(X, scaler_mean, scaler_scale)
    probs = _predict_proba(Xs, coefs, intercept)

    threshold_info: dict[str, Any] = {}
    if args.threshold is not None and args.threshold_method is not None:
        sys.exit("Provide either --threshold or --threshold-method, not both.")
    if args.threshold_method is not None:
        if args.calibration_features is None:
            sys.exit(
                "--threshold-method requires --calibration-features to derive it."
            )
        threshold, threshold_info = _derive_threshold(
            method=args.threshold_method,
            target_spec=args.target_specificity,
            cal_features=args.calibration_features,
            cal_outcome_col=args.calibration_outcome_col,
            sample_col=args.sample_col,
            feature_names=feature_names,
            scaler_mean=scaler_mean,
            scaler_scale=scaler_scale,
            coefs=coefs,
            intercept=intercept,
        )
    elif args.threshold is None:
        threshold = 0.5
        threshold_info = {"strategy": "default"}
    else:
        threshold = args.threshold
        threshold_info = {"strategy": "user-fixed"}

    predicted_class = (probs >= threshold).astype(int)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scores_df = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "predicted_prob": probs,
            "predicted_class": predicted_class,
        }
    )
    scores_path = args.output_dir / "panel_scores.csv"
    scores_df.to_csv(scores_path, index=False)

    summary = {
        "threshold": float(threshold),
        "threshold_info": threshold_info,
        "n_samples": int(len(sample_ids)),
        "n_positive_calls": int(predicted_class.sum()),
        "n_negative_calls": int((1 - predicted_class).sum()),
        "n_panel_features": len(feature_names),
        "n_features_missing": len(missing),
        "missing_features": missing,
    }
    summary_path = args.output_dir / "panel_score_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(
        "OK score_panel: %d samples scored, %d called positive at threshold %.3f"
        % (len(sample_ids), int(predicted_class.sum()), threshold)
    )
    print(f"  scores: {scores_path}")


if __name__ == "__main__":
    main()
