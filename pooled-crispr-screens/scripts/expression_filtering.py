"""Filter genes by group-wise mean expression.

Drops genes whose mean expression in every (or any) sgRNA / perturbation group
falls below a threshold — pure noise reduction, useful before DE testing.

Example:
    python expression_filtering.py --in normalized.h5ad --out expressed.h5ad \\
        --group-key sgRNA --min-mean 0.5 --require any
"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

try:
    import numpy as np
    import pandas as pd
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install numpy pandas anndata")


def grouped_obs_mean(
    adata: "ad.AnnData",
    group_key: str,
    layer: Optional[str] = None,
) -> "pd.DataFrame":
    """Return a (gene × group) matrix of within-group mean expression."""
    get_X = (lambda x: x.layers[layer]) if layer else (lambda x: x.X)
    grouped = adata.obs.groupby(group_key, observed=True)
    out = pd.DataFrame(
        np.zeros((adata.shape[1], grouped.ngroups), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names,
    )
    for group, idx in grouped.indices.items():
        sub_X = get_X(adata[idx])
        out[group] = np.ravel(sub_X.mean(axis=0, dtype=np.float64))
    return out


def filter_by_expression(
    adata: "ad.AnnData",
    group_key: str = "sgRNA",
    min_mean_expression: float = 0.5,
    require_expression_in: str = "any",
) -> "ad.AnnData":
    """Keep genes that pass the expression threshold in 'any' or 'all' groups."""
    n_before = adata.n_vars
    print(
        f"Filtering genes by group-wise mean expression "
        f"(group_key={group_key}, min={min_mean_expression}, require={require_expression_in})"
    )
    mean_exp = grouped_obs_mean(adata, group_key)

    if require_expression_in == "any":
        keep = (mean_exp > min_mean_expression).any(axis=1)
    elif require_expression_in == "all":
        keep = (mean_exp > min_mean_expression).all(axis=1)
    else:
        raise ValueError(f"require_expression_in must be 'any' or 'all', got {require_expression_in!r}")

    adata = adata[:, keep].copy()
    retention = adata.n_vars / n_before * 100
    print(f"✓ Expression filtering complete: {adata.n_vars}/{n_before} genes ({retention:.1f}%)")
    return adata


def summarize_expression(
    adata: "ad.AnnData", group_key: str = "sgRNA"
) -> "pd.DataFrame":
    """Per-gene summary across groups."""
    mean_exp = grouped_obs_mean(adata, group_key)
    return pd.DataFrame({
        "mean_across_groups": mean_exp.mean(axis=1),
        "max_across_groups": mean_exp.max(axis=1),
        "min_across_groups": mean_exp.min(axis=1),
        "expressed_in_n_groups": (mean_exp > 0).sum(axis=1),
        "expressed_in_frac_groups": (mean_exp > 0).sum(axis=1) / mean_exp.shape[1],
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--group-key", default="sgRNA")
    parser.add_argument("--min-mean", type=float, default=0.5)
    parser.add_argument("--require", choices=["any", "all"], default="any")
    args = parser.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    adata = sc.read_h5ad(args.in_path)
    adata = filter_by_expression(
        adata, group_key=args.group_key,
        min_mean_expression=args.min_mean, require_expression_in=args.require,
    )
    adata.write(args.out)
    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
