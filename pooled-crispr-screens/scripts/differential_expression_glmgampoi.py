"""Python wrapper for the glmGamPoi (R) batch-corrected DE tier.

Writes the relevant subset of the AnnData to a temporary h5ad, calls the R
script via subprocess, and reads back one CSV per perturbation. Requires R,
glmGamPoi, and zellkonverter on the R side.

Example:
    python differential_expression_glmgampoi.py --in raw_counts.h5ad \\
        --perturbations TREM2,CD86,STAT5A --donor-col batch \\
        --out-dir results/glmgampoi/
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

try:
    import numpy as np
    import pandas as pd
    import scanpy as sc
except ImportError:
    sys.exit("Install with: pip install scanpy pandas numpy")


def run_glmgampoi_batch(
    adata,
    perturbations: List[str],
    control_group: str = "non-targeting",
    donor_col: str = "donor",
    gene_col: str = "gene",
    min_expression: int = 100,
    output_dir: str = "results/glmgampoi/",
    r_script_path: Optional[str] = None,
    verbose: bool = True,
) -> Dict[str, "pd.DataFrame"]:
    """Run glmGamPoi for the listed perturbations via R subprocess."""
    if donor_col not in adata.obs.columns:
        raise ValueError(f"Donor/batch column '{donor_col}' not in adata.obs")

    Path(output_dir).mkdir(parents=True, exist_ok=True)
    if r_script_path is None:
        r_script_path = Path(__file__).parent / "run_glmgampoi.R"
    if not Path(r_script_path).exists():
        raise FileNotFoundError(f"R script not found: {r_script_path}")

    if verbose:
        print(f"glmGamPoi for {len(perturbations)} perturbations")
        print(f"  control='{control_group}' donor_col='{donor_col}' gene_col='{gene_col}'")

    keep = perturbations + [control_group]
    sub = adata[adata.obs[gene_col].isin(keep)].copy()
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = tmp.name
        if verbose:
            print(f"  writing subset h5ad → {tmp_path}")
        sub.write(tmp_path)

    try:
        cmd = [
            "Rscript", str(r_script_path), tmp_path,
            ",".join(perturbations), control_group, gene_col, donor_col,
            str(min_expression), output_dir,
        ]
        if verbose:
            print(f"  $ {' '.join(cmd)}")
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        if verbose and proc.stdout:
            print(proc.stdout)
        if proc.stderr:
            print(proc.stderr, file=sys.stderr)
    except subprocess.CalledProcessError as exc:
        print(f"R subprocess failed (returncode={exc.returncode})", file=sys.stderr)
        print(exc.stdout, file=sys.stderr)
        print(exc.stderr, file=sys.stderr)
        raise
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    results: Dict[str, pd.DataFrame] = {}
    for gene in perturbations:
        path = Path(output_dir) / f"{gene}_glmgampoi.csv"
        if path.exists():
            results[gene] = pd.read_csv(path)
            if verbose:
                n_sig = int((results[gene]["adj_pval"] < 0.05).sum())
                print(f"  {gene}: {n_sig} significant genes (FDR<0.05)")
        elif verbose:
            print(f"  {gene}: no result file (may have failed)")
    print(f"✓ glmGamPoi complete: {len(results)}/{len(perturbations)} perturbations")
    return results


def compare_ttest_vs_glmgampoi(
    ttest_results: Dict[str, "pd.DataFrame"],
    glmgampoi_results: Dict[str, "pd.DataFrame"],
    fdr_threshold: float = 0.05,
    output_file: str = "results/method_comparison.csv",
) -> "pd.DataFrame":
    """Per-perturbation overlap of significant-gene sets between the two tiers."""
    rows = []
    for gene, glm in glmgampoi_results.items():
        if gene not in ttest_results:
            continue
        tt = ttest_results[gene]
        sig_tt = set(tt[tt["qval"] < fdr_threshold]["gene"])
        sig_glm = set(glm[glm["adj_pval"] < fdr_threshold]["gene"])
        overlap = len(sig_tt & sig_glm)
        agree = overlap / max(len(sig_tt), len(sig_glm)) if (sig_tt or sig_glm) else np.nan
        rows.append({
            "perturbation": gene,
            "n_sig_ttest": len(sig_tt),
            "n_sig_glmgampoi": len(sig_glm),
            "n_overlap": overlap,
            "agreement": agree,
        })
    out = pd.DataFrame(rows)
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output_file, index=False)
    print(
        f"  median DE genes: t-test={out['n_sig_ttest'].median():.0f} "
        f"glmGamPoi={out['n_sig_glmgampoi'].median():.0f} "
        f"agreement={out['agreement'].mean():.1%}"
    )
    print(f"✓ Wrote {output_file}")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True,
                        help="h5ad with RAW counts (not log-normalized)")
    parser.add_argument("--perturbations", required=True,
                        help="Comma-separated list of perturbation names")
    parser.add_argument("--control", default="non-targeting")
    parser.add_argument("--donor-col", default="donor")
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--min-expression", type=int, default=100)
    parser.add_argument("--out-dir", default="results/glmgampoi/")
    parser.add_argument("--r-script", default=None)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.in_path)
    perturbations = [p.strip() for p in args.perturbations.split(",") if p.strip()]
    run_glmgampoi_batch(
        adata, perturbations,
        control_group=args.control,
        donor_col=args.donor_col,
        gene_col=args.gene_col,
        min_expression=args.min_expression,
        output_dir=args.out_dir,
        r_script_path=args.r_script,
    )


if __name__ == "__main__":
    main()
