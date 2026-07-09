"""Reconcile gene symbols between sgRNA library and expression matrix.

Detects symbols in ``obs['gene']`` that have no match in ``var_names`` (usually
because the library was designed against older HGNC symbols), suggests
corrections from a built-in alias map or optionally from mygene, and applies
them.

Example:
    python gene_name_corrections.py --in mapped.h5ad --out corrected.h5ad --use-synonyms
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, List

try:
    import anndata as ad
except ImportError:
    sys.exit("Install with: pip install anndata")


_KNOWN_ALIASES: Dict[str, str] = {
    "TMEM55A": "PIP4P2",
    "TMEM55B": "PIP4P1",
    "ATP5C1": "ATP5F1C",
    "ATP5H": "ATP5PD",
    "ATP5A1": "ATP5F1A",
    "ATP5B": "ATP5F1B",
    "ATP5L": "ATP5MG",
    "ATP5O": "ATP5PO",
    "C9orf72": "C9ORF72",
    "HIST1H2BC": "H2BC4",
    "HIST1H4C": "H4C3",
}


def detect_mismatches(adata: "ad.AnnData", gene_col: str = "gene") -> List[str]:
    """Return symbols present in obs[gene_col] but absent from var_names."""
    if gene_col not in adata.obs.columns:
        print(f"⚠ '{gene_col}' not found in adata.obs")
        return []
    in_obs = set(adata.obs[gene_col].unique())
    in_var = set(adata.var_names)
    missing = sorted(in_obs - in_var)
    if missing:
        print(f"Found {len(missing)} symbol mismatches:")
        for g in missing:
            print(f"  {g}")
    else:
        print("No symbol mismatches detected")
    return missing


def correct_gene_names(
    adata: "ad.AnnData",
    corrections: Dict[str, str],
    gene_col: str = "gene",
) -> "ad.AnnData":
    """Apply a mapping of old → new symbols to obs[gene_col]."""
    if gene_col not in adata.obs.columns:
        print(f"⚠ '{gene_col}' not found in adata.obs — no corrections applied")
        return adata
    if not corrections:
        print("No corrections requested")
        return adata

    print(f"Applying {len(corrections)} symbol corrections:")
    series = adata.obs[gene_col].copy()
    for old, new in corrections.items():
        n = int((series == old).sum())
        if n > 0:
            series = series.replace(old, new)
            print(f"  {old} → {new}  ({n} cells)")
        else:
            print(f"  {old} → {new}  (not present)")
    adata.obs[gene_col] = series

    remaining = detect_mismatches(adata, gene_col)
    if remaining:
        print(f"⚠ {len(remaining)} mismatches remain")
    else:
        print("✓ Gene name corrections applied")
    return adata


def suggest_corrections(
    adata: "ad.AnnData",
    gene_col: str = "gene",
    use_synonyms: bool = False,
) -> Dict[str, str]:
    """Suggest a corrections dict from the built-in alias map and (optionally) mygene."""
    missing = detect_mismatches(adata, gene_col)
    if not missing:
        return {}
    suggestions: Dict[str, str] = {}
    for sym in missing:
        if sym in _KNOWN_ALIASES and _KNOWN_ALIASES[sym] in adata.var_names:
            suggestions[sym] = _KNOWN_ALIASES[sym]
            print(f"  suggest {sym} → {_KNOWN_ALIASES[sym]}  (built-in)")

    if use_synonyms and len(suggestions) < len(missing):
        try:
            import mygene
        except ImportError:
            print("Install mygene for synonym lookup: pip install mygene")
            return suggestions
        mg = mygene.MyGeneInfo()
        unresolved = [g for g in missing if g not in suggestions]
        results = mg.querymany(
            unresolved, scopes="symbol,alias", fields="symbol", species="human"
        )
        for r in results:
            if "symbol" in r and r["symbol"] in adata.var_names:
                suggestions[r["query"]] = r["symbol"]
                print(f"  suggest {r['query']} → {r['symbol']}  (mygene)")
    return suggestions


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--in", dest="in_path", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--gene-col", default="gene")
    parser.add_argument("--use-synonyms", action="store_true",
                        help="Query mygene for additional synonyms (requires `pip install mygene`)")
    args = parser.parse_args()

    try:
        import scanpy as sc
    except ImportError:
        sys.exit("Install with: pip install scanpy")

    adata = sc.read_h5ad(args.in_path)
    suggested = suggest_corrections(adata, args.gene_col, args.use_synonyms)
    adata = correct_gene_names(adata, suggested, args.gene_col)
    adata.write(args.out)
    print(f"✓ Wrote {args.out}")


if __name__ == "__main__":
    main()
