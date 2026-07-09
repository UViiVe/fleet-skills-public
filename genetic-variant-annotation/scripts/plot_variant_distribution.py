#!/usr/bin/env python3
"""Five standard visualisations for a parsed variant-annotation table.

Public entry points (all save a PNG + an SVG):
  - ``plot_consequence_distribution(df, output_file)``
  - ``plot_impact_by_chromosome(df, output_file)``
  - ``plot_pathogenicity_scores(df, output_file)``
  - ``plot_allele_frequency(df, output_file)``
  - ``plot_gene_burden(gene_df, output_file, top_n=30)``

Example:
    python scripts/plot_variant_distribution.py \\
        --variants all_variants.csv --genes gene_summary.csv \\
        --output-dir figures/
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install required packages with: pip install matplotlib pandas numpy")


_IMPACT_ORDER  = ["HIGH", "MODERATE", "LOW", "MODIFIER"]
_IMPACT_COLORS = {"HIGH": "#c0392b", "MODERATE": "#f39c12",
                  "LOW": "#27ae60", "MODIFIER": "#7f8c8d"}


def _save(fig: plt.Figure, base_path: Path) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    png = base_path.with_suffix(".png")
    svg = base_path.with_suffix(".svg")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    fig.savefig(svg,          bbox_inches="tight")
    plt.close(fig)
    print(f"    -> {png}")
    print(f"    -> {svg}")


def plot_consequence_distribution(
    df: pd.DataFrame, output_file: str | Path,
) -> None:
    if df.empty or "consequence" not in df.columns:
        print("  consequence distribution: no data")
        return
    counts = df["consequence"].astype(str).value_counts().head(20)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(counts.index[::-1], counts.values[::-1], color="#3a7ca5")
    ax.set_xlabel("Variant count")
    ax.set_title("Top 20 variant consequences")
    _save(fig, Path(output_file).with_suffix(""))


def plot_impact_by_chromosome(
    df: pd.DataFrame, output_file: str | Path,
) -> None:
    if df.empty or "impact" not in df.columns or "chrom" not in df.columns:
        print("  impact-by-chromosome: no data")
        return
    chroms = sorted(df["chrom"].astype(str).unique(),
                    key=lambda s: (len(s), s))
    impact_levels = [l for l in _IMPACT_ORDER if l in df["impact"].unique()]
    counts = (df.groupby(["chrom", "impact"]).size()
                .unstack(fill_value=0)
                .reindex(index=chroms, columns=impact_levels, fill_value=0))
    fig, ax = plt.subplots(figsize=(12, 5))
    bottom = np.zeros(len(chroms))
    for level in impact_levels:
        ax.bar(chroms, counts[level], bottom=bottom,
               label=level, color=_IMPACT_COLORS.get(level, "#bdc3c7"))
        bottom += counts[level].values
    ax.set_xlabel("Chromosome"); ax.set_ylabel("Variant count")
    ax.set_title("Impact distribution by chromosome")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(title="Impact", loc="upper right")
    _save(fig, Path(output_file).with_suffix(""))


def plot_pathogenicity_scores(
    df: pd.DataFrame, output_file: str | Path,
) -> None:
    columns = [c for c in ("cadd_phred", "revel", "sift", "polyphen")
               if c in df.columns]
    if not columns:
        print("  pathogenicity scores: no plugin columns present")
        return
    fig, axes = plt.subplots(1, len(columns), figsize=(4 * len(columns), 4))
    if len(columns) == 1:
        axes = [axes]
    for ax, col in zip(axes, columns):
        values = pd.to_numeric(df[col], errors="coerce").dropna()
        if values.empty:
            ax.set_title(f"{col} (no data)")
            continue
        ax.hist(values, bins=30, color="#3a7ca5", edgecolor="white")
        ax.set_title(col)
        ax.set_xlabel(col); ax.set_ylabel("count")
    fig.suptitle("Pathogenicity score distributions")
    fig.tight_layout()
    _save(fig, Path(output_file).with_suffix(""))


def plot_allele_frequency(df: pd.DataFrame,
                          output_file: str | Path) -> None:
    col = "gnomad_af" if "gnomad_af" in df.columns else (
        "af" if "af" in df.columns else None)
    if col is None:
        print("  allele frequency: no AF column present")
        return
    af = pd.to_numeric(df[col], errors="coerce").dropna()
    if af.empty:
        print("  allele frequency: column is empty")
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    # Log-spaced bins, with a left-edge bin for AF == 0
    af_pos = af[af > 0]
    bins = np.logspace(np.log10(max(af_pos.min(), 1e-7)), 0, 30)
    ax.hist(af_pos, bins=bins, color="#3a7ca5", edgecolor="white")
    ax.set_xscale("log")
    ax.set_xlabel(f"{col} (log scale)"); ax.set_ylabel("Variant count")
    ax.set_title("Population allele-frequency distribution")
    ax.axvline(0.01, color="#c0392b", linestyle="--",
               label="Rare cutoff (AF=0.01)")
    ax.legend()
    _save(fig, Path(output_file).with_suffix(""))


def plot_gene_burden(
    gene_df: pd.DataFrame, *, output_file: str | Path,
    top_n: int = 30,
) -> None:
    if gene_df.empty:
        print("  gene burden: empty input")
        return
    top = gene_df.head(top_n).copy()
    fig, ax = plt.subplots(figsize=(10, max(4, 0.3 * len(top))))
    impact_cols = [c for c in ("n_HIGH", "n_MODERATE", "n_LOW", "n_MODIFIER")
                   if c in top.columns]
    bottom = np.zeros(len(top))
    for col in impact_cols:
        level = col[2:]
        ax.barh(top["gene_symbol"], top[col], left=bottom,
                label=level,
                color=_IMPACT_COLORS.get(level, "#bdc3c7"))
        bottom += top[col].values
    ax.invert_yaxis()
    ax.set_xlabel("Variant count")
    ax.set_title(f"Top {len(top)} genes by variant burden")
    ax.legend(title="Impact", loc="lower right")
    _save(fig, Path(output_file).with_suffix(""))


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n")[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--variants", type=Path, required=True)
    parser.add_argument("--genes",    type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    print("=== Variant plots ===")
    df = pd.read_csv(args.variants)
    plot_consequence_distribution(df, args.output_dir / "consequence_distribution.png")
    plot_impact_by_chromosome(df,    args.output_dir / "impact_by_chromosome.png")
    plot_pathogenicity_scores(df,    args.output_dir / "pathogenicity_scores.png")
    plot_allele_frequency(df,        args.output_dir / "allele_frequency.png")
    if args.genes is not None and args.genes.exists():
        genes = pd.read_csv(args.genes)
        plot_gene_burden(genes, output_file=args.output_dir / "gene_burden.png")
    print(f"✓ Variant plots completed: {args.output_dir}")


if __name__ == "__main__":
    main()
