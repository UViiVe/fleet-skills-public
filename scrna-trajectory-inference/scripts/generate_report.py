"""Generate a trajectory analysis report (PDF via reportlab, Markdown fallback).

Inputs:  AnnData + results dict + output dir with plots already generated.
Outputs: `trajectory_analysis_report.pdf` (if reportlab) and a `.md` fallback.
Example:
    python generate_report.py --input adata.h5ad --results results.pkl \\
        --output-dir trajectory_results
"""

from __future__ import annotations

import argparse
import pickle
import sys
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, TYPE_CHECKING

try:
    import numpy as np
    import pandas as pd
except ImportError:
    sys.exit("Install with: pip install numpy pandas")

if TYPE_CHECKING:
    from anndata import AnnData


def generate_report(
    adata: "AnnData",
    results: Dict[str, Any],
    output_dir: str = "trajectory_results",
) -> None:
    """Write a Markdown report and (if reportlab available) a PDF."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    md_path = out / "trajectory_analysis_report.md"
    md = _build_markdown(adata, results, out)
    md_path.write_text(md)
    print(f"✓ Wrote {md_path}")

    try:
        _build_pdf(adata, results, out, md)
    except Exception as exc:
        print(f"  PDF generation failed: {exc} (Markdown report still written)")


def _build_markdown(adata, results, out_dir: Path) -> str:
    params = results["parameters"]
    pt = results.get("pseudotime")
    finite = pt[np.isfinite(pt)] if isinstance(pt, pd.Series) else None
    lines: List[str] = [
        "# Single-Cell Trajectory Inference Report",
        "",
        f"Generated: {datetime.utcnow().isoformat()}Z",
        "",
        "## Methods",
        "",
        "- Diffusion pseudotime (Haghverdi et al. 2016) seeded at the root cell.",
        "- PAGA (Wolf et al. 2019) for cluster-level connectivity.",
        f"- Neighbor graph: n_neighbors={params['n_neighbors']}, n_pcs={params['n_pcs']}.",
        f"- Diffusion map: n_dcs={params['n_dcs']}.",
    ]
    if results.get("velocity_results"):
        lines.append(f"- RNA velocity ({results['velocity_results']['mode']} model, scVelo).")
    if results.get("cellrank_results"):
        lines.append("- Fate probabilities and driver genes via CellRank.")
    lines += [
        "",
        "## Dataset",
        f"- Cells: {adata.n_obs:,}",
        f"- Genes: {adata.n_vars:,}",
        f"- Root cell type: {params['root_cell_type']}",
        f"- Cluster column: {params['cluster_key']}",
    ]
    if finite is not None:
        lines += [
            f"- Pseudotime range: [{finite.min():.3f}, {finite.max():.3f}]",
            f"- Cells with finite pseudotime: {finite.shape[0]} / {adata.n_obs}",
        ]
    lines += ["", "## Top trajectory genes"]
    tg = results.get("trajectory_genes")
    if isinstance(tg, pd.DataFrame) and not tg.empty:
        top = tg.head(20)
        lines.append("")
        lines.append("| Gene | Correlation | FDR | Direction |")
        lines.append("| --- | --- | --- | --- |")
        for _, row in top.iterrows():
            lines.append(f"| {row['gene']} | {row['correlation']:.3f} | "
                          f"{row['fdr']:.2e} | {row['direction']} |")

    if results.get("velocity_results"):
        lines += ["", "## RNA velocity",
                  f"- Model: {results['velocity_results']['mode']}",
                  f"- Latent time available: {results['velocity_results']['has_latent_time']}"]
    if results.get("cellrank_results"):
        ts = results["cellrank_results"].get("terminal_states", [])
        lines += ["", "## CellRank fate mapping",
                  f"- Terminal states detected: {len(ts)}",
                  f"- Terminal state labels: {', '.join(map(str, ts))}"]

    lines += ["", "## Figures", ""]
    for fig in sorted(out_dir.glob("*.png")):
        lines.append(f"- `{fig.name}`")
    return "\n".join(lines) + "\n"


def _build_pdf(adata, results, out_dir: Path, markdown_text: str) -> None:
    try:
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak
    except ImportError:
        raise ImportError("reportlab")

    styles = getSampleStyleSheet()
    pdf_path = out_dir / "trajectory_analysis_report.pdf"
    doc = SimpleDocTemplate(str(pdf_path), pagesize=letter,
                              topMargin=36, bottomMargin=36,
                              leftMargin=36, rightMargin=36)
    story = []
    for line in markdown_text.splitlines():
        if line.startswith("# "):
            story.append(Paragraph(line[2:], styles["Title"]))
        elif line.startswith("## "):
            story.append(Spacer(1, 8))
            story.append(Paragraph(line[3:], styles["Heading2"]))
        elif line.startswith("- "):
            story.append(Paragraph("• " + line[2:], styles["BodyText"]))
        elif line.startswith("|"):
            story.append(Paragraph(line.replace("|", " | "), styles["Code"]))
        elif line.strip():
            story.append(Paragraph(line, styles["BodyText"]))

    for fig in sorted(out_dir.glob("*.png")):
        story.append(PageBreak())
        story.append(Paragraph(fig.name, styles["Heading3"]))
        story.append(Spacer(1, 4))
        story.append(Image(str(fig), width=520, height=400))

    doc.build(story)
    print(f"✓ Wrote {pdf_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, help=".h5ad input")
    parser.add_argument("--results", required=True, help="Pickle from run_trajectory_analysis")
    parser.add_argument("--output-dir", default="trajectory_results")
    args = parser.parse_args()

    import scanpy as sc
    adata = sc.read_h5ad(args.input)
    with open(args.results, "rb") as f:
        results = pickle.load(f)
    generate_report(adata, results, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
