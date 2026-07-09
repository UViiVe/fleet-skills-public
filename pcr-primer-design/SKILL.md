---
name: pcr-primer-design
id: "skill_992225c9dbb44ee88afaa48fe93d4309"
display-name: "PCR Primer Design"
short-description: "Design and validate primers for PCR, qPCR, TaqMan probes, multiplex, and sequencing assays with thermodynamic QC and MIQE-compliant reporting. Use when you need oligos for any nucleic-acid amplific..."
category: wet_lab
keywords: [PCR, primer, design, amplification, wet lab]
last-updated: "2026-07-07"
version: 1.0.1
author: "Fleet Skills contributors"
created: 2026-05-20
description: "Design and validate primers for PCR, qPCR, TaqMan probes, multiplex, and sequencing assays with thermodynamic QC and MIQE-compliant reporting. Use when you need oligos for any nucleic-acid amplification experiment."
source: https://github.com/UViiVe/fleet-skills-public
license: MIT
provenance:
  wraps: "Primer3 / NCBI Primer-BLAST"
  upstream_source: "https://github.com/UViiVe/fleet-skills-public"
  citation: "Untergasser A et al., Nucleic Acids Research 2012, doi:10.1093/nar/gks596; Ye J et al., BMC Bioinformatics 2012, doi:10.1186/1471-2105-13-134"
changelog:
  - version: 1.0.1
    date: 2026-07-07
    event: citation-provenance
    summary: "Added provenance block and References section with upstream tool citations."
status: active
promoted: true
---

# PCR Primer Design

End-to-end primer design that combines `primer3` for candidate generation with
thermodynamic QC (Tm, dimers, hairpins, 3' complementarity), in-silico
specificity checks, and publication-grade reports.

The pipeline is application-aware: it produces sane defaults for standard PCR,
qPCR (MIQE 2.0), TaqMan probes, multiplex sets, and Sanger sequencing — and
emits IDT-ready order sheets and a MIQE checklist alongside the usual CSV/Excel
output.

## When to Use This Skill

✅ Use this skill when you need:
- **qPCR primers** that satisfy MIQE 2.0 (70–140 bp amplicon, Tm matched within 2°C, in-silico specificity logged)
- **Standard PCR primers** for cloning, screening, or genotyping (100 bp – 5 kb amplicons)
- **TaqMan probe + primer trios** with the probe Tm set 5–10°C above the primers
- **Multiplex compatibility checks** — pairwise heterodimer screening across many primer pairs at once
- **Publication-ready documentation** — markdown/HTML reports, MIQE checklist, IDT order sheet

❌ Avoid this skill for:
- **CRISPR guide design** — use a guide-RNA scoring tool (off-target scoring is structurally different)
- **NGS adapter or index sequence selection** — vendor adapter sequences are fixed and not designed de novo
- **In-situ hybridization probes** — fluorescence/quenching constraints are different from PCR
- **Long-amplicon cloning (>10 kb)** — high-fidelity polymerase suppliers usually publish dedicated design recipes

**Application picker:**

| Application | Amplicon | Tm match | Strict 3' end |
|-------------|----------|----------|---------------|
| qPCR (SYBR) | 70–140 bp | ≤ 2°C | yes |
| qPCR (TaqMan) | 70–140 bp | ≤ 2°C | yes |
| Standard PCR | 100–1000 bp | ≤ 5°C | recommended |
| Multiplex | 100–600 bp | ≤ 3°C across all | yes |
| Sequencing | n/a (single primer) | n/a | recommended |

## Quick Start

The fastest path is to design qPCR primers against an in-memory target sequence.

```python
from scripts.design_qpcr_primers import design_qpcr_primers

# 700 bp human GAPDH exon 3-4 region (truncated for the example)
target = (
    "ATGGGGAAGGTGAAGGTCGGAGTCAACGGATTTGGTCGTATTGGGCGCCTGGTCACCAGGG"
    "CTGCTTTTAACTCTGGTAAAGTGGATATTGTTGCCATCAATGACCCCTTCATTGACCTCAA"
    "CTACATGGTTTACATGTTCCAATATGATTCCACCCATGGCAAATTCCATGGCACCGTCAAG"
    "GCTGAGAACGGGAAGCTTGTCATCAATGGAAATCCCATCACCATCTTCCAGGAGCGAGATC"
    "CCTCCAAAATCAAGTGGGGCGATGCTGGCGCTGAGTACGTCGTGGAGTCCACTGGCGTCTT"
    "CACCACCATGGAGAAGGCTGGGGCTCATTTGCAGGGGGGAGCCAAAAGGGTCATCATCTCT"
    "GCCCCCTCTGCTGATGCCCCCATGTTCGTCATGGGTGTGAACCATGAGAAGTATGACAACA"
    "GCCTCAAGATCATCAGCAATGCCTCCTGCACCACCAACTGCTTAGCACCCCTGGCCAAGGT"
    "CATCCATGACAACTTTGGTATCGTGGAAGGACTCATGACCACAGTCCATGCCATCACTGCC"
    "ACCCAGAAGACTGTGGATGGCCCCTCCGGGAAACTGTGGCGTGATGGCCGCGGGGCTCTCC"
    "AGAACATCATCCCTGCCTCTACTGGCGCTGCCAAGGCTGTGGGCAAGGTCATCCCTGAGCT"
    "GAACGGGAAGCTCACTGGCATGGCCTTCCGTGTCCCCACTGCCAACGTGTCAGTGGTG"
)

result = design_qpcr_primers(
    sequence=target,
    amplicon_size_range=(80, 120),
    num_return=5,
)
print(f"Returned {result['num_primers_found']} pairs, "
      f"{result['num_miqe_compliant']} MIQE-compliant")
```

Expected runtime: under 2 seconds on a laptop. You typically get 3–5
MIQE-compliant pairs.

## Installation

| Package | Version | License | Commercial use | Install |
|---------|---------|---------|----------------|---------|
| primer3-py | ≥ 2.0 | GPLv2 | yes (execution; redistribution restricted) | `pip install primer3-py` |
| biopython | ≥ 1.80 | BSD-3 | yes | `pip install biopython` |
| pandas | ≥ 1.5 | BSD-3 | yes | `pip install pandas` |
| openpyxl | ≥ 3.0 | MIT | yes | `pip install openpyxl` |
| matplotlib | ≥ 3.6 | matplotlib | yes | `pip install matplotlib` |
| requests | ≥ 2.28 | Apache-2.0 | yes | `pip install requests` |

```bash
python -m venv .venv && source .venv/bin/activate
pip install primer3-py biopython pandas openpyxl matplotlib requests
```

`primer3-py` ships compiled binaries via wheels, so no `cmake` is needed on
macOS or Linux.

NCBI Primer-BLAST (used for off-target screening) is rate-limited to 3
requests/sec without an API key and 10 requests/sec with one. The free key is
available at https://account.ncbi.nlm.nih.gov/.

## Inputs

**Required**
- `sequence` — the target DNA as one of:
  - a raw `str` of A/T/C/G/N (U is auto-converted to T)
  - a FASTA file path (read with `Bio.SeqIO.read`)
  - a GenBank/RefSeq accession (e.g. `NM_002046`), fetched via `Bio.Entrez`

**Optional**
- `target_region` — `(start, end)` to constrain primer placement around a SNP or exon
- `excluded_regions` — list of `(start, length)` blocks to avoid (repeats, paralogs, SNPs)
- `organism` — for downstream Primer-BLAST specificity check
- `application` — `"qpcr" | "standard" | "taqman" | "multiplex" | "sequencing"` (default: `qpcr`)
- Tm, GC, amplicon-size, and primer-length overrides — see [references/parameter_ranges.md](references/parameter_ranges.md)

**Minimum sequence length**: 150 bp for qPCR, 300 bp for standard PCR. Anything
shorter will refuse to design enough non-overlapping candidates.

## Outputs

**Primary**
- A `dict` per primer pair with: `forward_seq`, `reverse_seq`, `forward_tm`,
  `reverse_tm`, `tm_diff`, `forward_gc`, `reverse_gc`, `amplicon_size`,
  positions, primer3 penalty, and `miqe_compliant` flag.

**Files** (paths set by the caller; nothing is written unless you ask):
- `primers.csv` — flat table for spreadsheets
- `primers.xlsx` — Excel with `Primers`, `Parameters`, `Validation` sheets
- `primers.json` — full structured record for downstream tools
- `idt_order.tsv` — Name / Sequence / Scale / Purification, paste-ready
- `miqe_checklist.xlsx` — pre-filled MIQE 2.0 essentials (qPCR only)
- `report.md` or `report.html` — human-readable design + validation summary

**Plots** (300 DPI PNG or SVG):
- Primer binding-site map on target
- Tm distribution across pairs
- GC% and amplicon-size summary
- QC heatmap (Tm match / GC range / dimers / specificity)

Downstream consumers: ordering systems (IDT, Benchling), qPCR thermocyclers,
the `qpcr-data-analysis` skill (when Cq values come back), and any LIMS that
ingests the JSON.

## Clarification Questions

### Q1. ⚠️ CRITICAL: ASK THIS FIRST — What's the target sequence?

Pick one and provide it:
1. Raw sequence pasted in chat (≥ 150 bp)
2. FASTA file path
3. RefSeq/GenBank accession (e.g. `NM_001256799`) — also tell me the species
4. Gene symbol + organism — I'll fetch the canonical transcript from NCBI

If you only have a partial sequence, also tell me what region of the gene it
corresponds to (exon numbers, transcript variant).

### Q2. What's the assay format?

qPCR (SYBR), qPCR (TaqMan), standard PCR, multiplex PCR, Sanger sequencing,
SNP genotyping? Default is **qPCR (SYBR)** unless you say otherwise.

### Q3. Any region constraints?

- SNPs to avoid (give positions)
- Exon–exon junction the amplicon must span (qPCR best practice when gDNA could contaminate)
- Specific position the probe/primer must overlap

### Q4. How strict should specificity be?

- `basic` — local in-silico PCR on the supplied sequence (default; ~1 s)
- `standard` — basic + BLAST against representative RefSeq for the organism (~30 s with an API key)
- `publication` — standard + full Primer-BLAST submission with the organism database (~3–5 min)

### Q5. What outputs do you want?

CSV (default), Excel, JSON, IDT order sheet, MIQE checklist, markdown report,
HTML report, plots. Multiple are fine.

## Standard Workflow

1. **Load the target.** From file, raw string, or `Bio.Entrez.efetch`.
2. **Choose the application module** in `scripts/`:
   - `design_qpcr_primers.py` — MIQE-compliant qPCR (default)
   - `design_standard_primers.py` — broad-range PCR
   - `design_taqman_probes.py` — primers + internal probe
3. **Validate every returned pair:**
   - `check_dimers.analyze_dimers` — homodimer + heterodimer ΔG at the
     annealing temperature
   - `check_secondary_structures.analyze_secondary_structures` — hairpins and
     3'-self-complementarity
   - `validate_specificity.in_silico_pcr` — local exact-match scan; switch to
     Primer-BLAST for publication
4. **Render the report** with `generate_reports.generate_primer_report` and
   export with `export_results.export_primers`.
5. **Plot** binding sites and QC summary via `visualize_primers.*` if a figure
   is needed.

The scripts are intentionally single-purpose so they can be chained from a
notebook or a CLI driver.

## When Scripts Fail

Walk this escalation order:

1. **Read the error.** `primer3` raises `OSError` when no pairs satisfy the
   constraints — that's a parameter problem, not a bug. Loosen Tm or amplicon
   range and retry.
2. **Tweak parameters** via the script's keyword args. The `scripts/` files
   forward all primer3 parameters through `primer3_params` dicts, so any
   primer3 setting documented in
   [references/parameter_ranges.md](references/parameter_ranges.md) can be
   passed through.
3. **Fall back to a reference example.** [references/code_examples.md](references/code_examples.md)
   has working snippets for every application.
4. **Roll your own.** If you need degenerate primers or odd salt corrections,
   write a one-off script that calls `primer3.designPrimers` directly — the
   wrappers exist for ergonomics, not gatekeeping.

## Decision Points

| Choice | Default | When to change |
|--------|---------|----------------|
| Optimal Tm | 60°C | Move to 58°C for AT-rich targets, 62°C for GC-rich |
| `PRIMER_PAIR_MAX_DIFF_TM` | 2.0 (qPCR), 5.0 (PCR) | Drop to 1.0 if running multiple assays on a single block |
| `PRIMER_GC_CLAMP` | 1 | Set to 2 for difficult templates needing a stronger anchor |
| Amplicon size | 70–140 (qPCR), 100–1000 (PCR) | Shorten for FFPE/degraded DNA (target 60–90 bp) |
| `PRIMER_MAX_POLY_X` | 4 | Drop to 3 for repeats-prone targets |
| Specificity check depth | `basic` | Always `publication` if writing a manuscript |
| Number returned | 5 | Increase to 20 when designing for SNP-rich regions; you'll filter heavily |

## Common Issues

| Issue | Likely cause | Fix |
|-------|-------------|-----|
| `primer3` returns 0 pairs | Target too short or constraints too strict | Increase amplicon range; relax Tm window to ±5°C; widen GC to 35–65% |
| ΔG hetero-dimer < −9 kcal/mol on top pair | Complementary 3' ends | Drop the offending pair; re-rank by penalty |
| Multiple in-silico PCR products | Repeat element or paralog in target | Mask repeats with RepeatMasker first; design in unique flanking region |
| qPCR efficiency < 90% experimentally | Amplicon too long, secondary structure | Redesign with amplicon ≤ 120 bp; check hairpin ΔG > −2 kcal/mol |
| qPCR efficiency > 110% experimentally | Primer dimers, NTC contamination | Inspect melt curve sub-80°C peak; tighten heterodimer threshold |
| Tm mismatch between pairs > 2°C | Different GC content | Extend the cooler primer by 1–2 nt; use primer3 penalty weights |
| `ImportError: primer3` | Missing C extension wheel | `pip install --upgrade primer3-py`; on Apple Silicon use Python ≥ 3.10 |
| `Entrez.efetch` returns empty | Stale or wrong accession | Verify on https://www.ncbi.nlm.nih.gov/nuccore/; use accession **without** version suffix |
| MIQE compliance False but values look fine | Tm outside 58–62°C window | Tighten `primer_tm_range` upfront, or accept and document in methods |
| Primer-BLAST hangs | Rate-limited | Use the in-silico PCR fallback; add an NCBI API key |
| TaqMan probe starts with G | primer3 defaults | Re-run with `PRIMER_INTERNAL_FIRST_BASE_INDEX` shifted by 1 |
| Multiplex `is_compatible=False` | One pair dimerizes with another | Drop the worst pair; re-run pairwise check |

Detailed troubleshooting is in [references/troubleshooting_guide.md](references/troubleshooting_guide.md).

## Best Practices

- **Always check dimers at the annealing temperature you'll actually use** — `primer3.calcHeterodimer` is temperature-sensitive
- **Pin your primer3-py version in the project's pyproject/requirements** — the SantaLucia parameter set changed in 2.0
- **Save the JSON output, not just the CSV** — primer3 penalties and positional info are useful when an assay misbehaves months later
- **Treat the in-silico PCR pass as necessary but not sufficient** — confirm with Primer-BLAST or local BLAST before ordering
- **Record the exact Mg²⁺ concentration** used in Tm calculation; it changes Tm by ~1°C per mM

## Suggested Next Steps

After primers come back from the synthesizer:

1. **Wet-lab QC** — run a gradient PCR (Tm ± 3°C) to find the optimal annealing temperature
2. **For qPCR — generate a standard curve.** Compute efficiency (E = 10^(−1/slope) − 1, target 90–110%) and R² (≥ 0.98). Hand the Cq data to the **qpcr-data-analysis** skill.
3. **For TaqMan** — verify probe quenching by running a no-template control and confirming baseline fluorescence stays flat
4. **Archive the design** — keep the JSON output + the MIQE checklist in the experimental record
5. **For publication** — generate the markdown report with `include_miqe_checklist=True` and attach as supplementary

## Related Skills

- **qpcr-data-analysis** — fits standard curves, computes ΔΔCq, normalizes against reference genes
- **bulk-rnaseq-differential-expression** — produces the gene list you'll validate by qPCR
- **functional-enrichment-from-degs** — picks the candidate genes worth designing primers against

## References

### Local documentation
- [references/primer_design_best_practices.md](references/primer_design_best_practices.md) — design principles and challenging templates
- [references/parameter_ranges.md](references/parameter_ranges.md) — quick tables for Tm, GC, amplicon size by application
- [references/miqe_guidelines.md](references/miqe_guidelines.md) — MIQE 2.0 essentials, validation thresholds, checklist
- [references/code_examples.md](references/code_examples.md) — full working snippets for every script
- [references/troubleshooting_guide.md](references/troubleshooting_guide.md) — diagnostic flowcharts

### Primary citations
- Bustin SA et al. *Clinical Chemistry* 71(6):634–660 (2025) — MIQE 2.0 guidelines
- Untergasser A et al. *Nucleic Acids Research* 40(15):e115 (2012) — Primer3 algorithm
- Ye J et al. *BMC Bioinformatics* 13:134 (2012) — Primer-BLAST
- SantaLucia J. *PNAS* 95(4):1460–1465 (1998) — unified nearest-neighbor thermodynamic parameters

### External tools
- Primer3 web: https://primer3.ut.ee/
- NCBI Primer-BLAST: https://www.ncbi.nlm.nih.gov/tools/primer-blast/
- RDML / MIQE site: https://rdml.org/miqe.html

## References

- **[Primer3]**: Untergasser A, et al. "Primer3—new capabilities and interfaces." *Nucleic Acids Research* 2012;40:e115. doi:10.1093/nar/gks596 — [link](https://doi.org/10.1093/nar/gks596)
- **[Primer-BLAST]**: Ye J, et al. "Primer-BLAST: A tool to design target-specific primers for polymerase chain reaction." *BMC Bioinformatics* 2012;13:134. doi:10.1186/1471-2105-13-134 — [link](https://doi.org/10.1186/1471-2105-13-134)
- **[Upstream skill source]**: https://github.com/UViiVe/fleet-skills-public (MIT license)