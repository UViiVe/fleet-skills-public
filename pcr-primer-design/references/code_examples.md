# Code Examples

Working snippets for every wrapper script. Copy, adapt, run.

← Back to [SKILL.md](../SKILL.md)
Related: [Parameter Ranges](parameter_ranges.md) · [Best Practices](primer_design_best_practices.md) · [MIQE Guidelines](miqe_guidelines.md)

The examples assume you're working from the skill's root directory so that
`scripts/` and `references/` are direct subdirectories.

---

## 1. Loading the Target Sequence

### From a FASTA file

```python
from Bio import SeqIO

record = SeqIO.read("path/to/target.fasta", "fasta")
sequence = str(record.seq).upper()
print(f"Loaded {record.id}: {len(sequence)} bp")
```

If the file holds multiple records, use `SeqIO.parse` and pick by id:

```python
records = {r.id: str(r.seq).upper() for r in SeqIO.parse("multi.fasta", "fasta")}
sequence = records["NM_002046"]
```

### From a RefSeq/GenBank accession

```python
from Bio import Entrez, SeqIO

Entrez.email = "you@example.com"            # NCBI requires a contact email
Entrez.api_key = "YOUR_NCBI_API_KEY"        # optional but bumps rate limit to 10/sec

with Entrez.efetch(db="nucleotide", id="NM_002046",
                   rettype="fasta", retmode="text") as handle:
    record = SeqIO.read(handle, "fasta")

sequence = str(record.seq).upper()
print(f"Retrieved {record.id}: {len(sequence)} bp")
```

### From a gene symbol

```python
from Bio import Entrez

Entrez.email = "you@example.com"

with Entrez.esearch(db="gene", term="GAPDH[Gene] AND Homo sapiens[Organism]") as h:
    search = Entrez.read(h)

gene_id = search["IdList"][0]
print(f"Gene ID: {gene_id}")
# Follow up with Entrez.elink to find the canonical RefSeq, then efetch
```

---

## 2. Designing Primers

### Standard PCR

```python
from scripts.design_standard_primers import design_pcr_primers

result = design_pcr_primers(
    sequence=sequence,
    target_region=None,                  # or (start, end) to constrain
    primer_size_range=(18, 25),
    tm_range=(55, 65),
    gc_range=(40, 60),
    amplicon_size_range=(200, 600),
    num_return=5,
)

for i, pair in enumerate(result["primers"][:3], 1):
    print(f"Pair {i}: {pair['forward_seq']} / {pair['reverse_seq']}")
    print(f"   amplicon {pair['amplicon_size']} bp, "
          f"Tm {pair['forward_tm']:.1f} / {pair['reverse_tm']:.1f}, "
          f"penalty {pair['penalty']:.3f}")
```

### MIQE-compliant qPCR

```python
from scripts.design_qpcr_primers import (
    design_qpcr_primers,
    validate_miqe_compliance,
)

result = design_qpcr_primers(
    sequence=sequence,
    amplicon_size_range=(70, 140),
    tm_match_threshold=2.0,
    primer_tm_range=(58, 62),
    num_return=5,
)

print(f"MIQE-compliant pairs: {result['num_miqe_compliant']} of "
      f"{result['num_primers_found']}")

for pair in result["primers"]:
    audit = validate_miqe_compliance(pair)
    if not audit["miqe_compliant"]:
        print(f"  Pair {pair['pair_id']} failed: "
              f"{', '.join(audit['failed_checks'])}")
```

### TaqMan assay

```python
from scripts.design_taqman_probes import (
    design_taqman_assay,
    validate_taqman_assay,
)

assay = design_taqman_assay(
    sequence=sequence,
    probe_tm_offset=8.0,                 # probe Tm = primer Tm + 8°C
    primer_tm_range=(58, 62),
    amplicon_size_range=(70, 140),
)

best = assay["assays"][0]
print(f"F: {best['forward_seq']}  Tm {best['forward_tm']:.1f}°C")
print(f"R: {best['reverse_seq']}  Tm {best['reverse_tm']:.1f}°C")
print(f"P: {best['probe_seq']}    Tm {best['probe_tm']:.1f}°C")

audit = validate_taqman_assay(best)
for w in audit["warnings"]:
    print(f"  warning: {w}")
```

---

## 3. Validation

### Melting temperature

```python
from scripts.calculate_tm import (
    calculate_tm_comprehensive,
    calculate_tm_range,
)

tm = calculate_tm_comprehensive(
    "ATGCGTACGATCGATCG",
    method="all",
    salt_conc=50.0,
    mg_conc=3.0,                        # typical qPCR master-mix Mg
    dna_conc=200.0,
)
print(f"Nearest-neighbor Tm: {tm['tm_nn']:.1f}°C")
print(f"Salt-adjusted Wallace: {tm['tm_salt_adjusted']:.1f}°C")
print(f"%GC formula: {tm['tm_percent']:.1f}°C")

panel = calculate_tm_range([p["forward_seq"] for p in result["primers"]])
print(f"Tm spread across forward primers: {panel['tm_range']:.1f}°C")
```

### Dimer analysis

```python
from scripts.check_dimers import (
    analyze_dimers,
    check_3prime_complementarity,
)

pair = result["primers"][0]
dim = analyze_dimers(
    [pair["forward_seq"], pair["reverse_seq"]],
    temperature=60.0,
    dg_threshold=-5.0,
)
print(f"Problematic interactions: {dim['num_problematic']}")
for i in dim["interactions"]:
    flag = "FAIL" if i["problematic"] else "ok"
    print(f"  [{flag}] {i['type']:>12s}  ΔG = {i['dg']:>6.2f} kcal/mol")

three = check_3prime_complementarity(
    pair["forward_seq"], pair["reverse_seq"], window_size=5
)
print(f"3' complementary bases: {three['complementary_bases']}")
```

### Secondary structures

```python
from scripts.check_secondary_structures import (
    analyze_secondary_structures,
    comprehensive_qc_check,
    check_gc_clamp,
    check_poly_runs,
)

primer = result["primers"][0]["forward_seq"]

s = analyze_secondary_structures(primer, temperature=60.0)
print(f"Hairpin ΔG     {s['hairpin']['dg']:>6.2f}  problematic={s['hairpin']['problematic']}")
print(f"Self-dimer ΔG  {s['self_dimer']['dg']:>6.2f}  problematic={s['self_dimer']['problematic']}")

qc = comprehensive_qc_check(primer, temperature=60.0)
if not qc["passes_qc"]:
    for issue in qc["issues"]:
        print(f"  - {issue}")
```

### Multiplex compatibility

```python
from scripts.check_dimers import analyze_multiplex_compatibility

pairs = [
    {"forward_seq": "ATGCGTACGATCGATCG", "reverse_seq": "CGTAGCTAGCTAGCTAG",
     "forward_tm": 58.5, "reverse_tm": 59.2},
    {"forward_seq": "GCTATCGATCGTAGCTA", "reverse_seq": "TAGCTAGCTAGCGATCG",
     "forward_tm": 59.0, "reverse_tm": 58.8},
]

compat = analyze_multiplex_compatibility(pairs, temperature=60.0)
print(f"Multiplex compatible: {compat['is_compatible']}")
print(f"Tm range across all primers: {compat['tm_range']:.1f}°C")
for rec in compat["recommendations"]:
    print(f"  • {rec}")
```

### In-silico specificity

```python
from scripts.validate_specificity import in_silico_pcr

products = in_silico_pcr(
    forward_primer=pair["forward_seq"],
    reverse_primer=pair["reverse_seq"],
    sequence=sequence,
    max_product_size=5000,
)
print(f"In-silico products: {len(products)}")
for p in products:
    print(f"  {p['product_size']} bp at {p['forward_start']}–{p['reverse_start']}")
```

For Primer-BLAST against a public NCBI database, see the placeholder
function `check_primer_specificity` in `scripts/validate_specificity.py` —
production use requires submitting an asynchronous request and polling for
results.

---

## 4. Reports

```python
from scripts.generate_reports import (
    generate_primer_report,
    generate_summary_table,
)

# Pull everything together
validation = {
    "specificity": {"is_specific": True, "on_target_hits": 1, "off_target_hits": 0},
    "dimers": dim,
    "secondary_structures": s,
}

report_md = generate_primer_report(
    primers=result,
    validation_results=validation,
    output_format="markdown",
    include_miqe_checklist=True,
)
with open("primer_design_report.md", "w") as f:
    f.write(report_md)

# HTML
report_html = generate_primer_report(
    primers=result,
    validation_results=validation,
    output_format="html",
)
with open("primer_design_report.html", "w") as f:
    f.write(report_html)
```

Compare multiple designs:

```python
summary = generate_summary_table([result_v1, result_v2, result_v3])
print(summary)
```

---

## 5. Visualization

```python
from scripts.visualize_primers import (
    plot_primer_alignment,
    plot_tm_distribution,
    plot_amplicon_sizes,
)

plot_primer_alignment(
    sequence=sequence,
    pair=result["primers"][0],
    output_file="primer_alignment.png",
)

plot_tm_distribution(
    pairs=result["primers"],
    output_file="tm_distribution.png",
)

plot_amplicon_sizes(
    pairs=result["primers"],
    target_range=(70, 140),
    output_file="amplicon_sizes.png",
)
```

All plots default to 300 DPI raster PNG, sized for a US-Letter document at
~80% page width.

---

## 6. Export

```python
from scripts.export_results import (
    export_primers,
    export_for_benchling,
)

# CSV
export_primers(result, fmt="csv", output_file="primers.csv",
               include_validation=True, validation_results=validation)

# Excel with Primers / Parameters / Validation sheets
export_primers(result, fmt="excel", output_file="primers.xlsx",
               include_validation=True, validation_results=validation)

# JSON
export_primers(result, fmt="json", output_file="primers.json",
               include_validation=True, validation_results=validation)

# IDT order sheet
export_primers(result, fmt="idt_order", output_file="idt_order.tsv")

# MIQE checklist (qPCR only)
export_primers(result, fmt="miqe_checklist", output_file="miqe_checklist.xlsx")

# Benchling-compatible CSV
export_for_benchling(result, output_file="benchling.csv")
```

---

## 7. End-to-End Pipeline

A complete qPCR workflow in one script:

```python
from scripts.design_qpcr_primers import design_qpcr_primers
from scripts.check_dimers import analyze_dimers
from scripts.check_secondary_structures import analyze_secondary_structures
from scripts.validate_specificity import in_silico_pcr
from scripts.generate_reports import generate_primer_report
from scripts.export_results import export_primers

# 1. Design
designs = design_qpcr_primers(sequence=sequence, num_return=5)
best = designs["primers"][0]

# 2. Validate
dim = analyze_dimers([best["forward_seq"], best["reverse_seq"]], temperature=60.0)
sec_f = analyze_secondary_structures(best["forward_seq"], temperature=60.0)
sec_r = analyze_secondary_structures(best["reverse_seq"], temperature=60.0)
products = in_silico_pcr(best["forward_seq"], best["reverse_seq"], sequence)

validation = {
    "specificity": {"is_specific": len(products) == 1,
                    "on_target_hits": len(products),
                    "off_target_hits": 0},
    "dimers": dim,
    "secondary_structures": {"forward": sec_f, "reverse": sec_r},
}

# 3. Report and export
report = generate_primer_report(
    primers=designs, validation_results=validation,
    output_format="markdown", include_miqe_checklist=True,
)
with open("report.md", "w") as f:
    f.write(report)

export_primers(designs, fmt="csv", output_file="primers.csv")
export_primers(designs, fmt="idt_order", output_file="idt_order.tsv")
export_primers(designs, fmt="miqe_checklist", output_file="miqe.xlsx")
```

---

Last updated: 2026-05-20.
