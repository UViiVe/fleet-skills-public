# PCR Primer Troubleshooting Guide

Diagnostic recipes for the failure modes you'll actually hit. Each section
maps a symptom to its likely cause and the next experimental step.

← Back to [SKILL.md](../SKILL.md)
Related: [Parameter Ranges](parameter_ranges.md) · [Best Practices](primer_design_best_practices.md)

---

## How to Use This Guide

1. Identify the symptom (no amplification, weak signal, multiple products, etc.)
2. Run the **First Things to Try** at the top of the relevant section
3. If those don't resolve it, work through the **Cause and Action** table
4. The **Diagnostic Flow** at the end of each major section organizes the
   decision tree

---

## No Amplification

**Symptoms:** empty lane on a gel, qPCR amplification curve flat, all Cq
values undetected (or > 35 cycles)

**First things to try:**
1. Gradient PCR across Ta ± 5°C in 1–2°C steps
2. Re-quantify template; aim for 10–100 ng/reaction (gDNA) or equivalent
   cDNA
3. BLAST the primers; check the latest reference build for SNPs in the
   primer-binding sites

**Cause and action:**

| Cause | Check | Fix |
|-------|-------|-----|
| Primer 3' mismatch to template | BLAST against the species build the sample comes from | Move primer 5–10 bp upstream of any SNP |
| Degraded primers | Date on the synthesis order; freeze-thaw cycles | Reorder; aliquot working stocks to avoid > 5 freeze-thaws |
| Wrong primer concentration | Spectrophotometer of working stock | Remake at 10 µM; use 0.2–0.5 µM in the reaction |
| Template too dilute | Qubit reading | Increase input 5–10×; recheck for inhibitors |
| Template degraded | Run a 1% gel; check 260/280 (1.8–2.0) | Re-extract |
| Ta too high | Try Tm − 7°C instead of Tm − 5°C | Drop annealing temp 3–5°C |
| Insufficient Mg²⁺ | Buffer composition | Bump to 2.5–3.0 mM |
| Inhibitors carried over | Dilute template 1:10, retest | Re-extract with cleanup column; add 0.1 mg/mL BSA |
| Extension too short | Amplicon length | 1 min/kb for standard Taq, 30 s/kb for high-fidelity |
| Enzyme inactive | Positive control fails too | Use a fresh aliquot |

**Diagnostic flow:**

```
No amplification
├── Positive control works?
│    ├── No → reagent/instrument issue; replace polymerase or buffer
│    └── Yes → primer or template issue:
│        ├── BLAST primer → 3' mismatch found? Redesign
│        ├── Template gel → degraded? Re-extract
│        └── Gradient PCR → no clear optimum? Redesign primers
└── Check primer secondary structure (hairpin, self-dimer at Ta)
```

---

## Weak or Variable Amplification

**Symptoms:** faint band, Cq > 30 for a gene that should be abundant,
inconsistent replicates

**First things to try:**
1. Increase template input 2–10×
2. Lower Ta by 3°C
3. Bump Mg²⁺ to 2.5–3.0 mM

**Cause and action:**

| Cause | Check | Fix |
|-------|-------|-----|
| Low template | Qubit | Concentrate or use more input |
| Primer inefficient | Check Tm, GC, secondary structure | Redesign |
| Suboptimal Ta | Gradient PCR | Use the strongest-band temperature |
| Template secondary structure | Predict mfold around the binding sites | Add DMSO 2–10% or betaine 1–2 M |
| Primer degraded | Storage history | Reorder |
| Pipetting variance | Replicate CV | Use electronic pipettes; increase reaction volume to 20 µL |
| Cycle count too low | Number of cycles | qPCR runs 40, standard PCR 30–35 |

---

## Non-Specific Amplification

**Symptoms:** multiple bands on gel, multiple peaks on melt curve,
amplification in NTCs

**First things to try:**
1. Raise Ta by 3–5°C
2. Switch to a hot-start polymerase
3. Decrease primer concentration to 0.2 µM

**Cause and action:**

| Cause | Check | Fix |
|-------|-------|-----|
| Ta too low | Gradient PCR | Use a higher Ta |
| Primer too short / non-specific | BLAST hits with ≤ 2 mismatches | Extend primer or move site |
| Primer dimer | Heterodimer ΔG at Ta | Redesign one primer |
| Too many cycles | Cycle count | 25–30 cycles for standard PCR |
| Mg²⁺ too high | Reaction composition | Reduce to 1.5–2.0 mM |
| Hot-start absent | Polymerase product | Switch to hot-start Taq |
| gDNA in RT-qPCR | −RT control amplifies | Treat RNA with DNase; redesign primers to span exon junction |

**Touchdown PCR** can rescue moderate specificity problems without
redesigning. Start Ta 5°C above the calculated optimum and drop 0.5–1°C
per cycle for the first 10 cycles, then settle at Tm − 5°C for the rest.

---

## Primer Dimers

**Symptoms:** band under 100 bp on gel, low-temperature peak (< 80°C) on
melt curve, early amplification of NTCs

**First things to try:**
1. Raise Ta
2. Drop primer concentration to 0.2 µM
3. If neither helps, redesign — dimers usually mean the primers themselves
   are at fault

**Quick check from the design output:**

```python
from scripts.check_dimers import analyze_dimers

dimers = analyze_dimers(
    [pair["forward_seq"], pair["reverse_seq"]],
    temperature=60.0,
)
for d in dimers["interactions"]:
    if d["problematic"]:
        print(f"{d['type']}: ΔG = {d['dg']} kcal/mol — redesign")
```

A heterodimer with ΔG < −5 kcal/mol at the annealing temperature will form
in the tube. Below −9 kcal/mol, the dimer dominates over the on-target
product.

---

## qPCR-Specific Issues

### Efficiency Outside 0.90–1.10

| E | Slope | Likely cause | Action |
|---|-------|--------------|--------|
| < 0.90 | > −3.6 | Inhibitors, secondary structure, poor primer | Dilute template; redesign primers |
| 0.90–1.10 | −3.6 to −3.1 | OK | None |
| > 1.10 | < −3.1 | Primer dimers, contamination, pipetting variance | Inspect melt curve sub-80°C peak; clean workspace |

Run a fresh standard curve before redesigning. Pipetting variance shifts
slope without it being the primers' fault.

### High Cq Variation (CV > 5%)

Common causes are pipetting errors, bubbles, and template near the limit
of detection.

- Use electronic pipettes for volumes < 5 µL
- Spin plates briefly before cycling
- Avoid plate edges (evaporation)
- If Cq > 30, run technical quadruplicate or higher

### NTC Amplification

A no-template control that amplifies (Cq < 35) means contamination or
primer dimers.

- Separate pre- and post-PCR workspaces
- Use filter tips
- Aliquot every reagent (no pipetting from stock bottles)
- UV-treat bench and pipettes
- Inspect the melt curve — a peak under 80°C points to dimers, not
  contamination

### Multiple Melt-Curve Peaks

| Pattern | Cause | Action |
|---------|-------|--------|
| Main peak + low-T peak (< 80°C) | Primer dimer | Tighten Ta or redesign |
| Two similar-height peaks | Off-target amplification | Redesign for specificity |
| Broad single peak | Mixed products or SNP heterogeneity | Sequence the product; check for SNPs |

---

## Template-Specific Issues

### GC-Rich (> 65%)

- Add DMSO 5–10% or betaine 1–2 M
- Denaturation 98°C for 30–60 s
- Use a GC-optimized polymerase
- Design primers in AT-rich islands when possible

### AT-Rich (< 35%)

- Extend primers to 22–25 nt to raise Tm
- Lower Ta to match the lower primer Tm
- Hot-start polymerase reduces non-specific products

### Repetitive Regions

- Run RepeatMasker first
- Design primers in unique flanking sequence
- Consider nested PCR

### RNA Templates (RT-qPCR)

| Problem | Cause | Action |
|---------|-------|--------|
| Cq in −RT control | gDNA contamination | DNase treat; redesign primers to span exon junction |
| Poor RT efficiency | Low RIN, secondary structure | Verify RIN > 7; use random hexamer + oligo-dT mix |
| Variable Cq between days | RNA degrading during prep | Minimize freeze-thaw; store at −80°C; work on ice |

---

## Diagnostic Flow Reference

### Failure Cascade

```
PCR fails
├── Run positive control
│    ├── Fails → reagent or instrument problem
│    └── Works → primer/template problem:
│        ├── Gradient PCR finds optimum?
│        │    ├── Yes → use that Ta
│        │    └── No → redesign primer
│        ├── BLAST primers → mismatch found?
│        │    └── Yes → redesign at new site
│        └── Template gel → high MW band?
│             └── No → re-extract
```

### Non-Specific Cascade

```
Multiple products
├── Gel pattern:
│    ├── Band < 100 bp → primer dimer
│    ├── Expected + extra → off-target
│    └── Smear → degraded template or over-amplification
├── BLAST primers → off-targets with < 3 mismatches?
│    └── Yes → redesign
├── Raise Ta 3–5°C → resolves it?
│    ├── Yes → use that Ta
│    └── No → reduce primer concentration; redesign
```

---

## First-Step Recipes by Symptom

**No amplification:**
1. Gradient PCR
2. Verify template (Qubit + gel)
3. BLAST primers

**Weak signal:**
1. Increase template 5×
2. Lower Ta by 3°C
3. Bump Mg²⁺ to 2.5 mM

**Non-specific products:**
1. Raise Ta by 3–5°C
2. Lower primer concentration to 0.2 µM
3. Hot-start polymerase

**Primer dimers:**
1. Raise Ta
2. Lower primer concentration
3. Redesign

**Poor qPCR efficiency:**
1. Dilute template 1:10
2. Sweep primer concentration (0.2, 0.3, 0.5 µM)
3. Build a fresh standard curve before concluding the primers are at fault

---

Last updated: 2026-05-20.
