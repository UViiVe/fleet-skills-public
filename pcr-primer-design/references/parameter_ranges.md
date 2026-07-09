# Parameter Ranges for PCR Primer Design

A practical look-up table for the numbers you need when overriding the default
`design_*_primers` calls. Use this when the standard workflow won't return
enough candidates and you need to know which knob to turn first.

← Back to [SKILL.md](../SKILL.md)
Related: [Best Practices](primer_design_best_practices.md) · [MIQE Guidelines](miqe_guidelines.md) · [Troubleshooting](troubleshooting_guide.md)

---

## Primer Length

| Application | Sweet spot | Acceptable window | Why |
|-------------|-----------|-------------------|-----|
| qPCR | 18–22 nt | 18–24 nt | Short primers anneal faster; helps kinetics on a 40-cycle qPCR |
| Standard PCR | 20–22 nt | 18–25 nt | Balance specificity and yield |
| Sanger sequencing | 20–24 nt | 18–26 nt | Need clear initiation, so err longer |
| Multiplex | 20–22 nt | 18–24 nt | Keep similar across the panel |
| Allele-specific (SNP) | 18–25 nt | 18–28 nt | The allele-specific arm often needs an extra base for discrimination |

Going below 18 nt costs specificity. Going above 25 nt invites hairpins and
slows annealing.

## Melting Temperature (Tm)

The numbers below assume the SantaLucia 1998 nearest-neighbor model with
50 mM Na⁺ and no added Mg²⁺. Each 1 mM of Mg²⁺ raises Tm by roughly 1°C in the
Owczarzy correction.

| Application | Tm range | Forward/reverse ΔTm | Notes |
|-------------|----------|--------------------|-------|
| Standard PCR | 55–65°C | ≤ 5°C | Forgiving |
| qPCR | 58–62°C | ≤ 2°C | MIQE essential |
| qPCR (strict) | 59–61°C | ≤ 1°C | Best practice when running plates with multiple targets |
| Multiplex | 58–62°C | ≤ 3°C across all primers | A loose window costs specificity for some targets |
| TaqMan probe | 65–70°C | +5 to +10°C above primers | Probe must be bound when primers extend through it |
| Sequencing | 55–60°C | n/a | Single primer per reaction |

## GC Content

| GC% | Category | What it means |
|-----|----------|----------------|
| < 30% | Very AT-rich | Weak binding, low Tm — extend primer instead |
| 30–40% | AT-rich | Workable; raise primer length to compensate |
| 40–60% | Optimal | No special handling |
| 60–65% | GC-rich | Watch hairpin ΔG |
| > 65% | Very GC-rich | Redesign or add DMSO/betaine to the reaction |

## Amplicon Size

| Application | Optimal | Hard ceiling |
|-------------|---------|--------------|
| qPCR (SYBR) | 70–140 bp | 200 bp |
| qPCR (TaqMan) | 70–140 bp | 150 bp |
| Standard PCR | 200–600 bp | 3–5 kb with high-fidelity polymerase |
| Multiplex (gel) | 100–600 bp, ≥ 50 bp spacing between products | 1 kb |
| Multiplex (qPCR, multi-color) | 70–140 bp per pair | n/a — separated by fluorophore |
| Sanger sequencing | 400–800 bp from primer | 1 kb |

For FFPE-extracted DNA or RNA from blood collected on filter paper, shorten
qPCR amplicons to 60–90 bp — the template is fragmented and longer products
won't amplify.

## Thermodynamic Thresholds

### Primer dimers (ΔG at the annealing temperature)

| ΔG | Verdict |
|----|---------|
| > −5 kcal/mol | acceptable |
| −5 to −9 kcal/mol | marginal — monitor with melt curve |
| < −9 kcal/mol | redesign |

### Hairpins

| ΔG | Verdict |
|----|---------|
| > −2 kcal/mol | acceptable |
| −2 to −4 kcal/mol | marginal |
| < −4 kcal/mol | redesign |

### 3' end self-complementarity

| Complementary 3' bases | Verdict |
|------------------------|---------|
| < 8 | acceptable |
| 8–10 | marginal |
| > 10 | redesign |

### Total self-complementarity

| Total complementary bases | Verdict |
|--------------------------|---------|
| < 12 | acceptable |
| 12–16 | marginal |
| > 16 | redesign |

## Salt and Reagent Concentrations

| Component | Standard PCR | qPCR | High-fidelity PCR |
|-----------|-------------|------|-------------------|
| Monovalent (Na⁺ + K⁺) | 50 mM | 50 mM | 50 mM |
| Divalent (Mg²⁺) | 1.5–2.5 mM | 3.0–5.0 mM (in master mix) | 2.0–4.0 mM |
| dNTPs | 0.2–0.8 mM | 0.8 mM | 0.2 mM |
| Primer per reaction | 0.2–1.0 µM | 0.2–0.5 µM | 0.2–0.5 µM |
| TaqMan probe | n/a | 100–250 nM | n/a |

When you change Mg²⁺ between design and the actual reaction, re-calculate the
Tm. The `calculate_tm.calculate_tm_comprehensive` function accepts `mg_conc`
and `dntp_conc` and applies the Owczarzy correction automatically.

## Application Cheat Sheets

### Standard PCR

```
Primer length:     18–25 nt (opt 20–22)
Tm:                55–65°C
ΔTm:               ≤ 5°C
GC%:               40–60%
Amplicon:          100–1000 bp (opt 200–600)
Max poly-X:        4
GC clamp:          1–2 G/C in last 5 bases
```

Default cycling:
```
95°C 3 min initial denature
[ 95°C 30 s | annealing 30 s | 72°C 60 s/kb ] × 25–35 cycles
72°C 5 min final extension
```

### qPCR (MIQE)

```
Primer length:     18–22 nt (opt 20)
Tm:                58–62°C (opt 59–61°C)
ΔTm:               ≤ 2°C (strict: ≤ 1°C)
GC%:               40–60%
Amplicon:          70–140 bp
Max poly-X:        4
GC clamp:          1
3' end:            no mismatch tolerated
Exon-spanning:     yes, or intron > 1 kb
```

Validation pass/fail thresholds:
```
Efficiency E:      0.90–1.10
R²:                > 0.98
Linear range:      ≥ 5 logs
NTC Cq:            > 35 (or undetected)
−RT control Cq:    > target Cq + 5
Technical CV:      < 5%
```

### TaqMan

Use the qPCR primer parameters above and additionally:

```
Probe length:      18–30 nt (opt 20–25)
Probe Tm:          65–70°C
Probe GC%:         40–60%
No G at 5' end:    avoid (quenches the reporter)
Position:          between primers, closer to forward
```

### Multiplex

```
Primer length:     18–24 nt (all similar)
Tm:                58–62°C
ΔTm across panel:  ≤ 3°C
GC%:               40–60%
Amplicon spacing:  ≥ 50 bp for gel detection
Pairwise dimers:   ΔG > −5 kcal/mol for every primer pair, including cross-pairs
```

Validate by computing dimers for the *flat list* of every primer in the panel,
not just within pairs. `check_dimers.analyze_multiplex_compatibility` does
this automatically.

## Reaction Optimization Knobs

### Annealing temperature strategy

| Strategy | Ta | When |
|----------|----|------|
| Standard | Tm − 5°C | Default starting point |
| Stringent | Tm | Eliminating non-specific bands |
| Permissive | Tm − 7°C | No amplification at standard Ta |
| Gradient | Tm ± 5°C | Always run for new primer pairs |
| Touchdown | Tm + 5°C decreasing by 0.5°C/cycle | Multiple non-specific products |

### Mg²⁺ sweep

| [Mg²⁺] | Effect |
|--------|--------|
| 1.0–1.5 mM | Higher specificity |
| 1.5–2.5 mM | Standard |
| 2.5–4.0 mM | Higher yield, more secondary products |
| > 4.0 mM | Lower fidelity — not recommended |

### Cycle number

| Application | Cycles |
|-------------|--------|
| qPCR | 40–45 (stop when plateau reached) |
| Standard PCR | 25–35 |
| Colony PCR | 30–35 |
| Nested PCR | 20–25 + 20–30 |

## Tm Calculation Methods

Choose the model that matches the reaction you're actually running.

### Nearest-neighbor (default; SantaLucia 1998)

Most accurate. Used by `primer3.calcTm` and Biopython's `Tm_NN(saltcorr=7)`.

When to use: every qPCR design, anything published, anything you'll lock in.

### Salt-adjusted Wallace

```
Tm = 2(A + T) + 4(G + C) − 16.6·log₁₀[Na⁺] + 0.41·%GC
```

When to use: quick mental math, sanity-checking a screen.

### Basic Wallace

```
Tm = 2(A + T) + 4(G + C)
```

When to use: primers ≤ 13 nt where nearest-neighbor is overkill. Don't use for
qPCR.

### %GC formula

```
Tm = 64.9 + 41·(G + C − 16.4) / (A + T + G + C)
```

When to use: primers > 13 nt, quick estimate, no salt data available.

## Primer3 Parameter Reference

The wrapper scripts forward unknown parameters via `**kwargs`. Useful
parameters worth knowing:

```python
PRIMER_OPT_SIZE              = 20
PRIMER_MIN_SIZE              = 18
PRIMER_MAX_SIZE              = 25
PRIMER_OPT_TM                = 60.0
PRIMER_MIN_TM                = 58.0
PRIMER_MAX_TM                = 62.0
PRIMER_PAIR_MAX_DIFF_TM      = 2.0
PRIMER_MIN_GC                = 40.0
PRIMER_MAX_GC                = 60.0
PRIMER_PRODUCT_SIZE_RANGE    = [(70, 140)]
PRIMER_MAX_POLY_X            = 4
PRIMER_GC_CLAMP              = 1
PRIMER_MAX_SELF_ANY          = 8
PRIMER_MAX_SELF_END          = 3
PRIMER_PAIR_MAX_COMPL_ANY    = 8
PRIMER_PAIR_MAX_COMPL_END    = 3
PRIMER_MAX_HAIRPIN_TH        = 47.0   # max hairpin Tm in °C
PRIMER_MAX_NS_ACCEPTED       = 0      # disallow ambiguous bases
```

For stricter qPCR designs, tighten `PRIMER_PAIR_MAX_DIFF_TM` to 1.0 and
`PRIMER_MAX_SELF_END` to 2. For standard PCR on loose templates, widen
`PRIMER_PRODUCT_SIZE_RANGE` to `[(100, 1000)]` and bump
`PRIMER_PAIR_MAX_DIFF_TM` to 5.0.

## Decision Triggers

### Increase Ta (annealing temperature)
- Multiple bands on gel
- Multiple melt-curve peaks
- High NTC signal
- General non-specific background

### Decrease Ta
- No amplification with primers that BLAST perfectly
- Very weak signal
- Cq > 35 on a clearly expressed gene

### Redesign primers
- ΔTm > 2°C (or > 1°C for qPCR)
- Dimer ΔG < −5 kcal/mol
- Hairpin ΔG < −2 kcal/mol
- Off-target hits with ≤ 3 mismatches
- qPCR efficiency outside 90–110% after optimization
- Multiple peaks on melt curve after tightening Ta

---

Last updated: 2026-05-20.
