# Ensembl VEP — pipeline reference

This guide collects the Ensembl Variant Effect Predictor knobs that matter
in practice: installation, cache management, the plugins that change a
report's clinical usefulness, the most common command-line patterns, and
the schema of the annotations VEP writes into the VCF.

The skill's `run_vep.py` wrapper handles most of these for you. Use this
reference when you need to bend the wrapper to a non-standard situation
(custom plugin, offline cache, exotic species).

---

## Installation

The conda recipe is the path of least resistance. The Perl installer is
necessary only when conda is unavailable or you need to pin an exact
Ensembl release.

```bash
# Recommended: conda
conda create -n vep python=3.10 -y
conda activate vep
conda install -c bioconda ensembl-vep -y
vep --help                # sanity check
```

```bash
# Manual install (Perl)
git clone https://github.com/Ensembl/ensembl-vep.git
cd ensembl-vep
perl INSTALL.pl           # interactive; choose cache + plugin options
```

The wrapper script `install_tools.py` in this skill drives the conda path
non-interactively.

---

## Cache management

VEP is a thin wrapper around a per-species annotation cache. Without the
cache it falls back to live Ensembl REST queries, which are 10–100×
slower and have lower coverage. **Always run with `--cache`.**

```bash
# GRCh38 cache (current human reference)
vep_install -a c -s homo_sapiens -y GRCh38

# GRCh37 cache (legacy clinical pipelines still in flight)
vep_install -a c -s homo_sapiens -y GRCh37

# Non-human example
vep_install -a c -s mus_musculus -y GRCm39
```

Cache disk requirements: human ~15–20 GB per assembly, mouse ~8–10 GB,
zebrafish ~5 GB. The default location is `~/.vep/`; override with
`--dir_cache /custom/path` if you need to share the cache across users.

For offline clusters, build the cache on an internet-facing node and then
`rsync` the directory across — no environment variables are needed; VEP
finds it via the `--dir_cache` argument alone.

---

## Common command-line patterns

### Basic annotation

```bash
vep -i input.vcf \
    --cache --assembly GRCh38 \
    --vcf -o annotated.vcf \
    --force_overwrite
```

### Comprehensive (use this 90 % of the time)

```bash
vep -i input.vcf \
    --cache --assembly GRCh38 \
    --everything \
    --vcf -o annotated.vcf \
    --force_overwrite --fork 4
```

`--everything` is shorthand for a curated bundle: HGVS notation, protein
domains (Pfam + InterPro), regulatory features, existing dbSNP / COSMIC
IDs, population frequencies (gnomAD + 1000G + ESP), and SIFT / PolyPhen
predictions. It is the right default for any general-purpose analysis.

### Clinical-grade

```bash
vep -i input.vcf \
    --cache --assembly GRCh38 \
    --everything \
    --plugin CADD,/data/CADD.tsv.gz \
    --plugin REVEL,/data/revel_scores.tsv.gz \
    --plugin AlphaMissense,/data/alphamissense.tsv.gz \
    --custom ClinVar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
    --max_af gnomAD \
    --vcf -o annotated.vcf \
    --force_overwrite --fork 4
```

---

## The CSQ INFO field

VEP packs every annotation into the VCF `INFO/CSQ` field. The schema is
defined in the VCF header `##INFO=<ID=CSQ,...>` line, which lists the
ordered pipe-delimited fields. Always discover the schema from the header
of the file you are parsing — do not hard-code it.

Typical fields, in the order VEP emits them with `--everything`:

```
Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|
EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|
Amino_acids|Codons|Existing_variation|REF_ALLELE|UPLOADED_ALLELE|
DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|
MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRAIS|CCDS|ENSP|SWISSPROT|TREMBL|
UNIPARC|UNIPROT_ISOFORM|REFSEQ_MATCH|SOURCE|REFSEQ_OFFSET|GENE_PHENO|
SIFT|PolyPhen|DOMAINS|miRNA|HGVS_OFFSET|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|
SAS_AF|gnomADe_AF|gnomADe_AFR_AF|...|MAX_AF|MAX_AF_POPS|CLIN_SIG|
SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|
CADD_PHRED|CADD_RAW|REVEL|AlphaMissense_score|AlphaMissense_class
```

Multiple annotations per variant are separated by commas inside CSQ
because a single variant can hit multiple transcripts. Downstream tooling
must explode CSQ on `,` first, then on `|`.

---

## Performance knobs

| Knob | Effect | Recommendation |
|---|---|---|
| `--fork N` | N parallel workers | 4–8 for whole-exome; 8–16 for whole-genome on 32 GB+ RAM |
| `--buffer_size B` | Variants per read batch | Default 5000; reduce to 1000 on memory-constrained nodes |
| `--cache` | Local cache lookup | Mandatory; never skip |
| `--offline` | No Ensembl REST fallback | Add for clusters with no outbound network |
| `--quiet` | Suppress per-variant warnings | Use for very large VCFs |
| `--no_stats` | Skip summary HTML | Saves a few seconds on huge files |

Whole-exome (~25 K variants) with `--fork 4 --everything` plus CADD plus
ClinVar custom: ~10 min on a workstation. Whole-genome (~5 M variants):
~3–6 h. Splitting the VCF by chromosome with `bcftools view -r` and
fanning out across nodes is the standard scaling pattern.

---

## Recommended plugins

Plugins extend VEP with annotations that are not in the core cache.
Install once with `vep_install -a p -g PluginName` and download the
referenced data files separately.

| Plugin | What it adds | When to enable |
|---|---|---|
| CADD | Combined Annotation Dependent Depletion score | Always; mid-tier prioritisation |
| REVEL | Missense ensemble pathogenicity score | Clinical or rare disease |
| AlphaMissense | DeepMind missense pathogenicity classification | Recommended; trained on protein structures |
| dbNSFP | Bundles SIFT, PolyPhen, MutationTaster, etc. | When you need many predictors in one shot |
| SpliceAI | Splice-altering variants | RNA splicing-relevant phenotypes |
| MaxEntScan | Splice site scoring | Same as SpliceAI when you need explicit donor/acceptor scores |
| LoFtool | Gene-level LoF tolerance | Augments gnomAD constraint |
| LOEUF | gnomAD constraint metric | Same; preferred over LoFtool for human |
| GeneSplicer | Splice prediction | Cancer / rare disease |
| Mastermind | Curated literature-evidence | Clinical reporting (commercial-friendly licence to verify) |

---

## Output formats

- `--vcf` (default for downstream tools) — keeps the VCF as a VCF with
  annotations folded into `INFO/CSQ`. Recommended for all programmatic
  pipelines.
- `--tab` — human-readable tab-delimited table. Good for spot-checks.
- `--json` — one record per variant; convenient when downstream code is
  Python and you want to skip a CSQ parser.

If you need a flat table for spreadsheet review, run `--vcf` first and
let the skill's `parse_vep_output.py` produce the CSV — it explodes CSQ
correctly and emits a tidy per-variant per-transcript row.

---

## Common pitfalls

- **`chr1` vs `1` chromosome naming.** VEP errors out on mismatch. Convert
  with `bcftools annotate --rename-chrs` before annotating.
- **Mismatched assembly.** Annotating a GRCh37-coordinate VCF against the
  GRCh38 cache returns nearly everything as `intergenic_variant`. Confirm
  the assembly with the VCF `##reference=` header line.
- **API mode (`--database`)** is documented as a fallback but is unfit
  for production: 10–100× slower and partial annotation coverage. Use the
  cache.
- **Stale cache.** Annotation conventions evolve; rebuild the cache
  yearly. Pin the Ensembl release with `--cache_version N` for reproducible
  clinical reports.
- **Reference-allele mismatches.** Pass `--fasta /path/genome.fa` to make
  VEP verify the VCF REF column against the reference; mismatches almost
  always indicate an upstream pipeline bug.
