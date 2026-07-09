# Genome assembly conversion — GRCh37 ↔ GRCh38

Every annotation pipeline must agree on a single genome assembly. Mixing
GRCh37-coordinate VCFs with GRCh38 annotation caches is the single most
common source of "everything annotates as intergenic" failures. This
guide covers lift-over options, the difference between position-based
and name-based references, the chromosome-name harmonisation step, and
the situations where lift-over is unsafe.

---

## The two human reference assemblies in the wild

| Assembly | Released | Aliases | Notes |
|---|---|---|---|
| GRCh37 | 2009 | hg19, b37, NCBI build 37 | Legacy; some clinical pipelines still in flight; 1000 Genomes Phase 3 native |
| GRCh38 | 2013 | hg38, NCBI build 38 | Current standard; fixed many GRCh37 gaps and misassemblies; gnomAD v3+ and v4 are native |
| T2T-CHM13 | 2022 | T2T | Telomere-to-telomere reference; not yet widely supported by annotation databases |

GRCh38 is the right default for any new pipeline. Continue working in
GRCh37 only when the downstream consumer (a clinical reporting system,
a legacy database, a published reference panel) is fixed on it.

The differences between GRCh37 and GRCh38 are extensive in detail but
trivial in summary: GRCh38 fixed an estimated 100 K errors, filled gaps,
and updated centromeric assemblies. Functionally, the same genes are at
nearly the same positions — but "nearly" means most positions shift by a
few bases to a few kilobases.

---

## Lift-over tools

| Tool | Strengths | Notes |
|---|---|---|
| UCSC `liftOver` | Original; widely supported; chain files freely available | Works on BED / VCF / GFF; conservative |
| `picard LiftoverVcf` | VCF-native; preserves INFO/FORMAT; emits rejected records | Most common for VCFs |
| `bcftools +liftover` (plugin) | Lightweight; integrates with bcftools | Newer; check version availability |
| CrossMap | Python-native; supports many formats | Slower but handy in Python pipelines |

A standard `picard LiftoverVcf` invocation:

```bash
picard LiftoverVcf \
    I=input.GRCh37.vcf.gz \
    O=output.GRCh38.vcf.gz \
    REJECT=rejected.GRCh37.vcf.gz \
    CHAIN=hg19ToHg38.over.chain.gz \
    R=GRCh38.fa
```

Chain files (`hg19ToHg38.over.chain.gz` and reverse) are downloadable
from UCSC. The reference FASTA for the target assembly is mandatory
because lift-over verifies the REF allele matches the new assembly.

**Always inspect the `rejected.vcf.gz` output.** Even on a clean WES
dataset, 0.1–1 % of variants typically fail to lift over. Reasons:

- Variant lands on a region that was reassembled (e.g. a centromere) and
  has no unambiguous mapping.
- REF allele does not match the target reference at the mapped position
  (an actual error somewhere upstream, or a strand flip).
- The chain file's resolution is insufficient (rare; the available
  GRCh37 ↔ GRCh38 chains are well-developed).

---

## Position-based vs name-based references

A subtle but important distinction:

- **Position-based** references are coordinates: `chr1:155234567-155234567 A>G`.
  These move under lift-over. Two pipelines must agree on the assembly.
- **Name-based** references are identifiers: `rs7412`, `NM_000041.4:c.526C>T`,
  `BRCA1:p.Cys64Trp`. These do not move under lift-over — the rsID for
  a SNP is the same in GRCh37 and GRCh38.

Implications:

- **dbSNP rsIDs are assembly-agnostic** in their identity, but their
  resolved coordinate is assembly-specific. When you query dbSNP you
  must declare which assembly you want.
- **HGVS notation against a specific transcript** (e.g. `NM_000041.4:
  c.526C>T`) is assembly-agnostic provided the transcript itself is
  stable across assemblies. Most clinical HGVS strings are transcript-
  based for this reason.
- **HGVS notation against the genome** (e.g. `NC_000001.10:g.155234567A>G`
  with build-specific accession) is assembly-specific. The accession
  (`.10` vs `.11`) encodes the assembly.

For clinical and database work, prefer transcript-based HGVS over
position-based coordinates whenever the underlying call is on a coding
variant.

---

## Chromosome-name harmonisation

Two naming conventions are in active use:

- **UCSC style** — `chr1`, `chr2`, ..., `chrX`, `chrY`, `chrM`.
- **Ensembl / 1000 Genomes style** — `1`, `2`, ..., `X`, `Y`, `MT`.

The two conventions are not interchangeable. VEP errors out when the VCF
contig names do not match the cache's contig names; SnpEff silently
mis-annotates as intergenic. Convert with `bcftools annotate`:

```bash
# Build a rename map
cat > chr_rename.txt <<'TXT'
1   chr1
2   chr2
X   chrX
Y   chrY
MT  chrM
TXT

bcftools annotate --rename-chrs chr_rename.txt input.vcf.gz -Oz -o renamed.vcf.gz
tabix -p vcf renamed.vcf.gz
```

The skill's `validate_vcf.py` includes a chromosome-name check;
mismatches surface as a validation warning before annotation starts.

`chrM` vs `MT` is a separate gotcha: GRCh37 uses `MT` (16,569 bp,
NC_012920) while GRCh38 ships with both `chrM` (the same sequence) and
`chrMT` in some accessions. For mitochondrial variant calling, verify
the contig accession in the FASTA header before running.

---

## When lift-over is unsafe

Lift-over works well for the ~99 % of the genome with stable
coordinates. It is unsafe for:

- **Centromeric and telomeric repeats** — extensively reassembled in
  GRCh38. Many GRCh37 calls in these regions have no clean mapping.
- **Segmental duplications and paralogous regions** — multi-mapping
  reads in the original alignment may have collapsed to different
  positions in the two assemblies.
- **Indels of length > ~50 bp** — lift-over uses chain blocks of fixed
  resolution; large indels can fall across block boundaries and
  ambiguous positions.
- **Newly added or deleted regions** — GRCh38 added new alt contigs
  (HLA, KIR) and split a few previously-collapsed loci.

For these regions, prefer to **re-call variants on the target assembly**
rather than lift over. Lift-over is for SNVs and short indels in the
straightforward majority of the genome.

---

## Database-version pitfalls

Even after a successful lift-over, downstream annotations must use the
matching assembly:

- **gnomAD** — v2 is GRCh37; v4 is GRCh38. Mixing is the most common
  cause of "no frequency for a known variant".
- **ClinVar** — both assemblies supplied; pick the matching one.
- **CADD** — both supplied; the GRCh38 file is the current one.
- **REVEL** — both supplied.
- **AlphaMissense** — GRCh38 only at time of writing.
- **dbSNP** — both supplied; query the build-matched VCF.
- **TOPMed Bravo** — GRCh38 only.

Pin every database to the same assembly as the VCF before annotation
starts. The skill's installation step downloads matched-assembly caches
by default; if you override the cache path, manually verify the match.

---

## A clean workflow

1. **Declare the assembly upfront.** Read the VCF `##reference=` header
   line; if it is missing or ambiguous, ask the user.
2. **Validate chromosome naming.** Run `validate_vcf()` to check
   `chr` prefix consistency.
3. **If conversion is required**, choose between:
   - **Lift-over** (`picard LiftoverVcf` with the appropriate chain
     file) — fast, but verify rejected.vcf.gz and reject any pipeline
     where > 1 % of variants fail.
   - **Re-call** on the target assembly — slower but the right choice
     for unstable regions or large indels.
4. **Verify post-conversion.** Spot-check 5–10 known variants by their
   gene-position-allele triple against both assemblies (e.g.
   `BRCA1:c.5266dupC` should land at the same protein position).
5. **Match downstream annotations** — cache, plugin data files,
   custom-loaded ClinVar / gnomAD — to the target assembly.

The whole sequence takes 30 minutes on a well-prepared pipeline and
saves days of debugging downstream.
