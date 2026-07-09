# Software Requirements

Installation and license summary for every package referenced by this skill.

---

## R version

- Minimum: R ≥ 4.0.
- Recommended: R ≥ 4.2.
- RStudio is optional; the scripts work from any R interface.

---

## Required packages

### Bioconductor

| Package | Version | License | Purpose |
|---|---|---|---|
| RNASeqPower | ≥ 1.36 | GPL-2+ | Bulk RNA-seq power and sample size |
| ssizeRNA | ≥ 1.4 | GPL-2+ | Sample size with FDR control |
| DESeq2 | ≥ 1.40 | LGPL | Dispersion estimation from pilot data |
| powsimR | ≥ 1.2.3 | GPL-3 | scRNA-seq simulation-based power |
| OSAT | ≥ 1.38 | Artistic-2.0 | Optimal sample-to-batch assignment |
| IHW | ≥ 1.28 | Artistic-2.0 | Independent hypothesis weighting |
| qvalue | ≥ 2.32 | LGPL | Storey's q-value |
| sva | ≥ 3.48 | Artistic-2.0 | Surrogate variable analysis / ComBat |
| pasilla | ≥ 1.18 | Artistic-2.0 | Example RNA-seq data |

### CRAN

| Package | Version | License | Purpose |
|---|---|---|---|
| ggplot2 | ≥ 3.4 | MIT | Plotting |
| ggprism | ≥ 1.0.4 | GPL-3 | Publication-quality themes |
| pwr | ≥ 1.3 | GPL-3+ | General-purpose power calculators |
| jsonlite | ≥ 1.8 | MIT | JSON export |

### Optional advanced packages

| Package | License | When you'd want it |
|---|---|---|
| PROPER | GPL-2+ | Simulation-based RNA-seq power |
| ashr | GPL-3+ | Adaptive shrinkage for fold-change |
| locfdr | GPL-2+ | Local FDR estimates |
| blockTools | GPL-2+ | Incomplete block designs |
| SingleCellExperiment | GPL-3 | scRNA-seq pilot data structures |

---

## Quick install

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "RNASeqPower", "ssizeRNA", "DESeq2", "powsimR",
  "OSAT", "IHW", "qvalue", "sva", "pasilla"
), update = FALSE, ask = FALSE)

install.packages(c("ggplot2", "ggprism", "pwr", "jsonlite"),
                 repos = "https://cran.r-project.org")
```

Verify:

```r
pkgs <- c("RNASeqPower","ssizeRNA","DESeq2","OSAT","IHW","qvalue","ggplot2")
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
stopifnot(all(ok))
cat("✓ Core packages installed\n")
```

---

## License compliance

All bundled packages permit commercial use. Summary:

- **GPL/LGPL (GPL-2+, GPL-3, LGPL):** source-sharing obligations apply only if you redistribute binaries. Running them inside an AI agent does not create a redistribution event.
- **Artistic-2.0 (Bioconductor default):** permissive; commercial use OK.
- **MIT:** permissive with attribution.

Cite the primary papers when methods appear in publications or reports — see the SKILL.md References section for the canonical citations.

---

## Common installation issues

### BiocManager itself fails

```r
install.packages("BiocManager", repos = "https://cran.r-project.org")
```

### Bioconductor version mismatch

```r
BiocManager::install(version = "3.18")  # adjust to current release
```

### powsimR installation is slow

Expected — 10–30 minutes due to many dependencies. Run once and cache.

### macOS compilation errors

```bash
xcode-select --install
```

### Linux compilation errors

Debian / Ubuntu:

```bash
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
```

RHEL / CentOS:

```bash
sudo yum install libcurl-devel openssl-devel libxml2-devel
```

---

## System resources

| Resource | Minimum | Recommended |
|---|---|---|
| RAM | 4 GB | 8 GB (more for powsimR simulations) |
| Disk | 2 GB | 5 GB |
| OS | macOS / Linux / Windows 10+ | macOS / Linux |

---

## Getting help

- Package-specific: `?functionName` in R, or Bioconductor support at support.bioconductor.org.
- Workflow-specific: [troubleshooting_guide.md](troubleshooting_guide.md).
- Hard statistical questions: bring in a statistician at the design stage.
