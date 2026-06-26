# TCGA-LIHC Ketone Body Metabolism Analysis

R scripts supporting the exploratory TCGA-LIHC transcriptomic analysis associated with:

> Ketone body metabolism across the MASLD-to-HCC continuum: cell-type-specific mechanisms, evidence gaps, and a translational roadmap

## Scope and interpretation

TCGA-LIHC is used here only as mixed-etiology supportive / broad HCC exploratory evidence.
These analyses should not be interpreted as MASLD-HCC-specific proof, and no output from this repository should be interpreted as a clinical recommendation.

## Repository structure

```
TCGA-LIHC-ketone-metabolism/
├─ scripts/
│  └─ TCGA_LIHC_ketone_analysis.R
├─ output_pdf/
│  ├─ Fig5H_CTNNB1_OXCT1_vector.pdf
├─ .gitignore
├─ LICENSE
└─ README.md

Local-only directories required for rerunning the analysis:
├─ xena/
│  ├─ TCGA-LIHC.star_fpkm.tsv.gz
│  ├─ TCGA-LIHC.survival.tsv.gz
│  └─ TCGA-LIHC.clinical.tsv.gz
├─ GDCdata/
├─ data/
└─ output/
```

## Input data

Three UCSC Xena files are required and should be placed in the `xena/` directory:

| File | Purpose |
|---|---|
| `TCGA-LIHC.star_fpkm.tsv.gz` | RNA-seq expression matrix, log2(FPKM + 1) |
| `TCGA-LIHC.survival.tsv.gz` | overall survival data |
| `TCGA-LIHC.clinical.tsv.gz` | clinical covariates |

Somatic mutation data are accessed through `TCGAbiolinks` and cached locally under `GDCdata/`.

## Main outputs

| Output file | Figure | Description |
|---|---|---|
| `output_pdf/Fig5H_CTNNB1_OXCT1_vector.pdf` | Fig. 5H | OXCT1 expression by CTNNB1 mutation status |
## How to run

Open RStudio, set the working directory to the repository root, and run:

```r
# Replace this with your local clone of the repository
setwd("path/to/TCGA-LIHC-ketone-metabolism")

source("scripts/TCGA_LIHC_ketone_analysis.R", echo = TRUE)
```

The script uses project-relative paths, so users can replace the `setwd()` path with their own local repository path.

## Optional sections

Sections 7-9 of the script generate Fig. 5E-G and require ESTIMATE-compatible gene-symbol input.
By default, these optional sections are switched off:

```r
run_sections_7_to_9 <- FALSE
```

Set this value to `TRUE` only after preparing an ESTIMATE-compatible expression matrix.

## R environment

The analysis was prepared in R 4.5.2.

Required R packages include:

```r
# Bioconductor
BiocManager::install("TCGAbiolinks")

# CRAN
install.packages(c(
  "dplyr", "tidyr", "ggplot2", "ggpubr", "survminer",
  "survival", "maxstat", "estimate", "patchwork",
  "mice", "broom"
))
```

## License

MIT
