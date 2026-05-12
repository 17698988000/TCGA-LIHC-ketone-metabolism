# TCGA-LIHC Ketone Body Metabolism Analysis

R scripts for the TCGA-LIHC transcriptomic analysis associated with:

> Xu G, Tang B. **Dysregulation of ketone body metabolism in the progression from metabolic dysfunction-associated steatotic liver disease to hepatocellular carcinoma: mechanisms and translational opportunities.** (2026)

---

## What this repository contains

`TCGA_LIHC_ketone_analysis.R` — complete analysis pipeline including:

- Differential expression of HMGCS2, OXCT1, BDH1 (tumor vs. adjacent normal)
- Kaplan–Meier survival analysis and Cox regression (univariable + multivariable)
- Maxstat optimal cutpoint sensitivity analysis
- ESTIMATE immune infiltration scoring
- T-cell exhaustion marker analysis (FDR-corrected)
- CTNNB1 mutation subgroup analysis
- Multiple imputation sensitivity analysis for Cox regression (mice, m = 20)

---

## Requirements

**R version:** 4.5.2

**Packages:**
```r
# Bioconductor
BiocManager::install("TCGAbiolinks")

# CRAN
install.packages(c("dplyr", "tidyr", "ggplot2", "ggpubr", "survminer",
                   "survival", "maxstat", "estimate", "patchwork",
                   "mice", "broom"))
```

---

## Data

Three files must be downloaded manually from **UCSC Xena** before running the script:

| File | Source | Size |
|------|--------|------|
| `TCGA-LIHC.star_fpkm.tsv.gz` | [STAR-FPKM (n=424)](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.star_fpkm.tsv.gz) | ~150 MB |
| `TCGA-LIHC.survival.tsv.gz` | [Survival data (n=433)](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.survival.tsv.gz) | <1 MB |
| `TCGA-LIHC.clinical.tsv.gz` | [Clinical data (n=439)](https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.clinical.tsv.gz) | <1 MB |

All three are also accessible via the [UCSC Xena GDC TCGA-LIHC cohort page](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Liver%20Cancer%20(LIHC)).

Somatic mutation data (MAF) is downloaded automatically via `TCGAbiolinks::GDCdownload()` on first run and cached locally.

> **Note on expression data format:** The Xena STAR-FPKM file uses Ensembl gene IDs with version suffixes as row names (e.g. `ENSG00000134240.12`). The script extracts target genes by stable ID prefix using `grep()`. This replaces the gene-symbol-based extraction in v1, which failed on this file format.

---

## Usage

1. Download the three Xena files above and place them in the same directory (e.g. `D:/scRNA_project/`)
2. Open `TCGA_LIHC_ketone_analysis.R` and update the file paths at the top of each section to match your local directory
3. Create an output folder for figures (e.g. `D:/综述代码图/`)
4. Run the script sequentially from top to bottom

---

## Reproducibility note

TCGA-LIHC MAF data is retrieved via `TCGAbiolinks::GDCprepare()`. GDC periodically updates file metadata, which can cause `GDCprepare()` to fail if the local cache is outdated. If you encounter the error `"I couldn't find all the files from the query"`, re-run `GDCdownload()` to refresh the local cache before calling `GDCprepare()`.

Expression, survival, and clinical data are loaded directly from the Xena flat files and are not affected by GDC metadata updates.

---

## Output files

All figures are saved to the output directory specified in the script:

| File | Figure |
|------|--------|
| `ketone_genes_expression_v2.png` | Fig. 5A — Differential expression |
| `survival_OXCT1.png` | Fig. 5B — OXCT1 survival |
| `survival_HMGCS2.png` | Fig. 5C — HMGCS2 survival |
| `survival_combined.png` | Fig. 5D — Combined stratification |
| `OXCT1_ImmuneScore_correlation.png` | Fig. 5E — ESTIMATE immune score |
| `OXCT1_exhaustion_markers.png` | Fig. 5F — T-cell exhaustion |
| `HMGCS2_OXCT1_correlation.png` | Fig. 5G — HMGCS2/OXCT1 correlation |
| `Fig5H_CTNNB1_OXCT1.png` | Fig. 5H — CTNNB1 subgroup |
| `SuppFig_MI_sensitivity.png` | Supp. Fig. S1 — MI sensitivity |

---

## License

MIT
