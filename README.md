# Aging Proteomics Analysis

## Overview
This repository contains the R Markdown analysis pipeline for processing, analyzing, and visualizing aging-related proteomics data from middle-aged (MA) and young human muscle samples. The outputs include:
- Cleaned data tables
- Statistical summaries
- Pathway enrichment results
- Publication-ready figures in multiple formats:
  - PDF
  - Word
  - HTML
  - GitHub-renderable Markdown

---

## File Locations
- **Master Data**: `"13_MASTER_SEER proteome analysis (4-4-23).xlsx"`
  - Contains all proteomic information.
  - This code extracts Sarcoplasmic Protein Fraction.

---

## Repository Structure

```
├── README.md                # Main documentation file
├── data/                    # Raw data files
│   ├── master_seer_proteome.xlsx
├── scripts/                 # R scripts for analysis
│   ├── Aging-Proteomics.Rmd
│   ├── Clean_Long_t-test.R
│   ├── Clean_Long_t-test.Rmd
├── output/                  # Analysis outputs
│   ├── Clean_Long2.csv          # Cleaned long-format data
│   ├── Clean_Wide2.csv          # Cleaned wide-format data
│   ├── Group_Abundance.csv      # Grouped abundance data
│   ├── GO_Pre_vs_Young_Results.csv    # GO pathway enrichment results
│   ├── KEGG_Pre_vs_Young_Results.csv  # KEGG pathway enrichment results
│   ├── Reactome_Pre_vs_Young_Results.csv  # Reactome pathway enrichment results
│   ├── Aging-Proteomics.html     # HTML report
│   ├── Aging-Proteomics.md       # Markdown report
├── figures/                 # Figures and plots
│   └── zoom_in_plots/
```

---

## Software/Running Script

### Requirements
- **R** ≥ 4.1  
- **Pandoc** (installed with RStudio or standalone)  

### R Packages (Install before first use)
Run the following code in R to install the required packages:

```r
install.packages(c(
  "readxl", "dplyr", "stringr", "tidyr", "httr", "jsonlite",
  "purrr", "tibble", "broom", "ggplot2", "forcats", "fmsb", "BiocManager"
))

BiocManager::install(c(
  "ReactomePA", "clusterProfiler", "org.Hs.eg.db", "enrichplot"
))
```

---

## Notes
- All figures are rendered at 8×6″, 300 dpi, and scaled to 100% width via the setup chunk in `Aging-Proteomics.Rmd`.