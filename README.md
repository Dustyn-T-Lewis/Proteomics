# Aging Proteomics Analysis

## Overview
This repository contains the R Markdown analysis pipeline for processing, analyzing, and visualizing aging-related proteomics data from middle-aged (MA) and young human muscle samples. Outputs include cleaned data tables, statistical summaries, pathway enrichment results, and publication-ready figures in PDF, Word, HTML, and GitHub-renderable Markdown formats.
---

## File Locations
- **Master Data** "13_MASTER_SEER proteome analysis (4-4-23)""
    - Contains all proteomic information, this code extracts Sarcoplasmic Protein Fraction
    
---

## Software/Running Script
- **R** ≥ 4.1  
- **Pandoc** (installed with RStudio or standalone)  
- **R packages** (run before first use):
  ```r
  install.packages(c(
    "readxl","dplyr","stringr","tidyr","httr","jsonlite",
    "purrr","tibble","broom","ggplot2","forcats","fmsb","BiocManager"
  ))
  BiocManager::install(c(
    "ReactomePA","clusterProfiler","org.Hs.eg.db","enrichplot"
  ))
  ```
# Notes
- All figures are rendered at 8×6″, 300 dpi, and scaled to 100% width via the setup chunk in Aging-Proteomics.Rmd.
├── README.md
├── data/
│   └── 13_MASTER_SEER proteome analysis (4-4-23).xlsx
├── scripts/
│   └── Aging-Proteomics.Rmd
├── output/
│   ├── Aging-Proteomics.html
│   ├── Aging-Proteomics.md
│   ├── Clean_Long2.csv
│   ├── Clean_Wide2.csv
│   └── Group_Abundance.csv
└── figures/
    └── zoom_in_plots/
