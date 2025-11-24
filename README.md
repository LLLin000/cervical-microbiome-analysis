# cervical-microbiome-analysis
R scripts for ANCOM, LEfSe, and Random Forest analyses used in the study of vaginal and perianal microbiome shifts across cervical neoplasia progression and HPV infection.

## Included Scripts
- StdANCOM.R — ANCOM-BC differential abundance analysis
- StdLEfSe.R — LEfSe biomarker detection
- StdRF.R — Random Forest classification, variable importance, and ROC analysis

## Input Data
Each script reads OTU tables, taxonomy files, and metadata (CSV format).
Place the required CSV files in the same directory before running.

## Usage
Run each script directly in R:

```
source("StdANCOM.R")
source("StdLEfSe.R")
source("StdRF.R")
```
- Outputs include differential taxa tables, heatmaps, LEfSe plots, Random Forest performance, and top microbial predictors.

## Purpose
These scripts were used to support statistical analyses and result generation for the associated manuscript.
