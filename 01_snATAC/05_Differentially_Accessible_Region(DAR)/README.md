# Differentially Accessible Region (DAR) Analysis Using Nebula based on the whole-genome 100-kb count reads

This repository contains an R script for performing **Differentially Accessible Region (DAR)** analysis using the [Nebula](https://github.com/lhe17/nebula) package on cell-type specific peaks from single-cell ATAC-seq data. The analysis aims to identify regions of differential chromatin accessibility associated with specific pathologies (e.g., amyloid, Alzheimer's Disease Pathology), while controlling for key covariates such as sex, post-mortem interval (PMI), and library type (i.e. snATAC-seq data or snATAC data from 10x multiome-seq).

## Features

- **DAR Identification**: Detects differentially accessible regions (DARs) in relation to specific pathologies across different cell types.
- **Covariate Adjustment**: Adjusts for confounding factors such as age, sex, PMI, library type, and total counts.
- **Output**: Provides p-values, log fold changes, false discovery rate (FDR) adjusted p-values, overdispersion estimates, and convergence information.

## Dependencies

Before running the script, make sure to install the following R packages:

- **Nebula**: Used for DAR analysis.
- **Seurat**: Handles single-cell ATAC-seq data and converts Seurat objects into nebula-compatible formats.
