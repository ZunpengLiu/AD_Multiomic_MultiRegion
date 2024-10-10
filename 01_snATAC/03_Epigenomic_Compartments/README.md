# 01_HMMt_Epigenomic_Compartments.R:
## Epigenomic Compartment Analysis using Hidden Markov Models (HMM)

This script performs epigenomic compartment analysis on 100kb genomic bins, identifying active and repressive compartments across various brain regions using Hidden Markov Models (HMM). The analysis is performed on high-quality ATAC-seq data and annotated across different major cell types and different AD pathologies, i.e. nonAD, earlyAD, and lateAD.

## Overview

This pipeline uses a Hidden Markov Model (HMM) to classify genomic regions as either "Active" or "Repressive" based on scaled signal intensities.

## Features

- **HMM-based Classification**: The HMM identifies "Active" and "Repressive" genomic compartments.
- **Genomic Ranges Manipulation**: Uses `GenomicRanges` to efficiently handle genomic regions and overlaps.
- **Cell Type-Specific Analysis**: The script processes multiple cell types and generates compartment annotations for each.
- **High-quality 100kb Bins**: Uses high-resolution 100kb genomic bins for analysis.

## Dependencies

The following R packages are required to run the script:

- `devtools`
- `data.table`
- `GenomicRanges`
- `HMMt`
- `stringr`


# 02_Compartment_100kb_Embedding:
## 100-kb whole-genomic embeddings to visualize the Chromosomes, ATAC-seq Compartment, and the epigenomic features

This analysis is to process the transposed ATAC-seq data, perform dimensionality reduction, and visualize genomic compartments and chromosome correlations using UMAP and correlation matrices. The analysis involves the integration of epigenomic annotations, transposing matrices, and performing in silico chromosome correlation analysis.

## Description

The script performs the following tasks:

1. **Data Loading and Preprocessing**:
   - Loads ATAC-seq data stored in `.h5ad` format.
   - Transposes the matrix from cell x 100kb format to 100kb x cell format for further processing.

2. **UMAP Visualization**:
   - Uses UMAP to visualize epigenomic signatures such as LINE elements, SINE elements, and other repeat types.
   - Visualizes compartment groups using custom colors.
   
3. **Cell Class-Specific Matrix Generation**:
   - Generates transposed matrices for different cell types, storing results for each cell class separately.

4. **Chromosome Correlation Analysis**:
   - Computes in silico correlations between different chromosomes based on the ATAC-seq data from all cells.
   - Outputs a chromosome correlation matrix for further exploration.



# Fold_Enrichment_of_Features_across_Genomic_Bins.py
## Chromatin State Fold Enrichment Analysis

This Python script calculates fold enrichment for chromatin states across various genomic regions. The script processes multiple segment files, computes fold enrichment for specific chromatin states, and combines the results across different brain regions and chromatin marks such as Lamin B1 and SON.

## Description

The script performs the following tasks:

1. **Processing Segment Files**: For each input file ending with `_18_CALLS_segments.bed`, the script:
   - Intersects the regions with a high-quality hg38 100kb resolution BED file.
   - Calculates fold enrichment and log2 fold enrichment for various chromatin states.
   - Outputs the fold enrichment data in `.tsv` files for each segment file.

2. **Combining Data Across Samples**: The script combines the fold enrichment data across multiple brain regions (`BSS01125`, `BSS01126`, etc.) and computes the mean fold enrichment for specific chromatin states such as `EnhA1`, `TssA`, `ReprPC`, and more.

3. **Lamin B1 and SON Analysis**: The script processes additional files related to Lamin B1 and SON chromatin marks:
   - Computes log2 fold enrichment for these chromatin marks.
   - Saves the combined results in `.tsv` files.

4. **Final Dataset Combination**: The script merges the results of repeat elements, ChromHMM states, and Lamin B1/SON data into a final combined dataset, saved as `hg38_100kb.Repeats.ChromHMM_LaminB_SON.log2_fold_enrichment_all.tsv`.

## Dependencies

The script requires the following Python packages:
- `pandas`
- `numpy`
- `os`
- `sys`

Additionally, the script uses the `bedtools` command-line tool for region intersection. Ensure `bedtools` is installed and available in your system's `PATH`.