# run_LDSC.py: Linkage Disequilibrium Score Regression (LDSC) Pipeline

This Python script implements a pipeline for Linkage Disequilibrium Score Regression (LDSC) analysis. LDSC is commonly used to calculate the heritability enrichment and genetic correlations across traits using summary statistics and LD scores. The script performs tasks such as creating annotation files, computing LD scores, and calculating heritability enrichment based on GWAS data.

## Overview

This pipeline automates the process of LDSC analysis by following these steps:
1. **Create Annotation Files**: Using input data (e.g., BED files), it generates chromosome-specific annotation files.
2. **Compute LD Scores**: LDSC computes LD scores for each chromosome using the created annotation files.
3. **Calculate Heritability Enrichment**: The pipeline computes heritability enrichment for various GWAS traits.


## Script Workflow
- Step 1: Create Annotation Files
The script creates chromosome-specific annotation files using a BED file of genomic regions. It generates a shell script that executes this task for each chromosome (1 to 22).
Example output: SamplePrefix.chr1.annot.gz
- Step 2: Compute LD Scores
The LDSC software computes the linkage disequilibrium (LD) scores based on the annotation files from Step 1.
Example output: SamplePrefix.chr1.l2.ldscore.gz
- Step 3: Calculate Heritability Enrichment
The heritability enrichment analysis is performed using the computed LD scores and baseline annotations. The pipeline uses LDSC to calculate the heritability enrichment of each sample.
Example output: SamplePrefix.chr1.Enrichment.results





# run_SCANVAGE.R: SCAVENGE Chromatin Accessibility Analysis Pipeline

This pipeline performs chromatin accessibility analysis using **SCAVENGE**, **chromVAR**, and related tools. The main script `Run_SCAVENGE()` processes scATAC-seq data stored in a `SummarizedExperiment` object, performing trait enrichment analysis and identifying significant cell populations using KNN and random walk methods.

## Overview

This project analyzes chromatin accessibility signals in scATAC-seq data, integrating trait data to identify significant cell populations enriched for certain traits. It leverages weighted deviations and KNN-based random walk methods to identify significant cells. 

## Prerequisites

Before running this pipeline, install the following R packages:
- **SCAVENGE**
- **chromVAR**
- **gchromVAR**
- **BuenColors**
- **SummarizedExperiment**
- **data.table**
- **BiocParallel**
- **BSgenome.Hsapiens.UCSC.hg38**
- **dplyr**
- **igraph**
- **parallel**
- **Signac**
- **Seurat**
- **reshape2**


## Functionality

### Main Function: `Run_SCAVENGE(SE_Data, trait_file)`

#### Input Arguments:
- `SE_Data`: A `SummarizedExperiment` object containing peak-by-cell matrices from scATAC-seq.
- `trait_file`: A file containing trait information formatted as follows:  
  **chr**, **snp_loc**, **snp_loc**, **snp_id**, **posterior probability**, and **filename**.

#### Output:
- Returns a data frame (`mono_mat2`) of cells with trait relevance scores and the significant cells identified after permutation testing.

### Steps Performed:
1. **Data Preprocessing**: Adds GC bias information using UCSC's hg38 genome.
2. **Background Peaks**: Selects background peaks for normalization.
3. **Weighted Deviations**: Calculates trait enrichment scores using chromatin accessibility data.
4. **Random Walk on KNN Graph**: Propagates information across KNN graphs to compute neighborhood propagation scores.
5. **Permutation Testing**: Identifies significant cells by comparing observed and random walk-based scores through permutation testing.

## Usage

To run the SCAVENGE analysis, load your `SummarizedExperiment` object and pass it to the `Run_SCAVENGE()` function along with the trait file:

```r
SE_Data <- readRDS("path_to_SE_data.rds")
trait_file <- "path/to/trait_file.bed"
results <- Run_SCAVENGE(SE_Data, trait_file)
