# ChromHMM_enrichment.sh:
## ChromHMM OverlapEnrichment for Brain Tissues

This code contains a shell script to run ChromHMM's `OverlapEnrichment` tool on a set of brain tissue samples using Epimap data. The script loops through a list of predefined brain tissue samples and performs enrichment analysis to identify chromatin state overlaps in each bed file.

## Overview

This script automates the process of running ChromHMM's `OverlapEnrichment` command across multiple brain tissue samples. It handles the following tasks:

1. **Loops through each brain tissue sample**: Processes a list of predefined brain samples, each identified by a unique sample ID and tissue name.
2. **Runs ChromHMM's `OverlapEnrichment`**: For each sample, the script runs the `OverlapEnrichment` command using the provided segments files.
3. **Outputs results for each sample**: Results are saved in an output directory, with file names formatted as `sampleID_BRAIN_tissueName`.

## Prerequisites

Before running the script, ensure the following tools and data are available:

### Required Software
- **ChromHMM**: This script relies on `ChromHMM.jar` for running the `OverlapEnrichment` tool. You can download ChromHMM from its official repository: https://compbio.mit.edu/ChromHMM.

### Brain Tissue Samples

The script processes the following brain tissue samples from the Epimap dataset:

- **BSS00369**: Frontal Cortex (BRN.FTL.CTX)
- **BSS00371**: Frontal Cortex (BRN.FTL.CTX)
- **BSS01125**: Hippocampus (BRN.HPC)
- **BSS01126**: Hippocampus (BRN.HPC)
- **BSS01124**: Hippocampus (BRN.HPC)
- **BSS01271**: Middle Frontal Area (BRN.MID.FTL)
- **BSS01272**: Middle Frontal Area (BRN.MID.FTL)
- **BSS01273**: Middle Frontal Gyrus (BRN.MID.FTL.GYR)
- **BSS00077**: Angular Gyrus (BRN.AG)
- **BSS00078**: Angular Gyrus (BRN.AG)
- **BSS00089**: Astrocyte (BRN.AST)
- **BSS00090**: Astrocyte Cerebellum (BRN.AST.CRBLLM)
- **BSS00091**: Astrocyte Hippocampus (BRN.AST.HPC)


# GREAT_Functional_enrichment.R:
## rGREAT Annotation of Cell-Type Specific Peaks

This R script performs annotation of peaks from different cell subtypes or the peaks from the cCRE modules using the **rGREAT** package, with all peaks used as background. The script annotates peaks with Gene Ontology (GO) terms for molecular function (MF), biological processes (BP), and cellular components (CC), and saves the results in separate files.

## Overview

The script performs the following tasks:
1. **Reads the peak files** in BED format.
2. **Overlaps** the peaks with all the peaks.
3. **Runs rGREAT analysis** to annotate the peaks.
4. **Saves** the GO term enrichment results (MF, BP, CC) for each cell type in separate `.tsv` files.

## Requirements

The following R packages are required to run the script:

- `ggplot2`
- `rGREAT`
- `bedr`
- `stringr`