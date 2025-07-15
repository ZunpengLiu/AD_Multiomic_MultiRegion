# 02_ArchR_Proj_Reads_in_Chromatin_States.R: 
## Chromatin State Analysis using ArchR

This repository contains a script for analyzing chromatin states using the ArchR package. The script processes chromatin states from BED files, adds feature counts to an ArchR project, and writes out the results.

## Script Overview

- **Step 1: Load the ArchR Project**  
  The script begins by loading an existing ArchR project located in the directory `./03_ArchR/01_Overall`.

- **Step 2: Read and Process Chromatin States**  
  The script reads in chromatin states for a specific tissue and processes them for a predefined list of states, including `EnhA1`, `EnhA2`, `TssA`, etc. These BED files are filtered to include only certain chromosomes (chr1 to chr22, and chrX).

- **Step 3: Add Feature Counts**  
  The `addFeatureCounts` function from ArchR is used to add the processed chromatin states to the ArchR project and calculate the fraction of reads in each state.

- **Step 4: Write Output**  
  The final processed data for the chromatin state feature counts is written to the `./05_InformationScore/` directory.

## Output

The output of the script is a table of chromatin state feature counts per cell, which is saved in the `./05_InformationScore/` directory as `<tissue>.Count_Fractions.txt`.




# 03_Chromatin_state_erosion_score.py: 
## Chromatin State Erosion Score Analysis

This repository contains a Python script to compute the erosion scores for different chromatin states across various brain regions using Z-score normalization. The script reads in chromatin state ratio matrices, normalizes the data, and calculates erosion scores using a predefined sign matrix for each chromatin state.

## Script Overview
- **1.Load the Data**
The script reads in data from .txt files that contain chromatin state ratio matrices. It filters columns based on their ratio and the specific chromatin states being analyzed.


- **2. Normalize the Data**
The script normalizes the ratio matrix using Z-score normalization, which scales the data based on the mean and standard deviation of each chromatin state.


- **3. Calculate the Z-score epigenomic information Score**
Using the Z-score normalized matrix and a predefined sign_matrix, the epigenomic information score is calculated by summing the weighted Z-scores for each chromatin state.

- **4. Process Chromatin States in Brain Regions**
The script processes chromatin state data for several brain regions: frontal cortex, hippocampus, middle frontal area, and angular gyrus. It calculates erosion scores for each of these regions and stores the final result.

## Output
The final output is a tab-separated file Score.txt, which contains the Z-score epigenomic information scores for each chromatin state across all processed brain regions. The erosion score for each brain region is calculated and stored under the respective Z-score erosion score column. The overall erosion score is computed as the mean of these scores.

