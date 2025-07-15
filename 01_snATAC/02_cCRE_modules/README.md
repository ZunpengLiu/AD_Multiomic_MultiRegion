# Non-Negative Matrix Factorization (NMF) for Chromatin Accessibility Data Classification

## Overview

This Python script is used for performing **Non-Negative Matrix Factorization (NMF)** to classify over one million candidate cis-regulatory elements (cCREs) across 67 cell subtypes into distinct modules based on chromatin accessibility profiles. The analysis is performed using ATAC-seq data from 6 brain regions, enabling the identification of cCREs associated with specific cell subtypes and brain regions.


### Methodology

We employed NMF to classify cCREs across multiple cell subtypes based on their chromatin accessibility profiles. Given the sparse nature of single-nucleus ATAC-seq (snATAC-seq) data, chromatin accessibility signals were first aggregated across cCREs and cell subtypes from 6 brain regions, resulting in a matrix **V** with dimensions **N × M**:
- **N**: Represents the cCREs (~1 million).
- **M**: Represents the combined dimensions of the 67 cell subtypes and 6 brain regions.

Using the `sklearn.decomposition` package in Python, we factorized the non-negative matrix **V** into two smaller non-negative matrices:
- **W (N × R)**: The basis matrix, where each column represents a module, and each row represents the participation of each cCRE in the module.
- **H (R × M)**: The coefficient matrix, where each row corresponds to a module, and each column represents a cell subtype and brain region combination.

The factorization follows the formula:
V ≈ W × H


### Interpretation of the Matrices:
1. **Basis Matrix (W)**:
   - **W (N × R)** characterizes the module-related accessible cCREs. Each row represents a cCRE, and each column (R) represents a module, defining the extent to which a cCRE is associated with each module.
   - A "feature score" for each cCRE was calculated based on the W matrix, reflecting its specificity to a module, using the "kim" method from prior research. The feature score ranges from 0 to 1, where a higher score indicates greater specificity to a module.
   - **Module Assignment**: 
     - A cCRE was assigned to a module if its feature score was above the median plus one standard deviation. 
     - cCREs with feature scores below the median minus one standard deviation were aggregated into a non-specific module (M1). 
     - cCREs with intermediate feature scores were considered shared across multiple modules.

2. **Coefficient Matrix (H)**:
   - **H (R × M)** describes the contribution of each cell subtype and brain region to the identified modules. Each row corresponds to a module, and each column corresponds to a specific cell subtype and brain region.
   - The highest coefficient value for each cell subtype in matrix **H** was used to assign the primary association of each cell subtype to a particular module.


## Features

- **Matrix Factorization**: Factorizes the chromatin accessibility matrix into basis and coefficient matrices using NMF.
- **Feature Score Calculation**: Extracts feature scores for cCREs based on their specificity to modules.
- **Module Prediction**: Predicts and associates modules with cCREs and cell subtypes across brain regions based on the results of matrix factorization.
- **High Dimensionality Data**: Handles chromatin accessibility profiles from over one million cCREs and 67 cell subtypes across 6 brain regions.

## Dependencies

To run the script, you need the following Python packages:

- **numpy**: For numerical computations.
- **pandas**: For data manipulation and analysis.
- **sklearn.decomposition (NMF)**: To perform Non-Negative Matrix Factorization.
- **scipy**: For additional mathematical functions.
- **scikit-learn**: For matrix scaling and normalization.

### References:
- **Li, Yang Eric et al.** 2023. *A Comparative Atlas of Single-Cell Chromatin Accessibility in the Human Brain*. Science.
- **Zu, Songpeng et al.** 2023. *Single-Cell Analysis of Chromatin Accessibility in the Adult Mouse Brain*. Nature.
