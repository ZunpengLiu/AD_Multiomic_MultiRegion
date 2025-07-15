# CellType_Compositional_Changes_with_Pathology.R: 
## Cell Proportion Changes in snRNA and snATAC Across Classes, Subclasses, and Brain Regions

This script is for analyzing cell proportion changes in different cell classes, subclasses, and subtypes across multiple brain regions in snRNA and snATAC data. The script processes cell fraction data, performs statistical analysis using odds ratios, and writes the results to output files.

Script Overview
The script performs the following steps:

- **Step 1: Load the Dataset**
The script begins by loading the necessary metadata and assigning labels for brain regions and cell classes.

- **Step 2: Define Subclasses and Colors**
Subclasses and their respective colors for visualization are defined for plotting purposes.

- **Step 3: Calculate Cell Fraction Changes in classes, subclasses, and subtypes**
The script calculates cell fraction changes for each class, subclass, and subtype, and performs odds ratio analysis on cell fractions using a generalized linear model (glm).