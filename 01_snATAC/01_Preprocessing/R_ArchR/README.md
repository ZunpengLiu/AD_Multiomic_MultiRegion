# 01_ArchR_Generate_Arrows: ArchR Arrow Files Generation

This script generates Arrow files from fragment data using the ArchR(https://www.archrproject.com/bookdown/index.html) package. Arrow files are required for further analysis in single-cell ATAC-seq workflows in ArchR. The script also performs quality control and adds doublet scores to the Arrow files.

## Description

This script processes ATAC-seq fragment files to create Arrow files for downstream analysis with ArchR. The steps include:
1. Loading fragment data.
2. Generating Arrow files with quality control metrics such as transcription start site (TSS) enrichment and gene score matrices.
3. Adding doublet scores to identify potential doublets in the single-cell data.

## Dependencies

The following R packages are required to run the script:

- `GenomicRanges`
- `ArchR`
- `stringr`
- `dplyr`
- `BSgenome.Hsapiens.UCSC.hg38`



# 02_ArchR_Generate_Proj_QC_1st_UMAP: ArchR Project Creation and Dimensionality Reduction

This script is for generating an ArchR project from Arrow files, performing quality control, doublet filtering, and dimensionality reduction using IterativeLSI and UMAP. Batch effect correction is also applied using Harmony. The resulting project is ready for downstream analysis in single-cell ATAC-seq workflows.

## Description

This script processes Arrow files and applies the following steps:

1. **Reading Arrow Files**: Loads pre-generated Arrow files.
2. **Creating an ArchR Project**: Initializes an ArchR project from the Arrow files.
3. **Doublet Filtering**: Filters out potential doublets.
4. **Dimensionality Reduction**: Applies IterativeLSI, clustering, and UMAP for visualization.
5. **Batch Effect Correction**: Uses Harmony to correct batch effects.
6. **Clustering**: Performs clustering on the reduced dimensionality data.


# 03_ArchR_High-Quality_Cell_Subsetting_and_Dimensionality_Reduction: ArchR High-Quality Cell Processing and Dimensionality Reduction

This script processes an ArchR project by loading high-quality cells, performing dimensionality reduction using IterativeLSI and UMAP, and integrating UMAP embeddings from external sources (muon and Scanpy). The script also generates quality control plots and compares ArchR-generated UMAPs with external data.

## Description

The script performs the following steps:

1. **Loading an ArchR Project**: Loads an existing ArchR project for processing.
2. **Subsetting High-Quality Cells**: Filters high-quality cells based on muon and Scanpy results.
3. **Quality Control (QC) Plots**: Generates fragment size distribution and TSS enrichment plots for assessing data quality.
4. **Dimensionality Reduction**: Redoes dimensionality reduction with IterativeLSI, followed by clustering and UMAP embedding.
5. **UMAP from External Data**: Loads UMAP coordinates from external processing (muon and Scanpy) and integrates them into the ArchR project.
6. **UMAP Plotting**: Plots UMAP visualizations using cell subtype information.

## Output

- Subset ArchR project files containing only high-quality cells.
- Dimensionality reduction results, including UMAP embeddings from the ArchR.
- PDF files with QC plots (fragment size, TSS enrichment).
- UMAP plots comparing ArchR-generated and external UMAP data.


# ArchR Gene Score Analysis and Marker Gene Identification

This repository contains a script that processes an ArchR project to calculate gene scores, export gene score matrices for each cell subtype, and identify marker genes for major cell types. The script also generates UMAP plots for canonical marker genes and performs a marker gene analysis using the GeneScoreMatrix.

## Description

The script performs the following steps:

1. **Loading an ArchR Project**: Loads an existing ArchR project.
2. **Calculating Gene Scores**: Adds imputation weights and calculates gene scores using the GeneScoreMatrix.
3. **Exporting Gene Score Matrix**: Extracts the gene score matrix for each cell subtype and exports it to a TSV file.
4. **UMAP Plotting**: Plots UMAP embeddings for canonical marker genes using the GeneScoreMatrix.
5. **Marker Gene Identification**: Identifies marker genes for each cell subtype using the Wilcoxon test.

## Output

- **Gene Score Matrix**: A matrix of gene scores by cell subtype, saved as `GeneMat_by_Cell_Subtype.tsv`.
- **UMAP Plots**: UMAP visualizations for canonical marker genes, saved as `GeneScore.Key_Marker_genes_solarExtra.pdf`.
- **Marker Genes**: CSV files containing marker genes for each cell subtype, saved as `proj3.Cell_Subtype.Marker_genes.<Cell_Type>.csv`.

## Example Canonical Marker Genes

The script includes the following major cell type marker genes for UMAP plotting:

- **Excitatory neurons**: `SLC17A7`
- **Inhibitory neurons**: `GAD2`
- **Astrocytes**: `GFAP`
- **Oligodendrocytes**: `MOBP`
- **OPCs**: `VCAN`
- **Microglia**: `P2RY12`
