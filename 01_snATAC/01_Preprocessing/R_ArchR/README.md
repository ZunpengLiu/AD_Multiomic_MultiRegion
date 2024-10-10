# 01_ArchR_Generate_Arrows:
# ArchR Arrow Files Generation

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




# 02_ArchR_Generate_Proj_QC_1st_UMAP:
# ArchR Project Creation and Dimensionality Reduction

This script is for generating an ArchR project from Arrow files, performing quality control, doublet filtering, and dimensionality reduction using IterativeLSI and UMAP. Batch effect correction is also applied using Harmony. The resulting project is ready for downstream analysis in single-cell ATAC-seq workflows.

## Description

This script processes Arrow files and applies the following steps:

1. **Reading Arrow Files**: Loads pre-generated Arrow files.
2. **Creating an ArchR Project**: Initializes an ArchR project from the Arrow files.
3. **Doublet Filtering**: Filters out potential doublets.
4. **Dimensionality Reduction**: Applies IterativeLSI, clustering, and UMAP for visualization.
5. **Batch Effect Correction**: Uses Harmony to correct batch effects.
6. **Clustering**: Performs clustering on the reduced dimensionality data.



# 03_ArchR_High-Quality_Cell_Subsetting_and_Dimensionality_Reduction: 
# ArchR High-Quality Cell Processing and Dimensionality Reduction

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



# 04_ArchR_GeneScore_Marker_Genes:
# ArchR Gene Score Analysis and Marker Gene Identification

This script processes an ArchR project to calculate gene scores, export gene score matrices for each cell subtype, and identify marker genes for major cell types. The script also generates UMAP plots for canonical marker genes and performs a marker gene analysis using the GeneScoreMatrix.

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
- **Marker Genes**: CSV files containing marker genes for each cell subtype, saved as `proj3.Cell_Subtype.Marker_genes.Cell_Type.csv`.

## Example Canonical Marker Genes

The script includes the following major cell type marker genes for UMAP plotting:

- **Excitatory neurons**: `SLC17A7`
- **Inhibitory neurons**: `GAD2`
- **Astrocytes**: `GFAP`
- **Oligodendrocytes**: `MOBP`
- **OPCs**: `VCAN`
- **Microglia**: `P2RY12`



# 05_ArchR_BigWig_Generation:
# ArchR Group BigWig File Generation

This script generates group-specific BigWig files from an ArchR project. These BigWig files represent snATAC-seq signal tracks normalized by total fragment counts.

## Description

The script performs the following steps:

1. **Loading an ArchR Project**: Loads a pre-existing ArchR project.
2. **Generating Group BigWig Files**: Produces BigWig files for several cellular cell type and pathology groups:
   - `Cell_Class`
   - `Cell_Subclass`
   - `Cell_Subtype`
   - `Class_Pathology`
   - `Subclass_Pathology`

## Output

The output will consist of BigWig files for each specified group, saved to the output directory. These files can be used for downstream analysis or visualized using tools like IGV (Integrative Genomics Viewer).



# 06_ArchR_Peak_Calling： 
# ArchR Peak Calling and Export

This script performs peak calling for snATAC-seq data using the ArchR package. Peaks are grouped by cell subtype and exported in both `.tsv` and `.bed` formats for further analysis or downstream processing.

## Description

The script performs the following steps:

1. **Loading an ArchR Project**: Loads a pre-existing ArchR project containing high-quality cells.
2. **Peak Calling**: 
   - Grouping cells by `Cell_Subtype` and calling peaks using the MACS2 tool.
   - Generating reproducible peak sets for each cell subtype.
   - Exporting the union peak set as a `.tsv` file.
3. **Peak Matrix Generation**: Adds a peak matrix to the ArchR project after peak calling.
4. **Exporting Peaks**: Exports the reproducible peak sets in both `.tsv` and `.bed` formats for each cell subtype. Also exports a combined file with all reproducible peaks.
5. **Sorting `.bed` Files**: A system command is used to sort the `.bed` files by chromosome and start position for downstream compatibility.


## Output

The script generates the following outputs:
1. **Reproducible Peak Sets**: 
   - `.tsv` files for each cell subtype, containing the peak data.
   - `.bed` files for each cell subtype, containing the peak regions in BED format.
   - A combined `.tsv` file with all reproducible peaks across all cell subtypes.
   
2. **Sorted BED Files**: The `.bed` files are sorted by chromosome and start position for compatibility with downstream tools like IGV or BEDTools.

The final output will include files such as:
- `Cell_Subtype.tsv`
- `Cell_Subtype_macs2_peak.bed`
- `All_reproduciblePeaks.tsv`



# 07_ArchR_Motif_Deviation:
# ArchR Motif Deviations Analysis

This script processes an ArchR project to compute per-cell deviations for motif annotations, generate summarized motif z-scores across cell subtypes, and identify motifs with high variability across clusters. The analysis is performed using the `addMotifAnnotations()` and `addDeviationsMatrix()` functions from the ArchR package.

## Description

The script performs the following steps:

1. **Loading an ArchR Project**: Loads an existing ArchR project containing high-quality cells.
2. **Adding Motif Annotations**: Adds motif annotations from the **cisBP** motif database to the ArchR project.
3. **Adding Background Peaks**: Adds background peaks to the project for downstream analysis.
4. **Computing Deviations Matrix**: Computes per-cell deviations across all motif annotations.
5. **Identifying Variable Motifs**: Identifies motifs that show high variability in their z-scores across different cell subtypes.
6. **Exporting ChromVar Z-scores**: Generates a summarized z-score matrix for motifs across cell subtypes and exports it as a `.tsv` file.

## Output

The script produces the following key outputs:

1. **ChromVar Z-Score Matrix**:
   - The matrix is saved as `ChromVar_Zscore_Mat_by_Cell_subtype.tsv`.
   - Contains motif z-scores for each motif across cell subtypes.

2. **ArchR Project**:
   - The ArchR project is saved with all the calculated matrices and annotations for future use.