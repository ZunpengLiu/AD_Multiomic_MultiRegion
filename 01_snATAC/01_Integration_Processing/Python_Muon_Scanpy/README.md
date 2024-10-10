# Generate_TileMatrix.h5ad.py:
## snATAC-seq h5ad File Generation

This repository contains a Python script that processes single-nucleus ATAC-seq (snATAC-seq) data, renames fragment barcodes, counts fragments per peak, filters low-quality cells, and exports the data as an `.h5ad` file for use in **muon** and **scanpy**.

## Description

The script reads fragments from a `.tsv.gz` file, processes them by renaming barcodes, counts the fragments for each peak, and filters out cells based on total fragment counts. The resulting data is saved as an `.h5ad` file, which can be used for further analysis in **muon** or **scanpy**.

### Key Features:
1. **Renaming Fragments**: Barcodes are renamed by appending the sample ID to them.
2. **Counting Fragments**: The script counts the number of fragments overlapping predefined peaks.
3. **Quality Control (QC)**: Cells with total fragment counts outside a given range (1000–100,000) are filtered out.
4. **Output**: The processed data is saved as an `.h5ad` file with compressed storage.

## Dependencies

The following Python libraries are required:

- `pandas`
- `numpy`
- `anndata`
- `muon`
- `scanpy`
- `os`
- `sys`


## Usage
The script is designed to be run from the command line, with options provided to specify input directories, sample names, peak files, and output directories.

Command-line options:
- -i or --indir: Input directory containing the fragment .tsv.gz files.
- -s or --sample: Sample name used for renaming barcodes.
- -p or --peak: Path to the peak file in .bed format (default is ~/ref/hg38_5000bp.bed.gz).
- -o or --outdir: Output directory where the processed files will be saved.

### Example command:
1. For the generation of 5k Tilematrix
python /home/zunpeng/data/src/generate_h5ad_3.py \
-i ~/Fragments/ \
-o ~/outputs/ 
-p ~/hg38.5kb.bed \
-s Sample1

2. For the generation of the Cell-by-Peak count matrix
python /home/zunpeng/data/src/generate_h5ad_3.py \
-i ~/Fragments/ \
-o ~/outputs/ 
-p ~/Celltye_peakSet.tsv \
-s Sample1


# snATAC_Integration_processing.py
## snATAC-seq h5ad File Integration

This Python script for integrating multiple single-nucleus ATAC-seq (snATAC-seq) datasets stored in `.h5ad` format. The script performs batch correction, dimensionality reduction, clustering, and visualization using **muon**, **scanpy**, and optionally **Harmony** or **BBKNN**.

## Description

This script allows for the integration and analysis of snATAC-seq data using various techniques for batch correction, clustering, and visualization. The main steps in the workflow include:

1. **Reading the input `.h5ad` file**.
2. **Preprocessing and filtering cells and peaks** based on user-defined thresholds.
3. **Dimensionality reduction** using TF-IDF normalization and Latent Semantic Indexing (LSI).
4. **Batch correction** using either **Harmony** or **BBKNN** (optional).
5. **Clustering cells** using the Leiden algorithm.
6. **UMAP embedding** for visualization.
7. **Exporting results** including clustering results and UMAP plots.

## Dependencies

The following Python libraries are required to run the script:

- `pandas`
- `numpy`
- `muon`
- `scanpy`
- `harmony-pytorch` (if using Harmony)
- `bbknn` (if using BBKNN)

## Usage
The script is designed to be executed from the command line, with several input options to customize the analysis. The main input is the .h5ad file containing the snATAC-seq data.

## Example command:
python snATAC_Integration_processing \
    -i snATAC.AllRegion.h5ad \
    -o ./ --prefix snATAC.AllRegion.Harmony_LibType \
    --use_harmony \
    -b LibType