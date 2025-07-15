# Scanpy-based Integration Pipeline for snRNA-seq

This script performs integration, clustering, and visualization of single-nucleus RNA-seq (snRNA-seq) data using [Scanpy](https://scanpy.readthedocs.io), with optional batch correction via Harmony or BBKNN.

## Features

- Normalization and log1p transformation
- Highly variable gene (HVG) selection
- PCA and UMAP embedding
- Batch correction (optional): Harmony or BBKNN
- Leiden clustering
- QC and UMAP plots for metadata and quality metrics
- Saves final annotated `.h5ad` file

## Usage

```bash
python integrate_snRNAseq.py \
    -i input_file.h5ad \
    -o ./output/ \
    --prefix mydata \
    --batch LibraryType \
    --use_harmony \
    --leiden leiden_cluster \
    --resolution 1.0
```
