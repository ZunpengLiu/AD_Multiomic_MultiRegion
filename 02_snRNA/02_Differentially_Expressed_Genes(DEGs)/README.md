# NEBULA-based Differential Expression Analysis for Single-Cell RNA-seq

This project performs differential gene expression (DGE) analysis on single-cell RNA-seq data using the [NEBULA](https://cran.r-project.org/web/packages/nebula/index.html) framework, tailored for hierarchical modeling of cells nested within individuals or samples.

## Input

- A `Seurat` object in `.rds` format named:`atac2.All.<celltype>.rds`(e.g., `atac2.All.Exc.DG.granule.cells.rds`)
- The object must contain:

  - RNA assay
  - Metadata columns including:
    - Sample ID (`Sample`)
    - Covariates such as `age`, `msex`, `pmi`, `LibType`, `total_counts`
    - The pathology variable of interest (e.g., `Epigenetic_Information`)

## How to Run

```bash
Rscript run_nebula.R Exc.DG.granule.cells
```

## Dependencies

Before running the script, make sure to install the following R packages:

- **Nebula**: Used for DAR analysis.
- **Seurat**: Handles single-cell ATAC-seq data and converts Seurat objects into nebula-compatible formats.
