# Load required packages
library(nebula)  # For differential expression analysis using Negative Binomial GLMM
library(Seurat)  # For handling single-cell data

# Define input and output directories
indir <- "./indir/"
out_dir <- "./outdir/"

# Read command line arguments
argv <- commandArgs(T)
name <- argv[1]  # Cell type name (e.g., "Exc.DG.granule.cells")
sample.col <- "Sample"  # Column in metadata indicating sample identity

cat("Loading data...\n")

# Load the Seurat object corresponding to the cell type
Seurat_Obj <- readRDS(paste0(indir, "atac2.All.", name, ".rds"))

# Define the differential expression analysis function using nebula
deg.nebula <- function(Seurat_Obj, pathology, sample.col, 
                       offset = "total_counts", 
                       ncore = 12,
                       cpc = 0, reml = 1) {

    # Define covariates to include in the model
    covariates <- c("age", "msex", "pmi", "LibType", "total_counts")
    covariates <- covariates[covariates %in% colnames(Seurat_Obj@meta.data)]

    # Prepare input for nebula
    seuratdata <- scToNeb(obj = Seurat_Obj, 
                          assay = "RNA", 
                          id = sample.col, 
                          pred = c(covariates, pathology), 
                          offset = offset)

    # Construct model matrix
    design <- model.matrix(
        as.formula(paste0("~", paste(c(pathology, covariates), collapse = "+"))),
        data = Seurat_Obj@meta.data
    )

    # Run NEBULA differential expression analysis
    neb <- nebula(seuratdata$count,
                  seuratdata$id,
                  pred = design,
                  model = "NBGMM",
                  ncore = ncore,
                  offset = seuratdata$offset)

    # Parse NEBULA results
    neb$summary$FDR <- p.adjust(neb$summary[[paste0("p_", pathology)]], method = "fdr")
    neb$summary$log2FC <- neb$summary[[paste0("logFC_", pathology)]]

    neb_df <- neb$summary
    ovr <- neb$overdispersion
    colnames(ovr) <- paste0("overdispersion_", colnames(ovr))

    # Combine results
    neb_df <- cbind(neb_df, ovr)
    neb_df$convergence <- neb$convergence
    neb_df$algorithm <- neb$algorithm
    rownames(neb_df) <- neb_df$gene

    # Sort by p-value
    neb_df <- neb_df[order(neb_df[[paste0("p_", pathology)]]), ]

    # Write full result table
    write.table(
        neb_df,
        paste0(out_dir, name, ".", pathology, ".All.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )

    # Write simplified result table
    results <- neb_df[, c("gene", paste0("p_", pathology),
                          paste0("logFC_", pathology), "FDR", "log2FC")]
    write.table(
        results,
        paste0(out_dir, name, ".", pathology, ".Clean.tsv"),
        sep = "\t", row.names = FALSE, quote = FALSE
    )
}

# Run DEG analysis for the specified pathology
deg.nebula(Seurat_Obj, pathology, sample.col)
