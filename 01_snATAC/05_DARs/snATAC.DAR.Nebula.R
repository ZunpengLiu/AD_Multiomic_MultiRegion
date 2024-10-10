# Load necessary libraries
library(nebula)
library(Seurat)

# Set the input directory and output directory
indir <- "./05_Nebula_DAR/Subtype/"
out_dir <- "./Subtype_output/"

# Read arguments from the command line, assuming the first argument is the cell type name
argv <- commandArgs(T)
name <- argv[1]  # Example: "Exc.DG.granule.cells"

# Column name for sample information in Seurat metadata
sample.col <- "Sample"

# Print loading message
print("Loading data")

# Load the Seurat object for the specific cell type
Seurat_Obj <- readRDS(paste0(indir, "atac2.All.", name, ".rds"))

# Function to perform differential expression analysis using the nebula package
deg.nebula <- function(Seurat_Obj, pathology, sample.col, 
                       offset = "total_counts", 
                       ncore = 12,
                       cpc = 0, reml = 1) {
    
    # Define the covariates to be included in the model
    covariates <- c("msex", "pmi", "LibType", "total_counts")
    
    # Filter out covariates not present in the Seurat object metadata
    covariates <- covariates[covariates %in% colnames(Seurat_Obj@meta.data)]
    
    # Convert Seurat object to nebula-compatible format
    seuratdata <- scToNeb(obj = Seurat_Obj, 
                          assay = "RNA",         # Specify the assay (RNA in this case)
                          id = sample.col,       # ID column (usually the sample or batch column)
                          pred = c(covariates, pathology),  # Include covariates and the pathology of interest
                          offset = offset)       # Offset column (usually total counts)

    # Create a design matrix for the regression model
    design <- model.matrix(as.formula(paste0("~", paste0(c(pathology, covariates), collapse = "+"))), 
                           data = Seurat_Obj@meta.data)

    # Run the nebula model (Negative Binomial Generalized Linear Model with Mixed Effects)
    neb <- nebula(seuratdata$count,
                  seuratdata$id,
                  pred = design,
                  model = "NBGMM",      # Negative Binomial Gaussian Mixture Model
                  ncore = ncore,        # Number of cores for parallel processing
                  offset = seuratdata$offset)  # Offset for the model

    # Adjust p-values for multiple testing (False Discovery Rate)
    neb$summary$FDR <- p.adjust(neb$summary[[paste0("p_", pathology)]], "fdr")

    # Calculate log2 fold change
    neb$summary$log2FC <- neb$summary[[paste0("logFC_", pathology)]]

    # Extract summary results from nebula analysis
    neb_df <- neb$summary

    # Add overdispersion results (overdispersion estimates for each gene)
    ovr <- neb$overdispersion
    colnames(ovr) <- paste0("overdispersion_", colnames(ovr))
    neb_df <- cbind(neb_df, ovr)  # Merge overdispersion with the main results

    # Add convergence and algorithm information
    neb_df$convergence <- neb$convergence
    neb_df$algorithm <- neb$algorithm

    # Set gene names as row names
    rownames(neb_df) <- neb_df$gene

    # Sort the results by p-value for the pathology of interest
    neb_df <- neb_df[order(neb_df[[paste0("p_", pathology)]]),]

    # Write the full results to a file (including all columns)
    write.table(neb_df, paste0(out_dir, name, ".", pathology, ".All.tsv"), sep = "\t", row.names = F, quote = F)

    # Write a cleaner version of the results (select columns of interest)
    results <- neb_df[, c("gene", paste0("p_", pathology), paste0("logFC_", pathology), "FDR", "log2FC")]
    write.table(results, paste0(out_dir, name, ".", pathology, ".Clean.tsv"), sep = "\t", row.names = F, quote = F)
}

# Example function calls for amyloid pathology or global AD pathology
# deg.nebula(Seurat_Obj, "amyloid", sample.col)
# deg.nebula(Seurat_Obj, "Global_AD_Pathology", sample.col)
