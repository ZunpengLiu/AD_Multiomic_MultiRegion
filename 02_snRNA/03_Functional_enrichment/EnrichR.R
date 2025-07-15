########## Enrichr pathway enrichment for DEGs ##########

# Set working directory (if needed)
input_dir <- "./"        # Path to input DEG files
output_dir <- "./EnrichR" # Path to output results
dir.create(output_dir, showWarnings = FALSE)

# Load libraries
library(enrichR)
library(stringr)

# Set Enrichr base address (optional override)
options(base_address = "http://amp.pharm.mssm.edu/Enrichr/")

# Load database list once (for reference)
available_dbs <- listEnrichrDbs()

# Parameters
cutoff <- 0.05  # Enrichment p-value cutoff
input_files <- list.files(input_dir, pattern = "\\.tsv$", full.names = TRUE)

# Helper function to run enrichment and save result
run_enrichment <- function(file_path, db_list, file_tag = "Pathway") {
  cat("Processing:", file_path, "\n")

  # Read DEG result table
  deg_df <- read.table(file_path, header = TRUE)

  # Use all genes listed in the "Gene" column (no filtering by log2FC/FDR)
  if (!"Gene" %in% colnames(deg_df)) {
    warning("No 'Gene' column found in:", file_path)
    return(NULL)
  }

  gene_list <- unique(as.character(deg_df$Gene))
  if (length(gene_list) == 0) {
    warning("No genes found in:", file_path)
    return(NULL)
  }

  # Run Enrichr
  enriched <- enrichr(gene_list, db_list)

  # Extract sample name from file
  sample_name <- str_replace(basename(file_path), ".tsv$", "")
  output_file <- file.path(output_dir, paste0(file_tag, ".", sample_name, ".tsv"))

  # Combine all significant terms across databases
  combined_results <- data.frame()
  for (db in db_list) {
    if (!db %in% names(enriched)) next
    df <- enriched[[db]]
    df <- subset(df, P.value < cutoff)
    if (nrow(df) > 0) {
      df$Source <- db
      combined_results <- rbind(combined_results, df)
    }
  }

  # Save result if not empty
  if (nrow(combined_results) > 0) {
    combined_results <- combined_results[order(combined_results$P.value), ]
    write.table(combined_results, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

############ Enrichment Runs ############

# 1. GO & Pathways
db_set_pathway <- c(
  "GO_Molecular_Function_2023",
  "GO_Biological_Process_2023",
  "GO_Cellular_Component_2023",
  "WikiPathway_2023_Human",
  "KEGG_2021_Human",
  "Reactome_2022"
)

for (f in input_files) {
  run_enrichment(f, db_set_pathway, file_tag = "Pathway")
}

# 2. TF enrichment and histone modifications
db_set_tf <- c(
  "TRANSFAC_and_JASPAR_PWMs",
  "ChEA_2022",
  "ENCODE_TF_ChIP-seq_2015",
  "TRRUST_Transcription_Factors_2019",
  "Epigenomics_Roadmap_HM_ChIP-seq"
)

for (f in input_files) {
  run_enrichment(f, db_set_tf, file_tag = "TF_HistoneModifications")
}