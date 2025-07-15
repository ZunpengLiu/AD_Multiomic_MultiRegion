# Enrichr-Based Pathway and TF/Histone Modification Enrichment for DEGs

This R script performs **functional enrichment analysis** on gene lists using the [Enrichr](https://maayanlab.cloud/Enrichr/) API. It supports multiple gene set libraries, including:

- Gene Ontology (GO) terms
- KEGG, Reactome, WikiPathways
- Transcription factors (TFs) from ENCODE, ChEA, TRANSFAC
- Histone modification marks from the Epigenomics Roadmap

---

## Input

The script scans a specified directory (`input_dir`) for all `.tsv` files that contain **differentially expressed genes (DEGs)**. Each file must include a `Gene` column with gene symbols.

## Output

Results are saved under the folder `./EnrichR/`. Two types of output are generated for each input file:

- **Pathway enrichment**:

  - GO_Molecular_Function_2023
  - GO_Biological_Process_2023
  - GO_Cellular_Component_2023
  - WikiPathway_2023_Human
  - KEGG_2021_Human
  - Reactome_2022
- **TF & Histone modification enrichment**:

  - TRANSFAC_and_JASPAR_PWMs
  - ChEA_2022
  - ENCODE_TF_ChIP-seq_2015
  - TRRUST_Transcription_Factors_2019
  - Epigenomics_Roadmap_HM_ChIP-seq

Each output `.tsv` file contains enriched terms, their p-values, overlap statistics, and source database.

## Dependencies

- R >= 4.0
- [enrichR](https://cran.r-project.org/web/packages/enrichR/)
- stringr

### Install dependencies:

```r
install.packages(c("enrichR", "stringr"))
```
