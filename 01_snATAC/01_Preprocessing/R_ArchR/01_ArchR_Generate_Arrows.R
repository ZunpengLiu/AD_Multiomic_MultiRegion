library(GenomicRanges)
library(ArchR)
library(stringr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 60)
addArchRGenome("hg38")   ## use hg38

### ----------------------------------------------------
## 1. read in fragments data and generate Arrow files
### ----------------------------------------------------
setwd("./01_Fragments")

inputFiles <- list.files(".", pattern = ".gz$", full.names = TRUE)

names(inputFiles) <- str_remove(str_remove(inputFiles,".atac_fragments.tsv.gz"),"./")     ## extract names of the samples
samples <-  str_remove(str_remove(inputFiles,".atac_fragments.tsv.gz"),"./")

## create Arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 2,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  QCDir = "QualityControl",
  threads = getArchRThreads()
)

ArrowFiles <- list.files(".", pattern = ".arrow$", full.names = TRUE)

doubScores <- addDoubletScores(
    input = ArrowFiles,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1
)
