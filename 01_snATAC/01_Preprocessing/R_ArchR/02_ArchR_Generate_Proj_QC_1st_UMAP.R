library(GenomicRanges)
library(ArchR)
library(stringr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)


addArchRThreads(threads = 60)
addArchRGenome("hg38")   ## use hg38

### --------------------------
## 1. read in Arrow files
### --------------------------

setwd("./01_Fragments/all_arrow")

ArrowFiles <- list.files(".", pattern = ".arrow$", full.names = TRUE)

## extract names of the samples
names(inputFiles) <- str_remove(str_remove(ArrowFiles,".arrow"),"./") 
samples <-  str_remove(str_remove(ArrowFiles,".arrow"),"./")


### --------------------------
## 2. Creating an ArchRProject
### --------------------------
outd<-"./03_ArchR/"

proj1 <- ArchRProject(
                        ArrowFiles = ArrowFiles, 
                        outputDirectory = paste0(outd,"/01_TSS2"),
                        copyArrows = TRUE 
                        #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
saveArchRProject(ArchRProj = proj1, 
                 outputDirectory = workdir, 
                 load = TRUE)

### --------------------------
## 3. Filter Doublets
### --------------------------
proj2 <- filterDoublets(proj1)

saveArchRProject(ArchRProj = proj2, 
                 outputDirectory = paste0(outd,"/02_TSS2_filterDoublets"), 
                 load = TRUE)

write.table(proj2@cellColData,
            paste0(outd,"/02_TSS2_filterDoublets/proj2_TSS2_filterDoublets.cellColData.tsv"),
            sep="\t",
            quote=F,
            row.names=T)

### --------------------------
## 4. Initial Dimensionality Reduction, Clustering and Uniform Manifold Approximation and Projection (UMAP) 
### --------------------------
proj2 <- addIterativeLSI(
    ArchRProj = proj2,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 10000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

# Batch Effect Correction wtih Harmony
proj2 <- addHarmony(
    ArchRProj = proj2,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

# Clustering
proj2 <- addClusters(
    input = proj2,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

# Uniform Manifold Approximation and Projection (UMAP) 
proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
proj2 <- addUMAP(
    ArchRProj = proj2, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

saveArchRProject(ArchRProj = proj2, 
                 outputDirectory = paste0(outd,"/02_TSS2_filterDoublets"), 
                 load = TRUE)
