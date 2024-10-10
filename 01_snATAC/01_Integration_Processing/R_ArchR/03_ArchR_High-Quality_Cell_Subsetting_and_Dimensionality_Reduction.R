library(GenomicRanges)
library(ArchR)
library(stringr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 48)
addArchRGenome("hg38")   ## use hg38

### --------------------------
## 1. Loading an ArchRProject
### --------------------------
outd<-"./03_ArchR"
workdir<-paste0(outd,"/02_filterDoublets")

proj2<-loadArchRProject(workdir)

### --------------------------
## 2. get high quality cells processed from muon and scanpy
### --------------------------
indir="./02_h5ad/"
cells=read.csv(paste0(indir,"snATAC.AllBrainRegion.Harmony_LibType.obs.tsv.gz"),sep="\t")

# subset the cells
proj3=proj2[row.names(proj2@cellColData)%in%cells$X]
saveArchRProject(ArchRProj = proj3, 
                 outputDirectory = paste0(outd,"/03_High_Quality"), 
                 load = TRUE)
# write the cellColData
write.table(proj3@cellColData,
            paste0(outd,"/03_High_Quality/proj3_TSS2_filterDoublets.cellColData.tsv"),
            sep="\t",
            quote=F,
            row.names=T)

### --------------------------
## 3. TSS Enrichment and Fragment Size Distribution
### --------------------------
p1 <- plotFragmentSizes(ArchRProj = Subtype_proj3,  
                        groupBy = "Sample_Name")
p2 <- plotTSSEnrichment(ArchRProj = Subtype_proj3,  
                          groupBy = "Sample_Name",
                          flank = 3000)
plotPDF(p1,p2, 
    name = "QC-Sample-FragSizes-TSSProfile.pdf",
    ArchRProj = Subtype_proj3, addDOC = FALSE, width = 5, height = 5)

p1 <- plotFragmentSizes(ArchRProj = Subtype_proj3,  
                        groupBy = "CD45")
p2 <- plotTSSEnrichment(ArchRProj = Subtype_proj3,  
                          groupBy = "CD45",
                          flank = 3000)

plotPDF(p1,p2, 
    name = "QC-Sample-FragSizes-TSSProfile_CD45.pdf",
    ArchRProj = Subtype_proj3, addDOC = FALSE, width = 5, height = 5)


### --------------------------
## 4. Redo the Dimensionality Reduction in ArchR based on the high quality cells, Clustering, and UMAP.
## The cell type clustering from the ArchR is not used in the downstream analysis, 
## but it is worth to check the UMAP and compare it with the UMAP from the Munon and Scanpy processed data.
### --------------------------
proj3 <- addIterativeLSI(
    ArchRProj = proj3,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 2, 
    clusterParams = list( #See Seurat::FindClusters
        resolution = c(0.2), 
        sampleCells = 80000, 
        n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30
)

# Batch Effect Correction wtih Harmony
proj3 <- addHarmony(
    ArchRProj = proj3,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample"
)

# Clustering
proj3 <- addClusters(
    input = proj3,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.8
)

# Uniform Manifold Approximation and Projection (UMAP) 
proj3 <- addUMAP(
    ArchRProj = proj2_TSS4, 
    reducedDims = "proj3", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)
proj3 <- addUMAP(
    ArchRProj = proj3, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

saveArchRProject(ArchRProj = proj3, 
                 outputDirectory = paste0(outd,"/03_High_Quality"), 
                 load = TRUE)

### --------------------------
## 5. Get the UMAP from the Munon and Scanpy processed data.
### --------------------------
X_umap<-read.csv(paste0(indir,"snATAC.AllBrainRegion.Harmony_LibType.X_umap.tsv"),sep = "\t")
X_umap2<-X_umap[,-1]
rownames(X_umap2)<-X_umap[,1]
X_umap3<-X_umap2[rownames(proj3@cellColData),]
proj3@embeddings$X_umap$df<-X_umap3
colnames(proj3@embeddings$X_umap$df)<-c("IterativeLSI#UMAP_Dimension_1","IterativeLSI#UMAP_Dimension_2")

### --------------------------
## 6. Plot the finalized version of UMAP in ArchR for the high quality cells.
### --------------------------
p1 <- plotEmbedding(ArchRProj = proj3, 
                    colorBy = "cellColData", 
                    name = "Cell_Subtype", 
                    embedding = "X_umap")

plotPDF(p1, name = "Plot-X_umap-Cell_Subtype.pdf",
            ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)