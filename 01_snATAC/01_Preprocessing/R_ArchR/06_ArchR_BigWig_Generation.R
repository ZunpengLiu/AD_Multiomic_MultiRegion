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
workdir<-paste0(outd,"/02_TSS2_filterDoublets")

proj2<-loadArchRProject(workdir)

# get high quality cells processed from muon and scanpy
indir="./02_h5ad/"
cells=read.csv(paste0(indir,"snATAC.AllBrainRegion.Harmony_LibType.obs.tsv.gz"),sep="\t")

# subset the cells
proj3=proj2[row.names(proj2@cellColData)%in%cells$X]
saveArchRProject(ArchRProj = proj3, 
                 outputDirectory = paste0(outd,"/03_TSS2_High_Quality"), 
                 load = TRUE)
# write the cellColData
write.table(proj3@cellColData,
            paste0(outd,"/03_TSS2_High_Quality/proj3_TSS2_filterDoublets.cellColData.tsv"),
            sep="\t",
            quote=F,
            row.names=T)

### --------------------------
## 2. Dimensionality Reduction, Clustering and Uniform Manifold Approximation and Projection (UMAP) 
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
                 outputDirectory = paste0(outd,"/03_TSS2_High_Quality"), 
                 load = TRUE)


### --------------------------
## 2. Dimensionality Reduction, Clustering and Uniform Manifold Approximation and Projection (UMAP) 
### --------------------------
p1 <- plotEmbedding(ArchRProj = proj3, 
                    colorBy = "cellColData", 
                    name = "ATAC.Class.Mar5_2024", 
                    embedding = "X_umap")

p2 <- plotEmbedding(ArchRProj = proj3, 
                    colorBy = "cellColData", 
                    name = "ATAC.Subtype.Mar5_2024", 
                    embedding = "X_umap")

ggAlignPlots(p1, p2, type = "h")
plotPDF(p1,p2, name = "Plot-X_umap-Sample-Clusters.pdf",
            ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)


############ Plot QC figures
p1 <- plotFragmentSizes(ArchRProj = proj3,  
                        groupBy = "Sample_Name")
p2 <- plotTSSEnrichment(ArchRProj = proj3,  
                          groupBy = "Sample_Name",
                          flank = 3000)
plotPDF(p1,p2, 
    name = "QC-Sample-FragSizes-TSSProfile.pdf",
    ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)

p1 <- plotFragmentSizes(ArchRProj = proj3,  
                        groupBy = "CD45")
p2 <- plotTSSEnrichment(ArchRProj = proj3,  
                          groupBy = "CD45",
                          flank = 3000)

plotPDF(p1,p2, 
    name = "QC-Sample-FragSizes-TSSProfile_CD45.pdf",
    ArchRProj = proj3, addDOC = FALSE, width = 5, height = 5)


############
X_umap<-read.csv(paste0(indir,"snATAC6.AllBrainRegion.nCells_1225965.Mar5_2024.Harmony_LibType.X_umap.tsv"),sep = "\t")
X_umap2<-X_umap[,-1]
rownames(X_umap2)<-X_umap[,1]
X_umap3<-X_umap2[rownames(proj3@cellColData),]
proj3@embeddings$X_umap$df<-X_umap3
colnames(proj3@embeddings$X_umap$df)<-c("IterativeLSI#UMAP_Dimension_1","IterativeLSI#UMAP_Dimension_2")


############# Marker genes
proj3 <- addImputeWeights(proj3)


major_celltype_markergenes<-c("SLC17A7","GAD2","GFAP","MOBP","VCAN","P2RY12","CD8A")

p1 <- plotEmbedding(
    ArchRProj = proj3, 
    colorBy = "GeneScoreMatrix", 
    name = major_celltype_markergenes, 
    pal=ArchRPalettes$solarExtra,
    embedding = "X_umap",
    size = 0.1,
    quantCut = c(0.01, 0.99),
    imputeWeights = getImputeWeights(proj3)
)

plotPDF(p1, name = "GeneScore.Key_Marker_genes_solarExtra.pdf",
        ArchRProj = proj3,
        addDOC = FALSE,
        width = 5, height = 5)


p1 <- plotEmbedding(
    ArchRProj = proj3, 
    colorBy = "GeneScoreMatrix", 
    name = major_celltype_markergenes, 
    pal=ArchRPalettes$solarExtra,
    embedding = "X_umap",
    size = 0.1,
    quantCut = c(0.05, 0.95),
    imputeWeights = getImputeWeights(proj3)
)

plotPDF(p1, name = "GeneScore.Key_Marker_genes_solarExtra2.pdf",
        ArchRProj = proj3,
        addDOC = FALSE,
        width = 5, height = 5)




############################### marker genes for each cluster
markersGS <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "ATAC.Class.Mar5_2024",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

saveRDS(markersGS,"proj3.Class_Mar5_2024.Marker_genes.rds")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")

# save the marker genes to csv
for (i in 1:length(markerList)) {
    write.csv(markerList[[i]], file = paste0("proj3.Class_Feb12_2024.Marker_genes.",names(markerList)[i],".csv"),quote = FALSE,row.names = FALSE)
}


######
markersGS <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "ATAC.Class.Mar5_2024",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

saveRDS(markersGS,"proj3.Subclass_Feb17_2024.Marker_genes.rds")

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 0.25")

# save the marker genes to csv
for (i in 1:length(markerList)) {
    fname=names(markerList)[i]
    fname=stringr::str_replace(fname,"/",".")

    write.csv(markerList[[i]], file = paste0("proj3.Subclass_Feb17_2024.Marker_genes.",fname,".csv"),quote = FALSE,row.names = FALSE)
}



##------------------------ call peaks  ------------------------##

proj3 <- addGroupCoverages(
  ArchRProj = proj3,
  groupBy = "ATAC.Class_Region.Mar9_2024",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 10000,
  maxFragments = 25 * 10^6,
  minReplicates = 2,
  maxReplicates = 40, #The maximum number of cells to use during insertion coverage file generation. default 5.
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

pathToMacs2 <- "/net/bmc-lab5/data/kellis/users/zunpeng/04_Softwares/conda/anaconda3/envs/macs2/bin/macs2"

proj3 <- addReproduciblePeakSet(
  ArchRProj = proj3,
  maxPeaks = 500000,
  cutOff = 0.01,
  force = TRUE,
  groupBy = "ATAC.Class_Region.Mar9_2024",
  pathToMacs2 = pathToMacs2
)
####
peakSet_df<-data.frame(proj3@peakSet)
peakSet_df$Cell_type<-names(proj3@peakSet)
write.table(peakSet_df,"proj3.peakSet.ATAC.Class_Region.Mar9_2024.tsv",quote = F,sep = "\t",row.names = F)


proj3 <- addPeakMatrix(proj3)

saveArchRProject(ArchRProj = proj3, 
                outputDirectory = workdir,
                load = TRUE)

############### ---------------------------
############### Export peak as bed
############### ---------------------------
setwd(paste0(workdir,"/PeakCalls"))

samples<-list.files(pattern = "reproduciblePeaks.gr.rds")
library(stringr)
all=data.frame()
for (i in 1:length(samples)){
  tmp<-readRDS(samples[i])
  tmp2<-data.frame(tmp)
  Cell_type<-str_replace(samples[i],"-reproduciblePeaks.gr.rds","")
  tmp2$Cell_type<-str_replace(samples[i],"-reproduciblePeaks.gr.rds","")
  all<-rbind(all,tmp2)
  write.table(tmp2,paste0(Cell_type,".tsv"),quote = F,sep = "\t",row.names = F)
  write.table(tmp2[,c(1,2,3)],paste0(Cell_type,"_macs2_peak.bed"),quote = F,sep = "\t",row.names = F,col.names = F)

}
write.table(all,"ADMR_Class_Region_reproduciblePeaks.Mar7_2024.tsv",quote = F,sep = "\t",row.names = F)


system("for i in `ls -d *bed`; do sort -k1,1 -k2,2n $i > ${i%.*}".srt.bed";done")



















