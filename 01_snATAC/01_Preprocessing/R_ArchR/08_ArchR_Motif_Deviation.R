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













##### ---------------------------
############### Motif Enrichment
############### ---------------------------
Subtype_proj3 <- addMotifAnnotations(ArchRProj = Subtype_proj3, 
                                     motifSet = "cisbp", 
                                     name = "Motif")


############### ---------------------------
###############  Motif Deviations
############### ---------------------------
if("Motif" %ni% names(Subtype_proj3@peakAnnotation)){
    Subtype_proj3 <- addMotifAnnotations(ArchRProj = Subtype_proj3, motifSet = "cisbp", name = "Motif")
}

Subtype_proj3 <- addBgdPeaks(Subtype_proj3)
# We are now ready to compute per-cell deviations accross all of our motif annotations using the addDeviationsMatrix() function. This function has an optional parameter called matrixName that allows us to define the name of the deviations matrix that will be stored in the Arrow files. If we do not provide a value to this parameter, as in the example below, this function creates a matrix name by adding the word “Matrix” to the name of the peakAnnotation. The example below creates a deviations matrix in each of our Arrow files called “MotifMatrix”.

Subtype_proj3 <- addDeviationsMatrix(
  ArchRProj = Subtype_proj3, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(Subtype_proj3, name = "MotifMatrix", plot = TRUE)


outd<-"/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_ADMR_multiome_snRNA_snATAC/02_multiome_snATAC/03_ArchR/01_Overall"
workdir<-paste0(outd,"/03_TSS2_Subtype")
setwd(workdir)

saveArchRProject(ArchRProj = Subtype_proj3, 
                outputDirectory = workdir,
                load = TRUE)




###################### heatmap showing the motif deviation for each cell type ######################
ChromVar_Zscore_Mat <- getGroupSE(
  ArchRProj = Subtype_proj3,
  useMatrix = "MotifMatrix",
  groupBy = "ATAC.Subtype.Mar5_2024",
  divideN = TRUE,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)
saveRDS(ChromVar_Zscore_Mat,"ChromVar_Zscore_Mat.Mar11_2024.rds")

ChromVar_Zscore_Mat_filter <- ChromVar_Zscore_Mat[rowData(ChromVar_Zscore_Mat)$seqnames=="z",]

#Then we can identify the maximum delta in z-score between all clusters. This will be helpful in stratifying motifs based on the degree of variation observed across clusters.
rowData(ChromVar_Zscore_Mat_filter)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

mat.tmp<-as.data.frame(assay(ChromVar_Zscore_Mat_filter,"MotifMatrix"))
mat.tmp2<-mat.tmp
mat.tmp2$ID<-rownames(mat.tmp)

annotation<-as.data.frame(rowData(ChromVar_Zscore_Mat_filter))
annotation$ID<-rownames(annotation)

mat.tmp3<-merge(annotation,mat.tmp2,by="ID")

write.table(mat.tmp3,file = "ChromVar_Zscore_Mat_by_ATAC.Subtype.Mar5_2024.tsv",quote = FALSE, sep = "\t",row.names=FALSE)


####################### --------------------------------------- #######################
#######################     Get group Genescore matrix/peak score matrix
####################### --------------------------------------- #######################

####### impute weights 
Subtype_proj3 <- addImputeWeights(Subtype_proj3)



### Get genescore matirx
genemtx = getMatrixFromProject(
  ArchRProj = Subtype_proj3,
  useMatrix = "GeneScoreMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
saveRDS(genemtx,"Mic_Class.peak_mtx.Sep6_2023.rds")

# group expression matrix for major cell type and samples
geneMat <- getGroupSE(
      ArchRProj = Subtype_proj3,
      useMatrix = "GeneScoreMatrix",
      groupBy = "ATAC.Subtype.Mar5_2024",
      divideN = TRUE,
      verbose = TRUE,
      logFile = createLogFile("getGroupSE"))

mat.tmp <- (assay(geneMat,"GeneScoreMatrix"))
genelist <- getFeatures(Subtype_proj3)

rownames(mat.tmp)<-genelist
write.table(data.frame("gene"=rownames(mat.tmp),mat.tmp),
            file = "GeneMat_by_ATAC.Subtype.Mar5_2024.tsv",quote = FALSE, sep = "\t",row.names=TRUE)


#get expression matrix for major cell type and samples
PeakMatrix <- getGroupSE(
      ArchRProj = Subtype_proj3,
      useMatrix = "PeakMatrix",
      groupBy = "ATAC.Subtype.Mar5_2024",
      divideN = TRUE,
      verbose = TRUE,
      logFile = createLogFile("getGroupSE"))

Matrix.mat = assays(PeakMatrix)$PeakMatrix
row_anno<-rowData(PeakMatrix)
rownames(Matrix.mat) = paste0(row_anno$seqnames,":",row_anno$start,"-",row_anno$end)
write.table(Matrix.mat,"ATAC.Subtype.PeakScoreMatrix.Mar12.2024.tsv",quote = F,
            sep = "\t",row.names = T)

# peak matrix for each sample celltype
Subtype_proj3$Celltype_Sample=paste0(Subtype_proj3$ATAC.Subtype.Mar5_2024,".",Subtype_proj3$Sample)

PeakMatrix <- getGroupSE(
      ArchRProj = Subtype_proj3,
      useMatrix = "PeakMatrix",
      groupBy = "Celltype_Sample",
      divideN = TRUE,
      verbose = TRUE,
      logFile = createLogFile("getGroupSE"))

Matrix.mat = assays(PeakMatrix)$PeakMatrix
row_anno<-rowData(PeakMatrix)
rownames(Matrix.mat) = paste0(row_anno$seqnames,":",row_anno$start,"-",row_anno$end)
write.table(Matrix.mat,"ATAC.Subtype_Sample.PeakScoreMatrix.Mar12.2024.tsv",quote = F,
            sep = "\t",row.names = T)



# Peak matrix for each brain region celltype
Subtype_proj3$Celltype_Region=paste0(Subtype_proj3$ATAC.Subtype.Mar5_2024,".",Subtype_proj3$BrainRegion)

PeakMatrix <- getGroupSE(
      ArchRProj = Subtype_proj3,
      useMatrix = "PeakMatrix",
      groupBy = "Celltype_Region",
      divideN = TRUE,
      verbose = TRUE,
      logFile = createLogFile("getGroupSE"))

Matrix.mat = assays(PeakMatrix)$PeakMatrix
row_anno<-rowData(PeakMatrix)
rownames(Matrix.mat) = paste0(row_anno$seqnames,":",row_anno$start,"-",row_anno$end)
write.table(Matrix.mat,"ATAC.Subtype_BrainRegion.PeakScoreMatrix.Mar12.2024.tsv",quote = F,
            sep = "\t",row.names = T)






############################################################################################################
############################################################################################################
############################################################################################################
# Remove low quality cells
indir="/net/bmc-lab4/data/kellis/users/zunpeng/05_results/02_ADMR_multiome_snRNA_snATAC/02_multiome_snATAC/02_h5ad/01_Overall/"
cells=read.csv(paste0(indir,"snATAC7.AllBrainRegion.nCells_1223768.Mar12_2024.Harmony_LibType.obs.csv.gz"),sep="\t")

Subtype_proj4=Subtype_proj3[row.names(Subtype_proj3@cellColData)%in%cells$X]


Subtype_proj4@cellColData$ATAC.Class.Mar5_2024
saveArchRProject(ArchRProj = Subtype_proj4, 
                 outputDirectory = paste0(outd,"/03_TSS2_Subtype2"), 
                 load = TRUE)








####################### --------------------------------------- #######################
#######################     Get group bigwig
####################### --------------------------------------- #######################
getGroupBW(
  ArchRProj = proj4,
  groupBy = "Class_Pathology",
  normMethod = "nFrags",
  tileSize = 100000, #The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
  maxCells = 50000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)



####################### --------------------------------------- #######################
#######################         Repeats score matrix
####################### --------------------------------------- #######################

library(GenomicRanges)
## get group bw for each individual for each major cell type
library(ArchR)
library(stringr)
library(dplyr)
### 
repeats=read.table("/net/bmc-lab5/data/kellis/users/zunpeng/03_Database/00_genome/02_GRCh38/05_repeats/hg38.rmsk.custom.final.clean1.bed",sep = "\t",header = F)
chr.list<-paste0("chr",c(seq(1,22,1),"X"))
repeats=repeats[repeats$V1 %in% chr.list,]
repeats.GRanges=GRanges(
  seqnames = repeats$V1,
  ranges = IRanges(start = repeats$V2,end = repeats$V3 ),
  strand=repeats$V6,
  gene_id=paste0(repeats$V5,";",repeats$V4),
  symbol=repeats$V4) 

proj4<-addGeneScoreMatrix(
  input = proj4,
  genes = repeats.GRanges,
  geneModel = "exp(-abs(x)/5000) + exp(-1)",
  matrixName = "RepeatScoreMatrix_tileSize_10",
  extendUpstream = c(0, 0),
  extendDownstream = c(0, 0),
  geneUpstream = 0,
  geneDownstream = 0,
  useGeneBoundaries = TRUE,
  useTSS = FALSE,
  extendTSS = FALSE,
  tileSize = 10,
  ceiling = 4,
  geneScaleFactor = 5,
  scaleTo = 10000,
  excludeChr = c("chrY", "chrM"),
  blacklist = getBlacklist(proj4),
  threads = getArchRThreads(),
  parallelParam = NULL,
  subThreading = TRUE,
  force = FALSE,
  logFile = createLogFile("RepeatScoreMatrix_tileSize_10")
)



####
Matrix<- getGroupSE(ArchRProj = proj4, useMatrix = "RepeatScoreMatrix_tileSize_10", groupBy = "Subclass_Sep2")
Matrix.mat = assays(Matrix)$RepeatScoreMatrix_tileSize_10
rownames(Matrix.mat) = rowData(Matrix)$name
write.table(data.frame("gene"=rownames(Matrix.mat),Matrix.mat),
"RepeatScoreMatrix_tileSize_10/Subclass_Sep2.RepeatScoreMatrix.tsv",quote = F,sep = "\t",row.names = F)

####
proj4@cellColData$Region=ifelse(proj4@cellColData$region=="MFC","PFC",proj4@cellColData$region)
proj4@cellColData$Subclass_Region<-paste0(proj4@cellColData$Subclass_Sep2,":",proj4@cellColData$Region)

Matrix<- getGroupSE(ArchRProj = proj4, useMatrix = "RepeatScoreMatrix_tileSize_10", groupBy = "Subclass_Region")
Matrix.mat = assays(Matrix)$RepeatScoreMatrix_tileSize_10
rownames(Matrix.mat) = rowData(Matrix)$name
write.table(data.frame("gene"=rownames(Matrix.mat),Matrix.mat),
"RepeatScoreMatrix_tileSize_10/Subclass_Region.RepeatScoreMatrix.tsv",quote = F,sep = "\t",row.names = F)

####
Matrix<- getGroupSE(ArchRProj = proj4, useMatrix = "RepeatScoreMatrix_tileSize_10", groupBy = "Class_Aug4")
Matrix.mat = assays(Matrix)$RepeatScoreMatrix_tileSize_10
rownames(Matrix.mat) = rowData(Matrix)$name
write.table(data.frame("gene"=rownames(Matrix.mat),Matrix.mat),
"RepeatScoreMatrix_tileSize_10/Class_Aug4.RepeatScoreMatrix.tsv",quote = F,sep = "\t",row.names = F)



########## FRIP
  #Remove project from args
  args$ArchRProj <- NULL

  #Run With Parallel or lapply
  outList <- .batchlapply(args)

  readsInPeaks <- lapply(outList, function(x) x$RIP) %>% unlist
  FRIP <- lapply(outList, function(x) x$FRIP) %>% unlist
  ArchRProj <- addCellColData(ArchRProj, data = readsInPeaks, name = "ReadsInPeaks", names(readsInPeaks), force = force)
  ArchRProj <- addCellColData(ArchRProj, data = FRIP, name = "FRIP", names(readsInPeaks), force = force)
  
  .endLogging(logFile = logFile)

    out <- list(ArrowFile = ArrowFile, RIP = insertionsInPeaks, FRIP = insertionsInPeaks / totalInsertions)







