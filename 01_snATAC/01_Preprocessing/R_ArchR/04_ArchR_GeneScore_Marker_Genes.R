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
workdir<-paste0(outd,"/03_High_Quality")

proj3<-loadArchRProject(workdir)

### --------------------------
##  2. Calculating Gene Scores in ArchR
### --------------------------
proj3 <- addImputeWeights(proj3)

### --------------------------
##  3. Get Genescore matrix
### --------------------------
genemtx = getMatrixFromProject(
  ArchRProj = proj3,
  useMatrix = "GeneScoreMatrix",
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)
saveRDS(genemtx,"proj3.GeneScoreMatrix.rds")

### --------------------------
##  4. Export Group Summarized Experiment Genescore matrix
### --------------------------
geneMat <- getGroupSE(
      ArchRProj = proj3,
      useMatrix = "GeneScoreMatrix",
      groupBy = "Cell_Subtype",
      divideN = TRUE,
      verbose = TRUE,
      logFile = createLogFile("getGroupSE"))

mat.tmp <- (assay(geneMat,"GeneScoreMatrix"))
genelist <- getFeatures(proj3)

rownames(mat.tmp)<-genelist
write.table(data.frame("gene"=rownames(mat.tmp),mat.tmp),
            file = "GeneMat_by_Cell_Subtype.tsv",
            quote = FALSE, 
            sep = "\t",
            row.names=TRUE)

### --------------------------
##  5. UMAP Plotting of the canonical marker genes for major cell types
### --------------------------
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

### --------------------------
##  6. Identify marker genes for each cell subtypes
### --------------------------
markersGS <- getMarkerFeatures(
    ArchRProj = proj3, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Cell_Subtype",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

saveRDS(markersGS,"proj3.Cell_Subtype.Marker_genes.rds")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")

# save the marker genes to csv
for (i in 1:length(markerList)) {
    write.csv(markerList[[i]], file = paste0("proj3.Cell_Subtype.Marker_genes.",names(markerList)[i],".csv"),quote = FALSE,row.names = FALSE)
}