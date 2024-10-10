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
## 2. Motif Deviations, Adding motif annotations
### --------------------------
proj3 <- addMotifAnnotations(ArchRProj = proj3, 
                                     motifSet = "cisbp", 
                                     name = "Motif")

### --------------------------
## 3. Motif Deviations, Adding background peaks
### --------------------------
proj3 <- addBgdPeaks(proj3)

### --------------------------
## 4. Compute per-cell deviations accross all of our motif annotations using the addDeviationsMatrix() function. 
### --------------------------
proj3 <- addDeviationsMatrix(
  ArchRProj = proj3, 
  peakAnnotation = "Motif",
  force = TRUE
)

### --------------------------
# 5. Get Variable Deviations across cells in ArchRProject.
### --------------------------
plotVarDev <- getVarDeviations(proj3, name = "MotifMatrix", plot = TRUE)

### --------------------------
# 6. Get Group Summarized Experiment for MotifMatrix
### --------------------------
ChromVar_Zscore_Mat <- getGroupSE(
  ArchRProj = proj3,
  useMatrix = "MotifMatrix",
  groupBy = "Cell_subtype",
  divideN = TRUE,
  verbose = TRUE,
  logFile = createLogFile("getGroupSE")
)
saveRDS(ChromVar_Zscore_Mat,"ChromVar_Zscore_Mat.rds")

ChromVar_Zscore_Mat_filter <- ChromVar_Zscore_Mat[rowData(ChromVar_Zscore_Mat)$seqnames=="z",]

# Then we can identify the maximum delta in z-score between all clusters. 
# This will be helpful in stratifying motifs based on the degree of variation observed across clusters.
rowData(ChromVar_Zscore_Mat_filter)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

mat.tmp<-as.data.frame(assay(ChromVar_Zscore_Mat_filter,"MotifMatrix"))
mat.tmp2<-mat.tmp
mat.tmp2$ID<-rownames(mat.tmp)

annotation<-as.data.frame(rowData(ChromVar_Zscore_Mat_filter))
annotation$ID<-rownames(annotation)

mat.tmp3<-merge(annotation,mat.tmp2,by="ID")

write.table(mat.tmp3,
            file = "ChromVar_Zscore_Mat_by_Cell_subtype.tsv",
            quote = FALSE, 
            sep = "\t",
            row.names=FALSE)

### --------------------------
# 7. Save the ArchRProject
### --------------------------
saveArchRProject(ArchRProj = proj3, 
                outputDirectory = workdir,
                load = TRUE)