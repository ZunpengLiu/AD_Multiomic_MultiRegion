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
## 2. Peak Calling
### --------------------------
# Grouping cells by Cell_Subtype 
proj3 <- addGroupCoverages(
  ArchRProj = proj3,
  groupBy = "Cell_Subtype",
  useLabels = TRUE,
  minCells = 40,
  maxCells = 10000,
  maxFragments = 25 * 10^6,
  minReplicates = 2,
  maxReplicates = 40,
  sampleRatio = 0.8,
  kmerLength = 6,
  threads = getArchRThreads(),
  returnGroups = FALSE,
  parallelParam = NULL,
  force = TRUE,
  verbose = TRUE,
  logFile = createLogFile("addGroupCoverages")
)

pathToMacs2 <- "~/04_Softwares/conda/anaconda3/envs/macs2/bin/macs2"

proj3 <- addReproduciblePeakSet(
  ArchRProj = proj3,
  maxPeaks = 500000,
  cutOff = 0.01,
  force = TRUE,
  groupBy = "Cell_Subtype",
  pathToMacs2 = pathToMacs2
)

#### writing the union peak set to a tsv file
peakSet_df<-data.frame(proj3@peakSet)
peakSet_df$Cell_type<-names(proj3@peakSet)
write.table(peakSet_df,"proj3.peakSet.ATAC.Cell_Subtype.tsv",quote = F,sep = "\t",row.names = F)

#### Add peak matrix
proj3 <- addPeakMatrix(proj3)

saveArchRProject(ArchRProj = proj3, 
                outputDirectory = workdir,
                load = TRUE)

### --------------------------
## 3. Export peak as bed
### --------------------------
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
write.table(all,"All_reproduciblePeaks.tsv",quote = F,sep = "\t",row.names = F)

system("for i in `ls -d *bed`; do sort -k1,1 -k2,2n $i > ${i%.*}".srt.bed";done")
