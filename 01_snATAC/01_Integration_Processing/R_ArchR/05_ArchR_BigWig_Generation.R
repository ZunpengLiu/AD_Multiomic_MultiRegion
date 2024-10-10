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
## 2. Get group bigwig
### --------------------------
for (i in c("Cell_Class","Cell_Subclass","Cell_Subtype","Class_Pathology","Subclass_Pathology")){
  getGroupBW(
    ArchRProj = proj3,
    groupBy = i,
    normMethod = "nFrags",
    tileSize = 100000, #The numeric width of the tile/bin in basepairs for plotting ATAC-seq signal tracks. All insertions in a single bin will be summed.
    maxCells = 50000,
    ceiling = 4,
    verbose = TRUE,
    threads = getArchRThreads(),
    logFile = createLogFile("getGroupBW")
  )

}