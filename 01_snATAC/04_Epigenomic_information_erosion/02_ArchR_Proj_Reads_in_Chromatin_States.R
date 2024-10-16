library(GenomicRanges)
library(ArchR)
library(stringr)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)

addArchRThreads(threads = 60)
addArchRGenome("hg38")   ## use hg38

### --------------------------
## 1. load the ArchR project
### --------------------------

outd<-"./03_ArchR/01_Overall"
workdir<-paste0(outd,"/03_ArchRProject")
setwd(workdir)
proj3<-loadArchRProject(workdir)

### --------------------------
## 2. Read in the Chromatin States
### --------------------------
argv <- commandArgs(T)

tissue <- argv[1]

States<-c("EnhA1","EnhA2","EnhBiv","EnhG1","EnhG2",
          "EnhWk","Het","Quies","ReprPC","ReprPCWk",
          "TssA","TssBiv","TssFlnk","TssFlnkD","TssFlnkU",
          "Tx","TxWk","ZNF.Rpts")

outdir="./05_ErosionScore"

for (state in States){
  print(paste0("Processing ",tissue," and ",state))
  ChromatinState=read.table(paste0("epimap_Brain/",tissue,"/",state,".bed"),sep = "\t",header = F)
  chr.list<-paste0("chr",c(seq(1,22,1),"X"))
  ChromatinState=ChromatinState[ChromatinState$V1 %in% chr.list,]
  ChromatinState$V6<- "*"
  ChromatinState.GRanges=GRanges(
    seqnames = ChromatinState$V1,
    ranges = IRanges(start = ChromatinState$V2,end = ChromatinState$V3 ),
    strand=ChromatinState$V6,
    gene_id=paste0(ChromatinState$V1,":",ChromatinState$V2,"-",ChromatinState$V3),
    symbol=ChromatinState$V4)

# Add the Chromatin States to the ArchR project and calculate the fraction of reads in each state
proj3 <- addFeatureCounts(
  ArchRProj = proj3,
  features = ChromatinState.GRanges,
  name = paste0(tissue,".",state,"_"),
  addRatio = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("addFeatureCounts")
)
write.table(proj3@cellColData,file=paste0(outdir,"/",tissue,".Count_Fractions.txt"),sep="\t",quote=F)

}

### --------------------------
## 3. Run this script in the sbatch script
### --------------------------
# Example of the sbatch script
'''
#!/usr/bin/sh
#SBATCH -p kellis
#SBATCH -t 60:00:00
#SBATCH --mem=600G
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --job-name=BSS01125

cd ./05_ErosionScore
source activate R
Rscript ChromatinStates.R BSS01125
sbatch -p kellis --nodelist=b4 BSS01125
'''