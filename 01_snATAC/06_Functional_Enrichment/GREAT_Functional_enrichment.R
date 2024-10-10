## use rGREAT to annotate the cell-type specific peaks, use all the peaks as background
library(ggplot2)
library(rGREAT)
library("bedr")
library(stringr)

## set the working directory
setwd("./Module/Peaks")

## read the peak file
df = read.table("./Subtype.All_peakSet.tsv",head=T,sep="\t")

## get peak loc
df = df[,c(1,2,3)]
allpeaks = df
allpeaks.loc = paste0(allpeaks[,1],":",allpeaks[,2],"-",allpeaks[,3])
allpeaks.loc = bedr.sort.region(allpeaks.loc)

samples<-list.files("./",pattern = ".sorted.bed$",full.names = F)

sample_name<-str_remove(samples,".sorted.bed")

for(ct in sample_name){
  print(ct)
  df = read.csv(paste0(ct,".sorted.bed"),sep = "\t",header=F)
  colnames(df)<-c("seqnames","start","end")
  # do an overlap
  df.loc = paste0(df[,1],":",df[,2],"-",df[,3])
  df.loc = bedr.sort.region(df.loc)
  df.intersect = bedr.join.region(df.loc,allpeaks.loc)
  
  df = unique(df.intersect[,2:4])
  colnames(df) = colnames(allpeaks)
  df = df[df$seqnames != ".",]
  df$start = as.numeric(df$start)
  df$end = as.numeric(df$end)
  
  job = submitGreatJob(df,bg=allpeaks,species="hg38")
  tbl = getEnrichmentTables(job)
  write.table(tbl[1],paste0("../GREAT/GO_MF_",ct,".rGREAT.out.tsv"),sep = "\t",row.names = F,quote = F)
  write.table(tbl[2],paste0("../GREAT/GO_BP_",ct,".rGREAT.out.tsv"),sep = "\t",row.names = F,quote = F)
  write.table(tbl[3],paste0("../GREAT/GO_CC_",ct,".rGREAT.out.tsv"),sep = "\t",row.names = F,quote = F)
}