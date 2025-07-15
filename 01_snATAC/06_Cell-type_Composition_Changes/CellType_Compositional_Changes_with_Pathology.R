# Cell propotion changes in each Class Sub class and brain regions in snRNA and snATAC 
setwd("~/zunpeng@mit.edu - Google Drive/My Drive/01_Data/01_ADMR_snATAC_multiome_2024/01_Overall/CellFractionChanges")
library(dplyr)
library(data.table)
library(ggpubr)
library(ggplot2)
library(ggthemes)
#install.packages("oddsratio")
library(oddsratio)
library(reshape2)
library(pheatmap)
library(ComplexHeatmap)

###########################
# 1. Load the dataset
###########################
meta<-fread("./All.Cell.obs.tsv.gz",sep = "\t")
all<-meta
all$BrainRegion<-"All"
all<-rbind(all,meta2)
anno<-unique(meta[,c("Sample","projid","Pathology","BrainRegion","LibType","ATAC.Class","ATAC.Subclass","ATAC.Subtype")])

all$Class_f<-factor(all$ATAC.Class,levels = rev(c("Exc","Inh","Ast" ,"Oli" ,"OPC" ,"Mic_Immune","Vasc_Epithelia" )))
all$Region_f<-factor(all$BrainRegion,levels = rev(c("All","EC","HC" ,"TH" ,"AG" ,"MTC" ,"PFC" )))


#######################
# 2. Define Subclasses and colors for visualization
#######################
Subclasses<-c(
  "Exc L2-3 IT","Exc L3-4 IT","Exc L3-5 IT","Exc L4-5 IT-1","Exc L4-5 IT-2",
  "Exc L5-6 IT", "Exc L5/6 IT Car3","Exc L6 IT","Exc L5/6 NP","Exc L5 ET",
  "Exc L6 CT","Exc L6b",
  "Exc EC", #previously defined in EC
  "Exc CA pyramidal cells","Exc DG granule cells","Exc HC", # previously defined in HC_EC
  "Exc TH",
  "Inh VIP","Inh LAMP5","Inh SST","Inh PVALB","Inh PAX6","Inh MEIS2",
  "Ast","Oli","OPC",
  "Mic","T",
  "Per","End","SMC","Fib","CPEC","Epd")

Subclass_colors<-c("#69AED4","#B7B44B","#D96B7C","#69BE54","#939325",
                   "#E5E576","#127845","#761A55","#6A32A0","#77471F",
                   "#292774","#AA7947",
                   "#FEEE75",
                   "#61DDF2","#71ACD8","#BD80B8",
                   "#7472B5",
                   '#4B96D1', #Inh VIP
                   '#329CFC', # Inh LAMP5
                   "#8BCEA0", # Inh SST
                   "#D6569F",# Inh PVALB
                   "#FFB531", '#FDD3BA',# Inh PAX6
                   "#D62729","#FFBA7B","#B25828", # Ast Oli OPC 
                   "#CFA2D2","#9A6BBF", #Mic
                   "#E283B5","#F2C0D5","#6EC5A6","#A45627","#F38C64","#8BA0CC")

all$Subclass_f<-factor(all$ATAC.Subclass,levels = rev(Subclasses))

Subtypes<-c("Exc L2-3 CBLN2 LINC02306","Exc L3-4 RORB CUX2","Exc L3-5 RORB PLCH1","Exc L4-5 RORB GABRG1","Exc L4-5 RORB IL1RAPL2",
            "Exc L5-6 RORB LINC02196","Exc L5/6 IT Car3","Exc L6 THEMIS NFIA","Exc L5/6 NP","Exc L5 ET",
            "Exc L6 CT","Exc L6b",    # 1. Exc cortical layers
            "Exc RELN COL5A2","Exc RELN GPC5","Exc AGBL1 GPC5","Exc TOX3 GPC5","Exc TOX3 COL25A1",
            "Exc DLC1 SNTG2",    # 3. Exc EC layers
            "Exc CA1 pyramidal cells","Exc CA2, CA3 pyramidal cells","Exc DG granule cells","Exc COBLL1 UST",
            "Exc ZNF385D COL24A1","Exc GRIK1 CTXND1 (Subiculum)",     # 2. Exc HC layers
            "Exc NXPH1 RNF220",# Thalamus
            "Inh VIP ABI3BP","Inh RYR3 TSHZ2","Inh ALCAM TRPM3","Inh VIP THSD7B","Inh VIP TSHZ2",
            "Inh SORCS1 TTN","Inh SGCD PDE3A","Inh PTPRK FAM19A1","Inh VIP CLSTN2",
            "Inh LAMP5 NRG1 (Rosehip)","Inh L1-6 LAMP5 CA13","Inh LAMP5 RELN",
            "Inh L3-5 SST MAFB","Inh CUX2 MSR1","Inh ENOX2 SPHKAP","Inh FBN2 EPB41L4A",#SST
            "Inh PVALB SULF1","Inh PVALB HTR4","Inh PVALB CA8 (Chandelier)","Inh GPC5 RIT2",
            "Inh PAX6 RELN","Inh MEIS2 FOXP2",
            "Ast DPP10","Ast FOXB1","Ast GRM3", "Ast DCLK1",
            "Oli DPP10","Oli RASGRF1","Oli OPALIN",
            "OPC GPC5","OPC OLIG2","OPC TPST1",
            "Mic P2RY12","Mic TPT1",
            "CAMs","T",
            "Per","End","SMC","Fib","CPEC","Epd")

Subtype_colors<-c("#69AED4","#B7B44B","#D96B7C","#69BE54","#939325",
                  "#E5E576","#127845","#761A55","#6A32A0","#77471F",
                  "#292774","#AA7947",
                  "#DFAC7B","#7DC783","#0A8CED","#FEEE75","#E3AD28",
                  "#5576A6",
                  "#61DDF2","#4888C7","#71ACD8","#BD80B8",
                  "#34A07C","#AC4689",
                  "#7472B5",
                  '#FBFAD2','#CBE5C5',"#E83D36", "#62B8BF","#B0CFE0",
                  "#B6D88F","#FFDCAC","#DBDCDA","#F9DDE9", #Inh VIP
                  '#2FA147','#8466AF',"#EFE3CD", # Inh LAMP5
                  "#F8B6B4","#90D2D7","#B85E2C","#3581BA", # Inh SST
                  "#F5999D","#A27833","#FBBE6F","#BD96EA",# Inh PVALB
                  "#DFCCE4", '#FDD3BA',# Inh PAX6
                  "#FFDEDE","#FFB0BD","#D7517C","#991343", #Ast
                  "#F9CFA5","#f79c5e","#b55a00", # Oli
                  "#FFE2CF","#b38766","#75421f",# OPC
                  "#E1CAE5","#573B88", #Mic
                  "#7F0009","#6D0775", #Mac, T
                  "#E283B5","#F2C0D5","#6EC5A6","#A45627","#F38C64","#8BA0CC")

all$Subtype_f<-factor(all$ATAC.Subtype,levels = rev(Subtypes))

#######################
# 3. Cell fraction changes in each class
#######################
Major_count<-as.data.frame(table(all$Sample,all$ATAC.Class))

Major_count2<-Major_count %>% group_by(Var1) %>% mutate(per= prop.table(Freq) * 100)
colnames(Major_count2)<-c("Sample","Celltype","Count","Percentage")

Class_anno<-unique(meta[,c("Sample","Pathology","BrainRegion","Pathology")])

Major_count3<-merge(Major_count2,Class_anno,by="Sample")


Major_count3$Celltype_f<-factor(Major_count3$Celltype,levels = c("Exc","Inh",'Ast',"Oli","OPC",'Mic_Immune',"Vasc"))

Major_count3$Pathology_f<-factor(Major_count3$Pathology,levels = c('nonAD', 'earlyAD',"lateAD"))

Major_count4<-Major_count3
Major_count4$BrainRegion<-"All"

Major_count5<-rbind(Major_count4,Major_count3)
Major_count5$Region_f<-factor(Major_count5$BrainRegion,levels = c("All",'EC',"HC","TH","AG","MTC","PFC"))
Major_count5$log<-log(Major_count5$Percentage)


major_odds_p<-data.frame()


for (cell in c("Exc","Inh",'Ast',"Oli","OPC",'Mic_Immune',"Vasc_Epithelia") ){
  
  for (region in c("All",'EC',"HC","TH","AG","MTC","PFC")){
    tmp<-Major_count5[Major_count5$Celltype==cell & Major_count5$BrainRegion==region,]
    tmp$Percentage2<-tmp$Percentage/100
    
    fit_glm <- glm(Percentage2 ~ Pathology , 
                   data = tmp, 
                   family = "quasibinomial")
    
    summ<-summary(fit_glm)
    pvalues<-as.data.frame(summ$coefficients)$"Pr(>|t|)" [2]
    corrected_scores <- residuals(fit_glm, type = "response")
    
    coef_Pathology <- coef(fit_glm)["Pathology"]
    
    odds_ratio_Pathology <- exp(coef_Pathology)
    
    tmp_data<-data.frame(
      pathology=c('Pathology'),
      celltype=cell,
      region=region,
      Odds_ratio=odds_ratio_Pathology,
      pvalues=pvalues)
    
    major_odds_p<-rbind(major_odds_p,tmp_data)
    
  }
}

write.table(major_odds_p,"ATAC.Class.Pathology.OddsRatio.tsv",sep="\t",quote = F,row.names = F)


###################
# 4. Cell fraction changes in each subclass
###################
Subclass_count<-as.data.frame(table(all$Sample,all$ATAC.Subclass))

Subclass_count2<-Subclass_count %>% group_by(Var1) %>% mutate(per= prop.table(Freq) * 100)
colnames(Subclass_count2)<-c("Sample","Celltype","Count","Percentage")

Class_anno<-unique(meta[,c("Sample","Pathology","BrainRegion","Pathology")])

Subclass_count3<-merge(Subclass_count2,Class_anno,by="Sample")


Subclass_count3$Celltype_f<-factor(Subclass_count3$Celltype,levels = Subclasses)

Subclass_count3$Pathology_f<-factor(Subclass_count3$Pathology,levels = c('nonAD', 'earlyAD',"lateAD"))

Subclass_count4<-Subclass_count3
Subclass_count4$BrainRegion<-"Pathology"

Subclass_count5<-rbind(Subclass_count4,Subclass_count3)
Subclass_count5$Region_f<-factor(Subclass_count5$BrainRegion,levels = c("All",'EC',"HC","TH","AG","MTC","PFC"))
Subclass_count5$log<-log(Subclass_count5$Percentage)


Subclass_odds_p<-data.frame()


for (cell in Subclasses ){
  
  for (region in c("All",'EC',"HC","TH","AG","MTC","PFC")){
    print(cell)
    print(region)
    tmp<-Subclass_count5[Subclass_count5$Celltype==cell & Subclass_count5$BrainRegion==region,]
    tmp$Percentage2<-tmp$Percentage/100
    if (sum(tmp$Count>50)){
      
      fit_glm <- glm(Percentage2 ~ Pathology , 
                     data = tmp, 
                     family = "quasibinomial")
      
      summ<-summary(fit_glm)
      pvalues<-as.data.frame(summ$coefficients)$"Pr(>|t|)" [2]
      
      coef_Pathology <- coef(fit_glm)["Pathology"]
      
      odds_ratio_Pathology <- exp(coef_Pathology)
      
      tmp_data<-data.frame(
        pathology=c('Pathology'),
        celltype=cell,
        region=region,
        Odds_ratio=odds_ratio_Pathology,
        pvalues=pvalues)
      
      Subclass_odds_p<-rbind(Subclass_odds_p,tmp_data)
    }
  }
}

write.table(Subclass_odds_p,"ATAC.Subclass.Pathology.OddsRatio.tsv",sep="\t",quote = F,row.names = F)


###################
# 5. Cell fraction changes in each subtype
###################
Subtype_count<-as.data.frame(table(all$Sample,all$ATAC.Subtype))

Subtype_count2<-Subtype_count %>% group_by(Var1) %>% mutate(per= prop.table(Freq) * 100)
colnames(Subtype_count2)<-c("Sample","Celltype","Count","Percentage")

Class_anno<-unique(meta[,c("Sample","Pathology","BrainRegion","Pathology")])

Subtype_count3<-merge(Subtype_count2,Class_anno,by="Sample")


Subtype_count3$Celltype_f<-factor(Subtype_count3$Celltype,levels = Subtypes)

Subtype_count3$Pathology_f<-factor(Subtype_count3$Pathology,levels = c('nonAD', 'earlyAD',"lateAD"))

Subtype_count4<-Subtype_count3
Subtype_count4$BrainRegion<-"All"

Subtype_count5<-rbind(Subtype_count4,Subtype_count3)
Subtype_count5$Region_f<-factor(Subtype_count5$BrainRegion,levels = c("All",'EC',"HC","TH","AG","MTC","PFC"))
Subtype_count5$log<-log(Subtype_count5$Percentage)


Subtype_odds_p<-data.frame()


for (cell in Subtypes ){
  
  for (region in c("All",'EC',"HC","TH","AG","MTC","PFC")){
    print(cell)
    print(region)
    tmp<-Subtype_count5[Subtype_count5$Celltype==cell & Subtype_count5$BrainRegion==region,]
    tmp$Percentage2<-tmp$Percentage/100
    if (sum(tmp$Count>50)){
      
      fit_glm <- glm(Percentage2 ~ Pathology , 
                     data = tmp, 
                     family = "quasibinomial")
      
      summ<-summary(fit_glm)
      pvalues<-as.data.frame(summ$coefficients)$"Pr(>|t|)" [2]
      
      coef_Pathology <- coef(fit_glm)["Pathology"]
      
      odds_ratio_Pathology <- exp(coef_Pathology)
      
      tmp_data<-data.frame(
        pathology=c('Pathology'),
        celltype=cell,
        region=region,
        Odds_ratio=odds_ratio_Pathology,
        pvalues=pvalues)
      
      
      # Adding corrected scores to your data frame
      tmp_data$Corrected_Score <- corrected_scores
      
      # If you want to include the original and corrected scores in your output
      output_data <- data.frame(
        pathology = 'Pathology',
        celltype = cell,
        region = region,
        Odds_ratio = odds_ratio_Pathology,
        pvalues = pvalues,
        Corrected_Score = corrected_scores
      )
      
      Subtype_odds_p<-rbind(Subtype_odds_p,tmp_data)
    }
  }
}

write.table(Subtype_odds_p,"ATAC.Subtype.Pathology.OddsRatio.tsv",sep="\t",quote = F,row.names = F)

