#############
# Description: This script is used to generate the color scheme for the cell types

######################## 1. Group names and colors ########################
# Class - Major cell types
Class <- c("Exc", "Inh", "Ast", "Oli", "Opc","Mic_Immune","Vasc")


Class_colors<-c()

major_celltype_markergenes<-c("SLC17A7","GAD2","GFAP","MOBP","VCAN","P2RY12")

major_celltype_cols<-c("1-Exc" ="#34A047",
                       "2-Inh"="#74C476",
                       "3-Ast" = "#E11E26",
                       "4-Oli" ="#F47E1F",
                       "5-Opc" ="#B15928",
                       "6-Mic_Immune" = "#6A3E98",
                       "7-Vasc" = "#6A3E98")




# SubClass - Minor cell types



Subclass,levels=c("IT L2/3-1","IT L2/3-2","IT L3/4","IT L4/5-1","IT L4/5-2","IT L5/6","L5/6 NP","L6 CT","HCERC","ERC", # Exc
                  "VIP","SST-1","SST-2","LAMP5-1","LAMP5-2","PVALB-1","PVALB-2","MEIS2", # Inh
                  "Ast-1","Ast-2","Ast-3",
                  "Oli-1","Oli-2","Oli-3",
                  "OPC-1","OPC-2",
                  "Mic-1","Mic-2","Mic-3","Mic-4",
                 # Vasc 
                 )


Subclass_col= c("IT L2/3-1"="#3E6AB7",
"IT L2/3-2"= "#8496DF",
"IT L3/4"= "#C8CFED",
"IT L4/5-1"= "#AEEDE1",
"IT L4/5-2"= "#B3F2F2",
"IT L5/6"= "#71D1D3",
"L5/6 NP"= "#0CACC0",
"L6 CT"= "#00838F",
"HCERC"= "#006065",
"ERC"= "#324737",
"VIP"="#195E20",
"SST-1"="#23832C",
"SST-2"="#61AA7C",
"LAMP5-1"="#1ABC02",
"LAMP5-2"="#00CE81",
"PVALB-1"="#089992",
"PVALB-2"="#64D8CD",
"MEIS2-1"="#BAD8AB",
"MEIS2-2"="#99E099",
"Ast-1"= "#FFB2B2",
"Ast-2"= "#D7517C",
"Ast-3"= "#991343",
"Oli-1"="#C96145",
"Oli-2"="#EB8A47",
"Oli-3"="#FBC57F",
"OPC-1"= "#7B4F37",
"OPC-2"= "#D5A370",
"Mic-1"= "#CFA2D2",
"Mic-2"="#834BA0",
"Mic-3"="#573B88",
"Mic-4"= "#BA49B2")

subclass_colors<-c(
  '#3E6AB7', '#8496DF', '#C8CFED', '#AEEDE1', '#B3F2F2','#71D1D3','#0CACC0','#00838F','#006065','#324737', #Exc
  '#3E6AB7', '#63C66F', '#43D1C3', '#C93667', '#E69CAC','#DB874D','#F4CFB2','#D3D5ED','#6D79A0',# Inh
  '#FFB2B2', '#D7517C', '#991343', #Ast
  '#C96145','#EB8A47', '#FBC57F', #Oli
  '#7B4F37','#D5A370'#OPC
)



Inh_markers1<-c("PAX6","HMBOX1","CDH4","IL1RAPL2","DDR2",  # PAX6
                "LAMP5","FGF13","CACNA2D1","FBXL7","EYA4", # LAMP5
                "SST","KIF26B","EPB41L2","BCL11A","PDE1A", # SST
                "PVALB","ADCY1","DPP10","TMEM132C","CNTNAP3", # PVALB
                "VIP","CIT","HIPK2","ADAMTS6","CALB2",# VIP
                "MEIS2","FOXP1","CHRM2","GRM3","FOXP2") # VIP #MR_ATAC sub cell type


# Subtype



# Brain regions

brain_regions <- c(,,,,,)


brain_region_cols=c("1-AG"="#8DD3C7",
                    "2-MT"="#80B1D3",
                    "3-PFC"="#FDB462",
                    "4-EC"="#FFED6F",
                    "5-HC"="#BEBADA",
                    "6-TH"="#B3DE69")


marker_genes<-c("SLC17A7" ,"CUX2","RORB","NXPH2","VWA2","TOX3",
"GAD1","VIP","SST","PVALB","MEIS2",
"GFAP","MAG","PDGFRA",
"CX3CR1",
"P2RY12","MKI67")

major_celltype_markergenes<-c("SLC17A7","GAD2","GFAP","MOBP","VCAN","P2RY12")


######################## 2. Marker genes ########################

marker_genes<-c("SLC17A7" ,"CUX2","RORB","NXPH2","VWA2","TOX3",
"GAD1","VIP","SST","PVALB","MEIS2",
"GFAP","MAG","PDGFRA",
"CX3CR1",
"P2RY12","MKI67")

markerGenes1 <-c('AQP4','GFAP','ETNPPL','GLI3','SPI1',"AQP4", "APOE", "AQP1", "GRIA1", "PLCG1" # Ast
                 'SYT1','SNAP25','UCHL1','GABRA1', #Neurons
                 'SLC17A7',"SNAP25", "SATB2", "CUX2", "TSHZ2", "RORB", "FOXP2", "TLE4", # Exc
                 'GAD1','GAD2', #In
                 "VIP", "LAMP5","PVALB", "SST",'PAX6', "CPLX3", "CXCL14", "NDNF", "RELN", "LHX6", "NPY", "SOX6","ADARB2", # Inh
                 'MOBP','MOG','PLP1',"CNP", "MBP", "MAG", #Oligodendrocytes
                 'VCAN','PDGFRA','STK32A', #OPC
                 'P2RY12','SPI1','C3','C1QB','PTPRC','APBB1IP', #Microglia
                 'CAVIN2' #Vas
                 "VWF", "PECAM1", "CLDN5",'FLT1' # Endo
                 'PDGFRB',"ACTA2","MYH11","TAGLN", # Pericytes
                 "VWF","PECAM1","CLDN5", # Endo
                 "CEMIP"  #Fib
)

 

motifs <- c("NEUROD2","NEUROD4","NEUROD6","NEUROG1","NEUROG3","SPIB", "SPI1","ATOH7","ATOH8","BACH1","BACH2","BHLHA15","BHLHE22","BHLHE23","EPAS1","JUNB", "JUN","JUND", "FOS","FOSB","ID3","ID4","MEF2A","MEF2C","MEF2D","MESP1","MESP2","NFE2","OLIG1","OLIG2","OLIG3","SNAI1","TCF4","TCF23", # Neu
"NFIA","NFIB","NFIC","NFIX", #Ast
"CTCF","SOX4","SOX9","SOX10","SOX11","SOX12","SOX13","SOX17","SOX30", # Oli
"SPI1","SPIB","SPIC","BCL11A","BCL11B","CEBPB","EHF","ELF1","ELF2","ELF3","ELF4","ELF5","ETS2","ETV6","ETV7","IRF1","IRF2","IRF3","IRF4","STAT2",#MIC)
)



############### Ast
Ast_markerGene <- c("C4B","GFAP","DCLK1", "C3","SERPINA3","VIM", # 1.the most activated
                   "AQP4", #
                   "GRM3",
                   "HIF3A", "LSAMP", "GPC5","CXCL10","IL18",
                   "SLC1A2","DPP10","LUZP2","TRPM3","BRINP3") # xiaowei zhuang lab Ast activation score


Ast_markerGenes_2 <- c("C4B", "C3","SERPINA3", "CXCL10","GFAP","VIM","IL18","HIF3A","LSAMP","GPC5","SLC1A2","LUZP2","TRPM3","BRINP3","GFAP") # xiaowei zhuang lab Ast activation score

Ast_markerGenes_3 <- c("GRM3",
"AQP4",
"SLC1A2",
"SLC1A3",
"GFAP",
"CADM1",
'CACNB2',
'PDE3B',
'DPP10',
'SEMA6A',
"NHS",
"GRIA1",
"DCLK1",
"ST6GAL1",
"GRIN2A",
"CPAMD8",
"NRXN3",
"KCND2",
"SLC6A11",
'RALYL',
"LUZP2",
"CHI3L1",
"ARHGEF3")



####### Mic cell activation
Mic_markerGenes  <- c("TPT1", "P2RY12", "DUSP1", "MKI67",
"B2M", "TREM2", "CCL2", 'APOE', "AXL", "ITGAX", "CD9", "C1QA", "C1QC", "CTSS", # Xiaowei zhuang Cell activation score
"DPP10","MEG3","CHRM3","TRIM30A","IGTB2", # Xiaowei zhuang Cell marker genes
"CCL3","CSF3R","CX3CR1")




#### Oligodendrocyte
Oli_markerGenes_1  <- c("RASGRF1","OPALIN")
Oli_markerGenes_2  <- c("C4B","IL33","IL18") # Xiaowei zhuang lab for estimating the Oli activation score
Oli_markerGenes_3  <- c("ANK2","TRIM2","ROBO1","DGKI","SPOCK1", # OLIG1 to OLIG-3 decreased levels
"NEAT1",  #OLIG2 OLIG3 specific/ high
"IL33", # OLIG2 OLIG3 specific
"C4B", # OLIG3 specific
"IL18") # Xiaowei zhuang lab for estimating the Oli activation score shown in main figures




