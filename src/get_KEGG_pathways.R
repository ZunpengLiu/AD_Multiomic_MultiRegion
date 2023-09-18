library(KEGGREST)
library(org.Hs.eg.db)
library(tidyverse)


hsa_path_eg  <- keggLink("pathway", "hsa") %>% 
  tibble(pathway = ., eg = sub("hsa:", "", names(.)))


hsa_kegg_anno <- hsa_path_eg %>%
  mutate(
    symbol = mapIds(org.Hs.eg.db, eg, "SYMBOL", "ENTREZID"),
    ensembl = mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID")
  )


hsa_kegg_anno$pathway=str_replace_all(hsa_kegg_anno$pathway,"path:","")
write.table(hsa_kegg_anno,"hsa_KEGG_anno.All.Pathways.tsv",sep="\t",quote = F,row.names = F)

hsa05032=hsa_kegg_anno[hsa_kegg_anno$pathway=="hsa05032",]
write.table(hsa05032,"hsa_KEGG_anno.All.Pathways_hsa05032_Morphine_addiction.tsv",sep="\t",quote = F,row.names = F)

hsa05030=hsa_kegg_anno[hsa_kegg_anno$pathway=="hsa05030",]
write.table(hsa05030,"hsa_KEGG_anno.All.Pathways_hsa05030_Cocaine_addiction.tsv",sep="\t",quote = F,row.names = F)

Amphetamine
hsa05031=hsa_kegg_anno[hsa_kegg_anno$pathway=="hsa05031",]
write.table(hsa05031,"hsa_KEGG_anno.All.Pathways_hsa05031_Amphetamine_addiction.tsv",sep="\t",quote = F,row.names = F)

Nicotine
hsa05033=hsa_kegg_anno[hsa_kegg_anno$pathway=="hsa05033",]
write.table(hsa05033,"hsa_KEGG_anno.All.Pathways_hsa05033_Nicotine_addiction.tsv",sep="\t",quote = F,row.names = F)

Alcoholism
hsa05034=hsa_kegg_anno[hsa_kegg_anno$pathway=="hsa05034",]
write.table(hsa05034,"hsa_KEGG_anno.All.Pathways_hsa05033_Alcoholism.tsv",sep="\t",quote = F,row.names = F)

