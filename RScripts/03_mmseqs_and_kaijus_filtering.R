#taking mmseqs results in & cleaning up

library(tidyverse)
library(vroom)

#ran mmseqs2 easy-search and easy-taxonomy against NR protein database
#load in mmseqs2 results

#get taxonomy last common ancestory assignments
taxonomy.lca.all <- vroom("data/mmseqs_res/taxonomyResult_nr_GL_filtset3x5perc_lca.tsv", col_names=c("gene", "databaseID", "rank", "TaxaLCA", "FullLCA", "FullTopHit"), delim = "\t")

taxonomy.lca.all <- taxonomy.lca.all %>%
  mutate(FullLCA = gsub("unknown;", "", FullLCA)) %>%
  mutate(FullTopHit = gsub("-_cellular organisms;", "", FullTopHit)) %>%
  mutate(FullTopHit = gsub("-_cellular organisms", "unclassified", FullTopHit)) %>%
  replace(is.na(.), "unclassified") 


taxonomy.lca.all <- taxonomy.lca.all %>%
  mutate(Domain_Mmseqs = ifelse(str_detect(FullTopHit, "d_Bacteria"), "Bacteria", ifelse(str_detect(FullTopHit, "d_Archaea"), "Archaea", "Other"))) 


#write.table(taxonomy.lca.all, "results/gl_tax_filtset_genes_to_remove.tsv", quote=FALSE, sep = '\t', row.names = FALSE)

#import kaiju data
kaiju_res <- readRDS("../../../Box Sync/GrassLandMicrobes/GL_MetaG_Anvio/Kaiju/Sourcetracker_Meta/GL_kaiju_tax_ready_for_source_09132021.RDS")

#join the two taxonomy approaches and get combined results per gene
taxonomy.lca.all.join <- left_join(taxonomy.lca.all, kaiju_res, by=c("gene" = "ReadName")) 

taxonomy.lca.all.join <- taxonomy.lca.all.join %>%
  mutate(Joined_Domain = ifelse(Domain_Mmseqs == "Other" & Domain == "Unclassified", "Remove", ifelse(Domain == "Eukaryota" | Domain == "Viruses", "Remove" ,"Keep") )) %>%
  mutate(Joined_Domain = ifelse(str_detect(FullLCA, "Viridiplantae") | str_detect(FullLCA, "Eukaryota") | str_detect(FullTopHit, "Eukaryota") | str_detect(FullTopHit, "Viruses"), "Remove", as.character(Joined_Domain))) %>%
  mutate(Joined_Domain = ifelse(databaseID %in% c(0, 1, 131567, 2583588, 90371), "Remove", as.character(Joined_Domain)))
  

write.table(taxonomy.lca.all.join, "results/gl_tax_filtset_genes_to_remove_121422.tsv", quote=FALSE, sep = '\t', row.names = FALSE)

genes_to_remove <- taxonomy.lca.all.join$gene[taxonomy.lca.all.join$Joined_Domain == "Remove"]


