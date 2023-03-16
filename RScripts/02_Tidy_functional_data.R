#Tidy functional information from anvio
library(tidyverse)
library(vroom)

#exported functions from anvio
GL_fun <- vroom('GL_fun_v2.txt')

GL_fun <- GL_fun %>%
  filter(gene_callers_id < 8829098) %>%
  mutate(gene=paste0("gene_", gene_callers_id))

GL_fun <- GL_fun[-c(1,3,5)]

#spreading out by functional source to get gene id x each function DB prediction
GL_fun_table <- GL_fun %>% pivot_wider(names_from = source, values_from = `function`)

GL_fun_table <- apply(GL_fun_table,2,as.character)

write.csv(GL_fun_table, "GL_function_per_gene_update.csv", row.names=FALSE)

GL_fun_table <- vroom("GL_function_per_gene_update.csv")

length(unique(GL_fun_table$gene))
#5996849 total genes

#get filtered set of genes we want to extract functions for
norm_counts <- readRDS("decontam_norm_counts_w-controls-12082022.RDS")
norm_counts <- as.data.frame(norm_counts)
norm_counts$ReadName <- as.character(row.names(norm_counts))

norm_counts <- norm_counts %>%
  filter(ReadName < 8829098)

norm_counts_2 <- as.data.frame(norm_counts$ReadName)

#remove ribosomal genes and save
write_delim(norm_counts_2, "filt_gene_names.tsv")


row.names(norm_counts) <- paste0("gene_", norm_counts$ReadName)

filt_genes <- row.names(norm_counts)

length(unique(filt_genes))
#128022 genes

#filter functions and then assess annotation status
GL_fun_table_filt <- GL_fun_table[GL_fun_table$gene %in% filt_genes,]

GL_fun_table_filt.v2 <- GL_fun_table_filt %>%
  mutate(Annot = ifelse(KeggGhostKoala == "NULL", ifelse(Pfam == "NULL", ifelse(EGGNOG_BACT == "NULL", ifelse(COG_FUNCTION == "NULL", ifelse(COG_CATEGORY == "NULL", ifelse(GO_TERMS == "NULL",  ifelse(KEGG_PATHWAYS == "NULL", ifelse(BiGG_Reactions == "NULL", ifelse(COG20_FUNCTION == "NULL", ifelse(COG20_CATEGORY == "NULL", ifelse(COG20_PATHWAY == "NULL", "Unannotated", "Annotated"),"Annotated"), "Annotated"), "Annotated"), "Annotated"), "Annotated"), "Annotated"), "Annotated"), "Annotated"), "Annotated"), "Annotated"))

summary(as.factor(GL_fun_table_filt.v2$Annot))

#import taxonomy for some additional filtering steps
taxonomy.lca.all.join <- read.delim("/Users/Cassie/Desktop/Grass_cluster/mmseq2_results/results/gl_tax_filtset_genes_to_remove_121422.tsv")

genes_to_remove <- paste0("gene_", taxonomy.lca.all.join$gene[taxonomy.lca.all.join$Joined_Domain == "Remove"])

genes_to_remove
#10860

filt_genes_2 <- filt_genes[!filt_genes %in% genes_to_remove]
#117331

GL_fun_table_filt <- GL_fun_table_filt[!GL_fun_table_filt$gene %in% genes_to_remove,]
#110553

#also filter to only keep genes with COG functions 
GL_fun_table_filt.noCOG <- GL_fun_table_filt %>% 
  filter(!COG20_CATEGORY == "NULL") 
#104065 

#manual filtering out of some problematic host plant genes
GL_fun_table_filt.alt <- GL_fun_table_filt %>%
  filter(!str_detect(Pfam, "Retrotransposon gag protein"))%>%
  filter(!(str_detect(Pfam, "Reverse transcriptase") & COG20_CATEGORY == "NULL"))

GL_fun_table_filt
write.csv(GL_fun_table_filt.alt, "gl_fun_only_keep_alterate.csv", row.names=FALSE)


#% genes with annotations
(1 - ( length(unique(filt_genes_2)) - length(GL_fun_table_filt$gene) ) / length(unique(filt_genes_2)))*100


write.table(GL_fun_table_filt.noCOG, "gl_fun_only_keep_121422.tsv", quote=FALSE, sep = '\t', row.names = FALSE)
write.csv(GL_fun_table_filt.noCOG, "gl_fun_only_keep_121422.csv", row.names=FALSE)




