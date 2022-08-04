#Tidy functional information from anvio
library(tidyverse)
library(vroom)

#exported functions from anvio
GL_fun <- vroom('GL_fun_v2.txt')

GL_fun <- GL_fun %>%
  mutate(gene=paste0("gene_", gene_callers_id))

GL_fun <- GL_fun[-c(1,3,5)]

#spreading out by functional source to get gene id x each function DB prediction
GL_fun_table <- GL_fun %>% pivot_wider(names_from = source, values_from = `function`)

GL_fun_table <- apply(GL_fun_table,2,as.character)

write.csv(GL_fun_table, "GL_function_per_gene.csv", row.names=FALSE) # this is on dropbox
