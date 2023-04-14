### Script to analyze Metagenomic soil samples from McLaughlin ###

# Files needed (all are found within the R-Data-Files folder on box)
## 1. Kallisto file of raw gene counts
## 2. Functional pathways file
## 3. Joined MMSeq and Kaiju taxonomy file
## 4. Metadata file

rm(list=ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("microbiome", version = "3.14")

# Load libraries ####
library(DESeq2)
library(plyr)
library(dplyr)
library(tidyverse)
library(phyloseq)
library(decontam)
library(vegan)
library(Rmisc)
library(biobroom)
library(gridExtra)
library(microbiome)
library(ggordiplots)
library(ecole) #pairwise adonis2
library(pairwiseAdonis)
library(coin)
library(cowplot)
library(betapart)
library(ggpubr)
library(geosphere)
library(scales)

# Read in data ####
#count.gene <- readRDS("Data/Final-Data/txi.kallisto.08042021.RDS")
fun.gene <- read.csv("Data/Final-Data/gl_fun_only_keep_121422.csv")
tax.gene <- read.delim("Data/Final-Data/gl_tax_filtset_genes_to_remove_121422.tsv")
metadata <- read.csv("Data/Final-Data/McL_metag_metadata.csv")
latlong <- read.csv("Data/Lat-Long.csv")

# Prepare data ####

## Gene count prep ####
# colnames(count.gene$counts) <- sub("_S.*", "", colnames(count.gene$counts))
# colnames(count.gene$counts) <- sub("-", "_", colnames(count.gene$counts))
# colnames(count.gene$abundance) <- sub("_S.*", "", colnames(count.gene$abundance))
# colnames(count.gene$abundance) <- sub("-", "_", colnames(count.gene$abundance))

## Taxa prep ####
# Domain = Kaiju
tax.gene <- filter(tax.gene, Joined_Domain != "Remove")
row.names(tax.gene) <- tax.gene$gene
tax.gene <- tax.gene[order(tax.gene$gene),]
tax.gene <- tax.gene[,-1]
tax.gene <- tax.gene[,c(7:13)]

## Function prep ####
fun.gene$gene <- gsub("^.*\\_", "", fun.gene$gene)
row.names(fun.gene) <- fun.gene$gene
fun.gene <- fun.gene[order(fun.gene$gene),]
fun.gene <- fun.gene[,c(11,12,10)] # change order
fun.gene$COG20_CATEGORY <- gsub("!.*$", "", fun.gene$COG20_CATEGORY)
fun.gene$COG20_FUNCTION <- gsub("!.*$", "", fun.gene$COG20_FUNCTION)
fun.gene$COG20_PATHWAY <- gsub("!.*$", "", fun.gene$COG20_PATHWAY)

fun.gene <- fun.gene[rownames(fun.gene) %in% rownames(tax.gene),]
tax.gene <- tax.gene[rownames(tax.gene) %in% rownames(fun.gene),]
fun.gene <- as.matrix(fun.gene)
tax.gene <- as.matrix(tax.gene)

## Metadata prep ####
metadata <- merge(metadata, latlong, by = "Plot", all.x = T, all.y = F)

# removes controls we didn't sequence: GLM_0698, GLM_0696, GLM_0694, GLM_0693, GLM_0392, GLM_0391, GLM_0390
to.rm <- c("GLM_0698", "GLM_0696", "GLM_0694", "GLM_0693", "GLM_0392", "GLM_0391", "GLM_0390")
metadata <- metadata[!(metadata$SampleID %in% to.rm),] # sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID
metadata$plant.water <- paste(metadata$PlantTrt, metadata$WaterTrt, sep = ".")

metadata$PlotSubTreatmentName <- recode_factor(metadata$PlotSubTreatmentName, ST_forbs_X_grasses = "native mix 1 + invasives", SA_forbs_X_grasses = "native mix 2 + invasives", Stress_avoiding_forbs = "native mix 2", Stress_tolerant_forbs = "native mix 1")

#colnames(count.gene[[1]]) == metadata[,"SampleID"] # make sure everything is in the right order

controls <- c("GLM_0697", "GLM_0695", "GLM_0393")
metadata.nocontrols <- metadata[!(metadata$SampleID %in% controls),]

metadata.nocontrols$PlantTrt <- recode_factor(metadata.nocontrols$PlantTrt, forb = "natives",  grasses = "invasives")

metadata.nocontrols$PlantTrt <- factor(metadata.nocontrols$PlantTrt, levels = c("natives", "invasives", "mix"))

metadata.nocontrols$WaterTrt <- recode_factor(metadata.nocontrols$WaterTrt, Watering = "watered", Drought = "drought", Control = "control")

metadata.nocontrols$WaterTrt <- factor(metadata.nocontrols$WaterTrt, levels = c("drought", "control", "watered"))

metadata.nocontrols$PlotSubTreatmentName <- factor(metadata.nocontrols$PlotSubTreatmentName, levels = c("native mix 1", "native mix 2", "Invasive_grasses", "native mix 1 + invasives", "native mix 2 + invasives"))

# create table for supplement
# samp_data <- as_tibble(metadata.nocontrols) %>%
#   mutate(AllPlantSpeciesInCore = ifelse(str_detect(AllPlantSpeciesInCore, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInCore)),
#          Bromus = ifelse(str_detect(AllPlantSpeciesInCore, "Bromus"), 1, 0),
#          Avena = ifelse(str_detect(AllPlantSpeciesInCore, "Avena"), 1, 0),
#          Calycadenia = ifelse(str_detect(AllPlantSpeciesInCore, "Calycadenia"), 1, 0),
#          Hemizonia = ifelse(str_detect(AllPlantSpeciesInCore, "Hemizonia"), 1, 0),
#          Elymus = ifelse(str_detect(AllPlantSpeciesInCore, "Elymus"), 1, 0),
#          Plantago = ifelse(str_detect(AllPlantSpeciesInCore, "Plantago"), 1, 0),
#          Clarkia = ifelse(str_detect(AllPlantSpeciesInCore, "Clarkia"), 1, 0),
#          Agoseris = ifelse(str_detect(AllPlantSpeciesInCore, "Agoseris"), 1, 0),
#          Lasethenia = ifelse(str_detect(AllPlantSpeciesInCore, "Lasethenia"), 1, 0),
#          'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInCore, "LOMU"), 1, 0),
#          'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInCore, "VUMI"), 1, 0)) %>%
#   select(SampleID, Bromus, Avena, Calycadenia, Hemizonia, Elymus,Plantago,Clarkia,Agoseris,Lasethenia,'Festuca perennis','Festuca microstachys') %>%
#   arrange(SampleID)
# 
# write.csv(metadata.nocontrols, "Manuscripts/supplementarytable.csv", row.names = F)

## Filter data ####

# ddsTxi <- DESeqDataSetFromTximport(count.gene,
#                                    colData = metadata,
#                                    design = ~ plant.water)
# 
# saveRDS(ddsTxi, "GL_ddsTxi.deseq.obj.12082022.RDS")

# Remove rows with no counts and genes not seen more than 3 times in at least 5% of the samples. This protects against genes with small mean & trivially large C.V

# ddsTxi <- ddsTxi[rowSums(counts(ddsTxi)) > 1, ] # remove genes with no counts
# ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) > 3) >= (0.05*ncol(counts(ddsTxi))), ] 
# 
# 
# # saveRDS(ddsTxi, "GL_ddsTxi.deseq.obj.filt.3x.5perc.12082022.RDS")
# ddsTxi <- readRDS("Data/Final-Data/GL_ddsTxi.deseq.obj.filt.3x.5perc.12082022.RDS")

## Decontam Counts ####
# contam_counts <- counts(ddsTxi, normalized = F) # object for decontam
# 
# ps <- phyloseq(otu_table(contam_counts, taxa_are_rows = T), sample_data(metadata))
# 
# #first tell it which samples are the NC
# sample_data(ps)$is.neg <- sample_data(ps)$SampleType == "Kit"
# 
# #all sequences that are more prevalent in negative controls than in positive samples
# contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)
# table(contamdf.prev05$contaminant) # 143 contaminants 
# 
# #Make phyloseq object of presence-absence in negative controls
# ps.neg <- prune_samples(sample_data(ps)$is.neg == "TRUE", ps)
# ps.neg.presence <- transform_sample_counts(ps.neg, function(abund) 1*(abund>0))
# 
# #Make phyloseq object of presence-absence in true positive samples
# ps.pos <- prune_samples(sample_data(ps)$is.neg == "FALSE", ps)
# ps.pos.presence <- transform_sample_counts(ps.pos, function(abund) 1*(abund>0))
# 
# #Make data.frame of prevalence in positive and negative samples
# #using prev threshold = 0.5
# df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence),
#                       prevalence.neg = taxa_sums(ps.neg.presence),
#                       contam.prev = contamdf.prev05$contaminant)
# 
# ggplot(df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contam.prev)) +
#   geom_point()
# 
# contam.rows <- which(contamdf.prev05$contaminant)
# contaminants <- rownames(contamdf.prev05[contam.rows,])
# 
# #returns list of all taxa, except contaminants
# decontam_counts <- contam_counts[!(rownames(contam_counts) %in% contaminants),]
# 
# decontam_counts_nocontrols <- decontam_counts[, !(colnames(decontam_counts) %in% controls)] # get rid of controls
# 
# saveRDS(decontam_counts_nocontrols, "Data/Final-Data/decontam_non-norm_counts_nocontrols.RDS")
# 
# rm(ps, ps.neg, ps.neg.presence, ps.pos, ps.pos.presence, ddsTxi)

## Sourcetracker Prep ####
# dds <- DESeqDataSetFromMatrix(decontam_counts,
#                               colData = metadata,
#                               design = ~ plant.water)
# 
# 
# dds <- estimateSizeFactors(dds)
# 
# normalized_counts <- counts(dds, normalized = T) 
# 
# saveRDS(normalized_counts, "Data/Final-Data/decontam_norm_counts_w-controls.RDS") # output for Source Tracker

## Filter Counts ####
decontam_counts_nocontrols <- readRDS("Data/Final-Data/decontam_non-norm_counts_nocontrols.RDS") # not normalized

# filter out unwanted genes (ribosomes, plant associated etc)
decontam_counts_nocontrols <- decontam_counts_nocontrols[row.names(decontam_counts_nocontrols) < 8829098,]

decontam_counts_nocontrols <- decontam_counts_nocontrols[rownames(decontam_counts_nocontrols) %in% rownames(tax.gene),]

decontam_counts_nocontrols <- decontam_counts_nocontrols[rownames(decontam_counts_nocontrols) %in% rownames(fun.gene),]


## Set up contrasts ####
contrasts <- c("FW.FC", 
               "FD.FC", 
               "MD.FC", 
               "MW.FC", 
               "GD.GC",
               "GW.GC",
               "MD.GC",
               "MW.GC",
               "FC.GC",
               "FD.GD",
               "FW.GW",
               "MD.FD",
               "MW.FW",
               "MC.FC",
               "MC.GC")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(mix.Drought/forbs.Control) so if something is higher in mix drought, the number is greater than one, so the log fold change > 0, if something is higher in the control treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FW.FC = c("plant.water", "forb.Watering", "forb.Control"),
                      FD.FC = c("plant.water","forb.Drought","forb.Control"),
                      MD.FC = c("plant.water","mix.Drought","forb.Control"),
                      MW.FC = c("plant.water","mix.Watering","forb.Control"),
                      GD.GC = c("plant.water","grasses.Drought","grasses.Control"),
                      GW.GC = c("plant.water","grasses.Watering","grasses.Control"),
                      MD.GC = c("plant.water","mix.Drought","grasses.Control"),
                      MW.GC = c("plant.water","mix.Watering","grasses.Control"),
                      FC.GC = c("plant.water","forb.Control","grasses.Control"),
                      FD.GD = c("plant.water","forb.Drought","grasses.Drought"),
                      FW.GW = c("plant.water","forb.Watering","grasses.Watering"),
                      MD.FD = c("plant.water","mix.Drought","forb.Drought"),
                      MW.FW = c("plant.water", "mix.Watering", "forb.Watering"),
                      MC.FC = c("plant.water", "mix.Control", "forb.Control"), 
                      MC.GC = c("plant.water", "mix.Control", "grasses.Watering"))

plot.name.list <- list(FW.FC = "Native Watered v. Native Control",
                       FD.FC = "Native Drought v. Native Control",
                       MD.FC = "Mix Drought v. Native Control",
                       MW.FC = "Mix Watered v. Native Control",
                       GD.GC = "Invasive Drought v. Invasive Control",
                       GW.GC = "Invasive Watered v. Invasive Control",
                       MD.GC = "Mix Drought v. Invasive Control",
                       MW.GC  = "Mix Watered v. Invasive Control",
                       FC.GC = "Native Control v. Invasive Control",
                       FD.GD = "Native Drought v. Invasive Drought",
                       FW.GW = "Native Watered v. Invasive Watered",
                       MD.FD = "Mix Drought v. Native Drought",
                       MW.FW = "Mix Watered v. Native Watered",
                       MC.FC = "Mix Control v. Native Control",
                       MC.GC = "Mix Control v. Invasive Control")

# DESeq2 functions ####
## Run DESeq2 ####
ps.fun <- phyloseq(otu_table(decontam_counts_nocontrols, 
                         taxa_are_rows = T), 
               sample_data(metadata.nocontrols),
               tax_table(fun.gene))

ps.fun <- tax_glom(ps.fun, taxrank = "COG20_FUNCTION")

treat.fun <- phyloseq_to_deseq2(ps.fun, ~ plant.water) 

dds.pw.fun <- DESeq(treat.fun)

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw.fun, contrast = contrast.list[[i]], pAdjustMethod = "BH")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #tidy results 
  res.list[[i]] <- tidy(res.alpha)

 res.list[[i]]$gene <- factor(res.list[[i]]$gene, levels = res.list[[i]]$gene[order(res.list[[i]]$estimate)])
  
 if(nrow(res.list[[i]]) > 0) {
   p <- ggplot(res.list[[i]], aes(x = gene, y = estimate)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.title.x = element_blank(),
        ) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")))
  
   plot.list[[i]] = p
 }
   sig.genes <- c(as.numeric(as.character(res.list[[i]]$gene)), sig.genes)
}

grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], ncol = 2)

# results: 
length(unique(sig.genes)) #44

## Plot results: Function ####
test <- as.data.frame(cbind(fun.gene, gene = row.names(fun.gene)))

test <- test %>% mutate(COG20_FUN_clean = str_replace(COG20_FUNCTION, "!!!.*", "")) %>% 
  mutate(COG20_FUN_clean = str_replace(COG20_FUN_clean, " \\(PDB.*", "")) %>%
  mutate(COG20_FUN_clean = str_replace(COG20_FUN_clean, " \\(PUB.*", "")) %>%
  mutate(COG20_FUN_clean = ifelse(COG20_FUN_clean == "Uncharacterized conserved protein, contains Zn-finger domain of CDGSH type", "Uncharacterized Fe-S cluster protein YjdI (YjdI)", COG20_FUN_clean))

for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], test, by = "gene", all.y = F)
}

sig.dif.fun <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.fun[[i]] <- ddply(res.list[[i]], .(COG20_FUN_clean, COG20_CATEGORY), summarize, estimate = sum(estimate))
  }
}

plot.list2 <- list()

for(i in names(sig.dif.fun)){
  
  sig.dif.fun[[i]]$COG20_FUN_clean <- factor(sig.dif.fun[[i]]$COG20_FUN_clean, levels = sig.dif.fun[[i]]$COG20_FUN_clean[order(sig.dif.fun[[i]]$estimate)])
  
  p <- ggplot(sig.dif.fun[[i]], aes(x = COG20_FUN_clean, y = estimate, col = COG20_CATEGORY)) + 
      geom_point(size = 3) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)
          ) +
      ylim(-1.3,1.7) +
      labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change"))) +
      scale_x_discrete(labels = c("2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase MenH and related esterases, alpha/beta hydrolase fold (MenH)" = "2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase MenH and related\n esterases, alpha/beta hydrolase fold (MenH)")) +
      coord_flip() +
      scale_colour_manual(values = c(
        "Mobilome: prophages, transposons" = "dodgerblue2",                     
        "Cell wall/membrane/envelope biogenesis" = "#E31A1C",                      
        "Signal transduction mechanisms" = "green4",                               
        "Transcription" = "#6A3D9A",                                                
        "Intracellular trafficking, secretion, and vesicular transport" = "#FF7F00",
        "Carbohydrate transport and metabolism" = "black",                        
        "Energy production and conversion" = "gold1",                              
        "Inorganic ion transport and metabolism" = "skyblue2",                        
        "General function prediction only" = "palegreen2",                            
        "Coenzyme transport and metabolism" =  "#FDBF6F",                            
        "Function unknown" = "grey37",                                             
        "Defense mechanisms" = "#CAB2D6",                                           
        "Replication, recombination and repair" = "orchid1",                        
        "Amino acid transport and metabolism" = "darkturquoise",                      
        "Translation, ribosomal structure and biogenesis" = "darkorange4",            
        "Posttranslational modification, protein turnover, chaperones" = "yellow4" 
      )) 
   plot.list2[[i]] = p
}

legend <- get_legend(ggplot(sig.dif.fun[["MD.FC"]], aes(x = COG20_FUN_clean, y = estimate, col = COG20_CATEGORY)) + 
      geom_point(size = 3) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)
          ) +
      ylim(-1.3,1.7) +
      labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")), col = "COG20 Category") +
      scale_x_discrete(labels = c("2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase MenH and related esterases, alpha/beta hydrolase fold (MenH)" = "2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase MenH and related\n esterases, alpha/beta hydrolase fold (MenH)")) +
      coord_flip() +
      scale_colour_manual(values = c(
        "Mobilome: prophages, transposons" = "dodgerblue2",                     
        "Cell wall/membrane/envelope biogenesis" = "#E31A1C",                      
        "Signal transduction mechanisms" = "green4",                               
        "Transcription" = "#6A3D9A",                                                
        "Intracellular trafficking, secretion, and vesicular transport" = "#FF7F00",
        "Carbohydrate transport and metabolism" = "black",                        
        "Energy production and conversion" = "gold1",                              
        "Inorganic ion transport and metabolism" = "skyblue2",                        
        "General function prediction only" = "palegreen2",                            
        "Coenzyme transport and metabolism" =  "#FDBF6F",                            
        "Function unknown" = "grey37",                                             
        "Defense mechanisms" = "#CAB2D6",                                           
        "Replication, recombination and repair" = "orchid1",                        
        "Amino acid transport and metabolism" = "darkturquoise",                      
        "Translation, ribosomal structure and biogenesis" = "darkorange4",            
        "Posttranslational modification, protein turnover, chaperones" = "yellow4" 
      )) )

plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], plot.list2[[4]], align = "v")

fun.fig <- ggarrange(ggarrange(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], plot.list2[[4]], ncol = 1, labels = c("(a)", "(b)", "(c)", "(d)"), heights = c(1, 0.2, 0.15,0.17)), legend, ncol = 2, widths = c(1,0.5))

ggsave("Figures/deseq-fun2.jpeg", fun.fig, height = 18, width = 15, units = "in", dpi = 600)


## Plot results: Category ####

sig.dif.fun <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.fun[[i]] <- ddply(res.list[[i]], .(COG20_CATEGORY), summarize, estimate = sum(estimate))
  }
}

plot.list2 <- list()

for(i in names(sig.dif.fun)){
    
    sig.dif.fun[[i]]$COG20_CATEGORY <- factor(sig.dif.fun[[i]]$COG20_CATEGORY, levels = sig.dif.fun[[i]]$COG20_CATEGORY[order(sig.dif.fun[[i]]$estimate)])
    
    p <- ggplot(sig.dif.fun[[i]], aes(x = COG20_CATEGORY, y = estimate)) + 
      geom_point(size = 2) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 12.25),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)
          ) +
      ylim(-2,1.7) +
      labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change"))) +
      scale_x_discrete(labels = c("Posttranslational modification, protein turnover, chaperones" = "Posttranslational modification, protein\nturnover, chaperones", "Intracellular trafficking, secretion, and vesicular transport" = "Intracellular trafficking, secretion, and\nvesicular transport")) +
      coord_flip()
  
   plot.list2[[i]] = p
}

plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[3]], plot.list2[[2]], plot.list2[[4]], align = "v")

fun.fig <- ggarrange(plot.list2[[1]], ggarrange(plot.list2[[2]], plot.list2[[3]], plot.list2[[4]], ncol = 1, labels = c("(b)", "(c)", "(d)")), ncol = 2, labels = "(a)")

ggsave("Figures/deseq-fun.jpeg", fun.fig, height = 6, width = 11.75, units = "in", dpi = 600)

# DESeq2 taxa ####
## Run DESeq2 ####
ps.genus <- phyloseq(otu_table(decontam_counts_nocontrols, 
                         taxa_are_rows = T), 
               sample_data(metadata.nocontrols),
               tax_table(tax.gene))

ps.genus <- tax_glom(ps.genus, taxrank = "Genus")

treat.genus <- phyloseq_to_deseq2(ps.genus, ~ plant.water)

dds.pw.genus <- DESeq(treat.genus)

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw.genus, contrast = contrast.list[[i]], pAdjustMethod = "BH")
  #filter results by p-value
  res.alpha <- res[which(res$padj < alpha), ]
  #tidy results 
  res.list[[i]] <- tidy(res.alpha)

 res.list[[i]]$gene <- factor(res.list[[i]]$gene, levels = res.list[[i]]$gene[order(res.list[[i]]$estimate)])
  
 if(nrow(res.list[[i]]) > 0) {
   p <- ggplot(res.list[[i]], aes(x = gene, y = estimate)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.title.x = element_blank(),
        ) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")))
  
   plot.list[[i]] = p
 }
   sig.genes <- c(as.numeric(as.character(res.list[[i]]$gene)), sig.genes)
}

# results: 
length(unique(sig.genes)) # 41

## Plot results ####
grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]])


test2 <- as.data.frame(cbind(tax.gene, gene = row.names(tax.gene)))

for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], test2, by = "gene", all.y = F)
}

legend.prep <- rbind(res.list[[3]], res.list[[4]], res.list[[14]])

sig.dif.tax <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.tax[[i]] <- ddply(res.list[[i]], .(Family, Phylum), summarize, estimate = sum(estimate))
  }
}

plot.list2 <- list()


for(i in names(sig.dif.tax)){
    
  sig.dif.tax[[i]]$Family <- recode_factor(sig.dif.tax[[i]]$Family, 
                                         Burkholderiales = "Unclassified Burkholderiales",
                                         Chloroflexi = "Unclassified Chloroflexi",
                                         Deltaproteobacteria = "Unclassified Deltaproteobacteria",
                                         Coriobacteriia = "Unclassified Coriobacteriia",
                                         Myxococcales = "Unclassified Myxococcales",
                                         Thaumarchaeota = "Unclassified Thaumarchaeota"
                                         )
  
    sig.dif.tax[[i]]$Family <- factor(sig.dif.tax[[i]]$Family, levels = sig.dif.tax[[i]]$Family[order(sig.dif.tax[[i]]$estimate)])
    
    p <- ggplot(sig.dif.tax[[i]], aes(x = Family, y = estimate, col = Phylum)) + 
      geom_point(size = 2.5) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15)
          ) +
      ylim(-2.3, 5.5) +
      labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change"))) +
      scale_colour_manual(values = c(
        "Actinobacteria" = "dodgerblue2",                     
        "candidate division AD3" = "#E31A1C",                      
        "Candidatus Tectomicrobia" = "green4",                               
        "Chloroflexi" = "#6A3D9A",                                        
        "Firmicutes" = "#FF7F00",
        "Proteobacteria" = "gold1",
        "Bacteroidetes" = "skyblue2",
        "Acidobacteria" = "palegreen2",
        "Thaumarchaeota" =  "orchid1"
      )) +
      coord_flip()
  
   plot.list2[[i]] = p
}

legend <- get_legend(ggplot(legend.prep, aes(x = Family, y = estimate, col = Phylum)) + 
      geom_point(size = 2.5) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.text =  element_text(size = 15),
        legend.title = element_text(size = 17),
        legend.position = "right",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15)
          ) +
      ylim(-2.3, 5.5) +
      scale_colour_manual(values = c(
        "Actinobacteria" = "dodgerblue2",                     
        "candidate division AD3" = "#E31A1C",                      
        "Candidatus Tectomicrobia" = "green4",                               
        "Chloroflexi" = "#6A3D9A",                                        
        "Firmicutes" = "#FF7F00",
        "Proteobacteria" = "gold1",
        "Bacteroidetes" = "skyblue2",
        "Acidobacteria" = "palegreen2",
        "Thaumarchaeota" =  "orchid1"
      )) +
      coord_flip())

plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], align = "v")

tax.fig <- ggarrange(ggarrange(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], heights = c(1,0.4,0.72), labels = c("(a)", "(b)", "(c)"), ncol = 1), legend, ncol = 2, widths = c(1,0.5))

ggsave("Figures/deseq-tax.jpeg", tax.fig, height = 13, width = 10, units = "in", dpi = 600)

# Rarefy and agglomerate ####
ps <- phyloseq(otu_table(decontam_counts_nocontrols, 
                         taxa_are_rows = T), 
               sample_data(metadata.nocontrols))

## Look at library size
df <- as.data.frame(sample_data(ps)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

#plot library size with line at 25000 reads
ggplot(data=df, aes(x=Index, y=LibrarySize, color=plant.water)) + 
  #geom_point() + 
  geom_hline(yintercept = 25000, linetype="dashed") + 
  theme(text = element_text(size=18)) +
  geom_text(aes(label = SampleNumber))

summary(df$LibrarySize)

hist(df$LibrarySize, breaks = 50) #right skewed 

# dropping the smallest sample (library size of 35618) GLM 260
ps.rare668803 <- rarefy_even_depth(ps, sample.size = 66880, replace = FALSE, rngseed = 5311) 

ps.fun.rare66 <- ps.rare66880
tax_table(ps.fun.rare66) <- fun.gene
ps.fun.rare66 <- tax_glom(ps.fun.rare66, taxrank = "COG20_FUNCTION")

ps.genus.rare66 <- ps.rare66880
tax_table(ps.genus.rare66) <- tax.gene
ps.genus.rare66 <- tax_glom(ps.genus.rare66, taxrank = "Genus")

saveRDS(ps.rare66880, "Data/Data-Processing/ps-nocontrols-rare-66880.RDS")
saveRDS(ps.fun.rare66, "Data/Data-Processing/ps-nocontrols-rare-66880-fun.RDS")
saveRDS(ps.genus.rare66, "Data/Data-Processing/ps-nocontrols-rare-66880-tax.RDS")

ps.rare668802 <- readRDS("Data/Data-Processing/ps-nocontrols-rare-668802.RDS")
ps.fun.rare66 <- readRDS("Data/Data-Processing/ps-nocontrols-rare-66880-fun.RDS")
ps.genus.rare66 <- readRDS("Data/Data-Processing/ps-nocontrols-rare-66880-tax.RDS")

# Ordination - function ####
### Stats ####
PCoA_bray <- ordinate(physeq = ps.fun.rare66, method = "PCoA", distance = "bray")

Dist.fun <- phyloseq::distance(ps.fun.rare66, method = "bray", type = "samples")

set.seed(620)
adonis2(Dist.fun ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.fun.rare66), "data.frame"), by = "margin", permutations = 9999)  

set.seed(620)
adonis2(Dist.fun ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.fun.rare66), "data.frame"), by = "margin", permutations = 9999)

#### Contrasts ####
# nat.inv #
nat.inv <- prune_samples(sample_data(ps.fun.rare66)$PlantTrt != "mix", ps.fun.rare66) 

Dist.fun.nat.inv <- phyloseq::distance(nat.inv, method = "bray", type = "samples")

set.seed(620)
ad.nat.inv <- adonis2(Dist.fun.nat.inv ~ WaterTrt + PlantTrt + Plot, as(sample_data(nat.inv), "data.frame"), by = "margin", permutations = 9999)

# nat.mix #
nat.mix <- prune_samples(sample_data(ps.fun.rare66)$PlantTrt != "invasives", ps.fun.rare66) 

Dist.fun.nat.mix <- phyloseq::distance(nat.mix, method = "bray", type = "samples")

ad.nat.mix <- adonis2(Dist.fun.nat.mix ~ WaterTrt + PlantTrt + Plot, as(sample_data(nat.mix), "data.frame"), by = "margin", permutations = 9999)

# inv.mix #
inv.mix <- prune_samples(sample_data(ps.fun.rare66)$PlantTrt != "natives", ps.fun.rare66) 

Dist.fun.inv.mix <- phyloseq::distance(inv.mix, method = "bray", type = "samples")

set.seed(620)
ad.inv.mix <- adonis2(Dist.fun.inv.mix ~ WaterTrt + PlantTrt + Plot, as(sample_data(inv.mix), "data.frame"), by = "margin", permutations = 9999)

p.adjust(p = c(ad.inv.mix$`Pr(>F)`[2],
               ad.nat.mix$`Pr(>F)`[2],
               ad.nat.inv$`Pr(>F)`[2]
               ), "BH")

# w.d #
w.d <- prune_samples(sample_data(ps.fun.rare66)$WaterTrt != "control", ps.fun.rare66) 

Dist.fun.w.d <- phyloseq::distance(w.d, method = "bray", type = "samples")

set.seed(620)
ad.w.d <- adonis2(Dist.fun.w.d ~ WaterTrt + PlantTrt + Plot, as(sample_data(w.d), "data.frame"), by = "margin", permutations = 9999)

# w.c #
w.c <- prune_samples(sample_data(ps.fun.rare66)$WaterTrt != "drought", ps.fun.rare66) 

Dist.fun.w.c <- phyloseq::distance(w.c, method = "bray", type = "samples")

set.seed(620)
ad.w.c <- adonis2(Dist.fun.w.c ~ WaterTrt + PlantTrt + Plot, as(sample_data(w.c), "data.frame"), by = "margin", permutations = 9999)

# d.c #
d.c <- prune_samples(sample_data(ps.fun.rare66)$WaterTrt != "watered", ps.fun.rare66) 

Dist.fun.d.c <- phyloseq::distance(d.c, method = "bray", type = "samples")

set.seed(620)
ad.d.c <- adonis2(Dist.fun.d.c ~ WaterTrt + PlantTrt + Plot, as(sample_data(d.c), "data.frame"), by = "margin", permutations = 9999)

p.adjust(p = c(ad.w.d$`Pr(>F)`[1],
               ad.w.c$`Pr(>F)`[1],
               ad.d.c$`Pr(>F)`[1]
               ), "BH")

### Figure ####
PCoA_bray.fun <- ordinate(physeq = ps.fun.rare66, method = "PCoA", distance = "bray")

# differences between watering treatments
a <- plot_ordination(ps.fun.rare66, PCoA_bray.fun, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.87,0.9),
     legend.title = element_blank(),
     #legend.box.margin = margin(-1,-1,-1,-1),
     #legend.margin = margin(-3,-3,-3,-3),
     legend.text = element_text(size = 17),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18),
     plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  geom_point(size = 3) +
  lims(x = c(-0.13, 0.19)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (22.2%)", y = "PCoA2 (11.4%)")

# differences between plant treatments
b <- plot_ordination(ps.fun.rare66, PCoA_bray.fun, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.87,0.9),
    legend.title = element_blank(),
    #legend.box.margin = margin(-1,-1,-1,-1),
    #legend.margin = margin(-3,-3,-3,-3),
    legend.text = element_text(size = 17),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  lims(x = c(-0.13, 0.19)) +
  geom_point(size = 3) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (22.2%)", y = "PCoA2 (11.4%)")

# get centroids plots
C <- betadisper(Dist.fun, as(sample_data(ps.fun.rare66), "data.frame")$plant.water)

test <- data.frame(Category=rownames(C$centroids[,1:2]), C$centroids[,1:2])

forb <- test$Category[1:3]
grass <- test$Category[4:6]
drought <- test$Category[c(2,5,8)]
water <- test$Category[c(3,6,9)]

my.plot <- gg_ordiplot(C, groups = as(sample_data(ps.fun.rare66), "data.frame")$plant.water, pt.size = 0, kind = "se") 

my.plot$df_mean.ord$PlantTrt <- ifelse(my.plot$df_mean.ord$Group %in% forb, "natives", ifelse(my.plot$df_mean.ord$Group %in% grass, "invasives", "mix"))

my.plot$df_mean.ord$PlantTrt <- factor(my.plot$df_mean.ord$PlantTrt, levels = c("natives", "invasives", "mix"))

my.plot$df_mean.ord$WaterTrt <- ifelse(my.plot$df_mean.ord$Group %in% drought, "drought", ifelse(my.plot$df_mean.ord$Group %in% water, "watered", "control"))

my.plot$df_mean.ord$WaterTrt <- factor(my.plot$df_mean.ord$WaterTrt, levels = c("drought", "control", "watered"))

my.plot$df_ellipse$PlantTrt <- ifelse(my.plot$df_ellipse$Group %in% forb, "natives", ifelse(my.plot$df_ellipse$Group %in% grass, "invasives", "mix"))

my.plot$df_ellipse$PlantTrt <- factor(my.plot$df_ellipse$PlantTrt, levels = c("natives", "invasives", "mix"))

my.plot$df_ellipse$WaterTrt <- ifelse(my.plot$df_ellipse$Group %in% drought, "drought", ifelse(my.plot$df_ellipse$Group %in% water, "watered", "control"))

my.plot$df_ellipse$WaterTrt <- factor(my.plot$df_ellipse$WaterTrt, levels = c("drought", "control", "watered"))

c <- ggplot(my.plot$df_mean.ord, aes(x = x, y = y, col = WaterTrt, shape = PlantTrt)) +
  geom_point(data = my.plot$df_ord, mapping = aes(x = x, y = y), col = "lightgrey", size = 3, inherit.aes = F) +
  geom_point(size = 3.5) +
  geom_path(data = my.plot$df_ellipse, aes(x = x, y = y), show.legend = FALSE) +
  scale_color_manual(values = c("red", "black", "blue")) +
  theme_classic(base_size = 10) +
  theme(
     legend.box = "vertical",
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = "right",
     legend.title = element_blank(),
    # legend.box.margin = margin(-1,-1,-1,-1),
    # legend.margin = margin(-3,-3,-3,-3),
     legend.text = element_text(size = 17),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18),
    plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) +
  labs(x = "PCoA1 (22.2%)", y = "PCoA2 (11.4%)")
  #ylim(-0.075,0.075)

d <- ggarrange(ggarrange(a,b, ncol = 2, labels = c("(a)", "(b)"), font.label = list(size = 17)), c, labels = c("","(c)"), nrow = 2, heights = c(0.8, 1.2), font.label = list(size = 17))

ggsave("Figures/Ordination-fun2.jpeg", d, height = 10, width = 10, units = "in", dpi = 600)

# Ordination - Taxa ####
Dist.taxa <- phyloseq::distance(ps.genus.rare66, method = "bray", type = "samples")

### Stats ####
set.seed(620)
adonis2(Dist.taxa ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.genus.rare66), "data.frame"), by = "margin", permutations = 9999)

set.seed(620)
adonis2(Dist.taxa ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.genus.rare66), "data.frame"), by = "margin", permutations = 9999)

#### Contrasts ####
# nat.inv #
nat.inv <- prune_samples(sample_data(ps.genus.rare66)$PlantTrt != "mix", ps.genus.rare66) 

Dist.genus.nat.inv <- phyloseq::distance(nat.inv, method = "bray", type = "samples")

set.seed(620)
ad.nat.inv <- adonis2(Dist.genus.nat.inv ~ WaterTrt + PlantTrt + Plot, as(sample_data(nat.inv), "data.frame"), by = "margin", permutations = 9999)

# nat.mix #
nat.mix <- prune_samples(sample_data(ps.genus.rare66)$PlantTrt != "invasives", ps.genus.rare66) 

Dist.genus.nat.mix <- phyloseq::distance(nat.mix, method = "bray", type = "samples")

ad.nat.mix <- adonis2(Dist.genus.nat.mix ~ WaterTrt + PlantTrt + Plot, as(sample_data(nat.mix), "data.frame"), by = "margin", permutations = 9999)

# inv.mix #
inv.mix <- prune_samples(sample_data(ps.genus.rare66)$PlantTrt != "natives", ps.genus.rare66) 

Dist.genus.inv.mix <- phyloseq::distance(inv.mix, method = "bray", type = "samples")

set.seed(620)
ad.inv.mix <- adonis2(Dist.genus.inv.mix ~ WaterTrt + PlantTrt + Plot, as(sample_data(inv.mix), "data.frame"), by = "margin", permutations = 9999)

p.adjust(p = c(ad.inv.mix$`Pr(>F)`[2],
               ad.nat.mix$`Pr(>F)`[2],
               ad.nat.inv$`Pr(>F)`[2]
               ), "BH")

# w.d #
w.d <- prune_samples(sample_data(ps.genus.rare66)$WaterTrt != "control", ps.genus.rare66) 

Dist.genus.w.d <- phyloseq::distance(w.d, method = "bray", type = "samples")

set.seed(620)
ad.w.d <- adonis2(Dist.genus.w.d ~ WaterTrt + PlantTrt + Plot, as(sample_data(w.d), "data.frame"), by = "margin", permutations = 9999)

# w.c #
w.c <- prune_samples(sample_data(ps.genus.rare66)$WaterTrt != "drought", ps.genus.rare66) 

Dist.genus.w.c <- phyloseq::distance(w.c, method = "bray", type = "samples")

set.seed(620)
ad.w.c <- adonis2(Dist.genus.w.c ~ WaterTrt + PlantTrt + Plot, as(sample_data(w.c), "data.frame"), by = "margin", permutations = 9999)

# d.c #
d.c <- prune_samples(sample_data(ps.genus.rare66)$WaterTrt != "watered", ps.genus.rare66) 

Dist.genus.d.c <- phyloseq::distance(d.c, method = "bray", type = "samples")

set.seed(620)
ad.d.c <- adonis2(Dist.genus.d.c ~ WaterTrt + PlantTrt + Plot, as(sample_data(d.c), "data.frame"), by = "margin", permutations = 9999)

p.adjust(p = c(ad.d.c$`Pr(>F)`[1],
               ad.w.d$`Pr(>F)`[1],
               ad.w.c$`Pr(>F)`[1]
               ), "BH")

### Figure ####
PCoA_bray.taxa <- ordinate(physeq = ps.genus.rare66, method = "PCoA", distance = "bray")

# differences between watering treatments
a <- plot_ordination(ps.genus.rare66, PCoA_bray.taxa, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.87,0.9),
     legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3),
    legend.text = element_text(size = 14),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18)) +
  geom_point(size = 3) +
  lims(x = c(-0.23, 0.25)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (31.8%)", y = "PCoA2 (24.8%)")

# differences between plant treatments
b <- plot_ordination(ps.genus.rare66, PCoA_bray.taxa, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.15,0.9),
    legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
    legend.margin = margin(-3,-3,-3,-3),
    legend.text = element_text(size = 14),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18)) +
  lims(x = c(-0.23, 0.25)) +
  geom_point(size = 3) +
  stat_ellipse() +
  scale_y_continuous(breaks = c(-0.2,-0.1,0,0.1)) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (31.8%)", y = "PCoA2 (24.8%)")

# get centroids plots
C <- betadisper(Dist.taxa, as(sample_data(ps.genus.rare66), "data.frame")$plant.water)

test <- data.frame(Category=rownames(C$centroids[,1:2]), C$centroids[,1:2])

forb <- test$Category[1:3]
grass <- test$Category[4:6]
drought <- test$Category[c(2,5,8)]
water <- test$Category[c(3,6,9)]

my.plot <- gg_ordiplot(C, groups = as(sample_data(ps.genus.rare66), "data.frame")$plant.water, pt.size = 0, kind = "se") 

my.plot$df_mean.ord$PlantTrt <- ifelse(my.plot$df_mean.ord$Group %in% forb, "natives", ifelse(my.plot$df_mean.ord$Group %in% grass, "invasives", "mix"))

my.plot$df_mean.ord$PlantTrt <- factor(my.plot$df_mean.ord$PlantTrt, levels = c("natives", "invasives", "mix"))

my.plot$df_mean.ord$WaterTrt <- ifelse(my.plot$df_mean.ord$Group %in% drought, "drought", ifelse(my.plot$df_mean.ord$Group %in% water, "watered", "control"))

my.plot$df_mean.ord$WaterTrt <- factor(my.plot$df_mean.ord$WaterTrt, levels = c("drought", "control", "watered"))

my.plot$df_ellipse$PlantTrt <- ifelse(my.plot$df_ellipse$Group %in% forb, "natives", ifelse(my.plot$df_ellipse$Group %in% grass, "invasives", "mix"))

my.plot$df_ellipse$PlantTrt <- factor(my.plot$df_ellipse$PlantTrt, levels = c("natives", "invasives", "mix"))

my.plot$df_ellipse$WaterTrt <- ifelse(my.plot$df_ellipse$Group %in% drought, "drought", ifelse(my.plot$df_ellipse$Group %in% water, "watered", "control"))

my.plot$df_ellipse$WaterTrt <- factor(my.plot$df_ellipse$WaterTrt, levels = c("drought", "control", "watered"))

c <- ggplot(my.plot$df_mean.ord, aes(x = x, y = y, col = WaterTrt, shape = PlantTrt)) +
  geom_point(data = my.plot$df_ord, mapping = aes(x = x, y = y), col = "lightgrey", size = 3, inherit.aes = F) +
  geom_point(size = 3.5) +
  geom_path(data = my.plot$df_ellipse, aes(x = x, y = y), show.legend = FALSE) +
  scale_color_manual(values = c("red", "black", "blue")) +
  theme_classic(base_size = 10) +
  theme(
     legend.box = "vertical",
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = "right",
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3),
     legend.text = element_text(size = 14),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18)) +
  ylim(-0.16,0.14) +
  labs(x = "PCoA1 (31.8%)", y = "PCoA2 (24.8%)")

#d <- ggarrange(a,b,c, ncol = 3, labels = c("(a)", "(b)", "(c)"), widths = c(1,1,1.2))

d <- ggarrange(ggarrange(a,b, ncol = 2, labels = c("(a)", "(b)"), font.label = list(size = 17)), c, labels = c("","(c)"), nrow = 2, heights = c(0.8, 1.2), font.label = list(size = 17))

ggsave("Figures/Ordination-taxa.jpeg", d, height = 10, width = 10, units = "in", dpi = 600)

# Distance Relationships ####

# get lat and long
geo = data.frame(ps@sam_data$long, ps@sam_data$lat)

# defining geographic distance between two locations -
# haversine distance (accounts for spherical earth)
d.geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)

ps.hell <- transform(ps, "hellinger")

# getting community distances have to multiply by 1/sqrt(2)
# to get scaled 0-1 as current range is 0-sqrt(2)
comm.dist.hell = vegdist(1/sqrt(2) * as.matrix(t(ps.hell@otu_table)),
    method = "euclidean")

# test
mantel.test.hell = vegan::mantel(comm.dist.hell, dist.geo, method = "spearman",
    permutations = 9999, na.rm = TRUE)
mantel.test.hell

mantel.stats <- data.frame(label = paste("r = ", signif(mantel.test.hell$statistic,
    3), "\nP = ", signif(mantel.test.hell$signif, 3)))

# plot
dist_hell <- ggplot(mapping = aes(x = jitter(dist.geo, amount = 1),
    y = comm.dist.hell)) + 
  theme_bw() + 
  geom_point(shape = 16, size = 1, alpha = 0.1, color = "gray25") + 
  geom_smooth(method = "lm", color = "orange", se = F) + 
  labs(x = "Geographical Separation (m)", y = "Hellinger Distance") + 
  geom_text(data = mantel.stats, aes(x = 0.5, y = 0.5, label = label), hjust = 0, size = 3) +
  theme(text = element_text(size = 14))

dist_hell

# Supplement Table ####
samp_data <- metadata.nocontrols %>%
  mutate(AllPlantSpeciesInCore = ifelse(str_detect(AllPlantSpeciesInCore, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInCore)),
         'Bromus hordeaceus' = ifelse(str_detect(AllPlantSpeciesInCore, "Bromus"), 1, 0),
         'Avena fatua' = ifelse(str_detect(AllPlantSpeciesInCore, "Avena"), 1, 0),
         'Elymus caput-medusae' = ifelse(str_detect(AllPlantSpeciesInCore, "Elymus"), 1, 0),
         
         'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInCore, "Festuca perennis"), 1, 0),
         'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInCore, "Festuca microstachys"), 1, 0),
          'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInCore, "Vulpia microstachys"), 1, 0),
         'Calycadenia pauciflora' = ifelse(str_detect(AllPlantSpeciesInCore, "Calycadenia"), 1, 0),
         'Hemizonia congesta' = ifelse(str_detect(AllPlantSpeciesInCore, "Hemizonia"), 1, 0),
         'Plantago erecta' = ifelse(str_detect(AllPlantSpeciesInCore, "Plantago"), 1, 0),
         'Clarkia purpurea' = ifelse(str_detect(AllPlantSpeciesInCore, "Clarkia"), 1, 0),
         'Agoseris heterophylla' = ifelse(str_detect(AllPlantSpeciesInCore, "Agoseris"), 1, 0),
         'Lasthenia californica' = ifelse(str_detect(AllPlantSpeciesInCore, "Lasthenia"), 1, 0)) %>%
  select(WaterTrt, PlantTrt, SampleID, 'Bromus hordeaceus', 'Avena fatua', 'Elymus caput-medusae', 'Festuca perennis','Festuca microstachys', 'Calycadenia pauciflora', 'Hemizonia congesta', 'Plantago erecta', 'Clarkia purpurea', 'Agoseris heterophylla', 'Lasthenia californica') %>%
  arrange(WaterTrt, PlantTrt)

write.csv(samp_data, "Figures/Supp-Table.csv", row.names = F)
