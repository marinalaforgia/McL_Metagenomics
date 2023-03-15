### Script to analyze Metagenomic soil samples from McLaughlin ###

# Files needed (all are found within the R-Data-Files folder on box)
## 1. Kallisto file of raw gene counts
## 2. Functional pathways file
## 3. Joined MMSeq and Kaiju taxonomy file
## 4. Metadata file

## rarefy then agglomerate for ordination to Cassie for arrows
## non-rarefied agglomerate for deseq
## 
# Outline of analysis

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

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forb.Watering/forbs.Control) so if something is higher in forb watering, the number is greater than one, so the log fold change > 0, if something is higher in the control treatment then the number is less than one so the log2 fold change is negative


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

# DESeq2 genes ####
ps <- phyloseq(otu_table(decontam_counts_nocontrols, 
                         taxa_are_rows = T), 
               sample_data(metadata.nocontrols))

ps <- prune_samples(sample_data(ps)$SampleNumber != 241, ps)

treat <- phyloseq_to_deseq2(ps, ~ plant.water)

dds.pw <- DESeq(treat)

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
length(unique(sig.genes))  #0

# DESeq2 functions ####
ps.fun <- phyloseq(otu_table(decontam_counts_nocontrols, 
                         taxa_are_rows = T), 
               sample_data(metadata.nocontrols),
               tax_table(fun.gene))

ps.fun <- tax_glom(ps.fun, taxrank = "COG20_FUNCTION")

treat.fun <- phyloseq_to_deseq2(ps.fun, ~ plant.water) # but this imput shouldnt be rarefied omg im so annoyed 

dds.pw.fun <- DESeq(treat.fun)

#resultsNames(dds.pw.fun) #gives us comparisons for our contrasts

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


# results: 
length(unique(sig.genes)) #44

names(plot.list)

test <- rbind(data.frame(cbind(pair = "MD.FC", res.list[["MD.FC"]])), data.frame(cbind(pair = "MD.FD", res.list[["MD.FD"]])), data.frame(cbind(pair = "FMW.FC", res.list[["MW.FC"]])), data.frame(cbind(pair = "MC.FC", res.list[["MC.FC"]])))


write.csv(test, "Data/fun-gene-DEG.csv", row.names = F)

#plot results
grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], ncol = 2)

test2 <- as.data.frame(cbind(fun.gene, gene = row.names(fun.gene)))

for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], test2, by = "gene", all.y = F)
}

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

#fun.fig <- grid.arrange(plot.list2[[1]], arrangeGrob(plot.list2[[2]], plot.list2[[3]], plot.list2[[4]]), ncol = 2)

fun.fig <- ggarrange(plot.list2[[1]], ggarrange(plot.list2[[2]], plot.list2[[3]], plot.list2[[4]], ncol = 1, labels = c("(b)", "(c)", "(d)")), ncol = 2, labels = "(a)")

ggsave("Figures/deseq-fun.jpeg", fun.fig, height = 6, width = 11.75, units = "in", dpi = 600)

# DESeq2 taxa ####
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

#plot results
grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]])


test2 <- as.data.frame(cbind(tax.gene, gene = row.names(tax.gene)))

for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], test2, by = "gene", all.y = F)
}

sig.dif.tax <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.tax[[i]] <- ddply(res.list[[i]], .(Family), summarize, estimate = sum(estimate))
  }
}

plot.list2 <- list()

for(i in names(sig.dif.tax)){
    
    sig.dif.tax[[i]]$Family <- factor(sig.dif.tax[[i]]$Family, levels = sig.dif.tax[[i]]$Family[order(sig.dif.tax[[i]]$estimate)])
    
    p <- ggplot(sig.dif.tax[[i]], aes(x = Family, y = estimate)) + 
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
      coord_flip()
  
   plot.list2[[i]] = p
}

plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], align = "v")

tax.fig <- ggarrange(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], heights = c(1,0.4,0.72), labels = c("(a)", "(b)", "(c)"), ncol = 1)

ggsave("Figures/deseq-tax.jpeg", tax.fig, height = 12, width = 7, units = "in", dpi = 600)

# 
# df <- data.frame()
# 
# for(i in names(plot.list)){ #plot list is the one with all the significantly different ones
#   tmp <- cbind(pair = i, res.list[[i]])
#   df <- rbind(df, tmp)
# }
# 
# test2 <- as.data.frame(cbind(tax.gene, gene = row.names(tax.gene)))
# 
# df <- merge(df, test2, by = "gene", all.y = F)
# 
# write.csv(df, "Data/tax-gene-DEG.csv")

# Rarefy and agglomerate ####
ps <- phyloseq(otu_table(decontam_counts_nocontrols, 
                         taxa_are_rows = T), 
               sample_data(metadata.nocontrols))
## Look at library size
df <- as.data.frame(sample_data(ps.fun)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(ps.fun)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))

#plot library size with line at 25000 reads
ggplot(data=df, aes(x=Index, y=LibrarySize, color=plant.water)) + 
  #geom_point() + 
  geom_hline(yintercept = 25000, linetype="dashed") + 
  theme(text = element_text(size=18)) +
  geom_text(aes(label = SampleNumber))

summary(df$LibrarySize)

# Given that median is 10000, first qu. is ~4800, nice breaks around ~2000 or 4000 reads 
hist(df$LibrarySize, breaks = 50) #right skewed 

#which samples do we lose?
as_tibble(df) %>% filter(LibrarySize < 2000) %>%
  group_by(Location, Rel.Res, Wild.Captive) %>%
  tally() #lose a lot of captive sus

#log10(1) = 0, log10(10) = 1, log10(100) = 2, log10(1000) = 3, so rarefy at 1000 line would be:
ggplot(data=df, aes(x=log10(LibrarySize), fill=Wild.Captive)) + geom_histogram(binwidth  =1, position= "dodge", col="black") + theme(text = element_text(size=18))  + ylab("Frequency of samples") + xlab(expression(paste("lo", g[10]," transformed read counts")))+scale_x_continuous(breaks=c(0,1,2,3,4,5)) + geom_vline(xintercept = 3.5, linetype="solid", col="#EF7F4FFF", size =2)




# dropping the smallest sample (library size of 35618) GLM 260
ps.rare66880 <- rarefy_even_depth(ps.fun, sample.size = 66880, replace = FALSE, rngseed = 5311) 

ps.fun.rare66 <- ps.rare66880
tax_table(ps.fun.rare66) <- fun.gene
ps.fun.rare66 <- tax_glom(ps.fun.rare66, taxrank = "COG20_FUNCTION")

ps.genus.rare66 <- ps.rare66880
tax_table(ps.genus.rare66) <- tax.gene
ps.genus.rare66 <- tax_glom(ps.genus.rare66, taxrank = "Genus")

saveRDS(ps.rare66880, "Data/Data-Processing/ps-nocontrols-rare-66880.RDS")
saveRDS(ps.fun.rare66, "Data/Data-Processing/ps-nocontrols-rare-66880-fun.RDS")
saveRDS(ps.genus.rare66, "Data/Data-Processing/ps-nocontrols-rare-66880-tax.RDS")

# Ordination - genes ####
###Stats ####
#norm_counts <- counts(dds.pw, normalized = T)

#ps <- phyloseq(otu_table(norm_counts, taxa_are_rows = T), sample_data(metadata.nocontrols))

# ps.hell <- transform(ps.rare81000, "hellinger")
# 
# PCoA_hell <- ordinate(physeq = ps.hell, method = "PCoA", distance = "euclidean")
PCoA_bray <- ordinate(physeq = ps.rare81000, method = "PCoA", distance = "bray")

# Get distances for stats
#Dist.hell <- phyloseq::distance(ps.hell, method = "euclidean", type = "samples")
Dist.bray <- phyloseq::distance(ps.rare81000, method = "bray", type = "samples")

set.seed(620)
adonis2(Dist.bray ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.rare81000), "data.frame"), by = "margin", permutations = 9999)

set.seed(620)
adonis2(Dist.bray ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.rare81000), "data.frame"), by = "margin", permutations = 9999)

# #bpart <- beta.pair.abund(Dist.hell, index.family = "bray")
# 
# # set.seed(620)
# # adonis2(bpart$beta.bray.gra ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999)  #nestedness
# 
# # set.seed(620)
# # adonis2(bpart$beta.bray.bal ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) # turnover
# 
# set.seed(620)
# adonis2(Dist.hell ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999)
# 
# # set.seed(620)
# # adonis2(bpart$beta.bray.gra ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) # nestedness
# 
# # set.seed(620)
# # adonis2(bpart$beta.bray.bal ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) #turnover
# # # signficant turnover between watering treatments
# 
# set.seed(620)
# adonis2(Dist.hell ~ WaterTrt + PlotSubTreatment + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999)
# 
# C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$plant.water)
# 
# set.seed(50)
# permutest(C, permutations = 9999) # no sig dif
# 
# C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$PlantTrt)
# 
# set.seed(50)
# permutest(C, permutations = 9999) # no sig dif
# 
# C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$WaterTrt)
# 
# set.seed(50)
# permutest(C, permutations = 9999) # no sig dif
# 
# pairwise <- permanova_pairwise(
#   Dist.hell,
#   as(sample_data(ps.hell), "data.frame")$WaterTrt,
#   permutations = 9999,
#   padj = "BH"
# )
# 
# pairwise <- permanova_pairwise(
#   Dist.hell,
#   as(sample_data(ps.hell), "data.frame")$PlantTrt,
#   permutations = 9999,
#   padj = "BH"
# )

### Figure ####
# differences between watering treatments
a <- plot_ordination(ps.hell, PCoA_hell, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.8,0.15),
     legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (10.9%)", y = "PCoA2 (7.4%)")

plot_ordination(ps.rare81000, PCoA_bray, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.8,0.15),
     legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) 

# differences between plant treatments
b <- plot_ordination(ps.hell, PCoA_hell, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.85,0.15),
    legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
    legend.margin = margin(-3,-3,-3,-3)) +
  lims(x = c(-0.45, 0.55)) +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (10.9%)", y = "PCoA2 (7.4%)")


# differences between interactions
# ps.hell.sub <- subset_samples(ps.hell, plant.water == "forb.Control" | plant.water == "mix.Drought" )
# 
# plot_ordination(ps.hell.sub, PCoA_hell, color = "plant.water") +
#   theme_classic(base_size = 15) +
#   theme(
#      panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   geom_point(size = 2) +
#   stat_ellipse() +
#   scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))
# 
# 
# PCoA.df <- merge(data.frame(SampleID = rownames(PCoA_hell$vectors), PCoA_hell$vectors)[,1:3], metadata.nocontrols, by = "SampleID")
# 
# PCoA.df$plot <- ifelse(PCoA.df$plant.water == "forb.Control" | PCoA.df$plant.water == "mix.Drought", PCoA.df$plant.water, "none")
# 
# ggplot(PCoA.df[PCoA.df$plot != "none",], aes(x = Axis.1, y = Axis.2, col = plot)) +
#   theme_classic(base_size = 15) +
#   theme(
#      panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   geom_point(size = 2) +
#   stat_ellipse()
#   
# PCoA.df$plot <- factor(PCoA.df$plot, levels = c("none","forb.Control", "mix.Drought"))
# 
# #make grey more transparent, put gray points in the background
# ggplot(PCoA.df, aes(x = Axis.1, y = Axis.2, col = plot)) +
#   theme_classic(base_size = 15) +
#   theme(
#      panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
#   geom_point(size = 2) +
#   stat_ellipse() +
#   scale_color_manual(values = c("lightgrey", "black", "red"))

# get centroids plots
C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$plant.water)

test <- data.frame(Category=rownames(C$centroids[,1:2]), C$centroids[,1:2])

forb <- test$Category[1:3]
grass <- test$Category[4:6]
drought <- test$Category[c(2,5,8)]
water <- test$Category[c(3,6,9)]

test$PlantTrt <- ifelse(test$Category %in% forb, "forb", ifelse(test$Category %in% grass, "grass", "mix"))

test$WaterTrt <- ifelse(test$Category %in% drought, "shelter", ifelse(test$Category %in% water, "watered", "control"))

my.plot <- gg_ordiplot(C, groups = as(sample_data(ps.hell), "data.frame")$plant.water, pt.size = 0, kind = "se") 

my.plot$df_mean.ord$PlantTrt <- ifelse(my.plot$df_mean.ord$Group %in% forb, "forb", ifelse(my.plot$df_mean.ord$Group %in% grass, "grass", "mix"))

my.plot$df_mean.ord$WaterTrt <- ifelse(my.plot$df_mean.ord$Group %in% drought, "shelter", ifelse(my.plot$df_mean.ord$Group %in% water, "watered", "control"))

my.plot$df_ellipse$PlantTrt <- ifelse(my.plot$df_ellipse$Group %in% forb, "forb", ifelse(my.plot$df_ellipse$Group %in% grass, "grass", "mix"))

my.plot$df_ellipse$WaterTrt <- ifelse(my.plot$df_ellipse$Group %in% drought, "shelter", ifelse(my.plot$df_ellipse$Group %in% water, "watered", "control"))

c <- ggplot(my.plot$df_mean.ord, aes(x = x, y = y, col = WaterTrt, shape = PlantTrt)) +
  geom_point(data = my.plot$df_ord, mapping = aes(x = x, y = y), col = "lightgrey", size = 2, inherit.aes = F) +
  geom_point(size = 3) +
  geom_path(data = my.plot$df_ellipse, aes(x = x, y = y), show.legend = FALSE) +
  scale_color_manual(values = c("black", "red", "blue")) +
  theme_classic(base_size = 10) +
  theme(
     legend.box = "horizontal",
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.75,0.15),
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  labs(x = "PCoA1 (11.4%)", y = "PCoA2 (7.6%)")

d <- grid.arrange(a,b,c, ncol = 3)

ggsave("Figures/Ordination.jpeg", d, height = 4, width = 15, units = "in", dpi = 600)

# Ordination - function ####
### Stats ####

# ps.hell.fun <- transform(ps.fun.rare, "hellinger")
# 
# PCoA_hell <- ordinate(physeq = ps.hell.fun, method = "PCoA", distance = "euclidean")
PCoA_bray <- ordinate(physeq = ps.fun.rare66, method = "PCoA", distance = "bray")

#Dist.hell.fun <- phyloseq::distance(ps.hell.fun, method = "euclidean", type = "samples")
Dist.fun <- phyloseq::distance(ps.fun.rare66, method = "bray", type = "samples")

#bpart <- beta.pair.abund(Dist.hell.fun, index.family = "bray")

set.seed(620)
adonis2(Dist.fun ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.fun.rare66), "data.frame"), by = "margin", permutations = 9999)  


# set.seed(620)
# (adon.nest <- adonis2(bpart$beta.bray.gra ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999))  #nestedness
 
# set.seed(620)
# (adon.turn <- adonis2(bpart$beta.bray.bal ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999)) # turnover

set.seed(620)
adonis2(Dist.fun ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.fun.rare66), "data.frame"), by = "margin", permutations = 9999)

# set.seed(620)
# adonis2(bpart$beta.bray.gra ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) #nestedness
# 
# set.seed(620)
# adonis2(bpart$beta.bray.bal ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) # turnover
# # significant turnover in beta diversity between watering treatmnets, marginally signficant turnover between plant host treatments

# set.seed(620)
# pairwise <- permanova_pairwise(
#   Dist.hell.fun,
#   as(sample_data(ps.hell.fun), "data.frame")$WaterTrt,
#   permutations = 9999,
#   padj = "BH"
# )
# 
# set.seed(620)
# pairwise <- permanova_pairwise(
#   Dist.hell.fun,
#   as(sample_data(ps.hell.fun), "data.frame")$PlantTrt,
#   permutations = 9999,
#   padj = "BH"
# )

### Figure ####
PCoA_hell.fun <- ordinate(physeq = ps.hell.fun, method = "PCoA", distance = "euclidean")
PCoA_bray.fun <- ordinate(physeq = ps.fun.rare, method = "PCoA", distance = "bray")

plot_ordination(ps.fun.rare, PCoA_bray.fun, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
     panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_point(size = 2) +
  stat_ellipse() +
  facet_wrap(~WaterTrt) + 
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

# differences between watering treatments
a <- plot_ordination(ps.fun.rare, PCoA_bray.fun, color = "WaterTrt") +
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
  lims(x = c(-0.18, 0.2)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

# differences between plant treatments
b <- plot_ordination(ps.hell.fun, PCoA_hell.fun, color = "PlantTrt") +
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
  lims(x = c(-0.18, 0.2)) +
  geom_point(size = 3) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

# get centroids plots
C <- betadisper(Dist.hell.fun, as(sample_data(ps.hell.fun), "data.frame")$plant.water)

test <- data.frame(Category=rownames(C$centroids[,1:2]), C$centroids[,1:2])

forb <- test$Category[1:3]
grass <- test$Category[4:6]
drought <- test$Category[c(2,5,8)]
water <- test$Category[c(3,6,9)]
# 
# test$PlantTrt <- ifelse(test$Category %in% forb, "natives", ifelse(test$Category %in% grass, "invasives", "mix"))
# 
# test$WaterTrt <- ifelse(test$Category %in% drought, "drought", ifelse(test$Category %in% water, "watered", "control"))

my.plot <- gg_ordiplot(C, groups = as(sample_data(ps.hell.fun), "data.frame")$plant.water, pt.size = 0, kind = "se") 

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
     legend.box = "horizontal",
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.7,0.9),
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3),
     legend.text = element_text(size = 14),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18)) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

d <- ggarrange(a,b,c, ncol = 3, labels = c("(a)", "(b)", "(c)"))

ggsave("Figures/Ordination-fun.jpeg", d, height = 4, width = 15, units = "in", dpi = 600)

# Ordination - Taxa ####
# ps.hell.taxa <- transform(ps.genus.rare, "hellinger")
# 
# PCoA_hell <- ordinate(physeq = ps.hell.taxa, method = "PCoA", distance = "euclidean")
# PCoA_bray <- ordinate(physeq = ps.genus.rare, method = "PCoA", distance = "bray")
# 
# Dist.hell.taxa <- phyloseq::distance(ps.hell.taxa, method = "euclidean", type = "samples")
Dist.taxa <- phyloseq::distance(ps.genus.rare66, method = "bray", type = "samples")

### Stats ####

#bpart <- beta.pair.abund(Dist.hell.taxa, index.family = "bray")

set.seed(620)
adonis2(Dist.taxa ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.genus.rare66), "data.frame"), by = "margin", permutations = 9999)

# set.seed(620)
# adonis2(bpart$beta.bray.gra ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) #nestedness
# 
# set.seed(620)
# adonis2(bpart$beta.bray.bal ~ WaterTrt * PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) #turnover
# # marginally sig interaction effect on turnover between taxonomic differences

set.seed(620)
adonis2(Dist.taxa ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.genus.rare66), "data.frame"), by = "margin", permutations = 9999)

# set.seed(620)
# adonis2(bpart$beta.bray.gra ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) # nestedness
# 
# set.seed(620)
# adonis2(bpart$beta.bray.bal ~ WaterTrt + PlantTrt + Plot, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) #turnover
# 
# set.seed(620)
# pairwise <- permanova_pairwise(
#   bpart$beta.bray.bal,
#   as(sample_data(ps.hell.taxa), "data.frame")$plant.water,
#   permutations = 9999,
#   padj = "BH"
# )
# 
# set.seed(620)
# pairwise <- permanova_pairwise(
#   Dist.hell.taxa,
#   as(sample_data(ps.hell.taxa), "data.frame")$PlantTrt,
#   permutations = 9999,
#   padj = "BH"
# )

### Figure ####
#PCoA_hell.taxa <- ordinate(physeq = ps.hell.taxa, method = "PCoA", distance = "euclidean")
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
  lims(x = c(-0.3, 0.25)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

# differences between plant treatments
b <- plot_ordination(ps.hell.taxa, PCoA_hell.taxa, color = "PlantTrt") +
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
  lims(x = c(-0.3, 0.25)) +
  geom_point(size = 3) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

# get centroids plots
C <- betadisper(Dist.hell.taxa, as(sample_data(ps.hell.taxa), "data.frame")$plant.water)

test <- data.frame(Category=rownames(C$centroids[,1:2]), C$centroids[,1:2])

forb <- test$Category[1:3]
grass <- test$Category[4:6]
drought <- test$Category[c(2,5,8)]
water <- test$Category[c(3,6,9)]
# 
# test$PlantTrt <- ifelse(test$Category %in% forb, "natives", ifelse(test$Category %in% grass, "invasives", "mix"))
# 
# test$WaterTrt <- ifelse(test$Category %in% drought, "drought", ifelse(test$Category %in% water, "watered", "control"))

my.plot <- gg_ordiplot(C, groups = as(sample_data(ps.hell.taxa), "data.frame")$plant.water, pt.size = 0, kind = "se") 

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
     legend.box = "horizontal",
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.7,0.9),
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3),
     legend.text = element_text(size = 14),
     axis.text = element_text(size = 17),
     axis.title = element_text(size = 18)) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

d <- ggarrange(a,b,c, ncol = 3, labels = c("(a)", "(b)", "(c)"))


ggsave("Figures/Ordination-taxa.jpeg", d, height = 4, width = 15, units = "in", dpi = 600)

# Ordination - Forbs ####
### Forbs fun ####
ps.forb <- prune_samples(sample_data(ps.fun.rare66)$PlotSubTreatmentName == "native mix 1" | sample_data(ps.fun.rare66)$PlotSubTreatmentName == "native mix 2", ps.fun.rare66)

ps.forb <- prune_samples(sample_data(ps.forb)$SampleNumber !=  241, ps.forb)

ps.forb <- prune_samples(sample_data(ps.forb)$SampleNumber !=  345, ps.forb)

PCoA_forb <- ordinate(physeq = ps.forb, method = "PCoA", distance = "bray")

plot_ordination(ps.forb, PCoA_forb, color = "PlotSubTreatmentName", shape = "PlotSubTreatmentName", label = "SampleNumber") +
  theme_classic(base_size = 15) +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  geom_point(size = 2) +
  stat_ellipse()

 # Get distances for stats
Dist.forb <- phyloseq::distance(ps.forb, method = "bray", type = "samples")

set.seed(620)

adonis2(Dist.forb ~ WaterTrt + PlotSubTreatmentName + Plot, as(sample_data(ps.forb), "data.frame"), by = "margin", permutations = 9999) 

### Mix fun ####
ps.hell.mix <- prune_samples(sample_data(ps.hell.fun)$PlotSubTreatmentName == "native mix 1 + invasives" | sample_data(ps.hell.fun)$PlotSubTreatmentName == "native mix 2 + invasives", ps.hell.fun)

PCoA_hell_mix <- ordinate(physeq = ps.hell.mix, method = "PCoA", distance = "euclidean")

plot_ordination(ps.hell.mix, PCoA_hell_mix, color = "PlotSubTreatmentName", shape = "PlotSubTreatmentName") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title = element_blank()
  ) +
  geom_point(size = 2) +
  stat_ellipse()
  

# Get distances for stats
Dist.hell.mix <- phyloseq::distance(ps.hell.mix, method = "euclidean", type = "samples")

set.seed(620)
adonis2(Dist.hell.mix ~ WaterTrt + PlotSubTreatmentName + Plot, as(sample_data(ps.hell.mix), "data.frame"), by = "margin", permutations = 9999) 

### Forbs taxa ####
ps.hell.forb <- prune_samples(sample_data(ps.hell.taxa)$PlotSubTreatmentName == "native mix 1" | sample_data(ps.hell.taxa)$PlotSubTreatmentName == "native mix 2", ps.hell.taxa)

PCoA_hell_forb <- ordinate(physeq = ps.hell.forb, method = "PCoA", distance = "euclidean")

plot_ordination(ps.hell.forb, PCoA_hell_forb, color = "PlotSubTreatmentName", shape = "PlotSubTreatmentName") +
  theme_classic(base_size = 15) +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  geom_point(size = 2) +
  stat_ellipse()

 # Get distances for stats
Dist.hell.forb <- phyloseq::distance(ps.hell.forb, method = "euclidean", type = "samples")

set.seed(620)

adonis2(Dist.hell.forb ~ WaterTrt + PlotSubTreatmentName + Plot, as(sample_data(ps.hell.forb), "data.frame"), by = "margin", permutations = 9999) 

### Mix taxa ####
ps.hell.mix <- prune_samples(sample_data(ps.hell.taxa)$PlotSubTreatmentName == "native mix 1 + invasives" | sample_data(ps.hell.taxa)$PlotSubTreatmentName == "native mix 2 + invasives", ps.hell.taxa)

PCoA_hell_mix <- ordinate(physeq = ps.hell.mix, method = "PCoA", distance = "euclidean")

plot_ordination(ps.hell.mix, PCoA_hell_mix, color = "PlotSubTreatmentName", shape = "PlotSubTreatmentName") +
  theme_classic(base_size = 15) +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    legend.title = element_blank()
  ) +
  geom_point(size = 2) +
  stat_ellipse()
  

# Get distances for stats
Dist.hell.mix <- phyloseq::distance(ps.hell.mix, method = "euclidean", type = "samples")

set.seed(620)

adonis2(Dist.hell.mix ~ WaterTrt + PlotSubTreatmentName + Plot, as(sample_data(ps.hell.mix), "data.frame"), by = "margin", permutations = 9999) 


# Soil ####
soil <- read.csv("Data/Final-Data/Soil-CN.csv")
soil$perc.N <- (soil$Total.N..µg.*0.001)/soil$Sample.Weight..mg..from.Sample.List*100
soil$perc.C <- (soil$Total.C..µg.*0.001)/soil$Sample.Weight..mg..from.Sample.List*100
soil <- merge(soil, metadata, by.x = "Sample.ID", by.y = "PlotXSubTreatment")

soil$CN <- soil$Total.C..µg./soil$Total.N..µg.

soil <- filter(soil, Sample.ID != "27B") # notes that plot wasnt weeded

colnames(soil)[3] <- "C.ug"
colnames(soil)[6] <- "N.ug"

ggplot(soil, aes(x = PlantTrt, y = CN, group = interaction(WaterTrt, PlantTrt), fill = WaterTrt)) +
  geom_boxplot() 

ggplot(soil, aes(x = PlantTrt, y = CN)) +
  geom_boxplot() 

ggplot(soil, aes(x = WaterTrt, y = perc.N, group = interaction(WaterTrt, PlantTrt), fill = PlantTrt)) +
  geom_boxplot() 

ggplot(soil, aes(x = WaterTrt, y = perc.C, group = interaction(WaterTrt, PlantTrt), fill = PlantTrt)) +
  geom_boxplot() 

ggplot(soil, aes(x = perc.N, y = perc.C)) +
  geom_point() +
  geom_smooth(method = "lm")
         
ggplot(soil, aes(x = WaterTrt, y = perc.N/100, group = interaction(WaterTrt, PlantTrt), color = PlantTrt)) +
  geom_boxplot() +
  theme_classic() +
  theme(
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 12)
  ) +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(y = "Percent Nitrogen") +
  scale_y_continuous(labels = scales::percent)

ggsave("Figures/soil.jpeg", soil.fig, dpi = 600, units = "in", width = 6, height = 4)

library(lme4)
library(lmerTest)

# no interactive effects
hist(sqrt(soil$CN))
m1 <- lmer(CN ~ PlantTrt * WaterTrt + (1|Plot), data = soil)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1))
qqline(resid(m1), col = 2, lwd = 2, lty = 2)
anova(m1)

m1 <- lmer(CN ~ PlantTrt + WaterTrt + (1|Plot), data = soil)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1))
qqline(resid(m1), col = 2, lwd = 2, lty = 2)
summary(m1)
anova(m1) # effect of plants

pairs(emmeans(m1, ~PlantTrt))

m2 <- lmer(perc.C ~ PlantTrt * WaterTrt + (1|Plot), data = soil)
plot(fitted(m2), resid(m2))
qqnorm(resid(m2))
qqline(resid(m2), col = 2, lwd = 2, lty = 2)
anova(m2) # nothing

m2 <- lmer(perc.C ~ PlantTrt + WaterTrt + (1|Plot), data = soil)
plot(fitted(m2), resid(m2))
qqnorm(resid(m2))
qqline(resid(m2), col = 2, lwd = 2, lty = 2)
summary(m2)
anova(m2) # nothing

m3 <- lmer(asin(sqrt(perc.N)) ~ PlantTrt * WaterTrt + (1|Plot), data = soil)
plot(fitted(m3), resid(m3))
qqnorm(resid(m3))
qqline(resid(m3), col = 2, lwd = 2, lty = 2)
anova(m3)

m3 <- lmer(asin(sqrt(perc.N)) ~ PlantTrt + WaterTrt + (1|Plot), data = soil)
plot(fitted(m3), resid(m3))
qqnorm(resid(m3))
qqline(resid(m3), col = 2, lwd = 2, lty = 2)
summary(m3)
anova(m3)

# Biomass ####
calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

bio <- read.csv("Data/Final-Data/Plant-Biomass.csv")

#bio.w <- pivot_wider(bio, names_from = "Species", values_from = "Weight.g")

bio <- merge(bio, metadata.nocontrols[,c(1,6,7,10,12,13,15)], by.x = c("Plot", "Subplot"), by.y = c("Plot", "PlotSubTreatment"))

#write.csv(bio, "Data/bio-meta-comparison.csv", row.names = F)
bio.dups <- bio[duplicated(bio[,1:4]),] # 1 duplicate



ggplot(bio, aes(x = WaterTrt, y = Weight.g, fill = interaction(PlantTrt, FunGroup))) +
  geom_boxplot() 

bio.sum <- bio %>%
  dplyr::group_by(PlantTrt, WaterTrt, FunGroup, Plot, Subplot) %>%
  dplyr::summarize(sum.weight = sum(Weight.g))  %>%
  group_by(PlantTrt, WaterTrt, FunGroup) %>%
  dplyr::summarize(avg.weight = mean(sum.weight), se.weight = calcSE(sum.weight), n = length(unique(paste0(Plot,Subplot))))

ggplot(bio.sum, aes(x = WaterTrt, y = avg.weight, fill = FunGroup)) +
  geom_bar(stat = "identity", position = "dodge", col = "black") + 
  geom_errorbar(aes(ymin = avg.weight - se.weight, ymax = avg.weight + se.weight), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~PlantTrt)

bio.sum <- bio %>%
  dplyr::group_by(PlantTrt, WaterTrt, FunGroup, Plot, Subplot) %>%
  dplyr::summarize(sum.weight = sum(Weight.g))  %>%
  group_by(PlantTrt, WaterTrt) %>%
  dplyr::summarize(avg.weight = mean(sum.weight), se.weight = calcSE(sum.weight), n = length(sum.weight))

ggplot(bio.sum, aes(x = PlantTrt, y = avg.weight, fill = WaterTrt)) +
  geom_bar(stat = "identity", position = "dodge", col = "black") + 
  geom_errorbar(aes(ymin = avg.weight - se.weight, ymax = avg.weight + se.weight), width = 0.2, position = position_dodge(0.9)) 

ggplot(bio.sum, aes(x = WaterTrt, y = avg.weight, fill = PlantTrt)) +
  geom_bar(stat = "identity", position = "dodge", col = "black") + 
  geom_errorbar(aes(ymin = avg.weight - se.weight, ymax = avg.weight + se.weight), width = 0.2, position = position_dodge(0.9)) 

# Latlong BS ####
# distancy decay code
library(geosphere)
library(vegan)
library(microbiome)
library(phyloseq)

# phyloseq obj
ps

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

