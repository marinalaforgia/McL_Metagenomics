### Script to analyze Metagenomic soil samples from McLaughlin ###

# Files needed (all are found within the R-Data-Files folder on box)
## 1. Kallisto file of raw gene counts
## 2. Functional pathways file
## 3. Joined MMSeq and Kaiju taxonomy file
## 4. Metadata file

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


# Read in data ####
#count.gene <- readRDS("Data/Final-Data/txi.kallisto.08042021.RDS")
fun.gene <- read.csv("Data/Final-Data/gl_fun_only_keep_alterate.csv")
tax.gene <- read.delim("Data/Final-Data/gl_tax_filtset_genes_to_remove.tsv")
metadata <- read.csv("Data/Final-Data/McL_metag_metadata.csv")

# Prepare data ####

## Gene count prep ####
colnames(count.gene$counts) <- sub("_S.*", "", colnames(count.gene$counts))
colnames(count.gene$counts) <- sub("-", "_", colnames(count.gene$counts))
colnames(count.gene$abundance) <- sub("_S.*", "", colnames(count.gene$abundance))
colnames(count.gene$abundance) <- sub("-", "_", colnames(count.gene$abundance))

## Taxa prep ####
# domain.y = Kaiju, domain.x = mmseqs
colnames(tax.gene)[7] <- "Domain.mmseqs"
colnames(tax.gene)[8] <- "Domain.kaiju"
tax.gene <- filter(tax.gene, Joined_Domain != "Remove")
row.names(tax.gene) <- tax.gene$gene
tax.gene <- tax.gene[order(tax.gene$gene),]
tax.gene <- tax.gene[,-1]
tax.gene <- as.matrix(tax.gene)

## Function prep ####
fun.gene$gene <- gsub("^.*\\_", "", fun.gene$gene)
row.names(fun.gene) <- fun.gene$gene
fun.gene <- fun.gene[order(fun.gene$gene),]
fun.gene <- fun.gene[,-1]
fun.gene <- as.matrix(fun.gene)

## Metadata prep ####

# removes controls we didn't sequence: GLM_0698, GLM_0696, GLM_0694, GLM_0693, GLM_0392, GLM_0391, GLM_0390
to.rm <- c("GLM_0698", "GLM_0696", "GLM_0694", "GLM_0693", "GLM_0392", "GLM_0391", "GLM_0390")
metadata <- metadata[!(metadata$SampleID %in% to.rm),] # sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID
metadata$plant.water <- paste(metadata$PlantTrt, metadata$WaterTrt, sep = ".")

metadata$PlotSubTreatmentName <- recode_factor(metadata$PlotSubTreatmentName, ST_forbs_X_grasses = "forb mix 1 + grasses", SA_forbs_X_grasses = "forb mix 2 + grasses", Stress_avoiding_forbs = "forb mix 2", Stress_tolerant_forbs = "forb mix 1")

#colnames(count.gene[[1]]) == metadata[,"SampleID"] # make sure everything is in the right order

controls <- c("GLM_0697", "GLM_0695", "GLM_0393")
metadata.nocontrols <- metadata[!(metadata$SampleID %in% controls),]

metadata.nocontrols$PlantTrt <- recode_factor(metadata.nocontrols$PlantTrt, forb = "forbs",  mix = "forbs & grasses")
metadata.nocontrols$PlantTrt <- factor(metadata.nocontrols$PlantTrt, levels = c("forbs", "grasses", "forbs & grasses"))

metadata.nocontrols$WaterTrt <- factor(metadata.nocontrols$WaterTrt, levels = c("Drought", "Control", "Watering"))

metadata.nocontrols$PlotSubTreatmentName <- factor(metadata.nocontrols$PlotSubTreatmentName, levels = c("forb mix 1", "forb mix 2", "Invasive_grasses", "forb mix 1 + grasses", "forb mix 2 + grasses"))

## Filter data ####

ddsTxi <- DESeqDataSetFromTximport(count.gene,
                                   colData = metadata,
                                   design = ~ plant.water)

# Remove rows with no counts and genes not seen more than 3 times in at least 20% of the samples. This protects against genes with small mean & trivially large C.V
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi)) > 1, ] # remove rows with no counts
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) > 3) >= (0.2*ncol(counts(ddsTxi))), ] 

## Decontam Counts ####
contam_counts <- counts(ddsTxi, normalized = F) # object for decontam

ps <- phyloseq(otu_table(contam_counts, taxa_are_rows = T), sample_data(metadata))

#first tell it which samples are the NC
sample_data(ps)$is.neg <- sample_data(ps)$SampleType == "Kit"

#all sequences that are more prevalent in negative controls than in positive samples
contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev05$contaminant) # 53 contaminants 

#Make phyloseq object of presence-absence in negative controls
ps.neg <- prune_samples(sample_data(ps)$is.neg == "TRUE", ps)
ps.neg.presence <- transform_sample_counts(ps.neg, function(abund) 1*(abund>0))

#Make phyloseq object of presence-absence in true positive samples
ps.pos <- prune_samples(sample_data(ps)$is.neg == "FALSE", ps)
ps.pos.presence <- transform_sample_counts(ps.pos, function(abund) 1*(abund>0))

#Make data.frame of prevalence in positive and negative samples
#using prev threshold = 0.5
df.pres <- data.frame(prevalence.pos = taxa_sums(ps.pos.presence),
                      prevalence.neg = taxa_sums(ps.neg.presence),
                      contam.prev = contamdf.prev05$contaminant)

ggplot(df.pres, aes(x = prevalence.neg, y = prevalence.pos, color = contam.prev)) +
  geom_point()

contam.rows <- which(contamdf.prev05$contaminant)
contaminants <- rownames(contamdf.prev05[contam.rows,])

#returns list of all taxa, except contaminants
decontam_counts <- contam_counts[!(rownames(contam_counts) %in% contaminants),]

decontam_counts_nocontrols <- decontam_counts[, !(colnames(decontam_counts) %in% controls)] # get rid of controls

saveRDS(decontam_counts_nocontrols, "Data/Final-Data/decontam_non-norm_counts_nocontrols.RDS")

rm(ps, ps.neg, ps.neg.presence, ps.pos, ps.pos.presence, ddsTxi)

## Sourcetracker Prep ####
dds <- DESeqDataSetFromMatrix(decontam_counts,
                              colData = metadata,
                              design = ~ plant.water)


dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = T) 

saveRDS(normalized_counts, "Data/Final-Data/decontam_norm_counts_w-controls.RDS") # output for Source Tracker

## DESeq2 genes ####
decontam_counts_nocontrols <- readRDS("Data/Final-Data/decontam_non-norm_counts_nocontrols.RDS")

# filter out unwanted genes (ribosomes, plant associated etc)
decontam_counts_nocontrols <- decontam_counts_nocontrols[rownames(decontam_counts_nocontrols) %in% rownames(tax.gene),]

decontam_counts_nocontrols <- decontam_counts_nocontrols[rownames(decontam_counts_nocontrols) %in% rownames(fun.gene),]

dds <- DESeqDataSetFromMatrix(decontam_counts_nocontrols,
                                   colData = metadata.nocontrols,
                                   design = ~ plant.water)

dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

plotDispEsts(dds)

norm_counts <- counts(dds, normalized = T) 

## DESeq2 functions ####

ps <- phyloseq(otu_table(norm_counts, taxa_are_rows = T), sample_data(metadata.nocontrols), tax_table(fun.gene))

## DESeq2 Genus ####

# Ordination ####
## Genes ####
### Prep ####

ps <- phyloseq(otu_table(norm_counts, taxa_are_rows = T), sample_data(metadata.nocontrols), tax_table(tax.gene))

ps.hell <- transform(ps, "hellinger")

PCoA_hell <- ordinate(physeq = ps.hell, method = "PCoA", distance = "euclidean")

### Figure ####
plot_ordination(ps.hell, PCoA_hell, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
     panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_point(size = 2) +
  stat_ellipse() +
  facet_wrap(~WaterTrt) + 
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

plot_ordination(ps.hell, PCoA_hell, color = "WaterTrt") +
  theme_classic(base_size = 15) +
  theme(
     panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_point(size = 2) +
  stat_ellipse() +
  facet_wrap(~PlantTrt) + 
  scale_color_manual(values = c("red", "black", "blue"))


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

###Stats ####
# Get distances for stats
Dist.hell <- phyloseq::distance(ps.hell, method = "euclidean", type = "samples")

set.seed(50)

adonis2(Dist.hell ~ WaterTrt * PlantTrt, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999)  # interaction not significant

adonis2(Dist.hell ~ PlantTrt + WaterTrt, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) 

C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$plant.water)

set.seed(50)
permutest(C, permutations = 9999) # no sig dif

C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$PlantTrt)

set.seed(50)
permutest(C, permutations = 9999) # no sig dif

C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$WaterTrt)

set.seed(50)
permutest(C, permutations = 9999) # no sig dif

pairwise <- permanova_pairwise(
  Dist.hell,
  as(sample_data(ps.hell), "data.frame")$WaterTrt,
  permutations = 9999,
  padj = "BH"
)


## Taxa ####
ps.taxa <- tax_glom(ps, taxrank = "Genus")

ps.hell.taxa <- transform(ps.taxa, "hellinger")

### Figure ####
PCoA_hell.taxa <- ordinate(physeq = ps.hell.taxa, method = "PCoA", distance = "euclidean")

plot_ordination(ps.hell.taxa, PCoA_hell.taxa, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
     panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_point(size = 2) +
  stat_ellipse() +
  facet_wrap(~WaterTrt) + 
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

### Stats ####
Dist.hell.taxa <- phyloseq::distance(ps.hell.taxa, method = "euclidean", type = "samples")

set.seed(50)

adonis2(Dist.hell.taxa ~ WaterTrt * PlantTrt, as(sample_data(ps.hell.taxa), "data.frame"), by = "margin", permutations = 9999)  # interaction not significant

adonis2(Dist.hell.taxa ~ WaterTrt + PlantTrt, as(sample_data(ps.hell.taxa), "data.frame"), by = "margin", permutations = 9999) # both watering and plant treatment sig

pairwise <- permanova_pairwise(
  Dist.hell.taxa,
  as(sample_data(ps.hell.taxa), "data.frame")$PlantTrt,
  permutations = 9999,
  padj = "BH"
)

pairwise <- permanova_pairwise(
  Dist.hell.taxa,
  as(sample_data(ps.hell.taxa), "data.frame")$WaterTrt,
  permutations = 9999,
  padj = "BH"
)

## Forb Treatments ####
### Forbs ####
ps.hell.forb <- prune_samples(sample_data(ps.hell)$PlotSubTreatmentName == "forb mix 1" | sample_data(ps.hell)$PlotSubTreatmentName == "forb mix 2", ps.hell)

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

set.seed(50)

adonis2(Dist.hell.forb ~ PlotSubTreatmentName, as(sample_data(ps.hell.forb), "data.frame"), by = "margin", permutations = 9999) 

### Mix ####
ps.hell.mix <- prune_samples(sample_data(ps.hell)$PlotSubTreatmentName == "forb mix 1 + grasses" | sample_data(ps.hell)$PlotSubTreatmentName == "forb mix 2 + grasses", ps.hell)

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

set.seed(50)

adonis2(Dist.hell.mix ~ PlotSubTreatmentName, as(sample_data(ps.hell.mix), "data.frame"), by = "margin", permutations = 9999) 

# Genetic alpha diversity ####
GL_Alpha <- estimate_richness(ps.hell, measures = "Shannon")

GL_Alpha <- cbind(GL_Alpha, sample_data(ps.hell))

set.seed(50)
kruskal_test(Shannon ~ PlantTrt, distribution = approximate(nresample = 9999), data = GL_Alpha)

kruskal_test(Shannon ~ WaterTrt, distribution = approximate(nresample = 9999), data = GL_Alpha)

shannon.sum <- ddply(GL_Alpha, .(PlantTrt, WaterTrt), summarise, max = max(Shannon))

plot_richness(ps.hell, measures = "Shannon", x = "PlantTrt", color = "PlantTrt") + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  theme_bw(base_size = 15) +
  facet_wrap(~WaterTrt) + 
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

# Sig Diff Genes ####

contrasts <- c("FC.FW", 
               "FC.FD", 
               "FC.MD", 
               "FC.MW", 
               "GC.GD", 
               "GC.GW", 
               "GC.MD",
               "GC.MW",
               "GC.FC", 
               "GD.FD",
               "GW.FW",
               "FD.MD",
               "FW.MW")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forbs.Control/forb.Watering) so if something is higher in forb controls, the number is greater than one, so the log fold change > 0, if something is higher in the watering treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FC.FW = c("plant.water","forb.Control","forb.Watering"),
                      FC.FD = c("plant.water","forb.Control","forb.Drought"),
                      FC.MD = c("plant.water","forb.Control","mix.Drought"),
                      FC.MW = c("plant.water","forb.Control","mix.Watering"),
                      GC.GD = c("plant.water","grasses.Control","grasses.Drought"),
                      GC.GW = c("plant.water","grasses.Control","grasses.Watering"),
                      GC.MD = c("plant.water","grasses.Control","mix.Drought"),
                      GC.MW = c("plant.water","grasses.Control","mix.Watering"),
                      GC.FC = c("plant.water","grasses.Control","forb.Control"),
                      GD.FD = c("plant.water","grasses.Drought","forb.Drought"),
                      GW.FW = c("plant.water","grasses.Watering","forb.Watering"),
                      FD.MD = c("plant.water","forb.Drought","mix.Drought"),
                      FW.MW = c("plant.water", "forb.Watering", "mix.Watering"))

plot.name.list <- list(FC.FW = "Forb Control v Forb Watering",
                       FC.FD = "Forb Control v Forb Drought",
                       FC.MD = "Forb Control v Mix Drought",
                       FC.MW = "Forb Control v Mix Watering",
                       GC.GD = "Grass Control v Grass Drought",
                       GC.GW = "Grass Control v Grass Watering",
                       GC.MD = "Grass Control v Mix Drought",
                       GC.MW  = "Grass Control v Mix Watering",
                       GC.FC = "Grass Control v Forb Control",
                       GD.FD = "Grass Drought v Forb Drought",
                       GW.FW = "Grass Watering v Forb Watering",
                       FD.MD = "Forb Drought v Mix Drought",
                       FW.MW = "Forb Watering v Mix Watering")

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
length(unique(sig.genes))
## BH: 20 unique genes

#plot results
grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], ncol = 3)

df <- data.frame()

for(i in names(plot.list)){ #plot list is the one with all the significantly different ones
  tmp <- cbind(pair = i, res.list[[i]])
  df <- rbind(df, tmp)
}

nrow(df[df$pair == "FC.MD",]) #13
nrow(df[df$pair == "FC.MW",]) #1
nrow(df[df$pair == "FD.MD",]) #10

# Functional analysis ####
fun.gene.sig <- fun.gene[fun.gene$gene %in% sig.genes,] # 16 genes have functions!
no.fun <- unique(sig.genes[!sig.genes %in% fun.gene$gene]) #4 have no discernible function

write.csv(fun.gene.sig, "Data/Final-Data/20221205_fun-gene-cog-allsigs.csv", row.names = F)

#Edited to reflect more streamlined COG categories
fun.gene.cog <- read.csv("Data/Final-Data/20221205_fun-gene-cog-allsigs_ML-CE-Edit.csv")

length(unique(fun.gene.cog$gene))

unique(fun.gene.cog$Fun_Cat_CE_ML) # 13 unique functions yay!

for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], fun.gene.cog[,c(1,2,14)], by = "gene")
}

sig.dif.fun <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.fun[[i]] <- ddply(res.list[[i]], .(Fun_Cat_CE_ML), summarize, estimate = sum(estimate))
  }
}

## Summed by Function ####

names(sig.dif.fun)[[1]] <- "forbs (control) vs. mixes (drought)"
names(sig.dif.fun)[[2]] <- "forbs (control) vs. mixes (watered)"
names(sig.dif.fun)[[3]] <- "forbs (drought) vs. mixes (drought)"

plot.list2 <- list()

for(i in names(sig.dif.fun)){
    
    sig.dif.fun[[i]]$Fun_Cat_CE_ML <- factor(sig.dif.fun[[i]]$Fun_Cat_CE_ML, levels = sig.dif.fun[[i]]$Fun_Cat_CE_ML[order(sig.dif.fun[[i]]$estimate)])
    
    p <- ggplot(sig.dif.fun[[i]], aes(x = Fun_Cat_CE_ML, y = estimate)) + 
      geom_point(size = 2.5) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 12)
          ) +
      ylim(-26,10) +
      labs(title = i, y = expression(paste(log[2], " fold change"))) +
      scale_x_discrete(labels = c("Posttranslational modification, protein turnover, chaperones" = "Posttranslational modification, protein\nturnover, chaperones", "Translation, ribosomal structure & biogenesis" = "Translation, ribosomal structure\n& biogenesis")) +
      coord_flip()
  
   plot.list2[[i]] = p
}

plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[3]], plot.list2[[2]], align = "v")

p <- grid.arrange(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], ncol = 1, heights = c(0.6,0.5,0.3))

ggsave("Figures/change-genes2.jpeg", p, units = "in", width = 7, height = 10, dpi = 600)


sig.dif.fundf <- data.frame()

for(i in 1:length(names(sig.dif.fun))) {
  sig.dif.fundf <- rbind(cbind(plant.water = names(sig.dif.fun)[i],sig.dif.fun[[i]]), sig.dif.fundf)
}

## Colored by Gene ####
show_col(viridis_pal(option = "turbo")(12))
viridis_pal(option = "turbo")(12)

color_map <- c("Amino acid transport & metabolism" = "#30123BFF", 
               "Carbohydrate transport & metabolism" = "#4454C4FF", 
               "Cell wall/membrane/envelope biogenesis" = "#4490FEFF",
               "Defense mechanisms" = "#1FC8DEFF", 
               "Energy production & conversion" = "#29EFA2FF",
               "General function prediction only" = "#7DFF56FF",
               "Mobilome: prophages, transposons" = "#C1F334FF",
               "Posttranslational modification, protein turnover, chaperones" = "#F1CA3AFF",
               "Replication, recombination & repair" = "#FE922AFF",
               "Signal transduction mechanisms" = "#EA4F0DFF",
               "Symbiosis" = "#BE2102FF",
               "Translation, ribosomal structure & biogenesis" = "#7A0403FF")

sig.dif.fundf$Fun_Cat_CE_ML <- as.character(sig.dif.fundf$Fun_Cat_CE_ML)
sig.dif.fundf <- sig.dif.fundf[order(sig.dif.fundf$Fun_Cat_CE_ML),]

leg_plot <- ggplot(sig.dif.fundf, aes(x = Fun_Cat_CE_ML, y = estimate, col = Fun_Cat_CE_ML)) +
  theme_classic() +
  geom_point() +
  scale_color_manual(values = color_map)

legend <- g_legend(leg_plot)

plot.list3 <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
    
 p <- ggplot(res.list[[i]], aes(x = gene, y = estimate, col = Fun_Cat_CE_ML)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      #axis.text.y = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.title.x = element_blank(),
        ) +
    labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change"))) +
   ylim(-7,5) +
   scale_color_manual(values = color_map, drop = F) +
   coord_flip()
  
   plot.list3[[i]] = p
  }
}

p <- grid.arrange(plot.list3[[1]], 
             plot.list3[[3]], plot.list3[[2]], legend, ncol = 2)

ggsave("Figures/Gene_change.jpeg", p, height = 10, width = 10, units = "in", dpi = 600)

# Soil ####
soil <- read.csv("Data/Final-Data/Soil-CN.csv")
soil$perc.N <- (soil$Total.N..µg.*0.001)/soil$Sample.Weight..mg..from.Sample.List
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

soil.fig <- ggplot(soil, aes(x = WaterTrt, y = perc.N, group = interaction(WaterTrt, PlantTrt), color = PlantTrt)) +
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
m1 <- lmer(CN ~ PlantTrt + WaterTrt + (1|Plot), data = soil)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1))
qqline(resid(m1), col = 2, lwd = 2, lty = 2)
anova(m1)

m2 <- lmer(C.ug ~ plant.water + (1|Plot), data = soil)
plot(fitted(m2), resid(m2))
qqnorm(resid(m2))
qqline(resid(m2), col = 2, lwd = 2, lty = 2)
anova(m2)

m3 <- lm(perc.N ~ PlantTrt * WaterTrt, data = soil)
plot(fitted(m3), resid(m3))
qqnorm(resid(m3))
qqline(resid(m3), col = 2, lwd = 2, lty = 2)
anova(m3)

# Biomass ####
calcSE<-function(x){
  x2<-na.omit(x)
  sd(x2)/sqrt(length(x2))
}

bio <- read.csv("Data/Final-Data/Plant-Biomass.csv")
bio <- merge(bio, metadata.nocontrols[,c(6,7,10,12,13)], by.x = c("Plot", "Subplot"), by.y = c("Plot", "PlotSubTreatment"))

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

# what if we remove the contaminated treatments?
## Total of 15 contaminated treatments (out of 66 total)
## 23A, 23B, 27B, 29A, 29B, 49A, 55B, 87A, 91A, 97A, 30B, 98B, 93E, 23E, 29E
## 

# #### DESeq on Functions ####
# counts.fun <- data.frame(cbind(gene = row.names(decontam_counts_nocontrols), decontam_counts_nocontrols)) # need non norm counts for DESeq so go back through, maybe just make a gene to reduced function sheet... that might be easier 
# 
# # counts.fun <- left_join(norm.counts.fun, fun.gene, by = "gene")
# # 
# # norm.counts.fun$COG20_CATEGORY <- ifelse(norm.counts.fun$COG20_CATEGORY == "NULL", norm.counts.fun$COG_CATEGORY, norm.counts.fun$COG20_CATEGORY)
# # 
# # unique(norm.counts.fun[!nchar(norm.counts.fun$COG20_CATEGORY) > 5, ]$COG20_CATEGORY)
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "T"),]$COG20_CATEGORY <- "Signal transduction mechanisms"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "L"),]$COG20_CATEGORY <- "Replication, recombination and repair"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "S"),]$COG20_CATEGORY <- "Function unknown"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "E"),]$COG20_CATEGORY <- "Amino acid transport and metabolism"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "G"),]$COG20_CATEGORY <- "Carbohydrate transport and metabolism"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "C"),]$COG20_CATEGORY <- "Energy production and conversion"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "P"),]$COG20_CATEGORY <- "Inorganic ion transport and metabolism"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "M"),]$COG20_CATEGORY <- "Cell wall/membrane/envelope biogenesis"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "I"),]$COG20_CATEGORY <- "Lipid transport and metabolism"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "V"),]$COG20_CATEGORY <- "Defense mechanisms"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "N"),]$COG20_CATEGORY <- "Cell motility"
# # 
# # norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "Q"),]$COG20_CATEGORY <- "Secondary metabolites biosynthesis, transport and catabolism"
# # 
# # norm.counts.fun <- norm.counts.fun[,c(1:67,77)]
# # 
# # norm.counts.fun[,2:67] <- norm.counts.fun[,2:67] %>% mutate_all(as.numeric)
# # 
# # test <- separate(norm.counts.fun, COG20_CATEGORY, into = c("newcog1", "newcog2", "newcog3", "newcog4", "newcog5", "newcog6", "newcog7", "newcog8", "newcog9"), sep = "!!!")
# 
# #write.csv(test, "Data/norm-counts-fun.csv", row.names = F, na = "")
# fun.reduced <- read.csv("Data/function-reduced.csv", na = "")
# fun.reduced$gene <- as.factor(fun.reduced$gene)
# counts.fun <- left_join(counts.fun, fun.reduced, by = "gene")
# counts.fun[,2:67] <- counts.fun[,2:67] %>% mutate_all(as.numeric)
# 
# tmp2 <- data.frame()
# 
# for(i in 1:7) {
#   tmp <- counts.fun[,c(1:67,67+i)]
#   colnames(tmp)[68] <- "COG20_CATEGORY_reduced"
#   if(i %in% 2:7) {
#     tmp <- tmp[!is.na(tmp$COG20_CATEGORY_reduced),]
#   }
#   tmp2 <- rbind(tmp2, tmp)
# } # do we want to lump COG "Function Unknown" with our own category? what to do about when a gene has both COG Function unknown and other known functions? 
# tmp2[is.na(tmp2$COG20_CATEGORY_reduced),]$COG20_CATEGORY_reduced <- "NULL"
# 
# dds.fun <- tmp2 %>%
#   group_by(COG20_CATEGORY_reduced) %>%
#   dplyr::summarise(across(2:67, sum))
# 
# dds.fun <- as.data.frame(dds.fun)
# row.names(dds.fun) <- dds.fun$COG20_CATEGORY_reduced
# dds.fun <- as.matrix(dds.fun)[,-1]
# mode(dds.fun) <- "numeric"
# 
# dds.fun <- DESeqDataSetFromMatrix(dds.fun,
#                                    colData = metadata.nocontrols,
#                                    design = ~ plant.water)
# dds.fun <- estimateSizeFactors(dds.fun)
# 
# dds.fun <- DESeq(dds.fun)
# 
# plotDispEsts(dds.fun)
# 
# contrasts <- c("FC.FW", 
#                "FC.FD", 
#                "FC.MD", 
#                "FC.MW", 
#                "GC.GD", 
#                "GC.GW", 
#                "GC.MD",
#                "GC.MW",
#                "GC.FC", 
#                "GD.FD",
#                "GW.FW",
#                "FD.MD",
#                "FW.MW")
# 
# contrast.list <- list(FC.FW = c("plant.water","forb.Control","forb.Watering"),
#                       FC.FD = c("plant.water","forb.Control","forb.Drought"),
#                       FC.MD = c("plant.water","forb.Control","mix.Drought"),
#                       FC.MW = c("plant.water","forb.Control","mix.Watering"),
#                       GC.GD = c("plant.water","grasses.Control","grasses.Drought"),
#                       GC.GW = c("plant.water","grasses.Control","grasses.Watering"),
#                       GC.MD = c("plant.water","grasses.Control","mix.Drought"),
#                       GC.MW = c("plant.water","grasses.Control","mix.Watering"),
#                       GC.FC = c("plant.water","grasses.Control","forb.Control"),
#                       GD.FD = c("plant.water","grasses.Drought","forb.Drought"),
#                       GW.FW = c("plant.water","grasses.Watering","forb.Watering"),
#                       FD.MD = c("plant.water","forb.Drought","mix.Drought"),
#                       FW.MW = c("plant.water", "forb.Watering", "mix.Watering"))
# 
# plot.name.list <- list(FC.FW = "Forb Control v Forb Watering",
#                        FC.FD = "Forb Control v Forb Drought",
#                        FC.MD = "Forb Control v Mix Drought",
#                        FC.MW = "Forb Control v Mix Watering",
#                        GC.GD = "Grass Control v Grass Drought",
#                        GC.GW = "Grass Control v Grass Watering",
#                        GC.MD = "Grass Control v Mix Drought",
#                        GC.MW  = "Grass Control v Mix Watering",
#                        GC.FC = "Grass Control v Forb Control",
#                        GD.FD = "Grass Drought v Forb Drought",
#                        GW.FW = "Grass Watering v Forb Watering",
#                        FD.MD = "Forb Drought v Mix Drought",
#                        FW.MW = "Forb Watering v Mix Watering")
# 
# alpha = 0.05
# res.list <- list()
# plot.list <- list()
# sig.genes <- c()
# 
# for(i in contrasts) {
#   #get results for each contrast
#   res <- results(dds.fun, contrast = contrast.list[[i]], pAdjustMethod = "BH")
#   #filter results by p-value
#   res.alpha <- res[which(res$padj < alpha), ]
#   #tidy results 
#   res.list[[i]] <- tidy(res.alpha)
# 
#  res.list[[i]]$gene <- factor(res.list[[i]]$gene, levels = res.list[[i]]$gene[order(res.list[[i]]$estimate)])
#   
#  if(nrow(res.list[[i]]) > 0) {
#    p <- ggplot(res.list[[i]], aes(x = gene, y = estimate)) + 
#     geom_point(size = 3) +
#     theme_bw() +
#     theme(
#       axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
#       plot.title = element_text(hjust = 0.5, size = 12),
#       legend.position = "none",
#       axis.title.x = element_blank(),
#         ) +
#     labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change")))
#   
#    plot.list[[i]] = p
#  }
#    sig.genes <- c(as.character(res.list[[i]]$gene), sig.genes)
# }
# 
