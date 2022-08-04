### Script to analyze Metagenomic soil samples from McLaughlin ###

# Files needed (all are found within the R-Data-Files folder on box)
## 1. Kallisto file of raw gene counts
## 2. Functional pathways file
## 3. Kaiju taxonomy file
## 4. Metadata file

# Outline of analysis
## 1. Get normalized counts with negative controls to send to Cassie for SourceTracker
## 2. Create DESeq object, extract non-normalized counts to feed into decontam, save decontaminated object without controls of non-normalized counts
## 3. Run DESeq on non-normalized counts, use variance stabilizing transformation 
## 4. Generate list of differentially abundant genes
## 4. Gene-level ordinations
## 6. Investigate functions of differentially abundant genes
## 5. Taxonomy ordinations (Kaiju)

rm(list=ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("microbiome", version = "3.14")

#### Load Libraries ####
library(DESeq2)
library(tidyverse)
library(phyloseq)
library(decontam)
library(vegan)
library(Rmisc)
library(biobroom)
library(gridExtra)
library(microbiome)

source("Scripts/Adonis.pair.R")

#### Read in Data ####
count.gene <- readRDS("Data/Final-Data/txi.kallisto.08042021.RDS")
fun.gene <- read.csv("Data/Final-Data/GL_function_per_gene.csv")
tax.gene <- readRDS("Data/Final-Data/GL_kaiju_tax_ready_for_source_09132021.RDS")
tax2.fun <- read_tsv("Data/function.res.top1.function.tax.tsv") 
# remove everything but bacteria and compare to Cassie's first of sig dif genes, confirm that what I sent was all genes 
# second send cassie all unclassified and we can rerun it through cassie
# compare kaiju to rescued
# subset to bacteria and rescue uncalssified 
metadata <- read.csv("Data/Final-Data/McL_metag_metadata.csv")

soil <- read.csv("Data/Final-Data/Soil-CN.csv")

#### Prepare data ####

fun.gene$gene <- gsub("^.*\\_", "", fun.gene$gene)

## Gene count prep ##
colnames(count.gene$counts) <- sub("_S.*", "", colnames(count.gene$counts))
colnames(count.gene$counts) <- sub("-", "_", colnames(count.gene$counts))
colnames(count.gene$abundance) <- sub("_S.*", "", colnames(count.gene$abundance))
colnames(count.gene$abundance) <- sub("-", "_", colnames(count.gene$abundance))

## Metadata prep ##

# removes controls we didn't sequence: GLM_0698, GLM_0696, GLM_0694, GLM_0693, GLM_0392, GLM_0391, GLM_0390
to.rm <- c("GLM_0698", "GLM_0696", "GLM_0694", "GLM_0693", "GLM_0392", "GLM_0391", "GLM_0390")
metadata <- metadata[!(metadata$SampleID %in% to.rm),] # sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID
metadata$plant.water <- paste(metadata$PlantTrt, metadata$WaterTrt, sep = ".")

colnames(count.gene[[1]]) == metadata[,"SampleID"] # make sure everything is in the right order

controls <- c("GLM_0697", "GLM_0695", "GLM_0393")
metadata.nocontrols <- metadata[!(metadata$SampleID %in% controls),]

metadata.nocontrols$PlantTrt <- recode_factor(metadata.nocontrols$PlantTrt, forb = "forbs",  mix = "forbs & grasses")
metadata.nocontrols$PlantTrt <- factor(metadata.nocontrols$PlantTrt, levels = c("forbs", "grasses", "forbs & grasses"))

metadata.nocontrols$WaterTrt <- factor(metadata.nocontrols$WaterTrt, levels = c("Drought", "Control", "Watering"))

#### Filter and prep DESeq ####
ddsTxi <- DESeqDataSetFromTximport(count.gene,
                                   colData = metadata,
                                   design = ~ plant.water)

# Remove rows with no counts and genes not seen more than 3 times in at least 20% of the samples. This protects against genes with small mean & trivially large C.V
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi)) > 1, ] # remove rows with no counts
ddsTxi <- ddsTxi[rowSums(counts(ddsTxi) > 3) >= (0.2*ncol(counts(ddsTxi))), ] 

#### Decontam Counts ####
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

#### Normalize counts for ST ####
dds <- DESeqDataSetFromMatrix(decontam_counts,
                              colData = metadata,
                              design = ~ plant.water)


dds <- estimateSizeFactors(dds)

normalized_counts <- counts(dds, normalized = T) 

saveRDS(normalized_counts, "Data/Final-Data/decontam_norm_counts_w-controls.RDS") # output for Source Tracker

#### DESeq ####
decontam_counts_nocontrols <- readRDS("Data/Final-Data/decontam_non-norm_counts_nocontrols.RDS")

dds <- DESeqDataSetFromMatrix(decontam_counts_nocontrols,
                                   colData = metadata.nocontrols,
                                   design = ~ plant.water)
dds <- estimateSizeFactors(dds)

dds <- DESeq(dds)

plotDispEsts(dds)

#### Gene Ordination Set-up ####
norm_counts <- counts(dds, normalized = T) 

ps <- phyloseq(otu_table(norm_counts, taxa_are_rows = T), sample_data(metadata.nocontrols))

#https://rdrr.io/github/microbiome/microbiome/

# CLR Transformation
ps.clr <- transform(ps, "clr")

# CLR of counts (euclidean dist of CLR = Aitchison distance)
PCoA_clr <- ordinate(physeq = ps.clr, method = "PCoA", distance = "euclidean")

#### .___Main Figure ####
plot_ordination(ps.clr, PCoA_clr, color = "PlantTrt", shape = "PlantTrt") +
  theme_classic(base_size = 15) +
  geom_point(size = 2) +
  facet_wrap(~WaterTrt) +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  stat_ellipse() +
  scale_color_brewer(palette = "Dark2")

# Get distances for stats
Dist.clr <- phyloseq::distance(ps.clr, method = "euclidean", type = "samples")

set.seed(50)
# adonis(Dist.clr ~ plant.water, as(sample_data(ps), "data.frame"), permutations = 9999) # no sig interaction
# 
# BH.adjust.plant <- adonis.pair(Dist.clr, metadata.nocontrols$PlantTrt, nper = 9999, corr.method = "BH") 
# 
# BH.adjust.water <- adonis.pair(Dist.clr, metadata.nocontrols$WaterTrt, nper = 9999, corr.method = "BH") 

# but we really want to know...
adonis(Dist.clr ~ plant.water, as(sample_data(ps.clr), "data.frame"), permutations = 9999)

no.adjust <- adonis.pair(Dist.clr, metadata.nocontrols$plant.water, nper = 9999, corr.method = "none") 

# comparisons to test
## forb.Control <-> forb.Drought 
## forb.Control <-> mix.Drought 
## forb.Drought <-> mix.Drought

## forb.Control <-> forb.Watering
## forb.Control <-> mix.Watering
## forb.Watering <-> mix.Watering

tests <- c("forb.Control <-> forb.Drought", "forb.Control <-> forb.Watering", "forb.Control <-> mix.Drought", "forb.Control <-> mix.Watering", "forb.Watering <-> mix.Watering", "forb.Drought <-> mix.Drought")

df.p <- as.data.frame(cbind(tests = tests, p.val = p.adjust(no.adjust[no.adjust$combination %in% tests,]$P.value, "BH")))

# Hellinger Transformation
# ps.hell <- transform(ps, "hellinger")
# 
# PCoA_hell <- ordinate(physeq = ps.hell, method = "PCoA", distance = "euclidean")
# 
# plot_ordination(ps.hell, PCoA_hell, color = "PlantTrt", shape = "PlantTrt") +
#   theme_bw(base_size = 15) +
#   geom_point(size = 2) +
#   facet_wrap(~WaterTrt) +
#   stat_ellipse()
# 
# # Get distances for stats
# Dist.hell <- phyloseq::distance(ps.hell, method = "euclidean", type = "samples")
# 
# set.seed(50)
# adonis(Dist.hell ~ plant.water, as(sample_data(ps.hell), "data.frame"), permutations = 9999) # no sig interaction
# 
# no.adjust <- adonis.pair(Dist.hell, metadata.nocontrols$plant.water, nper = 9999, corr.method = "none")
# tests <- c("forb.Control <-> forb.Drought", "forb.Control <-> forb.Watering", "forb.Control <-> mix.Drought", "forb.Control <-> mix.Watering", "forb.Watering <-> mix.Watering", "forb.Drought <-> mix.Drought")
# 
# df.p2 <- as.data.frame(cbind(tests = tests, p.val = p.adjust(no.adjust[no.adjust$combination %in% tests,]$P.value, "BH")))

#### .___Forb treats ####
plot_ordination(ps.clr, PCoA_clr, color = "PlotSubTreatmentName", shape = "PlantTrt") +
  theme_classic(base_size = 15) +
  geom_point(size = 2) +
  facet_wrap(~WaterTrt) +
  theme(
    legend.title = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1)
  ) +
  #stat_ellipse() +
  scale_color_brewer(palette = "Dark2")

# Get distances for stats
Dist.clr <- phyloseq::distance(ps.clr, method = "euclidean", type = "samples")

set.seed(50)
# adonis(Dist.clr ~ plant.water, as(sample_data(ps), "data.frame"), permutations = 9999) # no sig interaction
# 
# BH.adjust.plant <- adonis.pair(Dist.clr, metadata.nocontrols$PlantTrt, nper = 9999, corr.method = "BH") 
# 
# BH.adjust.water <- adonis.pair(Dist.clr, metadata.nocontrols$WaterTrt, nper = 9999, corr.method = "BH") 

# but we really want to know...
adonis(Dist.clr ~ plant.water, as(sample_data(ps.clr), "data.frame"), permutations = 9999)

no.adjust <- adonis.pair(Dist.clr, metadata.nocontrols$plant.water, nper = 9999, corr.method = "none") 

# comparisons to test
## forb.Control <-> forb.Drought 
## forb.Control <-> mix.Drought 
## forb.Drought <-> mix.Drought

## forb.Control <-> forb.Watering
## forb.Control <-> mix.Watering
## forb.Watering <-> mix.Watering

tests <- c("forb.Control <-> forb.Drought", "forb.Control <-> forb.Watering", "forb.Control <-> mix.Drought", "forb.Control <-> mix.Watering", "forb.Watering <-> mix.Watering", "forb.Drought <-> mix.Drought")

df.p <- as.data.frame(cbind(tests = tests, p.val = p.adjust(no.adjust[no.adjust$combination %in% tests,]$P.value, "BH")))


#### Sig Diff Genes ####

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
## Bonferroni: 44-57 sig diff genes
## BH: ~457 sigdif genes 

#plot results
grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], plot.list[[5]], ncol = 3)

#### Functional analysis ####
fun.gene.sig <- fun.gene[fun.gene$gene %in% sig.genes,] # only 12 have associated functions with bonferroni; 211 with BH adjustment

write.csv(fun.gene.sig, "fun-gene-cog-allsigs.csv", row.names = F)

fun.gene.cog <- fun.gene.sig[fun.gene.sig$COG20_CATEGORY != "NULL",]



fun.gene.cog$new.cog <- fun.gene.cog$COG20_CATEGORY

fun.gene.cog[fun.gene.cog$new.cog == "Lipid transport and metabolism!!!Lipid transport and metabolism",]$new.cog <- "Lipid transport and metabolism"

fun.gene.cog[fun.gene.cog$new.cog == "Signal transduction mechanisms!!!Signal transduction mechanisms!!!Signal transduction mechanisms" ,]$new.cog <- "Signal transduction mechanisms"

fun.gene.cog[fun.gene.cog$new.cog == "Signal transduction mechanisms!!!Signal transduction mechanisms" ,]$new.cog <- "Signal transduction mechanisms"

fun.gene.cog[fun.gene.cog$new.cog == "Function unknown!!!Function unknown" ,]$new.cog <- "Function unknown" 

fun.gene.cog[fun.gene.cog$new.cog == "Nucleotide transport and metabolism!!!Nucleotide transport and metabolism" ,]$new.cog <- "Nucleotide transport and metabolism"

fun.gene.cog[fun.gene.cog$new.cog == "Replication, recombination and repair!!!Replication, recombination and repair" ,]$new.cog <- "Replication, recombination and repair"

write.csv(fun.gene.cog, "Data/fun-gene-cog.csv") # easier to do the rest by hand oops

fun.gene.cog <- read.csv("Data/fun-gene-cog.csv")

unique(fun.gene.cog$new.cog) # 21 unique functions yay!

# # comparisons worth checking out: FD.MD, FC.MD, FC.FW, FC.MW
# # Changes due to watering and grass competition
# FC.FW <- merge(res.list$FC.FW, fun.gene.sig, by = "gene", all.x = T, all.y = F) #BH: 18 out of 75 have functions
# 
# FC.MW <- merge(res.list$FC.MW, fun.gene.sig, by = "gene", all.x = T, all.y = F) # BH: kicked up nothing in any database
# 
# # Changes due to drought and grass competition
# FD.MD <- merge(res.list$FD.MD, fun.gene.sig, by = "gene", all.x = T, all.y = F) # BH: 95 of 172 have functions
# 
# FC.MD <- merge(res.list$FC.MD, fun.gene.sig, by = "gene", all.x = T, all.y = F) #BH: 157 of 268 have functions

# for(i in names(res.list)) {
#  res.list[[i]] <- merge(res.list[[i]], fun.gene.sig, by = "gene")
# }


for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], fun.gene.cog[,c(1,13)], by = "gene")
}

sig.dif.fun <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.fun[[i]] <- ddply(res.list[[i]], .(new.cog), summarize, estimate = sum(estimate))
  }
}

plot.list2 <- list()

for(i in names(sig.dif.fun)){
    
    sig.dif.fun[[i]]$new.cog <- factor(sig.dif.fun[[i]]$new.cog, levels = sig.dif.fun[[i]]$new.cog[order(sig.dif.fun[[i]]$estimate)])
    
    p <- ggplot(sig.dif.fun[[i]], aes(x = new.cog, y = estimate)) + 
      geom_point(size = 3) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 9),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.title.y = element_blank(),
          ) +
      labs(title = i, y = expression(paste(log[2], " fold change"))) +
      coord_flip()
  
   plot.list2[[i]] = p
}

p <- grid.arrange(plot.list2[[2]], plot.list2[[4]], plot.list2[[1]], plot.list2[[3]], ncol = 4)

ggsave("change-genes.pdf", p, units = "in", width = 10, height = 3)

sig.dif.fundf <- data.frame()
for(i in 1:length(names(sig.dif.fun))) {
  sig.dif.fundf <- rbind(cbind(plant.water = names(sig.dif.fun)[i],sig.dif.fun[[i]]), sig.dif.fundf)
}

test <- expand.grid(plant.water = unique(sig.dif.fundf$plant.water), new.cog = unique(sig.dif.fundf$new.cog))

sig.dif.fundf <- merge(sig.dif.fundf, test, by = c("plant.water", "new.cog"), all.y = T)


# test <- sig.dif.fundf[sig.dif.fundf$plant.water == "FC.MD",]
# test<- test[order(test$estimate), ]
# 
# sig.dif.fundf$new.cog <- factor(sig.dif.fundf$new.cog, levels = sig.dif.fundf[order(test$new.cog), ]

FC.MD <- sig.dif.fundf[sig.dif.fundf$plant.water == "FC.MD" , ]
FC.MD$new.cog <- factor(FC.MD$new.cog, levels = FC.MD[order(FC.MD$estimate, decreasing = T),]$new.cog)
FD.MD <- sig.dif.fundf[sig.dif.fundf$plant.water == "FD.MD" , ]
test <- rbind(FC.MD, FD.MD)

ggplot(test, aes(x = new.cog, y = estimate)) + 
      geom_point(size = 3) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 9),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.title.y = element_blank(),
          ) +
      labs(y = expression(paste(log[2], " fold change"))) +
      coord_flip() + 
      facet_wrap(~plant.water)
  


res.list[["FC.MD"]][res.list[["FC.MD"]]$gene == 3731342,] #iron transport protein decreases significantly in drought with grasses; this one didnt have a function in COG, but did in ghost koala and noticed it said iron so looked it up

#### DESeq on Functions ####
counts.fun <- data.frame(cbind(gene = row.names(decontam_counts_nocontrols), decontam_counts_nocontrols)) # need non norm counts for DESeq so go back through, maybe just make a gene to reduced function sheet... that might be easier 

# counts.fun <- left_join(norm.counts.fun, fun.gene, by = "gene")
# 
# norm.counts.fun$COG20_CATEGORY <- ifelse(norm.counts.fun$COG20_CATEGORY == "NULL", norm.counts.fun$COG_CATEGORY, norm.counts.fun$COG20_CATEGORY)
# 
# unique(norm.counts.fun[!nchar(norm.counts.fun$COG20_CATEGORY) > 5, ]$COG20_CATEGORY)
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "T"),]$COG20_CATEGORY <- "Signal transduction mechanisms"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "L"),]$COG20_CATEGORY <- "Replication, recombination and repair"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "S"),]$COG20_CATEGORY <- "Function unknown"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "E"),]$COG20_CATEGORY <- "Amino acid transport and metabolism"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "G"),]$COG20_CATEGORY <- "Carbohydrate transport and metabolism"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "C"),]$COG20_CATEGORY <- "Energy production and conversion"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "P"),]$COG20_CATEGORY <- "Inorganic ion transport and metabolism"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "M"),]$COG20_CATEGORY <- "Cell wall/membrane/envelope biogenesis"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "I"),]$COG20_CATEGORY <- "Lipid transport and metabolism"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "V"),]$COG20_CATEGORY <- "Defense mechanisms"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "N"),]$COG20_CATEGORY <- "Cell motility"
# 
# norm.counts.fun[which(norm.counts.fun$COG20_CATEGORY == "Q"),]$COG20_CATEGORY <- "Secondary metabolites biosynthesis, transport and catabolism"
# 
# norm.counts.fun <- norm.counts.fun[,c(1:67,77)]
# 
# norm.counts.fun[,2:67] <- norm.counts.fun[,2:67] %>% mutate_all(as.numeric)
# 
# test <- separate(norm.counts.fun, COG20_CATEGORY, into = c("newcog1", "newcog2", "newcog3", "newcog4", "newcog5", "newcog6", "newcog7", "newcog8", "newcog9"), sep = "!!!")

#write.csv(test, "Data/norm-counts-fun.csv", row.names = F, na = "")
fun.reduced <- read.csv("Data/function-reduced.csv", na = "")
fun.reduced$gene <- as.factor(fun.reduced$gene)
counts.fun <- left_join(counts.fun, fun.reduced, by = "gene")
counts.fun[,2:67] <- counts.fun[,2:67] %>% mutate_all(as.numeric)

tmp2 <- data.frame()

for(i in 1:7) {
  tmp <- counts.fun[,c(1:67,67+i)]
  colnames(tmp)[68] <- "COG20_CATEGORY_reduced"
  if(i %in% 2:7) {
    tmp <- tmp[!is.na(tmp$COG20_CATEGORY_reduced),]
  }
  tmp2 <- rbind(tmp2, tmp)
} # do we want to lump COG "Function Unknown" with our own category? what to do about when a gene has both COG Function unknown and other known functions? 
tmp2[is.na(tmp2$COG20_CATEGORY_reduced),]$COG20_CATEGORY_reduced <- "NULL"

dds.fun <- tmp2 %>%
  group_by(COG20_CATEGORY_reduced) %>%
  dplyr::summarise(across(2:67, sum))

dds.fun <- as.data.frame(dds.fun)
row.names(dds.fun) <- dds.fun$COG20_CATEGORY_reduced
dds.fun <- as.matrix(dds.fun)[,-1]
mode(dds.fun) <- "numeric"

dds.fun <- DESeqDataSetFromMatrix(dds.fun,
                                   colData = metadata.nocontrols,
                                   design = ~ plant.water)
dds.fun <- estimateSizeFactors(dds.fun)

dds.fun <- DESeq(dds.fun)

plotDispEsts(dds.fun)

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
  res <- results(dds.fun, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
   sig.genes <- c(as.character(res.list[[i]]$gene), sig.genes)
}


#### Taxonomic analysis ####
# Cassie: Source tracker and taxonomy 
# are these genes coming from the same taxa


#### Soil ####
soil <- merge(soil, metadata, by.x = "Sample.ID", by.y = "PlotXSubTreatment")

soil$CN <- soil$Total.C..µg./soil$Total.N..µg.

soil <- filter(soil, Sample.ID != "27B") # notes that plot wasnt weeded

colnames(soil)[3] <- "C.ug"
colnames(soil)[6] <- "N.ug"

ggplot(soil, aes(x = PlantTrt, y = CN, group = interaction(WaterTrt, PlantTrt), fill = WaterTrt)) +
  geom_boxplot() 

ggplot(soil, aes(x = PlantTrt, y = CN)) +
  geom_boxplot() 

ggplot(soil, aes(x = WaterTrt, y = C.ug, group = interaction(WaterTrt, PlantTrt), fill = PlantTrt)) +
  geom_boxplot() 

ggplot(soil, aes(x = WaterTrt, y = N.ug, group = interaction(WaterTrt, PlantTrt), fill = PlantTrt)) +
  geom_boxplot() 

library(lme4)
library(lmerTest)
m1 <- lmer(CN ~ PlantTrt + WaterTrt + (1|Plot), data = soil)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1))
qqline(resid(m1), col = 2, lwd = 2, lty = 2)
anova(m1)

m2 <- lmer(C.ug ~ PlantTrt + WaterTrt + (1|Plot), data = soil)
plot(fitted(m2), resid(m2))
qqnorm(resid(m2))
qqline(resid(m2), col = 2, lwd = 2, lty = 2)
anova(m2)

m3 <- lmer(N.ug ~ PlantTrt + WaterTrt + (1|Plot), data = soil)
plot(fitted(m3), resid(m3))
qqnorm(resid(m3))
qqline(resid(m3), col = 2, lwd = 2, lty = 2)
anova(m3)

# link microbial beta diversity as it relates to soil nutrient variables
# normalized expression of DEGs or FunGroups (cog pathways/categories)
# could also do families 
