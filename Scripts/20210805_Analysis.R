rm(list=ls())

#### Load Libraries ####
library(plyr)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(decontam)
library(microbiome)

#### Read in Data ####

# # DESeq data
# txi.kallisto <- readRDS("txi.kallisto.08042021.RDS")
# 
# #http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# 
# 
# colnames(txi.kallisto$counts) <- sub("_S.*", "", colnames(txi.kallisto$counts))
# colnames(txi.kallisto$counts) <- sub("-", "_", colnames(txi.kallisto$counts))
# colnames(txi.kallisto$abundance) <- sub("_S.*", "", colnames(txi.kallisto$abundance))
# colnames(txi.kallisto$abundance) <- sub("-", "_", colnames(txi.kallisto$abundance))
# 
metadata <- read.csv("Data/McL_metag_metadata.csv")
metadata <- metadata[-c(1,3,5,6,8:10),] # only sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID
# 
# 
# #need to make sampleTable (we probably have this) with metadata
# #condition -> formula for deseq where forumula based on sampleTable e.g. watering_status + plant_comm_treatment + watering_status*plant_comm_treatment

ddsTxi <- DESeqDataSetFromMatrix(norm.counts2,
                                   colData = metadata,
                                   design = ~ WaterTrt + PlantTrt + WaterTrt:PlantTrt)
txi.kallisto <- readRDS("txi.kallisto.08042021.RDS")

#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


colnames(txi.kallisto$counts) <- sub("_S.*", "", colnames(txi.kallisto$counts))
colnames(txi.kallisto$counts) <- sub("-", "_", colnames(txi.kallisto$counts))
colnames(txi.kallisto$abundance) <- sub("_S.*", "", colnames(txi.kallisto$abundance))
colnames(txi.kallisto$abundance) <- sub("-", "_", colnames(txi.kallisto$abundance))

metadata <- read.csv("McL_metag_metadata.csv")
metadata <- metadata[-c(1,3,5,6,8:10),] # only sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID

library(DESeq2)

#need to make sampleTable (we probably have this) with metadata
#condition -> formula for deseq where forumula based on sampleTable e.g. watering_status + plant_comm_treatment + watering_status*plant_comm_treatment
ddsTxi <- DESeqDataSetFromTximport(txi.kallisto,
                                   colData = metadata,
                                   design = ~ New_Treatment)

#remove any rows with no counts
dds <- ddsTxi[ rowSums(counts(ddsTxi)) > 10, ]

dds <- estimateSizeFactors(dds)

#run deseq diff abundance analysis
dds.test <- DESeq(dds)
saveRDS(dds.test, "GL_dds.deseq.obj.09132021.RDS")

#extract deseq normalized counts
normalized_counts <- counts(dds.test, normalized=TRUE)
saveRDS(normalized_counts, "GL_normcounts.09132021.RDS")
norm.counts <- readRDS("Data/DESeq/GL_normcounts.09132021.RDS")
norm.counts2 <- as.data.frame(norm.counts)

dds.old <- readRDS("Data/DESeq/GL_dds.deseq.obj.09132021.RDS")

# # Taxonomy data (Kaiju)
tax <- read.delim("Data/Post-processing/Salmon/Kaiju-taxonomy-table-for-anvio-2019.txt", header = F)
# 
# colnames(tax) <- c("classified", "gene_call", "NCBI_taxon_ID", "score", "taxon_ID", "acc", "match_seq", "taxonomy_string")
# 
# nrow(tax[tax$classified == "U",])/nrow(tax) #25% of gene calls could not be classified to a taxon
# 
# tax.c <- filter(tax, classified == "C")
# 
# tax.c.uni <- unique(tax.c[,c(1,8)])
# row.names(tax.c.uni) <- tax.c.uni$taxonomy_string
# 
# tax.c.uni <-
#   tax.c.uni %>% separate(taxonomy_string, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") # it thinks it is discarding information in the last "column" beacuse the string ends with a ";" but there is nothing in that column
# 
# for(i in 2:8){
#   tax.c.uni[,i] <- trimws(tax.c.uni[,i], which = "both")
# }
# 
# tax.c.bac <- filter(tax.c.uni, domain == "Bacteria") # filter to bacteria for now (98.5% of taxa are bacteria)
# 
# tax.c.bac[tax.c.bac == "NA"] <- NA
# tax_final <- unique(tax.c.bac[,-1])

#rm(tax)

# Functional data (Ghost Koala)
Fun <- read.delim("Data/Salmon/GL_KeggAnnotations_AnviImportable.txt")
colnames(Fun)[c(1,4)] <- c("gene_call", "Function")

# Metadata
metadata <- read.csv("Data/McL_metag_metadata.csv")
metadata <- metadata[-c(1,3,5,6,8:10),] # only sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID
#metadata.nocontrols <- metadata[-c(1:10),]


#### Prep Data ####

## Phyloseq object with taxonomy

# Gene call filtered for taxa
# gene.tax <- gene.call[row.names(gene.call) %in% tax.c$gene_call,]
# gene.tax <- cbind(gene_call = row.names(gene.tax), gene.tax)
# 
# gene.tax <- merge(gene.tax, tax[,c(2,8)], by = "gene_call", all.x = T, all.y = F)
# 
# gene.tax.sum <- gene.tax %>%
#   group_by(taxonomy_string) %>%
#   summarize(dplyr::across(GLM_0697:GLM_0241, sum)) %>%
#   ungroup()
# 
# gene.tax.sum <- column_to_rownames(gene.tax.sum, var = "taxonomy_string")
# 
tax_final <- as.matrix(tax_final)

otu_table <- otu_table(norm.counts, taxa_are_rows = T)
taxa_table <- tax_table(tax_final)
mapping_file <- sample_data(metadata)
ps <- phyloseq(otu_table, mapping_file, taxa_table)

#sample_data(ps)$is.neg <- sample_data(dds)$SampleType == "Kit"

#threshold=0.5, which will identify as contaminants 
#all sequences thare are more prevalent in negative controls than in positive samples

contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)
table(contamdf.prev05$contaminant) # 55 contaminants

#which ASV is contaminant
length(which(contamdf.prev05$contaminant))/nrow(contamdf.prev05) #1.2%

contam.rows <- which(contamdf.prev05$contaminant)
contaminants <- rownames(contamdf.prev05[contam.rows,])

allTaxa = taxa_names(ps)

#returns list of all taxa, except contaminants
allTaxa <- allTaxa[!(allTaxa %in% contaminants)]

#keeps taxa in allTaxa and gets rid of anthing else
ps_filt = prune_taxa(allTaxa, ps)

ps.nocontrols <- subset_samples(ps_filt, SampleSubType != "Control_ZymoMock")
ps.nocontrols <- subset_samples(ps.nocontrols, SampleSubType != "Control_H20")
ps.nocontrols <- subset_samples(ps.nocontrols, SampleSubType != "Control_Kit")

## Phyloseq object with function

#### Taxonomy ####

####.---Ordination: grass, forb, mix ####
vsd <- vst(dds.old, blind = F)
#rld <- rlog(dds, blind=FALSE)
bpcaData <- plotPCA(vsd, intgroup=c("New_Treatment", "PlotTreatment"), returnData=TRUE) #condition would have to be like watering_treatment or such 
percentVar <- round(100 * attr(bpcaData, "percentVar"))

ggplot(bpcaData[bpcaData$PlotTreatment != "None",], aes(PC1, PC2, color = New_Treatment)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    facet_wrap(~PlotTreatment) +
    stat_ellipse()


# ps.nocontrols_hell <- transform(ps.nocontrols, "hellinger")
# 
# tax_pcoa <- ordinate(
#   physeq = ps.nocontrols_hell, 
#   method = "PCoA", 
#   distance = "bray")
# 
# DistBC <- phyloseq::distance(ps.nocontrols, method = "bray", type = "samples")
# 
# set.seed(50)
# 
# adonis(DistBC ~ PlotTreatment, as(sample_data(ps.nocontrols), "data.frame"), permutations = 9999) 
# 
# C <- betadisper(DistBC, as(sample_data(ps.nocontrols), "data.frame")$PlotTreatment)
# set.seed(50)
# permutest(C, permutations = 9999) 
# plot(C, label = F)
# 
# ord.plot.a <- plot_ordination(ps.nocontrols_hell, tax_pcoa, color = "PlotTreatment", shape = "PlotTreatment") +
#   theme_bw(base_size = 15) +
#   geom_point(size = 2) +
#   stat_ellipse(aes(group = PlotTreatment)) +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "right",
#     legend.text = element_text(size = 10),
#     legend.margin = margin(c(.1,.1,.1,.1))
#   ) 
# 
# ord.plot.b <- plot_ordination(ps.nocontrols_hell, tax_pcoa, color = "New_Treatment", shape = "New_Treatment") +
#   theme_bw(base_size = 15) +
#   geom_point(size = 2) +
#   stat_ellipse(aes(group = New_Treatment)) +
#   theme(
#     legend.title = element_blank(),
#     legend.position = "right",
#     legend.text = element_text(size = 10),
#     legend.margin = margin(c(.1,.1,.1,.1))
#   ) 

####.---DESeq2 ####

#get results
#this is where we would need contrasts / loops
#if we had more than two groups to compare
resultsNames(dds) # why intercept? 

res <- results(dds) 

#get only results that are sig & have log fold change >2
res.sig <- subset(res, (padj < 0.05))
as.data.frame(res.sig) #only two? 3271634, 3964977

contrasts <- c("F.G", "FG.G")

contrast.list <- list(F.G = c("FunGroup", "Forb", "Grass"),
                      FG.G = c("FunGroup", "grass_x_forb", "Grass"))


plot.name.list <- list(F.G = "Forb vs. Grass",
                       FG.G = "Competition vs. Grass")

#### Function ####

####.---Ordination ####
## sum genes to function

####.---DESeq2 ####

# extract normalized abundances and redo ordinations and diversity metrics using deseq normalized abundances and look for differentially abundant genes 

# can make more complicated models

# can summarize by taxonomy and function later, do above first 

# Cassie's thoughts
# play with DA between treatments rerunning deseq
# 