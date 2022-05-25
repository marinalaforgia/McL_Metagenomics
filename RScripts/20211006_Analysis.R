rm(list=ls())

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("decontam", version = "3.14")

#### Load Libraries ####
library(DESeq2)
library(plyr)
library(tidyverse)
library(pheatmap)
library(phyloseq)
library(decontam)
library(vegan)
library(Rmisc)
library(broom)

#http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

#### Read in data ####

## Kallisto is abundance, counts, and length, raw GENE ABUNDANCE data fed into deseq
txi.kallisto <- readRDS("Data/Salmon/txi.kallisto.08042021.RDS")
colnames(txi.kallisto$counts) <- sub("_S.*", "", colnames(txi.kallisto$counts))
colnames(txi.kallisto$counts) <- sub("-", "_", colnames(txi.kallisto$counts))
colnames(txi.kallisto$abundance) <- sub("_S.*", "", colnames(txi.kallisto$abundance))
colnames(txi.kallisto$abundance) <- sub("-", "_", colnames(txi.kallisto$abundance))

## Metadata
metadata <- read.csv("Data/McL_metag_metadata.csv")
metadata <- metadata[-c(1,3,5,6,8:10),] # only sequenced one water, one kit, and one mock
#metadata <- filter(metadata, SampleSubType == "Rhizosphere_soil")
row.names(metadata) <- metadata$SampleID
metadata$plant.water <- paste(metadata$PlantTrt, metadata$WaterTrt, sep = ".") # make a new variable to avoid interactions

#need to make sampleTable (we probably have this) with metadata
#condition -> formula for deseq where forumula based on sampleTable e.g. watering_status + plant_comm_treatment + watering_status*plant_comm_treatment
txi.kallisto[[1]] <- txi.kallisto[[1]][,-c(1:3)]
txi.kallisto[[2]] <- txi.kallisto[[2]][,-c(1:3)]
txi.kallisto[[3]] <- txi.kallisto[[3]][,-c(1:3)]

# ddsTxi <- DESeqDataSetFromTximport(txi.kallisto,
#                                    colData = metadata,
#                                    design = ~ plant.water)
# 
# #remove any rows with no counts
# dds <- ddsTxi[ rowMeans(counts(ddsTxi)) > 0.5, ] # check on filtering of dds obj
# 
# dds <- estimateSizeFactors(dds)
# 
# #run deseq diff abundance analysis
# dds <- DESeq(dds)
#saveRDS(dds.test, "Data/GL_dds.deseq.obj.10062021.RDS") # Cassie's original
#dds <- readRDS("Data/DESeq/GL_dds.deseq.obj.09132021.RDS") # Controls removed
dds <- readRDS("Data/DESeq/GL_dds.deseq.obj.10062021.RDS") # controls and contaminants removed

#extract deseq normalized counts
#normalized_counts <- counts(dds, normalized = F) # raw counts to run through decontam
#normalized_counts <- counts(dds, normalized = T)
#saveRDS(normalized_counts, "Data/DESeq/GL_normcounts.10062021.RDS")
norm.counts <- readRDS("Data/GL_normcounts.10062021.RDS")
#vsd <- varianceStabilizingTransformation(dds)
vsd <- vst(dds, blind = T) # test blind = T for QA/QC
#rld <- rlog(dds, blind = F) #crashed my computer

bpcaData <- plotPCA(vsd, intgroup=c("SubPlotRichness"), returnData=TRUE) 
percentVar <- round(100 * attr(bpcaData, "percentVar"))

ggplot(bpcaData, aes(PC1, PC2, color = SubPlotRichness)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() 

bpcaData <- plotPCA(vsd, intgroup=c("WaterTrt", "PlantTrt"), returnData=TRUE) 
percentVar <- round(100 * attr(bpcaData, "percentVar"))

ggplot(bpcaData, aes(PC1, PC2, color = PlantTrt)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    facet_wrap(~WaterTrt) +
    stat_ellipse()

vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat))

df <- cbind(metadata, pca$x)
ggplot(df) + 
    geom_point(aes(x=PC3, y=PC4, color = PlantTrt)) +
    facet_wrap(~WaterTrt)

vsd_cor <- cor(vsd_mat) 
pheatmap(vsd_cor, annotation = metadata[,c(7,12)]) # no idea if this is telling me anything

plotDispEsts(dds)

### Old DDS (Cassie's run) 
# vsd <- vst(dds.old, blind = F)
# bpcaData <- plotPCA(vsd, intgroup=c("New_Treatment", "PlotTreatment"), returnData=TRUE) 
# percentVar <- round(100 * attr(bpcaData, "percentVar"))
# 
# ggplot(bpcaData[bpcaData$PlotTreatment != "None",], aes(PC1, PC2, color = New_Treatment)) +
#     geom_point(size=3) +
#     xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#     ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#     coord_fixed() +
#     facet_wrap(~PlotTreatment) +
#     stat_ellipse() ## why do these look different?
# 


### Trying to continue DESeq analysis
#DESeq2 doesn’t actually use normalized counts, rather it uses the raw counts and models the normalization inside the Generalized Linear Model (GLM). These normalized counts will be useful for downstream visualization of results, but cannot be used as input to DESeq2 or any other tools that perform differential expression analysis that use the negative binomial model.


#### ignore for now ####
# Specify contrast for comparison of interest
#contrasts we are interested in 
#f.c - f.d -> f.c-m.d?
#f.c - f.w -> f.c-m.w?
# 
# resultsNames(dds) # why intercept? 
# 
# res <- results(dds) 
# 
# contrasts <- c("F.G", "FG.G")
# 
# contrast.list <- list(F.G = c("WaterTrt.PlantTrt", "forb", "grasses"),
#                       FG.G = c("FunGroup", "grass_x_forb", "Grass"))
# 
# contrast <- c("condition", "level_to_compare", "base_level")
# 
# # Output results of Wald test for contrast
# res <- results(dds, 
# 	       contrast = contrast.list, 
# 	       alpha = 0.05)

# # Shrink the log2 fold changes to be more accurate
# res <- lfcShrink(dds, 
# 		 contrast = contrast, 
# 		 type = "normal")

## LRT
# dds.LRTwater <- DESeq(dds, test = "LRT", reduced = ~WaterTrt)
# dds.LRTplant <- DESeq(dds, test = "LRT", reduced = ~PlantTrt)
# dds.lrt <- nbinomLRT(dds, reduced = ~ WaterTrt + PlantTrt)

# ## Example 3: two conditions, three genotypes
# 
# res <- results(dds) 
# 
# #get only results that are sig & have log fold change >2
# res.sig <- subset(res, (padj < 0.1))
# as.data.frame(res.sig)
# 
# # ~~~ Using interaction terms ~~~
# dds <- makeExampleDESeqDataSet(n=100,m=18)
# dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
# design(dds) <- ~ genotype + condition + genotype:condition
# dds <- DESeq(dds)
# resultsNames(dds)
# # the condition effect for genotype I (the main effect)
# results(dds, contrast=c("PlantTrt","mix","forb")) # nothing
# results(dds, contrast=c("PlantTrt","mix","grasses")) # nothing
# results(dds, contrast=c("plant.water","forb.Control","forb.Watering")) # nothing

#results(dds, contrast=c("PlantTrt","forb","grasses")) # nothing
# the condition effect for genotype III.
# this is the main effect *plus* the interaction term
# (the extra condition effect in genotype III compared to genotype I).
# results(dds, contrast=list( c("WaterTrt_Drought_vs_Control","genotypeIII.conditionB") ))
# 
# # the interaction term for condition effect in genotype III vs genotype I.
# # this tests if the condition effect is different in III compared to I
# results(dds, name="WaterTrtDrought.PlantTrtmix")


# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
# results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))
# Note that a likelihood ratio could be used to test if there are any
# differences in the condition effect between the three genotypes.
# ~~~ Using a grouping variable ~~~
# This is a useful construction when users just want to compare
# specific groups which are combinations of variables.
# dds$group <- factor(paste0(dds$genotype, dds$condition))
# design(dds) <- ~ group
# dds <- DESeq(dds)
# resultsNames(dds)
# # the condition effect for genotypeIII
# results(dds, contrast=c("group", "IIIB", "IIIA"))
# 
# res <- results(dds)
# dds.LRTpw <- DESeq(dds, test = "LRT", reduced = ~ WaterTrt + PlantTrt)
# res_LRT <- results(dds_lrt)
# 
# # Create a tibble for LRT results
# res_LRT_tb <- res_LRT %>%
#   data.frame() %>%
#   rownames_to_column(var="gene") %>% 
#   as_tibble()
# 
# # Subset to return genes with padj < 0.05
# sigLRT_genes <- res_LRT_tb %>% 
#   filter(padj < padj.cutoff)
# 
# # Get number of significant genes
# nrow(sigLRT_genes)
# 
# # Compare to numbers we had from Wald test
# nrow(sigOE)
# nrow(sigKD)

# #### Decontam ####
# dds <- readRDS("Data/DESeq/GL_dds.deseq.obj.09132021.RDS") # cassie's data
# #extract deseq normalized counts
# normalized_counts <- counts(dds, normalized = F) # raw counts to run through decontam
# metadata <- read.csv("Data/McL_metag_metadata.csv")
# metadata <- metadata[-c(1,3,5,6,8:10),] # only sequenced one water, one kit, and one mock
# #metadata <- filter(metadata, SampleSubType == "Rhizosphere_soil")
# row.names(metadata) <- metadata$SampleID
# colnames(normalized_counts) <- metadata$SampleID
# 
# #first tell it which samples are the NC
# sample_data(ps)$is.neg <- sample_data(ps)$SampleType == "Kit"
# 
# #threshold=0.5, which will identify as contaminants 
# #all sequences thare are more prevalent in negative controls than in positive samples
# 
# contamdf.prev05 <- isContaminant(ps, method = "prevalence", neg = "is.neg", threshold = 0.5)
# table(contamdf.prev05$contaminant) # 4774 contaminants
# 
# # #which ASV is contaminant
# # which(contamdf.prev05$contaminant)
# #
# # #Make phyloseq object of presence-absence in negative controls
# # ps.neg <- prune_samples(sample_data(ps)$is.neg == "TRUE", ps)
# # ps.neg.presence <- transform_sample_counts(ps.neg, function(abund) 1*(abund>0))
# # 
# # #Make phyloseq object of presence-absence in true positive samples
# # ps.pos <- prune_samples(sample_data(ps)$is.neg == "FALSE", ps)
# # ps.pos.presence <- transform_sample_counts(ps.pos, function(abund) 1*(abund>0))
# # 
# # #Make data.frame of prevalence in positive and negative samples 
# # #using prev threshold = 0.5
# # df.pres <- data.frame(prevalence.pos=taxa_sums(ps.pos.presence),
# #                       prevalence.neg=taxa_sums(ps.neg.presence),
# #                       contam.prev=contamdf.prev05$contaminant)
# # ggplot(data=df.pres, aes(x=prevalence.neg, y=prevalence.pos, color=contam.prev)) + geom_point()
# 
# contam.rows <- which(contamdf.prev05$contaminant)
# contaminants <- rownames(contamdf.prev05[contam.rows,])
# rownames(normalized_counts)
# 
# #returns list of all taxa, except contaminants
# counts <- normalized_counts[!(rownames(normalized_counts) %in% contaminants),]
# counts <- counts[,-c(1:3)] # get rid of controls
# 
# counts <- readRDS("Data/DESeq/counts_decontam.RDS")
# 
# ps <- phyloseq(otu_table(counts, taxa_are_rows = T), sample_data(metadata))
# ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.2*length(x)), TRUE)
# 
# #run deseq diff abundance analysis
# 
# metadata <- read.csv("Data/McL_metag_metadata.csv")
# metadata <- filter(metadata, SampleSubType == "Rhizosphere_soil")
# row.names(metadata) <- metadata$SampleID
# metadata$plant.water <- paste(metadata$PlantTrt, metadata$WaterTrt, sep = ".")
# 
# dds <- phyloseq_to_deseq2(ps, ~ plant.water)
# 
# 
# # dds <- DESeqDataSetFromMatrix(counts,
# #                               colData = metadata,
# #                               design = ~ plant.water,
# #                               tidy = FALSE,
# #                               ignoreRank = FALSE)
# 
# #remove any rows with no counts
# 
# #dds <- dds[rowSums(counts(dds)) > 1, ]
# 
# #dds  <- dds[sum(counts(dds) > 3) > (0.2*length(counts(dds))),]
# dds <- estimateSizeFactors(dds)
# 
# #run deseq diff abundance analysis
# dds <- DESeq(dds)
# 
# saveRDS(dds, "GL_dds.deseq.obj.10212021.RDS")

dds <- readRDS("Data/DESeq/GL_dds.deseq.obj.10212021.RDS")

#### DESeq goodness of fit ####
dds <- readRDS("Data/DESeq/GL_dds.deseq.obj.10212021.RDS")
plotDispEsts(dds)

#### Ordination ####
vsd <- vst(dds, blind = T) # test blind = T for QA/QC
#rld <- rlog(dds, blind = F) #crashed my computer

bpcaData <- plotPCA(vsd, intgroup=c("WaterTrt", "PlantTrt"), returnData=TRUE) 
percentVar <- round(100 * attr(bpcaData, "percentVar"))

ggplot(bpcaData, aes(PC1, PC2, color = PlantTrt)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    coord_fixed() +
    facet_wrap(~WaterTrt) +
    stat_ellipse()

vsd_mat <- assay(vsd)
pca <- prcomp(t(vsd_mat))

df <- cbind(metadata, pca$x)
ggplot(df) + 
    geom_point(aes(x=PC3, y=PC4, color = PlantTrt)) +
    facet_wrap(~WaterTrt)

# Are they sig diff though?
dist.matrix <- t(data.frame(otu_table(ps.brack))) 
bray.not.na <- vegdist(dist.matrix, method = "bray")
set.seed(50)
adonis(bray.not.na ~ WaterTrt * PlantTrt, as(sample_data(ps.brack), "data.frame"), permutations = 9999)

#### Contrasts ####
res <- results(dds) 

res.sig <- subset(res, (padj < 0.05))
as.data.frame(res.sig)

res <- results(dds, contrast=c("plant.water","forb.Control","mix.Control")) # nothing
res.sig <- subset(res, (padj < 0.05))
as.data.frame(res.sig) 

res <- results(dds, contrast=c("plant.water","forb.Control","mix.Watering"), pAdjustMethod = "bonferroni")
res.sig <- subset(res, (padj < 0.05))
as.data.frame(res.sig)

res <- results(dds, contrast=c("plant.water","forb.Drought","mix.Drought")) 
res.sig <- subset(res, (padj < 0.05))
as.data.frame(res.sig)

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
               "FD.MD")

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
                      FD.MD = c("plant.water","forb.Drought","mix.Drought"))

plot.name.list <- list(FC.FW = "Forbs: Control v Watering",
                       FC.FD = "Forbs: Control v Drought",
                       FC.MD = "mix in drought vs forb control",
                       FC.MW = "mix in watering vs forb control",
                       GC.GD = "Grasses: Control v Drought",
                       GC.GW = "Grasses: Control v Watering",
                       GC.MD = "Mixed in drought vs grass control",
                       GC.MW  = "mix in water vs grass control",
                       GC.FC = "grass v forb control",
                       GD.FD = "grass v forb drought",
                       GW.FW = "grass v forb watering",
                       FD.MD = "forb drought v mix drought")

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
  res.list[[paste(i, sep = ".")]] <- tidy(res.alpha)

 res.list[[paste(i, sep = ".")]]$gene <- factor(res.list[[paste(i, sep = ".")]]$gene, levels = res.list[[paste(i, sep = ".")]]$gene[order(res.list[[paste(i, sep = ".")]]$estimate)])
  
   p <- ggplot(res.list[[paste(i, sep = ".")]], aes(x = gene, y = estimate)) + 
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
   sig.genes <- c(as.numeric(as.character(res.list[[paste(i, sep = ".")]]$gene)), sig.genes)
}

#plot results

plot.list$FC.FW + plot.list$FC.MW
plot.list$FC.FD + plot.list$FC.MD
plot.list$FD.MD
plot.list$GC.GD + plot.list$GC.MD
plot.list$GC.MW + plot.list$GC.GW
plot.list$GC.FC + plot.list$GD.FD

####
kaiju.tax <- read.delim("Data/Salmon/Kaiju-taxonomy-table-for-anvio-2019.txt", header = F) #taxonomy for gene calls
colnames(kaiju.tax) <- c("classified", "gene_call", "NCBI_taxon_ID", "score", "taxon_ID", "acc", "match_seq", "taxonomy_string")
KEGG.sum <- read.delim("Data/Salmon/GL_KeggAnnotations_AnviImportable.txt") #KEGG hits for our gene calls
KEGG <- read.delim("Data/Archive/KeggOrthology_Table1.txt", sep = ",") # KEGG pathways to match ours against

head(kaiju.tax)

kaiju.tax <- filter(kaiju.tax, gene_call %in% sig.genes, classified == "C")


kaiju.tax2 <-
  kaiju.tax %>% separate(taxonomy_string, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

for(i in ncol(kaiju.tax2)){
  kaiju.tax2[,12] <- trimws(kaiju.tax2[,12], which = "both")
}

KEGG.hits <- filter(KEGG.sum, gene_callers_id %in% sig.genes) # KEGG for our sig dif genes
KEGG3 <- filter(KEGG, accession %in% KEGG.hits$accession) # associated categories and pathways for those accessions

#### Groups to compare for functions ####
contrasts <- c("FC.FW", "FC.MW", "FC.MD", "GC.FC", "FD.MD")

res.df <- data.frame(contrast = character(),
                     gene = character(),
                     baseMean = numeric(),
                     estimate = numeric(),
                     stderror = numeric(),
                     statistic = numeric(),
                     p.value = numeric(),
                     p.adjusted = numeric()
                     )

for(i in contrasts){
  res.df <- rbind(res.df, cbind(contrast = i, res.list[[i]]))
}

res.df <- merge(res.df, KEGG.hits[,c(1,3)], by.x = "gene", by.y = "gene_callers_id", all.x = T) # so many unclassified

res.df <- merge(res.df, KEGG3, by = "accession")

test <- ddply(res.df, .(contrast, Category3), summarise, logchange = sum(estimate))


test$Category2 <- factor(test$Category2, levels = test$Category2[order(test$logchange)])
 
ggplot(test, aes(x = Category3, y = logchange)) +
  geom_hline(yintercept = 0) +
  geom_point() +
  coord_flip() +
  facet_wrap(~contrast)

#### Bracken ####
species <- read_tsv("Data/S_bracken_summary.tsv")
species <- filter(species, name != "Homo sapiens", )

species2 <- species %>% 
  dplyr::select(ends_with("frac"))
species2 <- rbind(species2, t(data.frame(unident = (1 - colSums(species2)))))

colnames(species2) <- sub("_S.*", "", colnames(species2))
colnames(species2) <- sub("-", "_", colnames(species2))

species2 <- as.data.frame(species2[,-c(67,68,69)])
rownames(species2) <- c(species$name, "unidentified")

ps.brack <- phyloseq(otu_table(species2, taxa_are_rows = T), sample_data(metadata2))

BK_pcoa <- ordinate(
  physeq = ps.brack, 
  method = "PCoA", 
  distance = "bray")
# library(devtools)
# install_version("vegan", version ="2.5-5", repos = "http://cran.us.r-project.org")


dist.matrix <- t(data.frame(otu_table(ps.brack)))
bray.not.na <- vegdist(dist.matrix, method = "bray")
DistBC <- phyloseq::distance(BK_pcoa, method = "bray", type = "samples") # not working, postentiall known issue with vegan/phyloseq

set.seed(50)
adonis(bray.not.na ~ WaterTrt * PlantTrt, as(sample_data(ps.brack), "data.frame"), permutations = 9999) # hm. plant treatments are different but not watering treatments? >_<

plot_ordination(ps.brack, BK_pcoa, color = "PlantTrt", shape = "PlantTrt") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  facet_wrap(~WaterTrt) +
  stat_ellipse()

## now at the family level
family <- read_tsv("Data/F_bracken_summary.tsv")
family <- filter(family, name != "Hominidae", )

family2 <- family %>% 
  dplyr::select(ends_with("frac"))
family2 <- rbind(family2, t(data.frame(unident = (1 - colSums(family2)))))

# family2 <- family %>% 
#   dplyr::select(ends_with("num"))

colnames(family2) <- sub("_S.*", "", colnames(family2))
colnames(family2) <- sub("-", "_", colnames(family2))

# unident <- read.csv("Data/kraken_total_reads.csv")
# unident$Sample <- sub("-", "_", unident$Sample)
# sums <- data.frame(names = colnames(family2), ident.reads  = colSums(family2))
# total <- merge(unident, sums, by.x = "Sample", by.y = "names")
# total$unident <- total$TotalReads - total$ident.reads
# test <- data.frame(unidentified = total[,4])
# row.names(test) <- total$Sample
# test <- t(test)
# family2 <- rbind(family2, test)

family2 <- as.data.frame(family2[,-c(67,68,69)])
rownames(family2) <- c(family$name, "unidentified")
metadata2 <- metadata[order(metadata$SampleNumber),]
metadata2$SampleID <- as.factor(as.character(metadata2$SampleID))
ps.brack <- phyloseq(otu_table(family2, taxa_are_rows = T), sample_data(metadata2))

BK_pcoa <- ordinate(
  physeq = ps.brack, 
  method = "PCoA", 
  distance = "bray")
# library(devtools)
# install_version("vegan", version ="2.5-5", repos = "http://cran.us.r-project.org")


dist.matrix <- t(data.frame(otu_table(ps.brack)))
bray.not.na <- vegdist(dist.matrix, method = "bray")
DistBC <- phyloseq::distance(BK_pcoa, method = "bray", type = "samples") # not working, postentiall known issue with vegan/phyloseq

set.seed(50)
adonis(bray.not.na ~ WaterTrt * PlantTrt, as(sample_data(ps.brack), "data.frame"), permutations = 9999) # hm. plant treatments are different but not watering treatments? >_<

plot_ordination(ps.brack, BK_pcoa, color = "PlantTrt", shape = "PlantTrt") +
  theme_bw(base_size = 15) +
  geom_point(size = 2) +
  facet_wrap(~WaterTrt) +
  stat_ellipse()

#### Taxonomic differences ####
# I guess start with the taxa I know from the other project?
# can also look at taxa that vary or have higher means
# what about taxa that the other drought paper found? Actinobacteria (though that is a higher taxonomic group, would have to have cassie gimme that shit)

#### SourceTracker ####
ST <- read.csv("Data/SourceTracker/sourcetracker_fam_bac.csv")
ST <- pivot_longer(ST, cols = c("forb", "grasses", "Unknown"), names_to = "PlantTrt", values_to = "prop")
ST.sum <- summarySE(ST, measurevar = "prop", groupvars = c("Treatment", "PlantTrt"))

ggplot(ST.sum, aes(y = prop, x = Treatment, fill = PlantTrt)) +
  geom_bar(stat = "identity") +
  #geom_errorbar(aes(ymin = prop - se, ymax = prop + se), width = 0.2) +
  scale_color_viridis_d()

#### Functional Analyses ####
fun <- read.csv("Data/GL_function_per_gene.csv")
head(fun)
colnames(fun)
