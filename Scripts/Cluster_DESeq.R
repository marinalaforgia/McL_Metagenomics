library(plyr)
library(tidyverse)
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(decontam)
library(microbiome)

gene.call <- readRDS("Data/Post-processing/Salmon/txi.kallisto.08042021.RDS")
colnames(gene.call$counts) <- sub("_S.*", "", colnames(gene.call$counts))
colnames(gene.call$counts) <- sub("-", "_", colnames(gene.call$counts))
colnames(gene.call$abundance) <- sub("_S.*", "", colnames(gene.call$abundance))
colnames(gene.call$abundance) <- sub("-", "_", colnames(gene.call$abundance))
head(gene.call$length) # not sure what this is
metadata <- read.csv("Data/Post-processing/McL_metag_metadata.csv")
metadata <- metadata[-c(1,3,5,6,8:10),] # only sequenced one water, one kit, and one mock
row.names(metadata) <- metadata$SampleID

ddsTxi <- DESeqDataSetFromTximport(gene.call,
                                    colData = metadata,
                                    design = ~ New_Treatment)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 1, ]


ddsTxi <- estimateSizeFactors(ddsTxi)

#run deseq diff abundance analysis
dds.test <- DESeq(ddsTxi)

#get results
#this is where we would need contrasts / loops
#if we had more than two groups to compare
res <- results(dds.test) 

#get only results that are sig & have log fold change >2
res.sig <- subset(res, (padj < 0.01 & abs(log2FoldChange) >2))

#write table with all differentially expressed genes (DEG) between Spo and Inf
write.table(res.sig, "Result_DEGs.tsv", 
            quote=F, 
            row.names=T, 
            sep="\t")
