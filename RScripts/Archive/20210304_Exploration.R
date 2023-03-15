# Data Exploration Metagenomics #


library(tidyverse)
library(lme4)
library(lmerTest)
library(decontam)

# metagenome is focus taxa approach
# table of kaiju and function can connect these and play with it; can also connect to coverage dataset; figure out how to connect them
# gene calls is the community wide approach 
# MetaB.mc <- read.delim("Data/Post-processing/mean_coverage.txt") # Metagenome bin mean coverage
# MetaB.RA <- read.delim("Data/Post-processing/relative_abundance.txt") # Metagenome bin relative abundance

cov.g <- read.delim("Data/Post-processing/GL_meta_merged_kaiju_2500bp_gene_coverage-GENE-COVERAGES.txt") # gene call coverage
det.g <- read.delim("Data/Post-processing/GL_meta_merged_kaiju_2500bp_gene_coverage-GENE-DETECTION.txt") # gene call detection

#kaiju.sum <- read.delim("Data/Post-processing/Kaiju_summary_2019.txt") # sum of gene sequences per taxon ID
kaiju.tax <- read.delim("Data/Post-processing/Kaiju_taxonomy_table_for_anvio_2019.txt", header = F) #taxonomy for gene calls
colnames(kaiju.tax) <- c("classified", "gene_call", "NCBI_taxon_ID", "score", "taxon_ID", "acc", "match_seq", "taxonomy_string")
KEGG.sum <- read.delim("Data/Post-processing/GL_KeggAnnotations_AnviImportable.txt") #KEGG hits for our gene calls
KEGG <- read.delim("Data/Post-processing/KeggOrthology_Table1.txt", sep = ",") # KEGG pathways to match ours against

#meta-data
metadata <- read.csv("Data/Post-processing/McL_metag_metadata.csv")

kaiju.tax2 <-
  kaiju.tax %>% separate(taxonomy_string, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";")

for(i in ncol(kaiju.tax2)){
  kaiju.tax2[,12] <- trimws(kaiju.tax2[,12], which = "both")
}

# let's start with something simple, methylophilaceae

kaiju.tax3 <- filter(kaiju.tax2, family == "Clostridiaceae")

cov.g2 <- filter(cov.g, key %in% kaiju.tax3$gene_call)
colnames(cov.g2) <- sub("_[^_]+$", "", colnames(cov.g2))
cov.g2 <- as.data.frame(t.data.frame(cov.g2))
cov.g2$cov.sum <- rowSums(cov.g2)
cov.g2$SampleID <- row.names(cov.g2)
cov.g2 <- cov.g2[-1, 98:99]

meth <- merge(metadata, cov.g2, by = "SampleID")
meth <- filter(meth, SampleType == "Soil")

ggplot(meth, aes(x = New_Treatment, y = cov.sum)) +
  geom_boxplot() +
  facet_wrap(~PlotTreatment)

hist(meth$cov.sum)
m1 <- lm(cov.sum ~ New_Treatment * PlotTreatment, data = meth)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
anova(m1)

# nothing going on with methylophilaceae in terms of gene coverage
# burkholderiaceae is lower in mixes than in forbs (maybe grasses too), in both control and drought but this difference disappears in watering


det.g2 <- filter(det.g, key %in% kaiju.tax3$gene_call)
colnames(det.g2) <- sub("_[^_]+$", "", colnames(det.g2))
det.g2 <- as.data.frame(t.data.frame(det.g2))
det.g2$det.sum <- rowSums(det.g2)
det.g2$SampleID <- row.names(det.g2)
det.g2 <- det.g2[-1, 515:516]

meth <- merge(metadata, det.g2, by = "SampleID")
meth <- filter(meth, SampleType == "Soil")

ggplot(meth, aes(x = New_Treatment, y = det.sum)) +
  geom_boxplot() +
  facet_wrap(~PlotTreatment)

hist(sqrt(meth$det.sum))
m1 <- lm(sqrt(det.sum) ~ New_Treatment * PlotTreatment, data = meth)
plot(fitted(m1), resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
anova(m1)
