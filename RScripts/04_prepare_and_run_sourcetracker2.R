#sourcetracker preparation in R

library(vroom)
library(tidyverse)
library(magrittr)

#read in Kaiju outputs
kaiju <-  read_tsv("../Kaiju.taxonomy.table.for.anvio.2019", col_names = c('Classified', 'ReadName', 'NCBI_ID', 'Score', 'TaxonID_Best_Match', 'Assession_Best_Match', 'Sequence_Match', 'Taxonomy'))

#split Taxonomy in kaiju files
kaiju_tax <- kaiju %>% separate('Taxonomy', c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep='; ')  %>% separate('Species', c('Species'), sep=';')

#saveRDS(kaiju_tax, "GL_kaiju_tax_09012021.RDS")

summary(as.factor(kaiju_tax$Classified))
# C       U 
# 6570442 2258838 

#read RDS in and then fix taxonomy
kaiju_tax <- readRDS("GL_kaiju_tax_09012021.RDS")

kaiju_tax %<>% mutate(Phylum = ifelse(Phylum == "NA" & Domain != "NA", paste("Unclassified", as.character(Domain)), as.character(Phylum)),
                               Class = ifelse(Class == "NA" & Phylum != "NA",  as.character(Phylum), as.character(Class)),
                               Order = ifelse(Order == "NA" & Class != "NA", as.character(Class), as.character(Order)),
                               Family = ifelse(Family == "NA" & Order != "NA",  as.character(Order), as.character(Family)),
                               Genus = ifelse(Genus == "NA" & Family != "NA",  as.character(Family), as.character(Genus)),
                               Species = ifelse(Species == "NA" & Genus != "NA",  as.character(Genus), as.character(Species)))

kaiju_tax %<>% mutate(Domain = ifelse(Domain == "NA", "Unclassified", as.character(Domain)),
                               Phylum = ifelse(Phylum == "NA", "Unclassified", as.character(Phylum)),
                               Class = ifelse(Class == "NA", "Unclassified", as.character(Class)),
                               Order = ifelse(Order == "NA", "Unclassified", as.character(Order)),
                               Family = ifelse(Family == "NA", "Unclassified", as.character(Family)),
                               Genus = ifelse(Genus == "NA", "Unclassified", as.character(Genus)),
                               Species = ifelse(Species == "NA", "Unclassified", as.character(Species)))


kaiju_tax %<>% mutate(Domain = fct_explicit_na(Domain, na_level = "Unclassified"), 
                      Phylum = fct_explicit_na(Phylum, na_level = "Unclassified"), 
                      Class = fct_explicit_na(Class, na_level = "Unclassified"), 
                      Order = fct_explicit_na(Order, na_level = "Unclassified"), 
                      Family = fct_explicit_na(Family, na_level = "Unclassified"), 
                      Genus = fct_explicit_na(Genus, na_level = "Unclassified"), 
                      Species = fct_explicit_na(Species, na_level = "Unclassified"))

kaiju_tax_file <- kaiju_tax[,-c(1,3,4,5,6,7)]
#saveRDS(kaiju_tax_file, "GL_kaiju_tax_ready_for_source_09132021.RDS")

## start here 

kaiju_tax_file <- readRDS("GL_kaiju_tax_ready_for_source_09132021.RDS")

kaiju_tax_file <- as.data.frame(kaiju_tax_file)
row.names(kaiju_tax_file) <- paste0("gene_", kaiju_tax_file$ReadName)

tax <- kaiju_tax_file[-c(1)]

#deseq normalized counts from raw data (so still have negative controls, etc)
norm_counts <- readRDS("decontam_norm_counts_w-controls-12082022.RDS")
norm_counts <- as.data.frame(norm_counts)
norm_counts$ReadName <- as.character(row.names(norm_counts))
row.names(norm_counts) <- paste0("gene_", norm_counts$ReadName)

#filter out genes we aren't including 
norm_counts <- norm_counts[!rownames(norm_counts) %in% genes_to_remove,]

norm_counts <- norm_counts %>%
  filter(ReadName < 8829098)

norm_counts <- norm_counts[rownames(norm_counts) %in% GL_fun_table_filt.noCOG$gene,]


#import into phyloseq
otu <- as.matrix(norm_counts[-c(70)]) 
tax <- as.matrix(kaiju_tax_file[-c(1)])

reads <- row.names(otu)

tax.sub <- subset(tax, rownames(tax) %in% reads)

#phyloseq
library("phyloseq")
OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax.sub)

physeq = phyloseq(OTU, TAX)

physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
physeq.forDEG <- physeq

physeq.bac <- subset_taxa(physeq, Domain == "Bacteria")
physeq.arc <- subset_taxa(physeq, Domain == "Archaea")


#whole dataset
physeq_tax_F <- tax_glom(physeq, taxrank = "Family") #650 total families

## family
asv.tab <- otu_table(physeq_tax_F)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq_tax_F))
rownames(asv.tab.2) <- tax.tab$Family
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL


write.table(x = rcb, file = "family_all_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



#bacteria only

physeq.bac_F <- tax_glom(physeq.bac, taxrank = "Family") #598 families

#generate asv table formatted for biom generation
asv.tab <- otu_table(physeq.bac_F)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq.bac_F))
rownames(asv.tab.2) <- tax.tab$Family
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

write.table(x = rcb, file = "family_bac_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#archea only
physeq.arc_F <- tax_glom(physeq.arc, taxrank = "Family") #51 family


#generate asv table formatted for biom generation
asv.tab <- otu_table(physeq.arc_F)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq.arc_F))
rownames(asv.tab.2) <- tax.tab$Family
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

write.table(x = rcb, file = "family_arcc_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



#no rarifing first
#sourcetracker commands python
#conda activate st2

# sourcetracker2 -i family_arcc_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_family_arc_12142022_rare600/ --jobs 5 --sink_rarefaction_depth 600 --source_rarefaction_depth 600 --per_sink_feature_assignments
# 
# sourcetracker2 -i family_all_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_family_all_12142022_rare25000/ --jobs 5 --sink_rarefaction_depth 25000 --source_rarefaction_depth 25000 --per_sink_feature_assignments
# 
# sourcetracker2 -i family_bac_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_family_bac_1214022_rare25000/ --jobs 5 --sink_rarefaction_depth 25000 --source_rarefaction_depth 25000 --per_sink_feature_assignments
# 




#function as taxonomy


func <- as.data.frame(GL_fun_table_filt.noCOG[c(1,11,10)])
rownames(func) <- func$gene
func <- func[-c(1)]


otu <- as.matrix(norm_counts[-c(70)]) 
func <- as.matrix(func)

reads <- row.names(otu)

func.sub <- subset(func, rownames(func) %in% reads)

#phyloseq
library("phyloseq")
OTU = otu_table(otu, taxa_are_rows = TRUE)
FUNC = tax_table(func.sub)

physeq.cog20 = phyloseq(OTU, FUNC)

keep = colnames(df[-c(1,2,3)])
physeq.cog20 <- prune_samples(keep, physeq.cog20)
physeq.cog20 <- prune_taxa(taxa_sums(physeq.cog20) > 0, physeq.cog20)


physeq.cog20_F <- tax_glom(physeq.cog20, taxrank = "COG20_FUNCTION") #4616 functions

#generate asv table formatted for biom generation
asv.tab <- otu_table(physeq.cog20_F)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq.cog20_F))
rownames(asv.tab.2) <- tax.tab$COG20_FUNCTION
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

write.table(x = rcb, file = "cogfucn_all_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#sourcetracker2 -i cogfucn_all_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_cog20_all_12212022_rare25000/ --jobs 5 --sink_rarefaction_depth 25000 --source_rarefaction_depth 25000 --per_sink_feature_assignments

