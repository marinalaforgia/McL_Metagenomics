#microbiome reads

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

#read RDS in
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
kaiju_tax_file <- readRDS("GL_kaiju_tax_ready_for_source_09132021.RDS")

kaiju_tax_file <- as.data.frame(kaiju_tax_file)
row.names(kaiju_tax_file) <- paste0("Read", kaiju_tax_file$ReadName)

tax <- kaiju_tax_file[-c(1)]

#deseq normalized counts from raw data (so still have negative controls, etc)
norm_counts <- readRDS('GL_normcounts.10142021.RDS')
norm_counts <- as.data.frame(norm_counts)
norm_counts$ReadName <- as.character(row.names(norm_counts))
row.names(norm_counts) <- paste0("Read", norm_counts$ReadName)

otu <- as.matrix(norm_counts[-c(70)]) 
tax <- as.matrix(kaiju_tax_file[-c(1)])

reads <- row.names(otu)

tax.sub <- subset(tax, rownames(tax) %in% reads)

OTU = otu_table(otu, taxa_are_rows = TRUE)
TAX = tax_table(tax.sub)

#phyloseq
library("phyloseq")
physeq = phyloseq(OTU, TAX)
# 
# 

df <- as.data.frame(otu_table(physeq))
keep = colnames(df[-c(1,2,3)])

physeq <- prune_samples(keep, physeq)

physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)

physeq.bac <- subset_taxa(physeq, Domain == "Bacteria")
physeq.arc <- subset_taxa(physeq, Domain == "Archaea")
physeq.euk <- subset_taxa(physeq, Domain == "Eukaryota")
physeq.vir <- subset_taxa(physeq, Domain == "Viruses")

physeq <- filter_taxa(physeq, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
physeq <- rarefy_even_depth(physeq, 20000, replace = FALSE, rngseed = 5311)

physeq.bac <- rarefy_even_depth(physeq.bac, 50000, replace = TRUE, rngseed = 5311)
physeq.arc <- rarefy_even_depth(physeq.arc, 800, replace = TRUE, rngseed = 5311)
physeq.euk <- rarefy_even_depth(physeq.euk, 300, replace = TRUE, rngseed = 5311)
#physeq.vir <- rarefy_even_depth(physeq.vir, 20000, replace = FALSE, rngseed = 5311)
#virus = too small

#whole dataset
physeq_tax_F <- tax_glom(physeq, taxrank = "Family")
physeq_tax_S <- tax_glom(physeq, taxrank = "Species")
physeq #Read level

#generate asv table formatted for biom generation
asv.tab <- otu_table(physeq_tax_S)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq_tax_S))
rownames(asv.tab.2) <- tax.tab$Species
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#read level
asv.tab <- otu_table(physeq)
asv.tab <- as.data.frame(asv.tab)
suppressWarnings(asv.tab <- as.matrix(asv.tab))
cb <- as.matrix(cbind(rownames(asv.tab), asv.tab))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#write tables for biom generation using biom
#asv.tab
write.table(x = rcb, file = "read_all_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "family_all_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "species_all_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# #generate tax table formatted for biom generation
# tax.tab$taxonomy <- tidyr::unite_(tax.tab, "out", c(colnames(tax.tab)), sep = ";")
# cbt <- as.matrix(cbind(rownames(tax.tab), tax.tab$taxonomy))
# rcbt <- as.matrix(rbind(c("#ASVID", "taxonomy"), cbt))
# rownames(cbt) <- NULL
# colnames(cbt) <- NULL
# 
# #tax.tab
# write.table(x = rcbt, file = "read_all_tax_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#bacteria only

physeq.bac_F <- tax_glom(physeq.bac, taxrank = "Family")
physeq.bac_S <- tax_glom(physeq.bac, taxrank = "Species")
physeq.bac #Read level

#generate asv table formatted for biom generation
asv.tab <- otu_table(physeq.bac_S)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq.bac_S))
rownames(asv.tab.2) <- tax.tab$Species
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#read level
asv.tab <- otu_table(physeq.bac)
asv.tab <- as.data.frame(asv.tab)
suppressWarnings(asv.tab <- as.matrix(asv.tab))
cb <- as.matrix(cbind(rownames(asv.tab), asv.tab))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#write tables for biom generation using biom
#asv.tab
write.table(x = rcb, file = "read_bac_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "family_bac_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "species_bac_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#archea only
physeq.arc_F <- tax_glom(physeq.arc, taxrank = "Family")
physeq.arc_S <- tax_glom(physeq.arc, taxrank = "Species")
physeq.arc #Read level

#generate asv table formatted for biom generation
asv.tab <- otu_table(physeq.arc_S)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq.arc_S))
rownames(asv.tab.2) <- tax.tab$Species
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#read level
asv.tab <- otu_table(physeq.arc)
asv.tab <- as.data.frame(asv.tab)
suppressWarnings(asv.tab <- as.matrix(asv.tab))
cb <- as.matrix(cbind(rownames(asv.tab), asv.tab))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#write tables for biom generation using biom
#asv.tab
write.table(x = rcb, file = "read_arc_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "family_arcc_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "species_arc_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")




#euks only
physeq.euk_F <- tax_glom(physeq.euk, taxrank = "Family")
physeq.euk_S <- tax_glom(physeq.euk, taxrank = "Species")
physeq.euk #Read level

#generate asv table formatted for biom generation
asv.tab <- otu_table(physeq.euk_S)
asv.tab.2 <- as.data.frame(asv.tab)
tax.tab <- as.data.frame(phyloseq::tax_table(physeq.euk_S))
rownames(asv.tab.2) <- tax.tab$Species
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
suppressWarnings(asv.tab.2 <- as.matrix(asv.tab.2))
cb <- as.matrix(cbind(rownames(asv.tab.2), asv.tab.2))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#read level
asv.tab <- otu_table(physeq.euk)
asv.tab <- as.data.frame(asv.tab)
suppressWarnings(asv.tab <- as.matrix(asv.tab))
cb <- as.matrix(cbind(rownames(asv.tab), asv.tab))
rcb <- as.matrix(rbind(colnames(cb), cb))
rcb[1,1] <- "#ASVID"
rownames(rcb) <- NULL
colnames(rcb) <- NULL

#write tables for biom generation using biom
#asv.tab
write.table(x = rcb, file = "read_euk_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "family_euk_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(x = rcb, file = "species_euk_asv_tab.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


#biom convert
#biom convert -i species_name_all_asv_tab.txt -o species_name_all_asv_tab.biom --table-type="OTU table" --to-hdf5

#sourcetracker commands python
# sourcetracker2 -i species_euk_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_species_euk_10292021_rare150/ --jobs 5 --sink_rarefaction_depth 150 --source_rarefaction_depth 150 --per_sink_feature_assignments 
# sourcetracker2 -i read_euk_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_read_euk_10292021_rare80/ --jobs 5 --sink_rarefaction_depth 80 --source_rarefaction_depth 80 --per_sink_feature_assignments 
# sourcetracker2 -i family_euk_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_family_euk_10292021_rare150/ --jobs 5 --sink_rarefaction_depth 150 --source_rarefaction_depth 150 --per_sink_feature_assignments 
# 
# sourcetracker2 -i species_arc_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_species_arc_10292021_rare800/ --jobs 5 --sink_rarefaction_depth 700 --source_rarefaction_depth 700 --per_sink_feature_assignments 
# sourcetracker2 -i read_arc_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_read_arc_10292021_rare80/ --jobs 5 --sink_rarefaction_depth 80 --source_rarefaction_depth 80 --per_sink_feature_assignments 
# sourcetracker2 -i family_arcc_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_family_arc_10292021_rare700/ --jobs 5 --sink_rarefaction_depth 700 --source_rarefaction_depth 700 --per_sink_feature_assignments 
# 
# 
# 
# sourcetracker2 -i family_all_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_family_all_10292021_rare15000/ --jobs 5 --sink_rarefaction_depth 15000 --source_rarefaction_depth 15000 --per_sink_feature_assignments 
# sourcetracker2 -i species_all_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_species_all_10292021_rare15000/ --jobs 5 --sink_rarefaction_depth 15000 --source_rarefaction_depth 15000 --per_sink_feature_assignments 
# sourcetracker2 -i read_all_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_read_all_10292021_rare1000/ --jobs 5 --sink_rarefaction_depth 1000 --source_rarefaction_depth 1000 --per_sink_feature_assignments 
# 
# sourcetracker2 -i family_bac_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_family_bac_10292021_rare25000/ --jobs 5 --sink_rarefaction_depth 25000 --source_rarefaction_depth 25000 --per_sink_feature_assignments 
# sourcetracker2 -i species_bac_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_species_bac_10292021_rare25000/ --jobs 5 --sink_rarefaction_depth 25000 --source_rarefaction_depth 25000 --per_sink_feature_assignments 
# sourcetracker2 -i read_bac_asv_tab.txt -m sourcetracker_meta.txt -o sourcetracker2_read_bac_10292021_rare7500/ --jobs 5 --sink_rarefaction_depth 7500 --source_rarefaction_depth 7500 --per_sink_feature_assignments 
# 



