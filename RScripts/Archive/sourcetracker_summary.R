#summarize sourcetracker results
library(vroom)
library(tidyverse)
library(magrittr)

meta <- vroom("McL_metag_metadata.csv")

sp_euk_asv_names <- t(read.delim("sourcetracker2_species_euk_06052022_rare150/mixing_proportions.txt", row.names = 1))
sp_euk_asv <- as_tibble(sp_euk_asv_names)
read_euk_asv <- as_tibble(t(read.delim("sourcetracker2_read_euk_06052022_rare145/mixing_proportions.txt", row.names = 1)))
fam_euk_asv <- as_tibble(t(read.delim("sourcetracker2_family_euk_06052022_rare150/mixing_proportions.txt", row.names = 1)))

sp_arc_asv_names <- t(read.delim("sourcetracker2_species_arc_06052022_rare350/mixing_proportions.txt", row.names = 1))
sp_arc_asv <-  as_tibble(sp_arc_asv_names)
read_arc_asv <- as_tibble(t(read.delim("sourcetracker2_read_arc_06052022_rare200/mixing_proportions.txt", row.names = 1)))
family_arc_asv <- as_tibble(t(read.delim("sourcetracker2_family_arc_06052022_rare400/mixing_proportions.txt", row.names = 1)))

sp_all_asv_names <- t(read.delim("sourcetracker2_species_all_06052022_rare20000/mixing_proportions.txt", row.names = 1))
sp_all_asv <-as_tibble(sp_all_asv_names)
fam_all_asv <- as_tibble(t(read.delim("sourcetracker2_family_all_06052022_rare20000/mixing_proportions.txt", row.names = 1)))
read_all_asv <- as_tibble(t(read.delim("sourcetracker2_read_all_06052022_rare18000/mixing_proportions.txt", row.names = 1)))

sp_bac_asv_names <- t(read.delim("sourcetracker2_species_bac_06052022_rare20000/mixing_proportions.txt", row.names = 1))
fam_bac_asv <- as_tibble(t(read.delim("sourcetracker2_family_bac_06052022_rare20000/mixing_proportions.txt", row.names = 1)))
sp_bac_asv <- as_tibble(t(read.delim("sourcetracker2_species_bac_06052022_rare20000/mixing_proportions.txt", row.names = 1)))
read_bac_asv <- as_tibble(t(read.delim("sourcetracker2_read_bac_06052022_rare15000/mixing_proportions.txt", row.names = 1)))

row_names <- rownames(sp_euk_asv_names)

sp_euk_asv$Sample <- row_names
read_euk_asv$Sample <- row_names
fam_euk_asv$Sample <- row_names

row_names <- rownames(sp_arc_asv_names)

sp_arc_asv$Sample <- row_names 
read_arc_asv$Sample <- row_names
family_arc_asv$Sample <- row_names

row_names <- rownames(sp_all_asv_names)

sp_all_asv$Sample <- row_names
fam_all_asv$Sample <- row_names
read_all_asv$Sample <- row_names

row_names <- rownames(sp_bac_asv_names)

fam_bac_asv$Sample <- row_names
sp_bac_asv$Sample <- row_names
read_bac_asv$Sample <- row_names

sp_euk_asv_met <- left_join(sp_euk_asv, meta, by = c("Sample" = "SampleID"))
#read_euk_asv_met <- left_join(read_euk_asv, meta, by = c("Sample" = "SampleID"))
fam_euk_asv_met <- left_join(fam_euk_asv, meta, by = c("Sample" = "SampleID"))
sp_arc_asv_met <- left_join(sp_arc_asv, meta, by = c("Sample" = "SampleID")) 
#read_arc_asv_met <- left_join(read_arc_asv, meta, by = c("Sample" = "SampleID"))
family_arc_asv_met <- left_join(family_arc_asv, meta, by = c("Sample" = "SampleID"))
sp_all_asv_met <- left_join(sp_all_asv, meta, by = c("Sample" = "SampleID"))
fam_all_asv_met <- left_join(fam_all_asv, meta, by = c("Sample" = "SampleID"))
#read_all_asv_met <- left_join(read_all_asv, meta, by = c("Sample" = "SampleID"))
fam_bac_asv_met <- left_join(fam_bac_asv, meta, by = c("Sample" = "SampleID"))
sp_bac_asv_met <- left_join(sp_bac_asv, meta, by = c("Sample" = "SampleID"))
#read_bac_asv_met <- left_join(read_bac_asv, meta, by = c("Sample" = "SampleID"))


### bacteria

grouped_res.fg <- group_by(fam_bac_asv_met, PlotTreatment)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res.fg.mean <- summarise(grouped_res.fg, 
                         grass=100*mean(grasses), 
                         forb=100*mean(forb), 
                         unknown=100*mean(Unknown))

avgs_res.fg.sd <- summarise(grouped_res.fg, 
                              grass=100*sd(grasses)/sqrt(length(grasses)),
                              forb=100*sd(forb)/sqrt(length(forb)),
                              unknown=100*sd(Unknown)/sqrt(length(Unknown)))


test <- avgs_res.fg.mean %>%
  gather(key = "source", value = "mean", -PlotTreatment)

test.2 <- avgs_res.fg.sd %>%
  gather(key="source", value="sd", -PlotTreatment)

test3 <- full_join(test,test.2)

wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"]) 

wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 

wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"]) 

sig.text <- as_data_frame(test3$source)
sig.text$cld <- c("a", "b", "a", # grass
                      "a", "a", "a", # forb
                      "b", "c", "b") # unknown
sig.text$max <- test3$mean + test3$sd


ggplot(test3, aes(x = source, y = mean, fill = source)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 14)) + 
  ylab("Mean % community derived from each source") + facet_wrap(~PlotTreatment)+
  scale_fill_manual(values = c("magenta4", "#1F968BFF", "#B8DE29FF")) + 
  xlab("") +
  geom_text(aes(x = sig.text$value, label = sig.text$cld, y = sig.text$max,
                  size = 8,
                  #vjust = 0.3),
                  vjust = -0.7),
            col = "black", show.legend = F) +
  ylim(NA, max(test3$mean + test3$sd) + 0.055*max(test3$mean))

  #+ theme(legend.position = "none")

ggsave("bacterial_family_sourcetracker2.png", dpi=300, device = "png", height = 5, width= 6, units = "in")



## ok but what about how forbs change with treatment?

#forbs, not sig
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"]) 

#grasses, not sig
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"]) 

#unknown, not sig
wilcox.test(grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 
wilcox.test(grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"]) 
wilcox.test(grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 



### archaea

grouped_res.fg <- group_by(family_arc_asv_met, PlotTreatment)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res.fg.mean <- summarise(grouped_res.fg, 
                              grass=100*mean(grasses), 
                              forb=100*mean(forb), 
                              unknown=100*mean(Unknown))

avgs_res.fg.sd <- summarise(grouped_res.fg, 
                            grass=100*sd(grasses)/sqrt(length(grasses)),
                            forb=100*sd(forb)/sqrt(length(forb)),
                            unknown=100*sd(Unknown)/sqrt(length(Unknown)))


test <- avgs_res.fg.mean %>%
  gather(key = "source", value = "mean", -PlotTreatment)

test.2 <- avgs_res.fg.sd %>%
  gather(key="source", value="sd", -PlotTreatment)

test3 <- full_join(test,test.2)

wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"]) 

wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 

wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"]) 
wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"]) 
wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"]) 

sig.text <- as_data_frame(test3$source)
sig.text$cld <- c("a", "a", "a", # grass
                  "a", "a", "a", # forb
                  "b", "b", "b") # unknown
sig.text$max <- test3$mean + test3$sd


ggplot(test3, aes(x = source, y = mean, fill = source)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
        text = element_text(size = 14)) + 
  ylab("Mean % community derived from each source") + facet_wrap(~PlotTreatment)+
  scale_fill_manual(values = c("magenta4", "#1F968BFF", "#B8DE29FF")) + 
  xlab("") +
  geom_text(aes(x = sig.text$value, label = sig.text$cld, y = sig.text$max,
                size = 8,
                #vjust = 0.3),
                vjust = -0.7),
            col = "black", show.legend = F) +
  ylim(NA, max(test3$mean + test3$sd) + 0.055*max(test3$mean))

#+ theme(legend.position = "none")

ggsave("archael_family_sourcetracker2.png", dpi=300, device = "png", height = 5, width= 6, units = "in")


## DEGs 

deg <- read_tsv("DEG.101722.txt")
#length(unique(deg$DEG)) 68 gene

deg$DEG <- paste0("gene_", deg$DEG)

physeq.deg <- prune_taxa(deg$DEG, physeq.forDEG)

physeq.deg.fam <- tax_glom(physeq.deg, "Family", NArm =FALSE)

physeq.deg.fam.melt <- psmelt(physeq.deg.fam)

physeq.deg.fam.melt_meta <- left_join(physeq.deg.fam.melt, only_treat, by=c("Sample"= "SampleID"))

grouped <- group_by(physeq.deg.fam.melt_meta, PlotTreatment, New_Treatment, Family)

avgs_grouped <- summarise(grouped, mean =  mean(Abundance), 
                           sd =  sd(Abundance))

plot = ggplot(avgs_grouped[avgs_grouped$New_Treatment == "mix",], aes(x = Family, y = (mean), fill = Family)) + 
  geom_bar(stat = "identity", position = "stack") + facet_wrap(~PlotTreatment)

plot = plot + theme(axis.text.x = element_text(angle = -70, hjust = 0, vjust = 0.5), text = element_text(size = 18)) + 
  ylab("Mean Relative Abundance") + xlab("Family") +
  scale_fill_viridis_d() #+  theme(legend.position = "none")









## what taxa ## 

file_dir <- "sourcetracker2_family_bac_06052022_rare20000"
files <- list.files(file_dir, full.names = TRUE, pattern = "feature_table.txt$")
combined_files <- bind_rows(lapply(files, read.delim), .id = "id")

samples <- str_replace(files, "sourcetracker2_family_bac_06052022_rare20000/", "")
sample_name <- str_replace(samples, ".feature_table.txt", "")
id <- as.character(1:29)

sample_df <- data.frame(sample_name, id)

combined_tax <- left_join(combined_files, sample_df)

#"Veillonellaceae", "Fibrobacteraceae", "Methylophilaceae" not in file 
#combined_tax_subset <- combined_tax[c("sample_name", "X", "Clostridiaceae", "Burkholderiaceae", "Rhodocyclaceae")] 

only_treat <- meta[c(1,7)]

#combined_tax_subset_meta <- left_join(combined_tax_subset, only_treat, by=c("sample_name"= "SampleID"))
combined_tax_subset_meta <- left_join(combined_tax, only_treat, by=c("sample_name"= "SampleID"))

#combined_tax_subset_meta_gather <- combined_tax_subset_meta %>%
#  gather(key="Family", value = "Counts", -sample_name, -X, -PlotTreatment)
combined_tax_subset_meta_gather <- combined_tax_subset_meta %>%
  gather(key="Family", value = "Counts", -id, -sample_name, -X, -PlotTreatment)


summed <- combined_tax_subset_meta_gather %>% 
  group_by(Family, PlotTreatment) %>%
  summarize(sum=sum(as.numeric(Counts)))

summed_2 <- combined_tax_subset_meta_gather %>% 
  group_by(Family, X, PlotTreatment) %>%
  summarize(RA = 100*Counts/20000, prop=Counts)

summed_3 <- left_join(summed_2, summed)

summed_4 <- summed_3 %>%
  group_by(Family, X, PlotTreatment) %>%
  summarize(mean_prop = sum(prop/sum))

fams <- summed_4 %>%
  filter(X == "forb") %>%
  filter(mean_prop > 0.7)

summed_5 <- summed_4 %>%
  filter(Family %in% fams$Family)

ggplot(summed_5, aes(y = mean_prop, x = Family, col = X, fill = X)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() + facet_grid(~PlotTreatment)+
  scale_fill_manual(values = c("magenta4", "#1F968BFF", "#B8DE29FF")) +
  labs(y = "Proportion of Families") +
  coord_flip()



## 

summed_4$X <- factor(summed_4$X, levels = c( "Unknown", "forb", "grasses"))

summed_5 <- summed_4 %>%
  filter(Family %in% physeq.deg.fam.melt$Family) %>%
  arrange(desc(mean_prop))


ggplot(summed_5, aes(y = mean_prop, x = Family, col = X, fill = X)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() + facet_grid(~PlotTreatment)+
  scale_fill_manual(name="source", values = c("#B8DE29FF","magenta4", "#1F968BFF"), labels = c("unknown", "forb", "grass")) +
  labs(y = "Proportion of Families") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "white") +
  theme_classic() +
  coord_flip()  

ggsave("bacterial_family_sourcetracker2_deg.png", dpi=300, device = "png", height = 6, width= 9, units = "in")

