#requires vegan, phyloseq, ggplot & R

library(ggplot2)
library(vegan)
library(phyloseq)
library(ggrepel)
library(patchwork)

#get data into phyloseq

#make sure any random things 
set.seed(5311)

# Taxonomy - Natives
#ps_gene_lvl <- readRDS("arrowplots/ps-nocontrols-rare-66880.RDS")
ps_txn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-tax.RDS")
#ps_fn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-fun.RDS")

ps_gene_lvl <- ps_txn_lvl
#ps_gene_lvl <- ps_fn_lvl

ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1", "native mix 2"))
#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1 + invasives", "native mix 2 + invasives"))

samp_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(AllPlantSpeciesInSubPlot = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInSubPlot)),
         Bromus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Bromus"), 1, 0),
         Avena = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Avena"), 1, 0),
         Calycadenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Calycadenia"), 1, 0),
         Hemizonia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Hemizonia"), 1, 0),
         Elymus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Elymus"), 1, 0),
         Plantago = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Plantago"), 1, 0),
         Clarkia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Clarkia"), 1, 0),
         Agoseris = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Agoseris"), 1, 0),
         Lasethenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Lasethenia"), 1, 0),
         'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU"), 1, 0),
         'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "VUMI"), 1, 0)) %>%
  select(SampleID, Bromus, Avena, Calycadenia, Hemizonia, Elymus,Plantago,Clarkia,Agoseris,Lasethenia,'Festuca perennis','Festuca microstachys') %>%
  arrange(SampleID)

samp_data <- samp_data[, colSums(samp_data != 0) > 0]

binding_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  select(SampleID, PlotXTreatment, PlotSubTreatmentName, WaterTrt) %>%
  arrange(SampleID)

#make metadata matrix; only has the environmental data you wish to correlate with and has no NA's or blanks
env=as.matrix(samp_data[-c(1)])
rownames(env) <- samp_data$SampleID        

#compute vegan dist & calculate pcoa values (must use bray curtis)

ps_gene_lvl.dist=vegdist(t(ps_gene_lvl@otu_table),method="bray")
ps_gene_lvl_pcoa_res =cmdscale(ps_gene_lvl.dist, eig = TRUE)
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa_res$points

var_exp <- round(ps_gene_lvl_pcoa_res$eig*100/sum(ps_gene_lvl_pcoa_res$eig),1)

adonis2(ps_gene_lvl.dist ~ ., data=as_data_frame(env), by="margin", permutations = 9999) #not sig

adonis2(ps_gene_lvl.dist ~ PlotSubTreatmentName + PlotXTreatment, data=binding_data, by="margin", permutations = 9999) #sig

#reorder the samples by SampleID b/c of coloring factors
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa[order(row.names(ps_gene_lvl_pcoa)), ]

#fit metadata to pcoa values
ps_gene_lvl_pcoa_env = envfit(ps_gene_lvl_pcoa,env)

#start buliding pcoa plot 

#makes a dataframe that will be used to add the pcoa sample points to plot
scrs <- as.data.frame(scores(ps_gene_lvl_pcoa, display = "sites"))
scrs$SampleID <- rownames(scrs)

#makes factors to color the plot by, these will color your samples, very important that the order matches
scrs <- left_join(scrs, binding_data)

#makes a dataframe that will be used to add metadata factors from envfit to plot as arrows
spp.scrs <- as.data.frame(scores(ps_gene_lvl_pcoa_env, display = "vectors"))

#have to name the arrows, in order they are in metadata file, will error if you don't name them all
spp.scrs$Species <- rownames(spp.scrs)

#add pvals 
spp.scrs$padj <- as_data_frame(p.adjust(ps_gene_lvl_pcoa_env$vectors$pvals, method='BH'))$value
spp.scrs$pvals <- as_data_frame(ps_gene_lvl_pcoa_env$vectors$pvals)$value


#make combined plot of pcoa + arrows
p_tax_natives <- ggplot(scrs) + 
  geom_point(mapping = aes(x = Dim1, y = Dim2, colour = PlotSubTreatmentName, shape=PlotSubTreatmentName), size=4) + 
  theme_bw() +
  geom_segment(data = spp.scrs, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.35, "cm")), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black")) + 
  geom_text_repel(data = spp.scrs, aes(x = Dim1, y = Dim2, label = Species), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black"), size = 4) +
  xlab(paste0("PCoA 1: ",var_exp[1],"% variance")) +
  ylab(paste0("PCoA 2: ",var_exp[2],"% variance")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("magenta4", "magenta"))

p_tax_natives

##stats

#calculate distance for metadata
env.dist = vegdist(env,method="euclidean")

#run mantel test to see if correlated with microbial distance (bray curtis must be used)
mantel(env.dist,ps_gene_lvl.dist, permutations = 9999) #not sig


## Taxonomy - Mixes

#ps_gene_lvl <- readRDS("arrowplots/ps-nocontrols-rare-66880.RDS")
ps_txn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-tax.RDS")
#ps_fn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-fun.RDS")

ps_gene_lvl <- ps_txn_lvl
#ps_gene_lvl <- ps_fn_lvl

#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1", "native mix 2"))
ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1 + invasives", "native mix 2 + invasives"))

samp_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(AllPlantSpeciesInSubPlot = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInSubPlot)),
         Bromus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Bromus"), 1, 0),
         Avena = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Avena"), 1, 0),
         Calycadenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Calycadenia"), 1, 0),
         Hemizonia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Hemizonia"), 1, 0),
         Elymus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Elymus"), 1, 0),
         Plantago = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Plantago"), 1, 0),
         Clarkia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Clarkia"), 1, 0),
         Agoseris = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Agoseris"), 1, 0),
         Lasethenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Lasethenia"), 1, 0),
         'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU"), 1, 0),
         'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "VUMI"), 1, 0)) %>%
  select(SampleID, Bromus, Avena, Calycadenia, Hemizonia, Elymus,Plantago,Clarkia,Agoseris,Lasethenia,'Festuca perennis','Festuca microstachys') %>%
  arrange(SampleID)

samp_data <- samp_data[, colSums(samp_data != 0) > 0]

binding_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  select(SampleID, PlotXTreatment, PlotSubTreatmentName, WaterTrt) %>%
  arrange(SampleID)

#make metadata matrix; only has the environmental data you wish to correlate with and has no NA's or blanks
env=as.matrix(samp_data[-c(1)])
rownames(env) <- samp_data$SampleID        

#compute vegan dist & calculate pcoa values (must use bray curtis)

ps_gene_lvl.dist=vegdist(t(ps_gene_lvl@otu_table),method="bray")
ps_gene_lvl_pcoa_res =cmdscale(ps_gene_lvl.dist, eig = TRUE)
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa_res$points

var_exp <- round(ps_gene_lvl_pcoa_res$eig*100/sum(ps_gene_lvl_pcoa_res$eig),1)

res <- adonis2(ps_gene_lvl.dist ~ ., data=as_data_frame(env), by="margin", permutations = 9999) #avena, caly, clarkia sig
p.adjust(res$`Pr(>F)`[1:length(colnames(env))], 'BH') < 0.05 #no sig after pval correction

adonis2(ps_gene_lvl.dist ~ PlotSubTreatmentName + PlotXTreatment, data=binding_data, by="margin", permutations = 9999) #sig

#reorder the samples by SampleID b/c of coloring factors
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa[order(row.names(ps_gene_lvl_pcoa)), ]

#fit metadata to pcoa values
ps_gene_lvl_pcoa_env = envfit(ps_gene_lvl_pcoa,env)

#start buliding pcoa plot 

#makes a dataframe that will be used to add the pcoa sample points to plot
scrs <- as.data.frame(scores(ps_gene_lvl_pcoa, display = "sites"))
scrs$SampleID <- rownames(scrs)

#makes factors to color the plot by, these will color your samples, very important that the order matches
scrs <- left_join(scrs, binding_data)

#makes a dataframe that will be used to add metadata factors from envfit to plot as arrows
spp.scrs <- as.data.frame(scores(ps_gene_lvl_pcoa_env, display = "vectors"))

#have to name the arrows, in order they are in metadata file, will error if you don't name them all
spp.scrs$Species <- rownames(spp.scrs)

#add pvals 
spp.scrs$padj <- as_data_frame(p.adjust(ps_gene_lvl_pcoa_env$vectors$pvals, method='BH'))$value
spp.scrs$pvals <- as_data_frame(ps_gene_lvl_pcoa_env$vectors$pvals)$value

#make combined plot of pcoa + arrows
p_tax_mixes <- ggplot(scrs) + 
  geom_point(mapping = aes(x = Dim1, y = Dim2, colour = PlotSubTreatmentName, shape=PlotSubTreatmentName), size=4) + 
  theme_bw() +
  geom_segment(data = spp.scrs, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.35, "cm")), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black")) + 
  geom_text_repel(data = spp.scrs, aes(x = Dim1, y = Dim2, label = Species), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black"), size = 4) +
  xlab(paste0("PCoA 1: ",var_exp[1],"% variance")) +
  ylab(paste0("PCoA 2: ",var_exp[2],"% variance")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("goldenrod3", "goldenrod4"))

p_tax_mixes

##stats

#calculate distance for metadata
env.dist = vegdist(env,method="euclidean")

#run mantel test to see if correlated with microbial distance (bray curtis must be used)
mantel(env.dist,ps_gene_lvl.dist, permutations = 9999) #not sig



## Function - Natives

#ps_gene_lvl <- readRDS("arrowplots/ps-nocontrols-rare-66880.RDS")
#ps_txn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-tax.RDS")
ps_fn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-fun.RDS")

#ps_gene_lvl <- ps_txn_lvl
ps_gene_lvl <- ps_fn_lvl

ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1", "native mix 2"))
#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1 + invasives", "native mix 2 + invasives"))

samp_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(AllPlantSpeciesInSubPlot = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInSubPlot)),
         Bromus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Bromus"), 1, 0),
         Avena = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Avena"), 1, 0),
         Calycadenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Calycadenia"), 1, 0),
         Hemizonia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Hemizonia"), 1, 0),
         Elymus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Elymus"), 1, 0),
         Plantago = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Plantago"), 1, 0),
         Clarkia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Clarkia"), 1, 0),
         Agoseris = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Agoseris"), 1, 0),
         Lasethenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Lasethenia"), 1, 0),
         'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU"), 1, 0),
         'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "VUMI"), 1, 0)) %>%
  select(SampleID, Bromus, Avena, Calycadenia, Hemizonia, Elymus,Plantago,Clarkia,Agoseris,Lasethenia,'Festuca perennis','Festuca microstachys') %>%
  arrange(SampleID)

samp_data <- samp_data[, colSums(samp_data != 0) > 0]

binding_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  select(SampleID, PlotXTreatment, PlotSubTreatmentName, WaterTrt) %>%
  arrange(SampleID)

#make metadata matrix; only has the environmental data you wish to correlate with and has no NA's or blanks
env=as.matrix(samp_data[-c(1)])
rownames(env) <- samp_data$SampleID        

#compute vegan dist & calculate pcoa values (must use bray curtis)

ps_gene_lvl.dist=vegdist(t(ps_gene_lvl@otu_table),method="bray")
ps_gene_lvl_pcoa_res =cmdscale(ps_gene_lvl.dist, eig = TRUE)
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa_res$points

var_exp <- round(ps_gene_lvl_pcoa_res$eig*100/sum(ps_gene_lvl_pcoa_res$eig),1)

res <- adonis2(ps_gene_lvl.dist ~ ., data=as_data_frame(env), by="margin", permutations = 9999) #avena, caly, clarkia sig
p.adjust(res$`Pr(>F)`[1:length(colnames(env))], 'BH') < 0.05 #no sig after pval correction

adonis2(ps_gene_lvl.dist ~ PlotSubTreatmentName + PlotXTreatment, data=binding_data, by="margin", permutations = 9999) #sig

#reorder the samples by SampleID b/c of coloring factors
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa[order(row.names(ps_gene_lvl_pcoa)), ]

#fit metadata to pcoa values
ps_gene_lvl_pcoa_env = envfit(ps_gene_lvl_pcoa,env)

#start buliding pcoa plot 

#makes a dataframe that will be used to add the pcoa sample points to plot
scrs <- as.data.frame(scores(ps_gene_lvl_pcoa, display = "sites"))
scrs$SampleID <- rownames(scrs)

#makes factors to color the plot by, these will color your samples, very important that the order matches
scrs <- left_join(scrs, binding_data)

#makes a dataframe that will be used to add metadata factors from envfit to plot as arrows
spp.scrs <- as.data.frame(scores(ps_gene_lvl_pcoa_env, display = "vectors"))

#have to name the arrows, in order they are in metadata file, will error if you don't name them all
spp.scrs$Species <- rownames(spp.scrs)

#add pvals 
spp.scrs$padj <- as_data_frame(p.adjust(ps_gene_lvl_pcoa_env$vectors$pvals, method='BH'))$value
spp.scrs$pvals <- as_data_frame(ps_gene_lvl_pcoa_env$vectors$pvals)$value

#make combined plot of pcoa + arrows
p_fun_natives <- ggplot(scrs) + 
  geom_point(mapping = aes(x = Dim1, y = Dim2, colour = PlotSubTreatmentName, shape=PlotSubTreatmentName), size=4) + 
  theme_bw() +
  geom_segment(data = spp.scrs, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.35, "cm")), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black")) + 
  geom_text_repel(data = spp.scrs, aes(x = Dim1, y = Dim2, label = Species), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black"), size = 4) +
  xlab(paste0("PCoA 1: ",var_exp[1],"% variance")) +
  ylab(paste0("PCoA 2: ",var_exp[2],"% variance")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("magenta4", "magenta"))

p_fun_natives

##stats

#calculate distance for metadata
env.dist = vegdist(env,method="euclidean")

#run mantel test to see if correlated with microbial distance (bray curtis must be used)
mantel(env.dist,ps_gene_lvl.dist, permutations = 9999) #not sig

## Function - Mixes

#ps_gene_lvl <- readRDS("arrowplots/ps-nocontrols-rare-66880.RDS")
#ps_txn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-tax.RDS")
ps_fn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-fun.RDS")

#ps_gene_lvl <- ps_txn_lvl
ps_gene_lvl <- ps_fn_lvl

#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1", "native mix 2"))
ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1 + invasives", "native mix 2 + invasives"))

samp_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(AllPlantSpeciesInSubPlot = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInSubPlot)),
         Bromus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Bromus"), 1, 0),
         Avena = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Avena"), 1, 0),
         Calycadenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Calycadenia"), 1, 0),
         Hemizonia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Hemizonia"), 1, 0),
         Elymus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Elymus"), 1, 0),
         Plantago = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Plantago"), 1, 0),
         Clarkia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Clarkia"), 1, 0),
         Agoseris = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Agoseris"), 1, 0),
         Lasethenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Lasethenia"), 1, 0),
         'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU"), 1, 0),
         'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "VUMI"), 1, 0)) %>%
  select(SampleID, Bromus, Avena, Calycadenia, Hemizonia, Elymus,Plantago,Clarkia,Agoseris,Lasethenia,'Festuca perennis','Festuca microstachys') %>%
  arrange(SampleID)

samp_data <- samp_data[, colSums(samp_data != 0) > 0]

binding_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  select(SampleID, PlotXTreatment, PlotSubTreatmentName, WaterTrt) %>%
  arrange(SampleID)

#make metadata matrix; only has the environmental data you wish to correlate with and has no NA's or blanks
env=as.matrix(samp_data[-c(1)])
rownames(env) <- samp_data$SampleID        

#compute vegan dist & calculate pcoa values (must use bray curtis)

ps_gene_lvl.dist=vegdist(t(ps_gene_lvl@otu_table),method="bray")
ps_gene_lvl_pcoa_res =cmdscale(ps_gene_lvl.dist, eig = TRUE)
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa_res$points

var_exp <- round(ps_gene_lvl_pcoa_res$eig*100/sum(ps_gene_lvl_pcoa_res$eig),1)

res <- adonis2(ps_gene_lvl.dist ~ ., data=as_data_frame(env), by="margin", permutations = 9999) #avena, caly sig
p.adjust(res$`Pr(>F)`[1:length(colnames(env))], 'BH') < 0.05 #one marg sig after pval correction; Calycadenia

adonis2(ps_gene_lvl.dist ~ PlotSubTreatmentName + PlotXTreatment, data=binding_data, by="margin", permutations = 9999) #sig

#reorder the samples by SampleID b/c of coloring factors
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa[order(row.names(ps_gene_lvl_pcoa)), ]

#fit metadata to pcoa values
ps_gene_lvl_pcoa_env = envfit(ps_gene_lvl_pcoa,env)

#start buliding pcoa plot 

#makes a dataframe that will be used to add the pcoa sample points to plot
scrs <- as.data.frame(scores(ps_gene_lvl_pcoa, display = "sites"))
scrs$SampleID <- rownames(scrs)

#makes factors to color the plot by, these will color your samples, very important that the order matches
scrs <- left_join(scrs, binding_data)

#makes a dataframe that will be used to add metadata factors from envfit to plot as arrows
spp.scrs <- as.data.frame(scores(ps_gene_lvl_pcoa_env, display = "vectors"))

#have to name the arrows, in order they are in metadata file, will error if you don't name them all
spp.scrs$Species <- rownames(spp.scrs)

#add pvals 
spp.scrs$padj <- as_data_frame(p.adjust(ps_gene_lvl_pcoa_env$vectors$pvals, method='BH'))$value
spp.scrs$pvals <- as_data_frame(ps_gene_lvl_pcoa_env$vectors$pvals)$value

#make combined plot of pcoa + arrows
p_fun_mixes <- ggplot(scrs) + 
  geom_point(mapping = aes(x = Dim1, y = Dim2, colour = PlotSubTreatmentName, shape=PlotSubTreatmentName), size=4) + 
  theme_bw() +
  geom_segment(data = spp.scrs, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.35, "cm")), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black")) + 
  geom_text_repel(data = spp.scrs, aes(x = Dim1, y = Dim2, label = Species), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black"), size = 4) +
  xlab(paste0("PCoA 1: ",var_exp[1],"% variance")) +
  ylab(paste0("PCoA 2: ",var_exp[2],"% variance")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values=c("goldenrod3", "goldenrod4"))


p_fun_mixes

##stats

#calculate distance for metadata
env.dist = vegdist(env,method="euclidean")

#run mantel test to see if correlated with microbial distance (bray curtis must be used)
mantel(env.dist,ps_gene_lvl.dist, permutations = 9999) #not sig



(p_tax_natives + p_tax_mixes) / (p_fun_natives + p_fun_mixes) + plot_annotation(tag_levels = 'A')


ggsave("envfit_ordinations.png", dpi=300, device = "png", height = 8, width= 12, units = "in")
ggsave("envfit_ordinations.pdf", dpi=300, device = "pdf", height = 8, width= 12, units = "in")



## Taxonomy - Grasses

#ps_gene_lvl <- readRDS("arrowplots/ps-nocontrols-rare-66880.RDS")
ps_txn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-tax.RDS")
#ps_fn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-fun.RDS")

ps_gene_lvl <- ps_txn_lvl
#ps_gene_lvl <- ps_fn_lvl

#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1", "native mix 2"))
#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1 + invasives", "native mix 2 + invasives"))
ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("Invasive_grasses"))

samp_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(AllPlantSpeciesInSubPlot = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInSubPlot)),
         Bromus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Bromus"), 1, 0),
         Avena = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Avena"), 1, 0),
         Calycadenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Calycadenia"), 1, 0),
         Hemizonia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Hemizonia"), 1, 0),
         Elymus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Elymus"), 1, 0),
         Plantago = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Plantago"), 1, 0),
         Clarkia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Clarkia"), 1, 0),
         Agoseris = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Agoseris"), 1, 0),
         Lasethenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Lasethenia"), 1, 0),
         'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU"), 1, 0),
         'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "VUMI"), 1, 0)) %>%
  select(SampleID, Bromus, Avena, Calycadenia, Hemizonia, Elymus,Plantago,Clarkia,Agoseris,Lasethenia,'Festuca perennis','Festuca microstachys') %>%
  arrange(SampleID)

samp_data <- samp_data[, colSums(samp_data != 0) > 0]

binding_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(PlotSubTreatmentName = ifelse(PlotSubTreatmentName == "Invasive_grasses", "invasive grasses", as.character(PlotSubTreatmentName))) %>%
  select(SampleID, PlotXTreatment, PlotSubTreatmentName, WaterTrt) %>%
  arrange(SampleID)

#make metadata matrix; only has the environmental data you wish to correlate with and has no NA's or blanks
env=as.matrix(samp_data[-c(1)])
rownames(env) <- samp_data$SampleID        

#compute vegan dist & calculate pcoa values (must use bray curtis)

ps_gene_lvl.dist=vegdist(t(ps_gene_lvl@otu_table),method="bray")
ps_gene_lvl_pcoa_res =cmdscale(ps_gene_lvl.dist, eig = TRUE)
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa_res$points

var_exp <- round(ps_gene_lvl_pcoa_res$eig*100/sum(ps_gene_lvl_pcoa_res$eig),1)

res <- adonis2(ps_gene_lvl.dist ~ ., data=as_data_frame(env), by="margin", permutations = 9999) #avena, caly sig
p.adjust(res$`Pr(>F)`[1:length(colnames(env))], 'BH') < 0.05 #one marg sig after pval correction; Calycadenia

#adonis2(ps_gene_lvl.dist ~ PlotSubTreatmentName + PlotXTreatment, data=binding_data, by="margin", permutations = 9999) #sig

#reorder the samples by SampleID b/c of coloring factors
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa[order(row.names(ps_gene_lvl_pcoa)), ]

#fit metadata to pcoa values
ps_gene_lvl_pcoa_env = envfit(ps_gene_lvl_pcoa,env)

#start buliding pcoa plot 

#makes a dataframe that will be used to add the pcoa sample points to plot
scrs <- as.data.frame(scores(ps_gene_lvl_pcoa, display = "sites"))
scrs$SampleID <- rownames(scrs)

#makes factors to color the plot by, these will color your samples, very important that the order matches
scrs <- left_join(scrs, binding_data)

#makes a dataframe that will be used to add metadata factors from envfit to plot as arrows
spp.scrs <- as.data.frame(scores(ps_gene_lvl_pcoa_env, display = "vectors"))

#have to name the arrows, in order they are in metadata file, will error if you don't name them all
spp.scrs$Species <- rownames(spp.scrs)

#add pvals 
spp.scrs$padj <- as_data_frame(p.adjust(ps_gene_lvl_pcoa_env$vectors$pvals, method='BH'))$value
spp.scrs$pvals <- as_data_frame(ps_gene_lvl_pcoa_env$vectors$pvals)$value

#make combined plot of pcoa + arrows
p_tax_grass <- ggplot(scrs) + 
  geom_point(mapping = aes(x = Dim1, y = Dim2, colour = PlotSubTreatmentName, shape=PlotSubTreatmentName), size=4) + 
  theme_bw() +
  geom_segment(data = spp.scrs, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.35, "cm")), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black")) + 
  geom_text_repel(data = spp.scrs, aes(x = Dim1, y = Dim2, label = Species), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black"), size = 4) +
  xlab(paste0("PCoA 1: ",var_exp[1],"% variance")) +
  ylab(paste0("PCoA 2: ",var_exp[2],"% variance")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("#1F968BFF"))


p_tax_grass

##stats

#calculate distance for metadata
env.dist = vegdist(env,method="euclidean")

#run mantel test to see if correlated with microbial distance (bray curtis must be used)
mantel(env.dist,ps_gene_lvl.dist, permutations = 9999) #not sig

## Function - Grasses

#ps_gene_lvl <- readRDS("arrowplots/ps-nocontrols-rare-66880.RDS")
#ps_txn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-tax.RDS")
ps_fn_lvl <-readRDS("arrowplots/ps-nocontrols-rare-66880-fun.RDS")

#ps_gene_lvl <- ps_txn_lvl
ps_gene_lvl <- ps_fn_lvl

#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1", "native mix 2"))
#ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("native mix 1 + invasives", "native mix 2 + invasives"))
ps_gene_lvl <- subset_samples(ps_gene_lvl, PlotSubTreatmentName %in% c("Invasive_grasses"))

samp_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(AllPlantSpeciesInSubPlot = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU/VUMI"), 'LOMU', as.character(AllPlantSpeciesInSubPlot)),
         Bromus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Bromus"), 1, 0),
         Avena = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Avena"), 1, 0),
         Calycadenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Calycadenia"), 1, 0),
         Hemizonia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Hemizonia"), 1, 0),
         Elymus = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Elymus"), 1, 0),
         Plantago = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Plantago"), 1, 0),
         Clarkia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Clarkia"), 1, 0),
         Agoseris = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Agoseris"), 1, 0),
         Lasethenia = ifelse(str_detect(AllPlantSpeciesInSubPlot, "Lasethenia"), 1, 0),
         'Festuca perennis' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "LOMU"), 1, 0),
         'Festuca microstachys' = ifelse(str_detect(AllPlantSpeciesInSubPlot, "VUMI"), 1, 0)) %>%
  select(SampleID, Bromus, Avena, Calycadenia, Hemizonia, Elymus,Plantago,Clarkia,Agoseris,Lasethenia,'Festuca perennis','Festuca microstachys') %>%
  arrange(SampleID)

samp_data <- samp_data[, colSums(samp_data != 0) > 0]

binding_data <- as_tibble(ps_gene_lvl@sam_data) %>%
  mutate(PlotSubTreatmentName = ifelse(PlotSubTreatmentName == "Invasive_grasses", "invasive grasses", as.character(PlotSubTreatmentName))) %>%
  select(SampleID, PlotXTreatment, PlotSubTreatmentName, WaterTrt) %>%
  arrange(SampleID)

#make metadata matrix; only has the environmental data you wish to correlate with and has no NA's or blanks
env=as.matrix(samp_data[-c(1)])
rownames(env) <- samp_data$SampleID        

#compute vegan dist & calculate pcoa values (must use bray curtis)

ps_gene_lvl.dist=vegdist(t(ps_gene_lvl@otu_table),method="bray")
ps_gene_lvl_pcoa_res =cmdscale(ps_gene_lvl.dist, eig = TRUE)
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa_res$points

var_exp <- round(ps_gene_lvl_pcoa_res$eig*100/sum(ps_gene_lvl_pcoa_res$eig),1)

res <- adonis2(ps_gene_lvl.dist ~ ., data=as_data_frame(env), by="margin", permutations = 9999) #avena, caly sig
p.adjust(res$`Pr(>F)`[1:length(colnames(env))], 'BH') < 0.05 #one marg sig after pval correction; Calycadenia

#adonis2(ps_gene_lvl.dist ~ PlotSubTreatmentName + PlotXTreatment, data=binding_data, by="margin", permutations = 9999) #sig

#reorder the samples by SampleID b/c of coloring factors
ps_gene_lvl_pcoa <- ps_gene_lvl_pcoa[order(row.names(ps_gene_lvl_pcoa)), ]

#fit metadata to pcoa values
ps_gene_lvl_pcoa_env = envfit(ps_gene_lvl_pcoa,env)

#start buliding pcoa plot 

#makes a dataframe that will be used to add the pcoa sample points to plot
scrs <- as.data.frame(scores(ps_gene_lvl_pcoa, display = "sites"))
scrs$SampleID <- rownames(scrs)

#makes factors to color the plot by, these will color your samples, very important that the order matches
scrs <- left_join(scrs, binding_data)

#makes a dataframe that will be used to add metadata factors from envfit to plot as arrows
spp.scrs <- as.data.frame(scores(ps_gene_lvl_pcoa_env, display = "vectors"))

#have to name the arrows, in order they are in metadata file, will error if you don't name them all
spp.scrs$Species <- rownames(spp.scrs)

#add pvals 
spp.scrs$padj <- as_data_frame(p.adjust(ps_gene_lvl_pcoa_env$vectors$pvals, method='BH'))$value
spp.scrs$pvals <- as_data_frame(ps_gene_lvl_pcoa_env$vectors$pvals)$value

#make combined plot of pcoa + arrows
p_fun_grass <- ggplot(scrs) + 
  geom_point(mapping = aes(x = Dim1, y = Dim2, colour = PlotSubTreatmentName, shape=PlotSubTreatmentName), size=4) + 
  theme_bw() +
  geom_segment(data = spp.scrs, aes(x = 0, xend = Dim1, y = 0, yend = Dim2), arrow = arrow(length = unit(0.35, "cm")), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black")) + 
  geom_text_repel(data = spp.scrs, aes(x = Dim1, y = Dim2, label = Species), colour = ifelse(spp.scrs$padj > 0.05, "grey", "black"), size = 4) +
  xlab(paste0("PCoA 1: ",var_exp[1],"% variance")) +
  ylab(paste0("PCoA 2: ",var_exp[2],"% variance")) +
  theme(legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_color_manual(values = c("#1F968BFF"))


p_fun_grass


##stats

#calculate distance for metadata
env.dist = vegdist(env,method="euclidean")

#run mantel test to see if correlated with microbial distance (bray curtis must be used)
mantel(env.dist,ps_gene_lvl.dist, permutations = 9999) #not sig


(p_tax_grass + p_fun_grass + plot_layout(guides = 'collect')) / (p_tax_natives + p_fun_natives + plot_layout(guides = 'collect')) / (p_tax_mixes + p_fun_mixes + plot_layout(guides = 'collect')) + 
  plot_annotation(tag_levels = list(c('(a)', '(b)','(c)', '(d)','(e)', '(f)')))

ggsave("envfit_ordinations_w_grass.png", dpi=300, device = "png", height = 12, width= 12, units = "in")
ggsave("envfit_ordinations_w_grass.pdf", dpi=300, device = "pdf", height = 12, width= 12, units = "in")

