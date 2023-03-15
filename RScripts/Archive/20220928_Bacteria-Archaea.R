# separate bacteria and archaea analysis

# tax.gene.arc <- filter(tax.gene, Joined_Domain != "Remove", Domain == "Archaea")
# tax.gene.bac <- filter(tax.gene, Joined_Domain != "Remove", Domain == "Bacteria")
# 
# row.names(tax.gene.arc) <- tax.gene.arc$gene
# tax.gene.arc <- tax.gene.arc[order(tax.gene.arc$gene),]
# tax.gene.arc <- tax.gene.arc[,-1]
# tax.gene.arc <- tax.gene.arc[,c(7:13)]
# 
# row.names(tax.gene.bac) <- tax.gene.bac$gene
# tax.gene.bac <- tax.gene.bac[order(tax.gene.bac$gene),]
# tax.gene.bac <- tax.gene.bac[,-1]
# tax.gene.bac <- tax.gene.bac[,c(7:13)]

# fun.gene.arc <- fun.gene[rownames(fun.gene) %in% rownames(tax.gene.arc),]
# tax.gene.arc <- tax.gene.arc[rownames(tax.gene.arc) %in% rownames(fun.gene.arc),]
# fun.gene.arc <- as.matrix(fun.gene.arc)
# tax.gene.arc <- as.matrix(tax.gene.arc)
# 
# fun.gene.bac <- fun.gene[rownames(fun.gene) %in% rownames(tax.gene.bac),]
# tax.gene.bac <- tax.gene.bac[rownames(tax.gene.bac) %in% rownames(fun.gene.bac),]
# fun.gene.bac <- as.matrix(fun.gene.bac)
# tax.gene.bac <- as.matrix(tax.gene.bac)

decontam_counts_nocontrols <- readRDS("Data/Final-Data/decontam_non-norm_counts_nocontrols.RDS")

# Bacteria
decontam_counts_nocontrols.bac <- decontam_counts_nocontrols[rownames(decontam_counts_nocontrols) %in% rownames(tax.gene.bac),]

decontam_counts_nocontrols.bac <- decontam_counts_nocontrols.bac[rownames(decontam_counts_nocontrols.bac) %in% rownames(fun.gene.bac),]

# Archaea
decontam_counts_nocontrols.arc <- decontam_counts_nocontrols[rownames(decontam_counts_nocontrols) %in% rownames(tax.gene.arc),]

decontam_counts_nocontrols.arc <- decontam_counts_nocontrols.arc[rownames(decontam_counts_nocontrols.arc) %in% rownames(fun.gene.arc),]

# Bacteria ####
## DESeq2 genes ####
ps <- phyloseq(otu_table(decontam_counts_nocontrols.bac, taxa_are_rows = T), sample_data(metadata.nocontrols))

treat <- phyloseq_to_deseq2(ps, ~ plant.water)

dds.pw <- DESeq(treat)

contrasts <- c("FW.FC", 
               "FD.FC", 
               "MD.FC", 
               "MW.FC", 
               "GD.GC",
               "GW.GC",
               "MD.GC",
               "MW.GC",
               "FC.GC",
               "FD.GD",
               "FW.GW",
               "MD.FD",
               "MW.FW",
               "MC.FC",
               "MC.GC")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forbs.Control/forb.Watering) so if something is higher in forb controls, the number is greater than one, so the log fold change > 0, if something is higher in the watering treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FW.FC = c("plant.water", "forb.Watering", "forb.Control"),
                      FD.FC = c("plant.water","forb.Drought","forb.Control"),
                      MD.FC = c("plant.water","mix.Drought","forb.Control"),
                      MW.FC = c("plant.water","mix.Watering","forb.Control"),
                      GD.GC = c("plant.water","grasses.Drought","grasses.Control"),
                      GW.GC = c("plant.water","grasses.Watering","grasses.Control"),
                      MD.GC = c("plant.water","mix.Drought","grasses.Control"),
                      MW.GC = c("plant.water","mix.Watering","grasses.Control"),
                      FC.GC = c("plant.water","forb.Control","grasses.Control"),
                      FD.GD = c("plant.water","forb.Drought","grasses.Drought"),
                      FW.GW = c("plant.water","forb.Watering","grasses.Watering"),
                      MD.FD = c("plant.water","mix.Drought","forb.Drought"),
                      MW.FW = c("plant.water", "mix.Watering", "forb.Watering"),
                      MC.FC = c("plant.water", "mix.Control", "forb.Control"), 
                      MC.GC = c("plant.water", "mix.Control", "grasses.Watering"))

plot.name.list <- list(FW.FC = "Forb Control v Forb Watering",
                       FD.FC = "Forb Control v Forb Drought",
                       MD.FC = "Forb Control v Mix Drought",
                       MW.FC = "Forb Control v Mix Watering",
                       GD.GC = "Grass Control v Grass Drought",
                       GW.GC = "Grass Control v Grass Watering",
                       MD.GC = "Grass Control v Mix Drought",
                       MW.GC  = "Grass Control v Mix Watering",
                       FC.GC = "Grass Control v Forb Control",
                       FD.GD = "Grass Drought v Forb Drought",
                       FW.GW = "Grass Watering v Forb Watering",
                       MD.FD = "Forb Drought v Mix Drought",
                       MW.FW = "Forb Watering v Mix Watering",
                       MC.FC = "Mix Control v Forb Control",
                       MC.GC = "Mix Control v Grass Control")

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
length(unique(sig.genes))  #0

## DESeq2 functions ####

ps.fun <- phyloseq(otu_table(decontam_counts_nocontrols.bac, taxa_are_rows = T), 
                   sample_data(metadata.nocontrols), 
                   tax_table(fun.gene.bac))

ps.fun <- tax_glom(ps.fun, taxrank = "COG20_FUNCTION")

treat.fun <- phyloseq_to_deseq2(ps.fun, ~ plant.water)

dds.pw.fun <- DESeq(treat.fun)

#resultsNames(dds.pw.fun) #gives us comparisons for our contrasts

contrasts <- c("FW.FC", 
               "FD.FC", 
               "MD.FC", 
               "MW.FC", 
               "GD.GC",
               "GW.GC",
               "MD.GC",
               "MW.GC",
               "FC.GC",
               "FD.GD",
               "FW.GW",
               "MD.FD",
               "MW.FW",
               "MC.FC",
               "MC.GC")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forbs.Control/forb.Watering) so if something is higher in forb controls, the number is greater than one, so the log fold change > 0, if something is higher in the watering treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FW.FC = c("plant.water", "forb.Watering", "forb.Control"),
                      FD.FC = c("plant.water","forb.Drought","forb.Control"),
                      MD.FC = c("plant.water","mix.Drought","forb.Control"),
                      MW.FC = c("plant.water","mix.Watering","forb.Control"),
                      GD.GC = c("plant.water","grasses.Drought","grasses.Control"),
                      GW.GC = c("plant.water","grasses.Watering","grasses.Control"),
                      MD.GC = c("plant.water","mix.Drought","grasses.Control"),
                      MW.GC = c("plant.water","mix.Watering","grasses.Control"),
                      FC.GC = c("plant.water","forb.Control","grasses.Control"),
                      FD.GD = c("plant.water","forb.Drought","grasses.Drought"),
                      FW.GW = c("plant.water","forb.Watering","grasses.Watering"),
                      MD.FD = c("plant.water","mix.Drought","forb.Drought"),
                      MW.FW = c("plant.water", "mix.Watering", "forb.Watering"),
                      MC.FC = c("plant.water", "mix.Control", "forb.Control"), 
                      MC.GC = c("plant.water", "mix.Control", "grasses.Watering"))

plot.name.list <- list(FW.FC = "Forb Control v Forb Watering",
                       FD.FC = "Forb Control v Forb Drought",
                       MD.FC = "Forb Control v Mix Drought",
                       MW.FC = "Forb Control v Mix Watering",
                       GD.GC = "Grass Control v Grass Drought",
                       GW.GC = "Grass Control v Grass Watering",
                       MD.GC = "Grass Control v Mix Drought",
                       MW.GC  = "Grass Control v Mix Watering",
                       FC.GC = "Grass Control v Forb Control",
                       FD.GD = "Grass Drought v Forb Drought",
                       FW.GW = "Grass Watering v Forb Watering",
                       MD.FD = "Forb Drought v Mix Drought",
                       MW.FW = "Forb Watering v Mix Watering",
                       MC.FC = "Mix Control v Forb Control",
                       MC.GC = "Mix Control v Grass Control")

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw.fun, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
length(unique(sig.genes)) #45

names(plot.list)

test <- rbind(data.frame(cbind(pair = "MD.FC", res.list[["MD.FC"]])), data.frame(cbind(pair = "MD.FD", res.list[["MD.FD"]])), data.frame(cbind(pair = "FW.GW", res.list[["FW.GW"]])), data.frame(cbind(pair = "MC.FC", res.list[["MC.FC"]])))


write.csv(test, "Data/fun-gene-DEG.csv", row.names = F)

#plot results
grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]], ncol = 2)

test2 <- as.data.frame(cbind(fun.gene.bac, gene = row.names(fun.gene.bac)))

for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], test2, by = "gene", all.y = F)
}

sig.dif.fun <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.fun[[i]] <- ddply(res.list[[i]], .(COG20_CATEGORY), summarize, estimate = sum(estimate))
  }
}

plot.list2 <- list()

for(i in names(sig.dif.fun)){
    
    sig.dif.fun[[i]]$COG20_CATEGORY <- factor(sig.dif.fun[[i]]$COG20_CATEGORY, levels = sig.dif.fun[[i]]$COG20_CATEGORY[order(sig.dif.fun[[i]]$estimate)])
    
    p <- ggplot(sig.dif.fun[[i]], aes(x = COG20_CATEGORY, y = estimate)) + 
      geom_point(size = 2.5) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 13),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15)
          ) +
      ylim(-2,1.7) +
      labs(title = i, y = expression(paste(log[2], " fold change"))) +
      scale_x_discrete(labels = c("Posttranslational modification, protein turnover, chaperones" = "Posttranslational modification, protein\nturnover, chaperones", "Translation, ribosomal structure and biogenesis" = "Translation, ribosomal structure\nand biogenesis", "Intracellular trafficking, secretion, and vesicular transport" = "Intracellular trafficking, secretion, and\nvesicular transport")) +
      coord_flip()
  
   plot.list2[[i]] = p
}

plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[3]], plot.list2[[2]], plot.list2[[4]], align = "v")

fun.fig <- grid.arrange(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], plot.list2[[4]], ncol = 2, labels = c("A", "B", "C"))

ggsave("Figures/deseq-fun.jpeg", fun.fig, height = 10, width = 18, units = "in", dpi = 600)

## DESeq2 taxa ####
ps <- phyloseq(otu_table(decontam_counts_nocontrols.bac, taxa_are_rows = T), sample_data(metadata.nocontrols), tax_table(tax.gene.bac))

ps.genus <- tax_glom(ps, taxrank = "Genus")
treat.genus <- phyloseq_to_deseq2(ps.genus, ~ plant.water)

dds.pw.genus <- DESeq(treat.genus)

contrasts <- c("FW.FC", 
               "FD.FC", 
               "MD.FC", 
               "MW.FC", 
               "GD.GC",
               "GW.GC",
               "MD.GC",
               "MW.GC",
               "FC.GC",
               "FD.GD",
               "FW.GW",
               "MD.FD",
               "MW.FW",
               "MC.FC",
               "MC.GC")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forbs.Control/forb.Watering) so if something is higher in forb controls, the number is greater than one, so the log fold change > 0, if something is higher in the watering treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FW.FC = c("plant.water", "forb.Watering", "forb.Control"),
                      FD.FC = c("plant.water","forb.Drought","forb.Control"),
                      MD.FC = c("plant.water","mix.Drought","forb.Control"),
                      MW.FC = c("plant.water","mix.Watering","forb.Control"),
                      GD.GC = c("plant.water","grasses.Drought","grasses.Control"),
                      GW.GC = c("plant.water","grasses.Watering","grasses.Control"),
                      MD.GC = c("plant.water","mix.Drought","grasses.Control"),
                      MW.GC = c("plant.water","mix.Watering","grasses.Control"),
                      FC.GC = c("plant.water","forb.Control","grasses.Control"),
                      FD.GD = c("plant.water","forb.Drought","grasses.Drought"),
                      FW.GW = c("plant.water","forb.Watering","grasses.Watering"),
                      MD.FD = c("plant.water","mix.Drought","forb.Drought"),
                      MW.FW = c("plant.water", "mix.Watering", "forb.Watering"),
                      MC.FC = c("plant.water", "mix.Control", "forb.Control"), 
                      MC.GC = c("plant.water", "mix.Control", "grasses.Watering"))

plot.name.list <- list(FW.FC = "Forb Control v Forb Watering",
                       FD.FC = "Forb Control v Forb Drought",
                       MD.FC = "Forb Control v Mix Drought",
                       MW.FC = "Forb Control v Mix Watering",
                       GD.GC = "Grass Control v Grass Drought",
                       GW.GC = "Grass Control v Grass Watering",
                       MD.GC = "Grass Control v Mix Drought",
                       MW.GC  = "Grass Control v Mix Watering",
                       FC.GC = "Grass Control v Forb Control",
                       FD.GD = "Grass Drought v Forb Drought",
                       FW.GW = "Grass Watering v Forb Watering",
                       MD.FD = "Forb Drought v Mix Drought",
                       MW.FW = "Forb Watering v Mix Watering",
                       MC.FC = "Mix Control v Forb Control",
                       MC.GC = "Mix Control v Grass Control")

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw.genus, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
length(unique(sig.genes)) # 49

#plot results
grid.arrange(plot.list[[1]], plot.list[[2]], plot.list[[3]], plot.list[[4]])


test2 <- as.data.frame(cbind(tax.gene.bac, gene = row.names(tax.gene.bac)))

for(i in names(res.list)) {
 res.list[[i]] <- merge(res.list[[i]], test2, by = "gene", all.y = F)
}

sig.dif.tax <- list()

for(i in names(res.list)) {
  if(nrow(res.list[[i]]) > 0) {
  sig.dif.tax[[i]] <- ddply(res.list[[i]], .(Family), summarize, estimate = sum(estimate))
  }
}

plot.list2 <- list()

for(i in names(sig.dif.tax)){
    
    sig.dif.tax[[i]]$Family <- factor(sig.dif.tax[[i]]$Family, levels = sig.dif.tax[[i]]$Family[order(sig.dif.tax[[i]]$estimate)])
    
    p <- ggplot(sig.dif.tax[[i]], aes(x = Family, y = estimate)) + 
      geom_point(size = 2.5) +
      geom_hline(yintercept = 0, col= "red") +
      theme_bw() +
      theme(
        axis.text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 15)
          ) +
      ylim(-4, 6) +
      labs(title = i, y = expression(paste(log[2], " fold change"))) +
      coord_flip()
  
   plot.list2[[i]] = p
}

plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], plot.list2[[4]], align = "v")

tax.fig <- grid.arrange(plot.list2[[1]], plot.list2[[2]], plot.list2[[3]], plot.list2[[4]], ncol = 2)

ggsave("Figures/deseq-tax.jpeg", tax.fig, height = 12, width = 9, units = "in", dpi = 600)







df <- data.frame()

for(i in names(plot.list)){ #plot list is the one with all the significantly different ones
  tmp <- cbind(pair = i, res.list[[i]])
  df <- rbind(df, tmp)
}

test2 <- as.data.frame(cbind(tax.gene, gene = row.names(tax.gene)))

df <- merge(df, test2, by = "gene", all.y = F)

write.csv(df, "Data/tax-gene-DEG.csv")

df <- ddply(df, .(Family, pair), summarize, logfoldchange = sum(estimate))

#df$Family <- factor(df$Family, levels = df$Family[order(df$logfoldchange)])

ggplot(df[df$pair == "FC.MD",], aes(x = Family, y = logfoldchange)) + 
    geom_point(size = 3) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 12),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none",
      axis.title.y = element_blank(),
        ) +
    labs(title = "Forb Control v Mix Drought", y = expression(paste(log[2], " fold change"))) +
  coord_flip()

## Ordination - genes ####

###Stats ####
dds <- DESeqDataSetFromMatrix(decontam_counts_nocontrols,
                                   colData = metadata.nocontrols,
                                   design = ~ plant.water)

dds <- DESeq(dds)

plotDispEsts(dds)

norm_counts <- counts(dds, normalized = T)
ps <- phyloseq(otu_table(norm_counts, taxa_are_rows = T), sample_data(metadata.nocontrols))


ps.hell <- transform(ps, "hellinger")

PCoA_hell <- ordinate(physeq = ps.hell, method = "PCoA", distance = "euclidean")

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

pairwise <- permanova_pairwise(
  Dist.hell,
  as(sample_data(ps.hell), "data.frame")$PlantTrt,
  permutations = 9999,
  padj = "BH"
)

### Figure ####
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




## Ordination - function ####
### Stats ####
ps.hell.fun <- transform(ps.fun, "hellinger")
#ps.hell.fun <- subset_samples(ps.hell.fun, WaterTrt != "Watering")

Dist.hell.fun <- phyloseq::distance(ps.hell.fun, method = "euclidean", type = "samples")

set.seed(620)
adonis2(Dist.hell.fun ~ WaterTrt * PlantTrt, as(sample_data(ps.hell.fun), "data.frame"), by = "margin", permutations = 9999)  

set.seed(620)
adonis2(Dist.hell.fun ~ WaterTrt + PlantTrt, as(sample_data(ps.hell.fun), "data.frame"), by = "margin", permutations = 9999)

set.seed(620)
pairwise <- permanova_pairwise(
  Dist.hell.fun,
  as(sample_data(ps.hell.fun), "data.frame")$WaterTrt,
  permutations = 9999,
  padj = "BH"
)

set.seed(620)
pairwise <- permanova_pairwise(
  Dist.hell.fun,
  as(sample_data(ps.hell.fun), "data.frame")$PlantTrt,
  permutations = 9999,
  padj = "BH"
)

### Figure ####
PCoA_hell.fun <- ordinate(physeq = ps.hell.fun, method = "PCoA", distance = "euclidean")

plot_ordination(ps.hell.fun, PCoA_hell.fun, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
     panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_point(size = 2) +
  stat_ellipse() +
  facet_wrap(~WaterTrt) + 
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

# differences between watering treatments
a <- plot_ordination(ps.hell.fun, PCoA_hell.fun, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.9,0.9),
     legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  geom_point(size = 2) +
  lims(x = c(-0.18, 0.2)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

# differences between plant treatments
b <- plot_ordination(ps.hell.fun, PCoA_hell.fun, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.9,0.85),
    legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
    legend.margin = margin(-3,-3,-3,-3)) +
  lims(x = c(-0.18, 0.2)) +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

# get centroids plots
C <- betadisper(Dist.hell.fun, as(sample_data(ps.hell.fun), "data.frame")$plant.water)

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
     legend.position = c(0.75,0.85),
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

d <- grid.arrange(a,b,c, ncol = 3)

ggsave("Figures/Ordination-fun.jpeg", d, height = 4, width = 15, units = "in", dpi = 600)

## Ordination - Taxa ####
ps.hell.taxa <- transform(ps.genus, "hellinger")

### Stats ####
Dist.hell.taxa <- phyloseq::distance(ps.hell.taxa, method = "euclidean", type = "samples")


set.seed(620)
adonis2(Dist.hell.taxa ~ WaterTrt * PlantTrt, as(sample_data(ps.hell.taxa), "data.frame"), by = "margin", permutations = 9999)  # interaction not significant

set.seed(620)
adonis2(Dist.hell.taxa ~ WaterTrt + PlantTrt, as(sample_data(ps.hell.taxa), "data.frame"), by = "margin", permutations = 9999) # watering sig

set.seed(620)
pairwise <- permanova_pairwise(
  Dist.hell.taxa,
  as(sample_data(ps.hell.taxa), "data.frame")$WaterTrt,
  permutations = 9999,
  padj = "BH"
)

set.seed(620)
pairwise <- permanova_pairwise(
  Dist.hell.taxa,
  as(sample_data(ps.hell.taxa), "data.frame")$PlantTrt,
  permutations = 9999,
  padj = "BH"
)

### Figure ####
PCoA_hell.taxa <- ordinate(physeq = ps.hell.taxa, method = "PCoA", distance = "euclidean")

# differences between watering treatments
a <- plot_ordination(ps.hell.taxa, PCoA_hell.taxa, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.9,0.9),
     legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  geom_point(size = 2) +
  lims(x = c(-0.3, 0.25)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

# differences between plant treatments
b <- plot_ordination(ps.hell.taxa, PCoA_hell.taxa, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.85,0.92),
    legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
    legend.margin = margin(-3,-3,-3,-3)) +
  lims(x = c(-0.28, 0.3)) +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

# get centroids plots
C <- betadisper(Dist.hell.taxa, as(sample_data(ps.hell.taxa), "data.frame")$plant.water)

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
     legend.position = c(0.8,0.9),
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

d <- grid.arrange(a,b,c, ncol = 3)

ggsave("Figures/Ordination-taxa.jpeg", d, height = 4, width = 15, units = "in", dpi = 600)

## Ordination - Forbs ####
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

# Archaea ####
## DESeq2 genes ####
ps <- phyloseq(otu_table(decontam_counts_nocontrols.arc, taxa_are_rows = T), sample_data(metadata.nocontrols))

treat <- phyloseq_to_deseq2(ps, ~ plant.water)

dds.pw <- DESeq(treat)

contrasts <- c("FW.FC", 
               "FD.FC", 
               "MD.FC", 
               "MW.FC", 
               "GD.GC",
               "GW.GC",
               "MD.GC",
               "MW.GC",
               "FC.GC",
               "FD.GD",
               "FW.GW",
               "MD.FD",
               "MW.FW",
               "MC.FC",
               "MC.GC")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forbs.Control/forb.Watering) so if something is higher in forb controls, the number is greater than one, so the log fold change > 0, if something is higher in the watering treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FW.FC = c("plant.water", "forb.Watering", "forb.Control"),
                      FD.FC = c("plant.water","forb.Drought","forb.Control"),
                      MD.FC = c("plant.water","mix.Drought","forb.Control"),
                      MW.FC = c("plant.water","mix.Watering","forb.Control"),
                      GD.GC = c("plant.water","grasses.Drought","grasses.Control"),
                      GW.GC = c("plant.water","grasses.Watering","grasses.Control"),
                      MD.GC = c("plant.water","mix.Drought","grasses.Control"),
                      MW.GC = c("plant.water","mix.Watering","grasses.Control"),
                      FC.GC = c("plant.water","forb.Control","grasses.Control"),
                      FD.GD = c("plant.water","forb.Drought","grasses.Drought"),
                      FW.GW = c("plant.water","forb.Watering","grasses.Watering"),
                      MD.FD = c("plant.water","mix.Drought","forb.Drought"),
                      MW.FW = c("plant.water", "mix.Watering", "forb.Watering"),
                      MC.FC = c("plant.water", "mix.Control", "forb.Control"), 
                      MC.GC = c("plant.water", "mix.Control", "grasses.Watering"))

plot.name.list <- list(FW.FC = "Forb Control v Forb Watering",
                       FD.FC = "Forb Control v Forb Drought",
                       MD.FC = "Forb Control v Mix Drought",
                       MW.FC = "Forb Control v Mix Watering",
                       GD.GC = "Grass Control v Grass Drought",
                       GW.GC = "Grass Control v Grass Watering",
                       MD.GC = "Grass Control v Mix Drought",
                       MW.GC  = "Grass Control v Mix Watering",
                       FC.GC = "Grass Control v Forb Control",
                       FD.GD = "Grass Drought v Forb Drought",
                       FW.GW = "Grass Watering v Forb Watering",
                       MD.FD = "Forb Drought v Mix Drought",
                       MW.FW = "Forb Watering v Mix Watering",
                       MC.FC = "Mix Control v Forb Control",
                       MC.GC = "Mix Control v Grass Control")

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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

## DESeq2 functions ####

ps.fun <- phyloseq(otu_table(decontam_counts_nocontrols.arc, taxa_are_rows = T), 
                   sample_data(metadata.nocontrols), 
                   tax_table(fun.gene.arc))

ps.fun <- tax_glom(ps.fun, taxrank = "COG20_FUNCTION")

treat.fun <- phyloseq_to_deseq2(ps.fun, ~ plant.water)

dds.pw.fun <- DESeq(treat.fun)

#resultsNames(dds.pw.fun) #gives us comparisons for our contrasts

contrasts <- c("FW.FC", 
               "FD.FC", 
               "MD.FC", 
               "MW.FC", 
               "GD.GC",
               "GW.GC",
               "MD.GC",
               "MW.GC",
               "FC.GC",
               "FD.GD",
               "FW.GW",
               "MD.FD",
               "MW.FW",
               "MC.FC",
               "MC.GC")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forbs.Control/forb.Watering) so if something is higher in forb controls, the number is greater than one, so the log fold change > 0, if something is higher in the watering treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FW.FC = c("plant.water", "forb.Watering", "forb.Control"),
                      FD.FC = c("plant.water","forb.Drought","forb.Control"),
                      MD.FC = c("plant.water","mix.Drought","forb.Control"),
                      MW.FC = c("plant.water","mix.Watering","forb.Control"),
                      GD.GC = c("plant.water","grasses.Drought","grasses.Control"),
                      GW.GC = c("plant.water","grasses.Watering","grasses.Control"),
                      MD.GC = c("plant.water","mix.Drought","grasses.Control"),
                      MW.GC = c("plant.water","mix.Watering","grasses.Control"),
                      FC.GC = c("plant.water","forb.Control","grasses.Control"),
                      FD.GD = c("plant.water","forb.Drought","grasses.Drought"),
                      FW.GW = c("plant.water","forb.Watering","grasses.Watering"),
                      MD.FD = c("plant.water","mix.Drought","forb.Drought"),
                      MW.FW = c("plant.water", "mix.Watering", "forb.Watering"),
                      MC.FC = c("plant.water", "mix.Control", "forb.Control"), 
                      MC.GC = c("plant.water", "mix.Control", "grasses.Watering"))

plot.name.list <- list(FW.FC = "Forb Control v Forb Watering",
                       FD.FC = "Forb Control v Forb Drought",
                       MD.FC = "Forb Control v Mix Drought",
                       MW.FC = "Forb Control v Mix Watering",
                       GD.GC = "Grass Control v Grass Drought",
                       GW.GC = "Grass Control v Grass Watering",
                       MD.GC = "Grass Control v Mix Drought",
                       MW.GC  = "Grass Control v Mix Watering",
                       FC.GC = "Grass Control v Forb Control",
                       FD.GD = "Grass Drought v Forb Drought",
                       FW.GW = "Grass Watering v Forb Watering",
                       MD.FD = "Forb Drought v Mix Drought",
                       MW.FW = "Forb Watering v Mix Watering",
                       MC.FC = "Mix Control v Forb Control",
                       MC.GC = "Mix Control v Grass Control")

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw.fun, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
length(unique(sig.genes)) # one gene function differed between forbs in watering and grasses in watering - 6047836


## DESeq2 taxa ####
ps <- phyloseq(otu_table(decontam_counts_nocontrols.arc, taxa_are_rows = T), sample_data(metadata.nocontrols), tax_table(tax.gene.arc))

ps.genus <- tax_glom(ps, taxrank = "Genus")
treat.genus <- phyloseq_to_deseq2(ps.genus, ~ plant.water)

dds.pw.genus <- DESeq(treat.genus)

contrasts <- c("FW.FC", 
               "FD.FC", 
               "MD.FC", 
               "MW.FC", 
               "GD.GC",
               "GW.GC",
               "MD.GC",
               "MW.GC",
               "FC.GC",
               "FD.GD",
               "FW.GW",
               "MD.FD",
               "MW.FW",
               "MC.FC",
               "MC.GC")

# contrasts: name of factor, numerator for fold change, denominator for fold change; i.e. log2(forbs.Control/forb.Watering) so if something is higher in forb controls, the number is greater than one, so the log fold change > 0, if something is higher in the watering treatment then the number is less than one so the log2 fold change is negative


contrast.list <- list(FW.FC = c("plant.water", "forb.Watering", "forb.Control"),
                      FD.FC = c("plant.water","forb.Drought","forb.Control"),
                      MD.FC = c("plant.water","mix.Drought","forb.Control"),
                      MW.FC = c("plant.water","mix.Watering","forb.Control"),
                      GD.GC = c("plant.water","grasses.Drought","grasses.Control"),
                      GW.GC = c("plant.water","grasses.Watering","grasses.Control"),
                      MD.GC = c("plant.water","mix.Drought","grasses.Control"),
                      MW.GC = c("plant.water","mix.Watering","grasses.Control"),
                      FC.GC = c("plant.water","forb.Control","grasses.Control"),
                      FD.GD = c("plant.water","forb.Drought","grasses.Drought"),
                      FW.GW = c("plant.water","forb.Watering","grasses.Watering"),
                      MD.FD = c("plant.water","mix.Drought","forb.Drought"),
                      MW.FW = c("plant.water", "mix.Watering", "forb.Watering"),
                      MC.FC = c("plant.water", "mix.Control", "forb.Control"), 
                      MC.GC = c("plant.water", "mix.Control", "grasses.Watering"))

plot.name.list <- list(FW.FC = "Forb Control v Forb Watering",
                       FD.FC = "Forb Control v Forb Drought",
                       MD.FC = "Forb Control v Mix Drought",
                       MW.FC = "Forb Control v Mix Watering",
                       GD.GC = "Grass Control v Grass Drought",
                       GW.GC = "Grass Control v Grass Watering",
                       MD.GC = "Grass Control v Mix Drought",
                       MW.GC  = "Grass Control v Mix Watering",
                       FC.GC = "Grass Control v Forb Control",
                       FD.GD = "Grass Drought v Forb Drought",
                       FW.GW = "Grass Watering v Forb Watering",
                       MD.FD = "Forb Drought v Mix Drought",
                       MW.FW = "Forb Watering v Mix Watering",
                       MC.FC = "Mix Control v Forb Control",
                       MC.GC = "Mix Control v Grass Control")

alpha = 0.05
res.list <- list()
plot.list <- list()
sig.genes <- c()

for(i in contrasts) {
  #get results for each contrast
  res <- results(dds.pw.genus, contrast = contrast.list[[i]], pAdjustMethod = "BH")
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
length(unique(sig.genes)) # one taxa dif abundant 2594854 (MC.FC)


## Ordination - genes ####
###Stats ####
dds <- DESeqDataSetFromMatrix(decontam_counts_nocontrols.arc,
                                   colData = metadata.nocontrols,
                                   design = ~ plant.water)

dds <- DESeq(dds)

plotDispEsts(dds)

norm_counts <- counts(dds, normalized = T)
ps <- phyloseq(otu_table(norm_counts, taxa_are_rows = T), sample_data(metadata.nocontrols))

ps.hell <- transform(ps, "hellinger")

PCoA_hell <- ordinate(physeq = ps.hell, method = "PCoA", distance = "euclidean")

# Get distances for stats
Dist.hell <- phyloseq::distance(ps.hell, method = "euclidean", type = "samples")

set.seed(620)
adonis2(Dist.hell ~ WaterTrt * PlantTrt, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999)  # interaction not significant

set.seed(620)
adonis2(Dist.hell ~ PlantTrt + WaterTrt, as(sample_data(ps.hell), "data.frame"), by = "margin", permutations = 9999) # watering treatment sig differs

C <- betadisper(Dist.hell, as(sample_data(ps.hell), "data.frame")$WaterTrt)

set.seed(620)
permutest(C, permutations = 9999) # no sig dif

pairwise <- permanova_pairwise(
  Dist.hell,
  as(sample_data(ps.hell), "data.frame")$WaterTrt,
  permutations = 9999,
  padj = "BH"
)

### Figure ####
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




## Ordination - function ####
### Stats ####
ps.hell.fun <- transform(ps.fun, "hellinger")

Dist.hell.fun <- phyloseq::distance(ps.hell.fun, method = "euclidean", type = "samples")

set.seed(620)
adonis2(Dist.hell.fun ~ WaterTrt * PlantTrt, as(sample_data(ps.hell.fun), "data.frame"), by = "margin", permutations = 9999)  

set.seed(620)
adonis2(Dist.hell.fun ~ WaterTrt + PlantTrt, as(sample_data(ps.hell.fun), "data.frame"), by = "margin", permutations = 9999)

set.seed(620)
pairwise <- permanova_pairwise(
  Dist.hell.fun,
  as(sample_data(ps.hell.fun), "data.frame")$WaterTrt,
  permutations = 9999,
  padj = "BH"
)

### Figure ####
PCoA_hell.fun <- ordinate(physeq = ps.hell.fun, method = "PCoA", distance = "euclidean")

plot_ordination(ps.hell.fun, PCoA_hell.fun, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
     panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
  geom_point(size = 2) +
  stat_ellipse() +
  facet_wrap(~WaterTrt) + 
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3"))

# differences between watering treatments
a <- plot_ordination(ps.hell.fun, PCoA_hell.fun, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.9,0.9),
     legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  geom_point(size = 2) +
  lims(x = c(-0.18, 0.2)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

# differences between plant treatments
b <- plot_ordination(ps.hell.fun, PCoA_hell.fun, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.9,0.85),
    legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
    legend.margin = margin(-3,-3,-3,-3)) +
  lims(x = c(-0.18, 0.2)) +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

# get centroids plots
C <- betadisper(Dist.hell.fun, as(sample_data(ps.hell.fun), "data.frame")$plant.water)

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
     legend.position = c(0.75,0.85),
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  labs(x = "PCoA1 (21.2%)", y = "PCoA2 (11.4%)")

d <- grid.arrange(a,b,c, ncol = 3)

ggsave("Figures/Ordination-fun.jpeg", d, height = 4, width = 15, units = "in", dpi = 600)

## Ordination - Taxa ####
ps.hell.taxa <- transform(ps.genus, "hellinger")

### Stats ####
Dist.hell.taxa <- phyloseq::distance(ps.hell.taxa, method = "euclidean", type = "samples")


set.seed(620)
adonis2(Dist.hell.taxa ~ WaterTrt * PlantTrt, as(sample_data(ps.hell.taxa), "data.frame"), by = "margin", permutations = 9999)  # interaction not significant

set.seed(620)
adonis2(Dist.hell.taxa ~ WaterTrt + PlantTrt, as(sample_data(ps.hell.taxa), "data.frame"), by = "margin", permutations = 9999) # watering sig

set.seed(620)
pairwise <- permanova_pairwise(
  Dist.hell.taxa,
  as(sample_data(ps.hell.taxa), "data.frame")$WaterTrt,
  permutations = 9999,
  padj = "BH"
)

set.seed(620)
pairwise <- permanova_pairwise(
  Dist.hell.taxa,
  as(sample_data(ps.hell.taxa), "data.frame")$PlantTrt,
  permutations = 9999,
  padj = "BH"
)

### Figure ####
PCoA_hell.taxa <- ordinate(physeq = ps.hell.taxa, method = "PCoA", distance = "euclidean")

# differences between watering treatments
a <- plot_ordination(ps.hell.taxa, PCoA_hell.taxa, color = "WaterTrt") +
  theme_classic(base_size = 10) +
  theme(
     #panel.border = element_rect(colour = "black", fill = NA, size = 1),
     legend.position = c(0.9,0.9),
     legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  geom_point(size = 2) +
  lims(x = c(-0.3, 0.25)) +
  stat_ellipse() +
  scale_color_manual(values = c("red", "black", "blue")) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

# differences between plant treatments
b <- plot_ordination(ps.hell.taxa, PCoA_hell.taxa, color = "PlantTrt") +
  theme_classic(base_size = 10) +
  theme(
    #panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.position = c(0.85,0.92),
    legend.title = element_blank(),
    legend.box.margin = margin(-1,-1,-1,-1),
    legend.margin = margin(-3,-3,-3,-3)) +
  lims(x = c(-0.28, 0.3)) +
  geom_point(size = 2) +
  stat_ellipse() +
  scale_color_manual(values = c("magenta4", "#1F968BFF", "goldenrod3")) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

# get centroids plots
C <- betadisper(Dist.hell.taxa, as(sample_data(ps.hell.taxa), "data.frame")$plant.water)

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
     legend.position = c(0.8,0.9),
     legend.title = element_blank(),
     legend.box.margin = margin(-1,-1,-1,-1),
     legend.margin = margin(-3,-3,-3,-3)) +
  labs(x = "PCoA1 (28.7%)", y = "PCoA2 (22.6%)")

d <- grid.arrange(a,b,c, ncol = 3)

ggsave("Figures/Ordination-taxa.jpeg", d, height = 4, width = 15, units = "in", dpi = 600)

## Ordination - Forbs ####
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

##  Alpha diversity ####
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



#plot results
grid.arrange(plot.list[[1]])

df <- data.frame()

for(i in names(plot.list)){ #plot list is the one with all the significantly different ones
  tmp <- cbind(pair = i, res.list[[i]])
  df <- rbind(df, tmp)
}

## OLD Functional analysis ####
# fun.gene.sig <- fun.gene[fun.gene$gene %in% sig.genes,] # 16 genes have functions!
# # fun.gene.sig <- filter(fun.gene.sig, !is.null(fun.gene.sig$COG20_CATEGORY))
# # no.fun <- unique(sig.genes[!sig.genes %in% fun.gene$gene]) #4 have no discernible function
# 
# #write.csv(fun.gene.sig, "Data/Final-Data/20221205_fun-gene-cog-allsigs.csv", row.names = F)
# 
# #Edited to reflect more streamlined COG categories
# #fun.gene.cog <- read.csv("Data/Final-Data/20221205_fun-gene-cog-allsigs_ML-CE-Edit.csv")
# 
# length(unique(fun.gene.cog$gene))
# 
# unique(fun.gene.cog$Fun_Cat_CE_ML) # 13 unique functions yay!
# 
# for(i in names(res.list)) {
#  res.list[[i]] <- merge(res.list[[i]], fun.gene[,c(1,11)], by = "gene", all.y = F)
# }
# 
# sig.dif.fun <- list()
# 
# for(i in names(res.list)) {
#   if(nrow(res.list[[i]]) > 0) {
#   sig.dif.fun[[i]] <- ddply(res.list[[i]], .(COG20_CATEGORY), summarize, estimate = sum(estimate))
#   }
# }
# 
## OLD Summed by Function ####
# 
# names(sig.dif.fun)[[1]] <- "forbs (control) vs. forb (watered)"
# names(sig.dif.fun)[[2]] <- "forbs (control) vs. mixes (drought)"
# names(sig.dif.fun)[[3]] <- "forbs (drought) vs. mixes (drought)"
# 
# plot.list2 <- list()
# 
# for(i in names(sig.dif.fun)){
#     
#     sig.dif.fun[[i]]$COG20_CATEGORY <- factor(sig.dif.fun[[i]]$COG20_CATEGORY, levels = sig.dif.fun[[i]]$COG20_CATEGORY[order(sig.dif.fun[[i]]$estimate)])
#     
#     p <- ggplot(sig.dif.fun[[i]], aes(x = COG20_CATEGORY, y = estimate)) + 
#       geom_point(size = 2.5) +
#       geom_hline(yintercept = 0, col= "red") +
#       theme_bw() +
#       theme(
#         axis.text = element_text(size = 12),
#         plot.title = element_text(hjust = 0.5, size = 13),
#         legend.position = "none",
#         axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 12)
#           ) +
#       #ylim(-26,10) +
#       labs(title = i, y = expression(paste(log[2], " fold change"))) +
#       scale_x_discrete(labels = c("Posttranslational modification, protein turnover, chaperones" = "Posttranslational modification, protein\nturnover, chaperones", "Translation, ribosomal structure & biogenesis" = "Translation, ribosomal structure\n& biogenesis")) +
#       coord_flip()
#   
#    plot.list2[[i]] = p
# }
# 
# #plot.list2 <- align_plots(plot.list2[[1]], plot.list2[[3]], plot.list2[[2]], align = "v")
# 
# p <- grid.arrange(plot.list2[[1]], ncol = 1)
# 
# ggsave("Figures/change-genes2.jpeg", p, units = "in", width = 7, height = 10, dpi = 600)
# 
# 
# sig.dif.fundf <- data.frame()
# 
# for(i in 1:length(names(sig.dif.fun))) {
#   sig.dif.fundf <- rbind(cbind(plant.water = names(sig.dif.fun)[i],sig.dif.fun[[i]]), sig.dif.fundf)
# }
# 
# ## OLD Colored by Gene ####
# show_col(viridis_pal(option = "turbo")(12))
# viridis_pal(option = "turbo")(12)
# 
# color_map <- c("Amino acid transport & metabolism" = "#30123BFF", 
#                "Carbohydrate transport & metabolism" = "#4454C4FF", 
#                "Cell wall/membrane/envelope biogenesis" = "#4490FEFF",
#                "Defense mechanisms" = "#1FC8DEFF", 
#                "Energy production & conversion" = "#29EFA2FF",
#                "General function prediction only" = "#7DFF56FF",
#                "Mobilome: prophages, transposons" = "#C1F334FF",
#                "Posttranslational modification, protein turnover, chaperones" = "#F1CA3AFF",
#                "Replication, recombination & repair" = "#FE922AFF",
#                "Signal transduction mechanisms" = "#EA4F0DFF",
#                "Symbiosis" = "#BE2102FF",
#                "Translation, ribosomal structure & biogenesis" = "#7A0403FF")
# 
# sig.dif.fundf$Fun_Cat_CE_ML <- as.character(sig.dif.fundf$Fun_Cat_CE_ML)
# sig.dif.fundf <- sig.dif.fundf[order(sig.dif.fundf$Fun_Cat_CE_ML),]
# 
# leg_plot <- ggplot(sig.dif.fundf, aes(x = Fun_Cat_CE_ML, y = estimate, col = Fun_Cat_CE_ML)) +
#   theme_classic() +
#   geom_point() +
#   scale_color_manual(values = color_map)
# 
# legend <- g_legend(leg_plot)
# 
# plot.list3 <- list()
# 
# for(i in names(res.list)) {
#   if(nrow(res.list[[i]]) > 0) {
#     
#  p <- ggplot(res.list[[i]], aes(x = gene, y = estimate, col = Fun_Cat_CE_ML)) + 
#     geom_point(size = 3) +
#     theme_bw() +
#     theme(
#       #axis.text.y = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 12),
#       plot.title = element_text(hjust = 0.5, size = 12),
#       legend.position = "none",
#       axis.title.x = element_blank(),
#         ) +
#     labs(title = plot.name.list[[i]], y = expression(paste(log[2], " fold change"))) +
#    ylim(-7,5) +
#    scale_color_manual(values = color_map, drop = F) +
#    coord_flip()
#   
#    plot.list3[[i]] = p
#   }
# }
# 
# p <- grid.arrange(plot.list3[[1]], 
#              plot.list3[[3]], plot.list3[[2]], legend, ncol = 2)
# 
# ggsave("Figures/Gene_change.jpeg", p, height = 10, width = 10, units = "in", dpi = 600)
