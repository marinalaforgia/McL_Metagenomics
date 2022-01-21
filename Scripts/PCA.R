rm(list=ls())

library(plyr)
library(dplyr)
trait <- read.csv("/Users/Marina/Google Drive/02_McLaughlin_80Sites_Organized/Modified_CombinedFiles/McL_80SitesSpeciesTraits_012615.csv")
#trait <- filter(trait, Annual.Perennial == "Annual", Native.Exotic == "Native" & Grass.Forb.Shrub == "Forb" | Native.Exotic == "Exotic" & Grass.Forb.Shrub == "Grass", Family != "Fabaceae")
trait <- filter(trait, Annual.Perennial == "Annual", Native.Exotic == "Native", Grass.Forb.Shrub == "Forb", Family != "Fabaceae")
trait <- trait[,c(1:8,18)]
trait <- trait[complete.cases(trait),]
colnames(trait)[4] <- "Height"
colnames(trait)[3] <- "Seed Mass"
colnames(trait)[6] <- "Leaf Water Content"
colnames(trait)[5] <- "SLA"
colnames(trait)[7] <- "C:N"
trait.pca <- prcomp(trait[,3:7], scale = T)
biplot(trait.pca)

# 
plot(trait.pca$x[,1], trait.pca$x[,2] )
subset <- c("Clarkia purpurea", "Agoseris heterophylla", "Plantago erecta", "Calycadenia pauciflora", "Hemizonia congesta", "Lasthenia californica")
trait$row <- rownames(trait)
sub <- filter(trait, trait$Species_Name %in% subset)
new <- data.frame(row = rownames(trait.pca$x), PC1 = trait.pca$x[,1], PC2 = trait.pca$x[,2])
tog <- merge(sub, new, by = "row", all.y = F)

plot(trait.pca$x[,1], trait.pca$x[,2] )
text(tog$PC1, tog$PC2, tog$Species_Name, pos= 3 )
# 
# 
# new.all <- merge(trait, new, by = "row", all.y = F)
# plot(trait.pca$x[,1], trait.pca$x[,2] )
# text(new.all$PC1, new.all$PC2, new.all$Species_Name, pos= 3 )
# text(tog$PC1, tog$PC2, tog$Species_Name, pos= 3 )

#install.packages("ggfortify")
library(ggfortify)

#trait.pca <- prcomp(trait[,3:7], scale = T)
#arrow_ends <- layer_data(g, 2)[,c(2,4)]

autoplot(trait.pca, data = trait[trait$Grass.Forb.Shrub == "Forb",], loadings = T, loadings.label = T, label = F, loadings.label.vjust = -1.5, loadings.label.hjust = .6) +
  theme_classic() +
  ylim(-0.45, .4)

# new <- as.data.frame(trait.pca$x)
# new$rows <- row.names(trait.pca$x)
# trait$rows <- row.names(trait)
# trait <- merge(trait, new, by = "rows")
subset <- c("Clarkia purpurea", "Agoseris heterophylla", "Plantago erecta", "Calycadenia pauciflora", "Hemizonia congesta", "Lasthenia californica")
# subset.pca <- trait[trait$Species_Name %in% subset,]
# subset.pca$PC1 <- scale(subset.pca$PC1)
# subset.pca$PC2 <- scale(subset.pca$PC2)
# my.sp.s <- c("AGHE","CLPU","HECO","LACA","PLER", "CAPA")

sub <- as.data.frame(cbind(Species = c("AGHE","CLPU","HECO","LACA","PLER", "CAPA"), 
                              x = c(0.194, -0.08, -0.217, 0.08, -0.07, -.165), 
                              y = c(-0.2, -0.01, 0.15, 0, -0.15, -0.15)))
sub$x <- as.numeric(as.character(sub$x))
sub$y <- as.numeric(as.character(sub$y))



trait$yes <- ifelse(trait$Species_Name %in% subset, "Y", "N")

autoplot(trait.pca, data = trait, loadings = T, loadings.label = T, label = F, loadings.label.vjust = -1.8, loadings.label.hjust = 0.8, scaling = "columns", col = "yes") +
  # geom_text(aes(x = x, y = y, label = Species), data = subset) + 
  geom_label(aes(x = x, y = y, label = Species), data = sub, size = 3) +
  theme_classic() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("black", "blue")) + 
  ylim(-0.45, 0.4)

autoplot(trait.pca, data = trait, loadings = T, loadings.label = T, label = T, loadings.label.vjust = -1.8, loadings.label.hjust = 0.8, scaling = "none", col = "yes") +
  theme_classic() +
  ylim(-0.45, 0.4)
