#summarize sourcetracker results
library(vroom)
library(tidyverse)
library(magrittr)

meta <- vroom("McL_metag_metadata.csv")

sp_arc_asv_names <- t(read.delim("sourcetracker2_family_arc_12142022_rare600/mixing_proportions.txt", row.names = 1))
family_arc_asv <- as_tibble(t(read.delim("sourcetracker2_family_arc_12142022_rare600/mixing_proportions.txt", row.names = 1)))

sp_all_asv_names <- t(read.delim("sourcetracker2_family_all_12142022_rare25000/mixing_proportions.txt", row.names = 1))
fam_all_asv <- as_tibble(t(read.delim("sourcetracker2_family_all_12142022_rare25000/mixing_proportions.txt", row.names = 1)))

sp_bac_asv_names <- t(read.delim("sourcetracker2_family_bac_1214022_rare25000/mixing_proportions.txt", row.names = 1))
fam_bac_asv <- as_tibble(t(read.delim("sourcetracker2_family_bac_1214022_rare25000/mixing_proportions.txt", row.names = 1)))

func_all_asv_names <- t(read.delim("sourcetracker2_cog20_all_12212022_rare25000/mixing_proportions.txt", row.names = 1))
func_fam_all_asv <- as_tibble(t(read.delim("sourcetracker2_cog20_all_12212022_rare25000/mixing_proportions.txt", row.names = 1)))


row_names <- rownames(func_all_asv_names)
func_fam_all_asv$Sample <- row_names

row_names <- rownames(sp_arc_asv_names)
family_arc_asv$Sample <- row_names

row_names <- rownames(sp_all_asv_names)
fam_all_asv$Sample <- row_names

row_names <- rownames(sp_bac_asv_names)
fam_bac_asv$Sample <- row_names

family_arc_asv_met <- left_join(family_arc_asv, meta, by = c("Sample" = "SampleID"))
fam_all_asv_met <- left_join(fam_all_asv, meta, by = c("Sample" = "SampleID"))
fam_bac_asv_met <- left_join(fam_bac_asv, meta, by = c("Sample" = "SampleID"))
func_fam_all_asv_met <- left_join(func_fam_all_asv, meta, by=c("Sample" = "SampleID"))

### bacteria + archaea combined

grouped_res.fg <- group_by(fam_all_asv_met, PlotTreatment)
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
sig.text$cld <- c("b", "b", "b", # grass
                      "a", "a", "a", # forb
                      "c", "c", "c") # unknown
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

ggsave("all_family_sourcetracker2_25000set.png", dpi=300, device = "png", height = 5, width= 6, units = "in")



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


## function

grouped_res.fg <- group_by(func_fam_all_asv_met, PlotTreatment)
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
sig.text$cld <- c("b", "b", "b", # grass
                  "a", "a", "a", # forb
                  "c", "c", "c") # unknown
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

ggsave("cog20_family_sourcetracker2_25000set.png", dpi=300, device = "png", height = 5, width= 6, units = "in")



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





# ### archaea
# 
# grouped_res.fg <- group_by(family_arc_asv_met, PlotTreatment)
# #calculate mean of means + sd (sd(x)/sqrt(length(x)))
# avgs_res.fg.mean <- summarise(grouped_res.fg, 
#                               grass=100*mean(grasses), 
#                               forb=100*mean(forb), 
#                               unknown=100*mean(Unknown))
# 
# avgs_res.fg.sd <- summarise(grouped_res.fg, 
#                             grass=100*sd(grasses)/sqrt(length(grasses)),
#                             forb=100*sd(forb)/sqrt(length(forb)),
#                             unknown=100*sd(Unknown)/sqrt(length(Unknown)))
# 
# 
# test <- avgs_res.fg.mean %>%
#   gather(key = "source", value = "mean", -PlotTreatment)
# 
# test.2 <- avgs_res.fg.sd %>%
#   gather(key="source", value="sd", -PlotTreatment)
# 
# test3 <- full_join(test,test.2)
# 
# wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"]) 
# wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"]) 
# wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Control"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Control"]) 
# 
# wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"]) 
# wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 
# wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Drought"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Drought"]) 
# 
# wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"]) 
# wilcox.test(grouped_res.fg$forb[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"]) 
# wilcox.test(grouped_res.fg$grasses[grouped_res.fg$PlotTreatment=="Watering"], grouped_res.fg$Unknown[grouped_res.fg$PlotTreatment=="Watering"]) 
# 
# sig.text <- as_data_frame(test3$source)
# sig.text$cld <- c("a", "a", "b", # grass
#                   "a", "a", "a", # forb
#                   "b", "b", "c") # unknown
# sig.text$max <- test3$mean + test3$sd
# 
# 
# ggplot(test3, aes(x = source, y = mean, fill = source)) + 
#   geom_bar(stat = "identity", position = position_dodge()) +
#   geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
#                 width = .4, position = position_dodge(.9)) + 
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
#         text = element_text(size = 14)) + 
#   ylab("Mean % community derived from each source") + facet_wrap(~PlotTreatment)+
#   scale_fill_manual(values = c("magenta4", "#1F968BFF", "#B8DE29FF")) + 
#   xlab("") +
#   geom_text(aes(x = sig.text$value, label = sig.text$cld, y = sig.text$max,
#                 size = 8,
#                 #vjust = 0.3),
#                 vjust = -0.7),
#             col = "black", show.legend = F) +
#   ylim(NA, max(test3$mean + test3$sd) + 0.055*max(test3$mean))
# 
# #+ theme(legend.position = "none")
# 
# ggsave("archael_family_sourcetracker2_25000set.png", dpi=300, device = "png", height = 5, width= 6, units = "in")


## All, remove unknowns

grouped_res.fg <- group_by(fam_all_asv_met, PlotTreatment)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res.fg.mean <- summarise(grouped_res.fg, 
                              grass=100*mean(grasses), 
                              forb=100*mean(forb))

avgs_res.fg.sd <- summarise(grouped_res.fg, 
                            grass=100*sd(grasses)/sqrt(length(grasses)),
                            forb=100*sd(forb)/sqrt(length(forb)))


test <- avgs_res.fg.mean %>%
  gather(key = "source", value = "mean", -PlotTreatment)

test.2 <- avgs_res.fg.sd %>%
  gather(key="source", value="sd", -PlotTreatment)

test3 <- full_join(test,test.2)


test3$PlotTreatment <- factor(test3$PlotTreatment, levels = c("Drought", "Control", "Watering"), labels=c("Drought", "Control", "Watered"))

test3$Source <- factor(test3$source, labels = c("forb" = "Native", "grass" = "Invasive"))

sig.text <- as_data_frame(test3$Source)
sig.text$cld <- c("b", "b", "b", # grass
                  "a", "a", "a" # forb
) 
sig.text$max <- test3$mean + test3$sd


ggplot(test3, aes(x = Source, y = mean, fill = Source)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = .5),
        text = element_text(size = 14)) + 
  ylab("Mean % community derived from each source") + facet_wrap(~PlotTreatment)+
  scale_fill_manual(values = c("magenta4", "#1F968BFF", "#B8DE29FF")) + 
  xlab("Source") +
  geom_text(aes(x = sig.text$value, label = sig.text$cld, y = sig.text$max,
                size = 8,
                #vjust = 0.3),
                vjust = -0.7),
            col = "black", show.legend = F) +
  ylim(NA, max(test3$mean + test3$sd) + 0.055*max(test3$mean)) + theme(legend.position = "none")

ggsave("all_family_sourcetracker2_nounk_25000est.png", dpi=300, device = "png", height = 5, width= 6, units = "in")
ggsave("all_family_sourcetracker2_nounk_25000est.pdf", dpi=300, device = "pdf", height = 5, width= 6, units = "in")


## function, remove unknowns

grouped_res.fg <- group_by(func_fam_all_asv_met, PlotTreatment)
#calculate mean of means + sd (sd(x)/sqrt(length(x)))
avgs_res.fg.mean <- summarise(grouped_res.fg, 
                              grass=100*mean(grasses), 
                              forb=100*mean(forb))

avgs_res.fg.sd <- summarise(grouped_res.fg, 
                            grass=100*sd(grasses)/sqrt(length(grasses)),
                            forb=100*sd(forb)/sqrt(length(forb)))


test <- avgs_res.fg.mean %>%
  gather(key = "source", value = "mean", -PlotTreatment)

test.2 <- avgs_res.fg.sd %>%
  gather(key="source", value="sd", -PlotTreatment)

test3 <- full_join(test,test.2)


test3$Source <- factor(test3$source, labels = c("forb" = "Native", "grass" = "Invasive"))
test3$PlotTreatment <- factor(test3$PlotTreatment, levels = c("Drought", "Control", "Watering"), labels=c("Drought", "Control", "Watered"))

sig.text <- as_data_frame(test3$Source)
sig.text$cld <- c("b", "b", "b", # grass
                  "a", "a", "a" # forb
) 
sig.text$max <- test3$mean + test3$sd


ggplot(test3, aes(x = Source, y = mean, fill = Source)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
                width = .4, position = position_dodge(.9)) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = .5),
        text = element_text(size = 14)) + 
  ylab("Mean % community derived from each source") + facet_wrap(~PlotTreatment)+
  scale_fill_manual(values = c("magenta4", "#1F968BFF", "#B8DE29FF")) + 
  xlab("Source") +
  geom_text(aes(x = sig.text$value, label = sig.text$cld, y = sig.text$max,
                size = 8,
                #vjust = 0.3),
                vjust = -0.7),
            col = "black", show.legend = F) +
  ylim(NA, max(test3$mean + test3$sd) + 0.055*max(test3$mean)) + theme(legend.position = "none")

ggsave("all_cog20_sourcetracker2_nounk_25000est.png", dpi=300, device = "png", height = 5, width= 6, units = "in")
ggsave("all_cog20_sourcetracker2_nounk_25000est.pdf", dpi=300, device = "pdf", height = 5, width= 6, units = "in")


# ## Archaeal, remove unknowns
# 
# grouped_res.fg <- group_by(family_arc_asv_met, PlotTreatment)
# #calculate mean of means + sd (sd(x)/sqrt(length(x)))
# avgs_res.fg.mean <- summarise(grouped_res.fg, 
#                               grass=100*mean(grasses), 
#                               forb=100*mean(forb))
# 
# avgs_res.fg.sd <- summarise(grouped_res.fg, 
#                             grass=100*sd(grasses)/sqrt(length(grasses)),
#                             forb=100*sd(forb)/sqrt(length(forb)))
# 
# 
# test <- avgs_res.fg.mean %>%
#   gather(key = "source", value = "mean", -PlotTreatment)
# 
# test.2 <- avgs_res.fg.sd %>%
#   gather(key="source", value="sd", -PlotTreatment)
# 
# test3 <- full_join(test,test.2)
# 
# sig.text <- as_data_frame(test3$source)
# sig.text$cld <- c("a", "a", "b", # grass
#                   "a", "a", "a" # forb
# ) 
# sig.text$max <- test3$mean + test3$sd
# 
# test3$PlotTreatment <- factor(test3$PlotTreatment, levels = c("Drought", "Control", "Watering"), labels=c("Drought", "Control", "Watered"))
# 
# ggplot(test3, aes(x = source, y = mean, fill = source)) + 
#   geom_bar(stat = "identity", position = position_dodge()) +
#   geom_errorbar(aes(ymin = (mean - sd), ymax = (mean + sd)), 
#                 width = .4, position = position_dodge(.9)) + 
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = .5),
#         text = element_text(size = 14)) + 
#   ylab("Mean % community derived from each source") + facet_wrap(~PlotTreatment)+
#   scale_fill_manual(values = c("magenta4", "#1F968BFF", "#B8DE29FF")) + 
#   xlab("") +
#   geom_text(aes(x = sig.text$value, label = sig.text$cld, y = sig.text$max,
#                 size = 8,
#                 #vjust = 0.3),
#                 vjust = -0.7),
#             col = "black", show.legend = F) +
#   ylim(NA, max(test3$mean + test3$sd) + 0.055*max(test3$mean))
# 
# #+ theme(legend.position = "none")
# 
# ggsave("archaeal_family_sourcetracker2_nounk_25000et.png", dpi=300, device = "png", height = 5, width= 6, units = "in")
# ggsave("archaeal_family_sourcetracker2_nounk_25000set.pdf", dpi=300, device = "pdf", height = 5, width= 6, units = "in")











## Plotting SourceTracker2 results of specific families/functions found with DESeq 

deg.tax <- read_csv("tax-gene-DEG_121922.csv")
#length(unique(deg$DEG)) 20 gene

deg.tax$gene <- paste0("gene_", deg.tax$gene)

#so get all the genes that have the same function as the DEG (since collapsed)
deg_tax <- kaiju_tax_file$Family[rownames(kaiju_tax_file) %in% deg.tax$gene]

genes_in_fam <- rownames(kaiju_tax_file[kaiju_tax_file$Family %in% deg_tax,])

physeq.deg <- prune_taxa(genes_in_fam, physeq.forDEG)

physeq.deg.fam <- tax_glom(physeq.deg, "Family", NArm =FALSE)

physeq.deg.fam.melt <- psmelt(physeq.deg.fam)

only_treat <- meta[c(1,7,12)]
physeq.deg.fam.melt_meta <- left_join(physeq.deg.fam.melt, only_treat, by=c("Sample"= "SampleID"))



## what taxa ## 

file_dir <- "sourcetracker2_family_all_12142022_rare25000"
files <- list.files(file_dir, full.names = TRUE, pattern = "feature_table.txt$")
combined_files <- bind_rows(lapply(files, read.delim), .id = "id")

samples <- str_replace(files, "sourcetracker2_family_all_12142022_rare25000/", "")
sample_name <- str_replace(samples, ".feature_table.txt", "")
id <- as.character(1:29)

sample_df <- data.frame(sample_name, id)

combined_tax <- left_join(combined_files, sample_df)

only_treat <- meta[c(1,7)]

combined_tax_subset_meta <- left_join(combined_tax, only_treat, by=c("sample_name"= "SampleID"))

combined_tax_subset_meta_gather <- combined_tax_subset_meta %>%
  gather(key="Family", value = "Counts", -id, -sample_name, -X, -PlotTreatment)


summed <- combined_tax_subset_meta_gather %>% 
  filter(!X == "Unknown") %>%
  group_by(Family, PlotTreatment) %>%
  summarize(sum=sum(as.numeric(Counts)))

summed_2 <- combined_tax_subset_meta_gather %>% 
  filter(!X == "Unknown") %>%
  group_by(Family, X, PlotTreatment) %>%
  summarize(prop=Counts)

summed_3 <- left_join(summed_2, summed)

summed_4 <- summed_3 %>%
  group_by(Family, X, PlotTreatment) %>%
  summarize(mean_prop = sum(prop/sum))

summed_4$X <- factor(summed_4$X, levels = c( "forb", "grasses"))


summed_5 <- summed_4 %>%
  filter(Family %in% physeq.deg.fam.melt$Family) %>%
  arrange(desc(mean_prop)) 

summed_5$PlotTreatment <- factor(summed_5$PlotTreatment, levels = c("Drought", "Control", "Watering"), labels=c("Drought", "Control", "Watered"))

ggplot(summed_5, aes(y = mean_prop, x = Family, col = X, fill = X)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() + facet_grid(~PlotTreatment)+
  scale_fill_manual(name="Source", values = c("magenta4", "#1F968BFF"), labels = c("Native", "Invasive")) +
  labs(y = "Proportion of Families") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "white") +
  theme_classic() +
  coord_flip()  



ggsave("all_family_sourcetracker2_deg_nounk_121922.png", dpi=300, device = "png", height = 6, width= 9, units = "in")
ggsave("all_family_sourcetracker2_deg_nounk_121922.pdf", dpi=300, device = "pdf", height = 6, width= 9, units = "in")



## FUNCTION ###


## what taxa ## 

file_dir <- "sourcetracker2_cog20_all_12212022_rare25000"
files <- list.files(file_dir, full.names = TRUE, pattern = "feature_table.txt$")
combined_files <- bind_rows(lapply(files, read.delim), .id = "id")

samples <- str_replace(files, "sourcetracker2_cog20_all_12212022_rare25000/", "")
sample_name <- str_replace(samples, ".feature_table.txt", "")
id <- as.character(1:29)

sample_df <- data.frame(sample_name, id)

combined_tax <- left_join(combined_files, sample_df)

only_treat <- meta[c(1,7)]

combined_tax_subset_meta <- left_join(combined_tax, only_treat, by=c("sample_name"= "SampleID"))

combined_tax_subset_meta_gather <- combined_tax_subset_meta %>%
  gather(key="COG_FUN", value = "Counts", -id, -sample_name, -X, -PlotTreatment)


## DEGs 

#cog function
deg <- vroom("fun-gene-DEG_121922.csv")

deg$gene <- paste0("gene_", deg$gene)

deg_fun <- left_join(deg, select(GL_fun_table_filt.noCOG, gene, COG20_FUNCTION, COG20_CATEGORY))

deg_fun1 <- deg_fun %>% 
  mutate(COG20_FUNCTION2 = str_replace_all(COG20_FUNCTION, " ", ".")) %>% 
  mutate(COG20_FUNCTION2 = str_replace_all(COG20_FUNCTION2, "-", ".")) %>% 
  mutate(COG20_FUNCTION2 = str_replace_all(COG20_FUNCTION2, "[(]", ".")) %>%
  mutate(COG20_FUNCTION2 = str_replace_all(COG20_FUNCTION2, "[)]", ".")) %>%
  mutate(COG20_FUNCTION2 = str_replace_all(COG20_FUNCTION2, ":", ".")) %>%
  mutate(COG20_FUNCTION2 = str_replace_all(COG20_FUNCTION2, "!", ".")) %>%
  mutate(COG20_FUNCTION2 = str_replace_all(COG20_FUNCTION2, ",", "."))

#so get all the genes that have the same function as the DEG (since collapsed)
deg_cogs <- GL_fun_table_filt.noCOG$COG20_FUNCTION[GL_fun_table_filt.noCOG$gene %in% deg$gene]
deg_cogs2 <- unique(deg_fun1$COG20_FUNCTION2[deg_fun1$gene %in% deg$gene])

genes_with_cogs <- GL_fun_table_filt.noCOG$gene[GL_fun_table_filt.noCOG$COG20_FUNCTION %in% deg_cogs]

physeq.cog20.deg <- prune_taxa(genes_with_cogs, physeq.cog20)

physeq.cog20.deg.F <- tax_glom(physeq.cog20.deg, "COG20_FUNCTION", NArm =FALSE)

physeq.cog20.deg.F.melt <- psmelt(physeq.cog20.deg.F)

summed <- combined_tax_subset_meta_gather %>% 
  filter(!X == "Unknown") %>%
  group_by(COG_FUN, PlotTreatment) %>%
  summarize(sum=sum(as.numeric(Counts)))

summed_2 <- combined_tax_subset_meta_gather %>% 
  filter(!X == "Unknown") %>%
  group_by(COG_FUN, X, PlotTreatment) %>%
  summarize(prop=Counts)

summed_3 <- left_join(summed_2, summed)

summed_4 <- summed_3 %>%
  group_by(COG_FUN, X, PlotTreatment) %>%
  summarize(mean_prop = sum(prop/sum))

summed_4$X <- factor(summed_4$X, levels = c( "forb", "grasses"))

summed_5 <- summed_4 %>%
  filter(COG_FUN %in% deg_cogs2) %>%
  left_join(select(deg_fun1, COG20_FUNCTION, COG20_CATEGORY, COG20_FUNCTION2), by=c("COG_FUN" = "COG20_FUNCTION2"))%>%
  arrange(desc(mean_prop)) 

fam_to_rem <- unique(summed_5$COG_FUN[summed_5$mean_prop == "NaN"])

summed_5$PlotTreatment <- factor(summed_5$PlotTreatment, levels = c("Drought", "Control", "Watering"), labels=c("Drought", "Control", "Watered"))

summed_6 <- summed_5 %>%
  filter(!COG_FUN %in% fam_to_rem)

summed_6_plus <- distinct(summed_6) %>%
  mutate(COG20_FUN_clean = str_replace(COG20_FUNCTION, "!!!.*", "")) %>% 
  mutate(COG20_FUN_clean = str_replace(COG20_FUN_clean, " \\(PDB.*", "")) %>%
  mutate(COG20_FUN_clean = str_replace(COG20_FUN_clean, " \\(PUB.*", "")) %>%
  mutate(COG20_FUN_clean = ifelse(COG20_FUN_clean == "Uncharacterized conserved protein, contains Zn-finger domain of CDGSH type", "Uncharacterized Fe-S cluster protein YjdI (YjdI)", COG20_FUN_clean))

ggplot(summed_6_plus, aes(y = mean_prop, x = COG20_FUN_clean, col = X, fill = X)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() + facet_grid(~PlotTreatment)+
  scale_fill_manual(name="Source", values = c("magenta4", "#1F968BFF"), labels = c("Native", "Invasive")) +
  labs(y = "Proportion of COG functions") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "white") +
  theme_classic() +
  coord_flip() + xlab("")

ggsave("cogfunction_sourcetracker2.png", dpi=300, device = "png", height = 6, width= 12, units = "in")




summed_7 <- distinct(summed_6) %>%
  mutate(COG20_CATEGORY_clean = str_replace(COG20_CATEGORY, "!!!.*", "")) %>%
  group_by(COG20_CATEGORY_clean, X, PlotTreatment) %>%
  summarize(mean_prop2 = mean(mean_prop))

ggplot(summed_7, aes(y = mean_prop2, x = COG20_CATEGORY_clean, col = X, fill = X)) +
  geom_bar(position="stack", stat="identity", col = "black") + 
  theme_bw() + facet_grid(~PlotTreatment)+
  scale_fill_manual(name="Source", values = c("magenta4", "#1F968BFF"), labels = c("Native", "Invasive")) +
  labs(y = "Proportion of COG categories") +
  geom_hline(yintercept=0.5, linetype="dashed", color = "white") +
  theme_classic() +
  coord_flip() + xlab("")

ggsave("cogfunction_meanofcat_sourcetracker2.png", dpi=300, device = "png", height = 6, width= 12, units = "in")










