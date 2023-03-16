library(vroom)


#obtain read counts per sample sequenced 

#import data
read_counts <- vroom("read_counts_per_sample.txt") 
meta <- vroom("McL_metag_metadata.csv")
lat_long <-vroom("Lat-Long.csv")
meta$Plot <- as.character(meta$Plot)


meta <- left_join(meta,lat_long)
write.csv(meta, "meta.for.ncbi.csv")

read_counts.sum <- read_counts %>% 
  group_by(Sample) %>%
  summarize(count = sum(Reads),
            count_paired = count/2) %>%
  mutate(Sample = str_replace(Sample, '-', '_'))

            

read_counts.meta <- left_join(read_counts.sum, meta, by=c("Sample" = "SampleID")) %>%
  mutate(New_Treatment = ifelse(New_Treatment == "grasses", "invasives", ifelse(New_Treatment == "forbs", "natives", "mix")),
         PlotTreatment = ifelse(PlotTreatment == "Watering", "Watered", as.character(PlotTreatment))) %>%
  mutate(TotalTreat = paste0(PlotTreatment, ' x ', New_Treatment)) %>%
  mutate(TotalTreat = ifelse(TotalTreat == "None x mix", paste0('Experimental control x ', SampleType), as.character(TotalTreat))) %>%
  select(Sample, count, TotalTreat)


write.csv(read_counts.meta, "supp.table1.csv")




