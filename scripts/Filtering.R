#Filtering BIOtime (collated pairs)
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping

#load data
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/tidy/collated_pairs.RData")#collated pairs of overlapping studies
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/prep_biotime/bio_pairs_10km.RData") #metadata of overlapping studies

#subseted<-collated.pairs[1:10000,]

working<-collated.pairs %>%
  group_by(ID,YEAR,LATITUDE,LONGITUDE,SAMPLE_DESC,GENUS) %>%
  summarize(mean_abun=mean(sum.allrawdata.ABUNDANCE,na.rm=T),
            median_abun=median(sum.allrawdata.ABUNDANCE,na.rm=T),
            min_abun=min(sum.allrawdata.ABUNDANCE),
            max_abun=max(sum.allrawdata.ABUNDANCE),
            sd_abun=sd(sum.allrawdata.ABUNDANCE,na.rm=T),
            mean_bio=mean(sum.allrawdata.BIOMASS,na.rm=T),
            median_bio=median(sum.allrawdata.BIOMASS,na.rm=T),
            min_bio=min(sum.allrawdata.BIOMASS),
            max_bio=max(sum.allrawdata.BIOMASS),
            sd_bio=sd(sum.allrawdata.BIOMASS,na.rm=T))
            
#sampling effort
working <- collated.pairs %>%
  group_by(ID, YEAR) %>%
  summarize(EFFORT.YEAR = n_distinct(PLOT)) %>%
  left_join(working, ., by = join_by(ID, YEAR))

hist(working$EFFORT.YEAR)
hist(working$mean_bio)

hist(working$mean_abun)

table(working$ID[which(working$EFFORT.YEAR>600)])












