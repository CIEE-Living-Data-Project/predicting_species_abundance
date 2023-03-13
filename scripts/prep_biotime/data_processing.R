# Aims:
# 1. Filter collated BioTIME dataset to 10km pairs
# 2. Exclude biomass-only studies
# 3. Convert fake 0s to NA

# Author: Nathalie Chardon
# Date created: 16 Feb 2023
# Date updated: 16 Feb 2023 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)


rm(list=ls()) 


# # INPUT FILES # #
biotime.raw <- read.csv('data/prep_biotime/BioTIMEQuery_24_06_2021.csv')

# load('bio_pairs_10km.RData') #unique 10 km geographic-overlap years-taxa pairs (explore_biotime.R)
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use in lit review (explore_biotime.R)


# # OUTPUT FILES # #
load('data/tidy/collated_pairs.RData') #collated data filtered for 10km & overlapping pairs with fixed NAs (data_processing.R)




####################################################################################################

# # FILTER COLLATED DATA TO 10 KM PAIRS # # 

####################################################################################################

ii <- unique(meta.pairs$STUDY_ID) #unique study IDs for 10 km pairs to include


# Select matched studies from metadata file containing only pairs within 10km and temporally overlapping

collated.pairs <- biotime.raw %>% 
  filter(STUDY_ID %in% ii) #keep only study IDs from identified pairs

length(unique(collated.pairs$STUDY_ID)) == length(unique(meta.pairs$STUDY_ID)) #check filter




####################################################################################################

# # EXCLUDE BIOMASS-ONLY STUDIES # # 

####################################################################################################

# Create an exclusion dataframe

ee <- c(492, 334) ##IN PROGRESS: need to add any studies to exclude from HB

exclude <-  collated.pairs %>% 
  filter(STUDY_ID %in% ee) #filter for studies to exclude


# Remove biomass-only studies from collated.pairs

collated.pairs <- anti_join(collated.pairs, exclude) #keep only study IDs from identified pairs

length(unique(collated.pairs$STUDY_ID)) #should be original N - N of studies to exclude

collated.pairs %>% 
  filter(STUDY_ID %in% ee) #should be 0 rows




####################################################################################################

# # CONVERT FAKE 0s TO NA # # 

####################################################################################################

# Problem: authors report they don’t measure abundance or biomass, but then put 0’s in that column 
# Solution: change to NA with for loop

summary(collated.pairs$sum.allrawdata.ABUNDANCE) #n = 0 NAs
summary(collated.pairs$sum.allrawdata.BIOMASS) #n = 27141 NAs


# Add relevant metadata to collated pairs

x <- select(meta.pairs, STUDY_ID, ABUNDANCE_TYPE, BIOMASS_TYPE) #retain only relevant metadata 
collated.pairsx <- left_join(collated.pairs, x) #add to collated pairs


# Make abundance and biomass = NA if listed as NA in methods

collated.pairs <- mutate(collated.pairsx, sum.allrawdata.ABUNDANCE = #change abundance column
                         if_else(is.na(ABUNDANCE_TYPE), NA_real_, sum.allrawdata.ABUNDANCE)) %>%
                     mutate(sum.allrawdata.BIOMASS = #change biomass column
                              if_else(is.na(BIOMASS_TYPE), NA_real_, sum.allrawdata.BIOMASS))
                
         
# Save updated dataframe
save(collated.pairs, file = 'data/tidy/collated_pairs.RData')



