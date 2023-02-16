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


# # WORKING DIRECTORIES # #
raw_dat <- '~/Desktop/Code/predicting_species_abundance_offline/' #WD for NC (biotime query 1.2 GB)
figs <- '~/Desktop/Code/predicting_species_abundance/figures/explore_biotime_Dec2022/'
tidy_dat <- '~/Desktop/Code/predicting_species_abundance/data/prep_biotime/'

# # INPUT FILES # #
setwd(raw_dat)
biotime.raw <- read.csv('BioTIMEQuery_24_06_2021.csv')

setwd(tidy_dat)
# load('bio_pairs_10km.RData') #unique 10 km geographic-overlap years-taxa pairs (explore_biotime.R)
load('meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use in lit review (explore_biotime.R)


# # OUTPUT FILES # #




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

ee <- c(334, 274) ##IN PROGRESS: need to add any studies to exclude from NC & HB & RL

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






