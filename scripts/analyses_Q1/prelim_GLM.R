# Aims:
# 1. Establish best random effect structure
# 2. Exclude biomass-only studies
# 3. Convert fake 0s to NA

# Author: Nathalie Chardon
# Date created: 16 Feb 2023
# Date updated: 15 Mar 2023 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)


rm(list=ls()) 


# # INPUT FILES # #
biotime.raw <- read.csv('data/prep_biotime/BioTIMEQuery_24_06_2021.csv')
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use in lit review (explore_biotime.R)


# # OUTPUT FILES # #
#load('data/tidy/collated_pairs.RData') #collated data filtered for 10km & overlapping pairs with fixed NAs (data_processing.R)




####################################################################################################

# # FILTER COLLATED DATA TO 10 KM PAIRS # # 

####################################################################################################
