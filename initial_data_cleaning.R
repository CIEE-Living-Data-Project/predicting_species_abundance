# Aims:
# 1. Visualize BioTIME data in space and complete the following tasks sprinkled throughout code below
# 2. TASK 1 - Figure out what to do with site visits and plots 
# 3. TASK 2 - Aggregate to the genus level 

# Authors: Emily Black
# Date created: 25 Apr 2023
# Date updated: 

# # LIBRARIES # #
library(tidyverse)
library("dplyr")
library("tibble")
library("readr") 
library("ggplot2")
library("magrittr")
library('ggmap')

rm(list=ls()) 


####################################################################################################

# # VISUALIZE DATA # # 

####################################################################################################

## Load Data
#these are the datasets in bioTime that overlap both in time (>=1 year) and space (within a 10km distance) 
load("data/tidy/collated_pairs.RData")
load('data/prep_biotime/bio_pairs_10km.RData') #metadata of overlapping studies
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use


#let's look at the data to get a sense of their structure 
head(bio.pairs)
head(collated.pairs)
months <- collated.pairs %>%
  group_by(ID)%>%
  filter(!is.na(MONTH))
plots_days <- collated.pairs %>%
  filter(!is.na(PLOT), !is.na(DAY))%>%
  nrow()
days <- collated.pairs %>%
  filter(!is.na(DAY))%>%
  nrow()
abundance <- collated.pairs %>%
  filter(!is.na(sum.allrawdata.ABUNDANCE)) %>%
  nrow()
biomass <- collated.pairs %>%
  filter(!is.na(sum.allrawdata.BIOMASS)) %>%
  nrow()
species <- collated.pairs %>%
  filter(is.na(GENUS_SPECIES)) %>%
  nrow()

#First pass: year and no plot 
#Sampling effort column
#Go coarse with the year 
#Summary statistics: mean, min, max, median, stdev, CoV
#abundance and biodiversity
#Scan through for unidentified

# #Find studies in southern hemisphere
# southern_hemisphere <- collated.pairs %>%
#   filter(LATITUDE<=0)
# #note: all southern hemisphere has months
# 
# southern_hemisphere <- southern_hemisphere %>%
#   mutate(aug_jul_year = if_else(as.numeric(MONTH) >= 8, as.numeric(YEAR), as.numeric(YEAR) - 1))
# 
# northern_hemisphere <- collated.pairs %>%
#   filter(LATITUDE>0)

collated_pairs_sampling_effort <- collated.pairs %>%
  group_by(ID, YEAR) %>%
  summarize(sampling_effort = n_distinct(PLOT))

