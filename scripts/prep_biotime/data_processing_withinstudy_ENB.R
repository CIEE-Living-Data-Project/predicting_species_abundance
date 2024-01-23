# Aims:
# 1. Filter to terrestrial-only studies 
# 2. Exclude biomass-only studies
# 3. Convert fake 0s to NA

# Author: Nathalie Chardon, Emily Black
# Original version created: 16 Feb 2023
# Updated version created: 23 Jan 2024
# Date updated:  Jan 2024 (ENB)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)


rm(list=ls()) 


# # INPUT FILES # #
biotime.raw <- read.csv("~/Documents/Work and Career/LDP/Working Group/BioTIMEQuery_24_06_2021.csv")
metadata <- read.csv('data/prep_biotime/BioTIMEMetadata_24_06_2021.csv') #biotime metadata for 10 km pairs to use in lit review (explore_biotime.R)


# # OUTPUT FILES # #
#load('data/tidy/collated_within_studies.RData') #collated data filtered for 10km & overlapping pairs with fixed NAs (data_processing.R)


####################################################################################################

# # FILTER TO ONLY STUDIES IN THE TERRESTRIAL REALM # # 

####################################################################################################


x <- select(metadata, STUDY_ID, ABUNDANCE_TYPE, BIOMASS_TYPE, REALM) #retain only relevant metadata 
biotimex <- left_join(biotime.raw, x) #add to collated pairs

biotimex <- biotimex %>%
  filter(REALM=="Terrestrial")

#Check how many studies remain 
unique(biotimex$STUDY_ID)




####################################################################################################

# # CONVERT FAKE 0s TO NA # # 

####################################################################################################

# Problem: authors report they don’t measure abundance or biomass, but then put 0’s in that column 
# Solution: change to NA with for loop
# Add relevant metadata to collated pairs

# Make abundance and biomass = NA if listed as NA in methods
biotime <- mutate(biotimex, sum.allrawdata.ABUNDANCE = #change abundance column
                           if_else(is.na(ABUNDANCE_TYPE), NA_real_, sum.allrawdata.ABUNDANCE)) %>%
  mutate(sum.allrawdata.BIOMASS = #change biomass column
           if_else(is.na(BIOMASS_TYPE), NA_real_, sum.allrawdata.BIOMASS))


# Save updated dataframe
biotime<-rename(biotime, ID=STUDY_ID) #rename for clarity CGC
head(biotime)

save(biotime, file = 'data/prep_biotime/collated.biotime.nopairs.RData')



