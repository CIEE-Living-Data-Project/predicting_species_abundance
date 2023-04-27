# Aims:
# 1. Explore approaches to answer Q2

# Author: Nathalie Chardon
# Date created: 27 April 2023
# Date updated: 26 April 2023 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(lattice)


rm(list=ls()) 


# # INPUT FILES # #
dat <- readRDS('data/log.prop.change.with.meta.w.taxa.RDS')


# # OUTPUT FILES # #




####################################################################################################

# # DATA EXPLORATION # # 

####################################################################################################

# Look at data
dat <- readRDS('data/log.prop.change.with.meta.w.taxa.RDS')

head(dat)

unique(dat$REALM1)
unique(dat$CLIMATE1)


# Make sure CLIMATE is always listed as the same within pair
foo <- dat %>% 
  
  mutate(mismatch = if_else(CLIMATE1 != CLIMATE2, 1, 0))

unique(foo$mismatch) #CLIMATE always the same


# Look at climate data
xyplot(dat$Log.prop.change.abun.Gn1 ~ dat$Log.prop.change.abun.Gn2 | dat$CLIMATE1)

# # Not much tropical data => climate a valuable fixed effect? 


# Filter to temperate-terrestrial (where almost all data is)
tt <- dat %>% 
  filter(REALM1 == 'Terrestrial' | REALM2 == 'Terrestrial') %>% 
  filter(CLIMATE1 == 'Temperate') 


# Look at taxa pairs
unique(tt$RESOLVED.TAXA1)
unique(tt$RESOLVED.TAXA2)

tt <- mutate(tt, taxa.pair = paste(RESOLVED.TAXA1, RESOLVED.TAXA2, sep = '_')) #create taxa pairs
unique(tt$taxa.pair)

# Taxa pairs for abundance
abund.tt <- filter(tt, is.na(tt$Log.prop.change.abun.Gn1)!=TRUE & #filter out non-abund studies
                     is.na(tt$Log.prop.change.abun.Gn2)!=TRUE)

xyplot(abund.tt$Log.prop.change.abun.Gn1 ~ abund.tt$Log.prop.change.abun.Gn2 | abund.tt$taxa.pair, 
       xlab = 'Abundance Change Genus 2', ylab = 'Abundance Change Genus 1')

# Taxa pairs for biomass
bio.tt <- filter(tt, is.na(tt$Log.prop.change.bio.Gn1)!=TRUE & #filter out non-biomass studies
                     is.na(tt$Log.prop.change.bio.Gn2)!=TRUE)

xyplot(bio.tt$Log.prop.change.bio.Gn1 ~ bio.tt$Log.prop.change.bio.Gn2 | bio.tt$taxa.pair, 
       xlab = 'Biomass Change Genus 2', ylab = 'Biomass Change Genus 1')


