#Repair NA taxa names 
# Created by: ENB
# Created: 28 Sept 2023 (ENB)
# Last Modified: 


#Part 0: Read in neccesary packages
rm(list=ls()) 

# # LIBRARIES # #
library(tidyverse)
library("dplyr")
library("tibble")
library("readr") 
library("ggplot2")
library("magrittr")
library(progress)




#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 1. Read in data with resolved taxa pairs and disturbances

dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23
head(dat)
#read in the disturbance data
load("data/prep_biotime/meta_pairs_10km.RData")
meta.pairs$STUDY_ID <- as.character(meta.pairs$STUDY_ID)
dat <- left_join(dat, meta.pairs[, c(1,46)],
                 by=c("ID1" = "STUDY_ID"))

#Make resolved taxa pairs column
sorted_words <- apply(dat[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
dat$resolved_taxa_pair <- sorted_words

#Get NA rows from the dataset
dat_na <- dat %>%
  filter(grepl("\\.NA", resolved.taxa.pair))

