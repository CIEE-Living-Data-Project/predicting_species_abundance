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
library(taxize)




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
  filter(grepl("\\.NA|NA\\.", resolved_taxa_pair)) %>%
  filter(Type=="Within") %>%
  filter(REALM1=="Terrestrial")

#Check what group we are working with
unique(dat_na$ORGANISMS1)
unique(dat_na$ORGANISMS2)

#Great! These are easy enough to fix, 
dat_na$RESOLVED.TAXA1 <- ifelse(
  is.na(dat_na$RESOLVED.TAXA1),
  ifelse(
    dat_na$ORGANISMS1 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat_na$ORGANISMS1 == "birds", "Aves", 
           ifelse(dat_na$ORGANISMS1 == "rodents", "Mammalia", NA)
    )
  ),
  dat_na$RESOLVED.TAXA1
)



dat_na$RESOLVED.TAXA2 <- ifelse(
  is.na(dat_na$RESOLVED.TAXA2),
  ifelse(
    dat_na$ORGANISMS2 %in% c("insects", "grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat_na$ORGANISMS2 == "birds", "Aves", 
           ifelse(dat_na$ORGANISMS2 == "rodents", "Mammalia", NA)
    )
  ),
  dat_na$RESOLVED.TAXA2
)


#Now run on the ful dataset

dat_filtered <-   dat %>%
  filter(Type=="Within") %>%
  filter(REALM1=="Terrestrial")

dat_filtered$RESOLVED.TAXA1 <- ifelse(
  is.na(dat_filtered$RESOLVED.TAXA1),
  ifelse(
    dat_filtered$ORGANISMS1 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat_filtered$ORGANISMS1 == "birds", "Aves", 
           ifelse(dat_filtered$ORGANISMS1 == "rodents", "Mammalia", NA)
    )
  ),
  dat_filtered$RESOLVED.TAXA1
)



dat_filtered$RESOLVED.TAXA2 <- ifelse(
  is.na(dat_filtered$RESOLVED.TAXA2),
  ifelse(
    dat_filtered$ORGANISMS2 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat_filtered$ORGANISMS2 == "birds", "Aves", 
           ifelse(dat_filtered$ORGANISMS2 == "rodents", "Mammalia", NA)
    )
  ),
  dat_filtered$RESOLVED.TAXA2
)

sorted_words <- apply(dat_filtered[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
dat_filtered$resolved_taxa_pair <- sorted_words
unique(dat_filtered$resolved_taxa_pair)
table(dat_filtered$resolved_taxa_pair)

#check what didnt work
#Get NA rows from the dataset
dat_na_check <- dat_filtered %>%
  filter(grepl("\\.NA", resolved_taxa_pair)) %>%
  filter(Type=="Within") %>%
  filter(REALM1=="Terrestrial")

#Great! It worked! I'll add this to the figure generation to check it out 
