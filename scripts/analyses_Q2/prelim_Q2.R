# Aims:
# 1. Preliminary code to answer Q2: 
# prediction accuracy ~ climate + interaction type + realm + treatment + distance between pairs

# Authors: Nathalie Chardon
# Date created: 12 May 2023
# Date updated: 12 May 2023 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)


rm(list=ls()) 


# # INPUT FILES # #
dat <- readRDS('data/log.prop.change.with.meta.w.taxa.RDS')


# # OUTPUT FILES # #




####################################################################################################

# # DATA EXPLORATION # # 

####################################################################################################
