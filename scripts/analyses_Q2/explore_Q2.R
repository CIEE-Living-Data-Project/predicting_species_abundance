# Aims:
# 1. Visualize data to explore Q2

# Author: Gavia Lertzman-Lepofsky
# Date created: 27 April 2023
# Date updated: 27 April 2023 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(ggplot2)
library(ggthemes)

rm(list=ls()) 


# # INPUT FILES # #
dat <- readRDS('data/log.prop.change.with.meta.w.taxa.RDS')
int <- read.csv("data/preprocessing/genus_interaction_list.csv")



####################################################################################################

# # DATA EXPLORATION # # 

####################################################################################################

# climate1 and climate2 are the same
table(dat$CLIMATE1) # only temperate, tropical
table(dat$CLIMATE2) # only temperate, tropical

# studies compare between different realms
table(dat$REALM1)
table(dat$REALM2)
length(unique(dat$PairID)) # 63 pairs 

# is there marine tropical? 
# nope
subset(dat, CLIMATE1=="tropical" & REALM2=="Marine")

# how many types of interactions?
table(int$interaction)

# visualize pairs across realm and climate
ggplot(dat, 
       aes(Log.prop.change.abun.Gn1, Log.prop.change.abun.Gn2, colour=UNIQUE.PAIR.ID)) +
  geom_point(pch=1) +
  facet_wrap(REALM1~CLIMATE1) +
  theme_base() +
  theme(legend.position = "none")

# how many data points per panel?
dat %>% 
  group_by(REALM1, CLIMATE1) %>% 
  summarize(length(UNIQUE.PAIR.ID))

# plot abundance change through time only terrestrial temperate, coloured by something
# Year.T on x axis
dat %>% 
  filter(REALM1=="Terrestrial" & REALM2=="Terrestrial") %>% 
  mutate(taxa.pair=paste(TAXA1, TAXA2, sep="_")) %>% 
ggplot(., 
       aes(YEAR.T, Log.prop.change.abun.Gn1, colour=taxa.pair)) +
  geom_point(pch=1) +
  theme_base() +
  theme(legend.position = "none")


# a priori yes/no interaction classification in model
# network centrality in model
# to me, doesn't make sense to include type of trophic interaction--type pred/prey shouldn't
# can't do BioTime climate classification bc not enough variation, could maybe download worldclim or abs lat, but what is hypothesis? 
# are there certain types of taxon pairs that have better predictive accuracy
# eg. Check if prediction accuracy varies for different taxon pairs (e.g., mammal-bird)





