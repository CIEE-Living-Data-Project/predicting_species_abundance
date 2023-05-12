# Aims:
# 1. Preliminary code to answer Q2: 
# prediction accuracy ~ climate + interaction type + realm + treatment + distance between pairs

# Authors: Nathalie Chardon
# Date created: 12 May 2023
# Date updated: 12 May 2023 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(brms)
library(priorsense)
library(bayesplot)


rm(list=ls()) 


# # INPUT FILES # #
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')


# # OUTPUT FILES # #


# # FUNCTIONS # # 

# Skew function (from: https://towardsdatascience.com/evaluating-bayesian-mixed-models-in-r-python-27d344a03016_)
skew <- function(y){ # Fisher-Pearson Skew function based on NIST definition
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}




####################################################################################################

# # DATA EXPLORATION (revised from explore_Q2.R) # # 

####################################################################################################

dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')

# CLIMATE
table(dat$CLIMATE1) # climate1 and climate2 are the same
table(dat$CLIMATE2) 


# REALM
table(dat$REALM1) #studies compare between different realms
table(dat$REALM2)

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


# INTERACTION
table(dat$interaction_found) 

# NOTE: interaction types are listed as separate columns but should be in 1 column as different numbers
# to use as explanatory variable (EB working on this)


# TREATMENT

# NOTE: need to create treatment variable (SS working on this)


# DISTANCE BETWEEN PAIRS

# NOTE: need to include this variable in the main dataset


# ENOUGH DATA TO SEPARATE BY TAXON-TAXON PAIRS?

# Create unique taxa-taxa pairs & preserving direction of pair interaction
dat <- dat %>% 
  mutate(taxa.taxa = paste(TAXA1, TAXA2, sep = '-')) 

# All pairs have enough data to potentially separate pairs in analyses
table(dat$taxa.taxa)


# IS TAXON-TAXON PAIRS CORRELATED WITH INTERACTION TYPE?

# NOTE: need updated interaction type data




####################################################################################################

# # FRAMEWORK FOR DATA ANALYSIS # # 

####################################################################################################

# # Prep data
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')

# generate random predictive accuracy data
dat$fake.pred.acc <- rnorm(nrow(dat))
hist(dat$fake.pred.acc, breaks = 50)

# fixed effects as factor
dat <- dat %>% 
  mutate(CLIMATE1 = factor(CLIMATE1)) %>% 
  mutate(REALM1 = factor(REALM1)) %>% 
  mutate(interaction_found = factor(interaction_found)) #add treatment and distance once ready

# NOTE: assumes that REALM1 is always the realm of the response variable, need to check if this 
# is true for analyses Gn 2 ~ Gn 1 (checking with IE & CC)

# check structure of variables
str(dat[, c('fake.pred.acc', 'CLIMATE1', 'REALM1', 'interaction_found')])

         
# Randomly select 1% of data for trail run
dat <- dat %>% 
  sample_frac(0.01)


# # Model Framework for all taxa
# prediction accuracy ~ climate + interaction type + realm + treatment + distance between pairs

FAM <- gaussian(link = 'identity')

MODFORM <- bf(fake.pred.acc ~ CLIMATE1 + REALM1 + interaction_found #intercept + fixed effect
                
               )  #don't need random slopes because included in Q1 models

mod <-brm(MODFORM, data = dat, family = FAM, seed = 042023, #set seed
              
              # control = list(adapt_delta=0.99, max_treedepth = 12), 
              
              chains = 3, iter = 5000, warmup = 1000, cores = 4) #fitting information         


# # Posterior Distribution

# Summary of posterior distribution for fixed and random effects:
# Estimate: mean of posterior distribution
# Est. Error: error associated with mean (standard error)
# CI intervals: if CI intervals encompass 0 => can't be sure effect isn't 0
summary(mod)

plot(conditional_effects(mod), ask = FALSE) #fitted parameters and their CI


# # Prior distribution

# What about priors? brms will try to pick reasonable defaults, which you can see:
prior_summary(mod) #can define priors in brm(priors = ...)


# Do priors overwhelm likelihood?
ps <- powerscale_sensitivity(mod) #look at 'diagnosis' column to see if prior isn't appropriate
unique(ps$sensitivity$diagnosis)


# # Model Fit

# sample size (Bulk_ESS & Tail_ESS) should > 1000 & rhat < 1.1
summary(mod)

plot(mod) #model convergence (L: does distribution mean match estimate? R: did all values get explored?)

# Posterior predictive check: Does data match model? (could be wrong distribution, not all effects modeled)
pp_check(mod, ndraws = 100) #posterior predictive checks - are predicted values similar to posterior distribution?

# Pairs plots to diagnose sampling problems (should show Gaussian blobs)
pairs(mod)

# Skewness: observed (black line) and simulated (grey distribution) SKEW metric (1000 simulated datasets)
ppc_stat(y = dat$fake.pred.acc,
         yrep = posterior_predict(mod, ndraws = 1000),
         stat = "skew")

# Leave-one-out Crossvalidation (LOO) for marginal (pointwise) posterior predictive checks
model_loo <- loo(mod, save_psis = TRUE, cores=4) #higher elpd => better fit

w <- weights(model_loo$psis_object)

# Probability integral transform (PIT) uses CDF properties to convert continuous random
# variables into uniformly distributed ones. However, this only holds exactly if the distribution
# used to generate the continuous random variables is the true distribution.
# -> If LOO-PIT values concentrated at 0/1 => possibly under-dispersed model
# -> If LOO-PIT values concentrated near 0.5 => possibly over-dispersed model
ppc_loo_pit_overlay(dat$fake.pred.acc,
                    posterior_predict(mod),
                    lw = w)

# # CONCLUSION: model fits very well on all counts except for skew. 


# # NEXT STEPS: 

# - split dataframe into unique taxon-taxon pairs and run model for each (i.e. importance of each of 
# these factors could be different dependent on taxon-taxon pair)

# - could also include taxon-taxon pair as fixed or random effect, depending on what information we 
# want to get out of it (fixed = want to know the difference; random = want to simply account for 
# different responses)

# - generate fake.pred.acc with different distribution to more realistically simulate what these 
# values are likely to be



