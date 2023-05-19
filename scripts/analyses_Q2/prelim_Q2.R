# Aims:
# 1. Preliminary code to answer Q2: 
# prediction accuracy ~ climate + interaction type + realm + treatment + distance between pairs

# Authors: Nathalie Chardon, GLL
# Date created: 12 May 2023
# Date updated: 15 May 2023 (NC)


#### LIBRARIES ####
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(brms)
library(priorsense)
library(bayesplot)


rm(list=ls()) 


#### INPUT FILES ####
dat <- readRDS('data/data_processing/log.prop.change.full.data.RDS')


# # OUTPUT FILES # #


#### FUNCTIONS ####

# Skew function (from: https://towardsdatascience.com/evaluating-bayesian-mixed-models-in-r-python-27d344a03016_)
skew <- function(y){ # Fisher-Pearson Skew function based on NIST definition
  n <- length(y)
  dif <- y - mean(y)
  skew_stat <- (sqrt(n-1)/(n-2))*n *(sum(dif^3)/(sum(dif^2)^1.5))
  return(skew_stat)
}




####################################################################################################

#### DATA EXPLORATION ####
# (revised from explore_Q2.R) # # 

####################################################################################################

# Mtg notes 10 May 2023 describe dataset Isaac made at end of working group:
# https://docs.google.com/document/d/1XvBeeiNZJZmDtivf1hJX9zN3PND9AgM6GxQeh-LpVY4/edit#heading=h.ebszuu9ld6v2

dat <- readRDS('data/data_processing/log.prop.change.full.data.RDS')
dat.int <- readRDS('data/data_processing/log.prop.change.interactions.RDS')

# CLIMATE
table(dat$CLIMATE1) # climate1 and climate2 are the same
table(dat$CLIMATE2) 


# REALM
table(dat$REALM1) #studies compare between different realms
table(dat$REALM2)

# visualize pairs across realm and climate
ggplot(dat, 
       aes(Prop.Change.Gn1, Prop.Change.Gn2, colour=UNIQUE.PAIR.ID)) +
  geom_point(pch=1) +
  facet_wrap(REALM1~CLIMATE1) +
  theme_base() +
  theme(legend.position = "none")

# how many data points per panel?
dat %>% 
  group_by(REALM1, CLIMATE1) %>% 
  summarize(length(UNIQUE.PAIR.ID))


# INTERACTION
table(dat.int$interaction_found) 

# NOTE: interaction types are listed as separate columns in different dataset 
# (log.prop.change.interactions.RDS) but should be in 1 column in this dataset as different numbers
# to use as explanatory variable (EB working on this)


# TREATMENT

# NOTE: need to create treatment variable (SS working on this)


# DISTANCE BETWEEN PAIRS
hh <- hist(dat$dist, breaks = 8)
hh$counts #number of pairs per 1km interval

# is distance between pairs related to between vs within study classification?
# yes, within study classifications only at 0 distances
# I suspect this will cause the model problems with convergence if no variation in distance across within category
plot(as.factor(dat$Type), as.numeric(dat$dist))
sum(subset(dat, Type=="Within")$dist) # all within distances are 0

# Create unique taxa-taxa pairs & preserving direction of pair interaction
dat <- dat %>% 
  mutate(taxa.taxa = paste(TAXA1, TAXA2, sep = '-')) %>% 
  mutate(resolved.taxa.taxa = paste(RESOLVED.TAXA1, RESOLVED.TAXA2, sep = '-')) 

# All original taxa pairs have enough data to potentially separate pairs in analyses
# note: what is ALL category?
table(dat$taxa.taxa)

# enough variation in climate and realm to use as explanatory variables?
summary.taxa <- dat %>% 
  group_by(taxa.taxa) %>% 
  summarize(n_distinct(CLIMATE1), n_distinct(REALM1))
View(summary.taxa)
# no there are not, withing taxa.taxa pair, mostly only 1 climate or 1 realm covered
# probably because climate is such a coarse binning

# But much fewer pairs in resolved taxa names
table(dat$resolved.taxa.taxa)


# IS TAXON-TAXON PAIRS CORRELATED WITH INTERACTION TYPE?

# NOTE: need updated interaction type data




####################################################################################################

#### FRAMEWORK FOR DATA ANALYSIS ####

####################################################################################################

# # Prep data

# generate random predictive accuracy data
# bounded between 0 and 1 using truncated normal distribution
# x <- rnorm(mean=.2, sd=.1, n=nrow(dat)+30000) # plus some buffer
# x <- x[x < 1 & x >= 0]
# x <- x[1:50000]
# hist(x)

# predictions: 
# higher for shorter distances between pairs
# higher predictive ability for terrestrial than aquatic, aquatic also has more variation
# takes about 2 min
out <- c()
for( i in 1:nrow(dat)){
  print(i)
  # higher prediction accuracy for terrestrial
  if(dat$REALM1[i]=="Terrestrial") {
    repeat{
      # 0 distances give mean of .8
      if(dat$dist[i]==0){
        x <- rnorm(mean=0.65, sd=.1, n=1) 
      }
      
      # as sites get further apart, give proportionally lower prediction accuracy
      if(dat$dist[i]>0){
        x <- rnorm(mean=(.65/(dat$dist[i]*10000))*2.5, sd=.1, n=1) # multiply km by scalar so no numbers less than 1
      }
      
      if (x < 1 & x >= 0){ # change these numbers to truncate normal at different points
        break
      }
    }
    out[i] <- x
    
  } else {
    # lower base prediction accuracy for aquatic and higher sd
    repeat{
      # 0 distances give mean of .8
      if(dat$dist[i]==0){
        x <- rnorm(mean=0.5, sd=.15, n=1) 
      }
      
      # as sites get further apart, give proportionally lower prediction accuracy
      if(dat$dist[i]>0){
        x <- rnorm(mean=(.5/(dat$dist[i]*10000))*2.5, sd=.15, n=1) # convert km to m, scales by log
      }
      
      if (x < 1 & x >= 0){ # change these numbers to truncate normal at different points
        break
      }
    }
    out[i] <- x
  }
 
}

hist(out)

dat$fake.pred.acc <- out

# looks okay
ggplot(dat, aes(fake.pred.acc)) + geom_histogram() + facet_wrap(~REALM1)

# fixed effects as factor
dat <- dat %>% 
  mutate(CLIMATE1 = factor(CLIMATE1)) %>% 
  mutate(REALM1 = factor(REALM1)) # add treatment and interaction once ready

# NOTE: assumes that REALM1 is always the realm of the response variable, and currently analyses only 
# done for Gn1 ~ Gn2. Need to switch when reverse direction is analyzed.

# check structure of variables
str(dat[, c('fake.pred.acc', 'CLIMATE1', 'REALM1', 'dist')])


# Randomly select 1% of data for trial run
data_model <- dat %>% 
  sample_frac(0.01)
data_model <- dat

# # Model Framework for all taxa
# prediction accuracy ~ climate + interaction type + realm + treatment + distance between pairs

#FAM <- gaussian(link = 'identity') # update based on requirements of y data
FAM <- Beta() # for moddelling proportiond data between 0 and 1

MODFORM <- bf(fake.pred.acc ~ CLIMATE1 + REALM1 + dist # intercept + fixed effect
                
               )  #don't need random slopes because included in Q1 models

mod <-brm(MODFORM, data = data_model, family = FAM, seed = 042023, #set seed
          
          chains = 3, iter = 5000, warmup = 1000, cores = 4#, #fitting information
          
          #file = 'outputs/brms_May2023/fake-pred-acc_mod.rds'
          )

mod <- readRDS(file="outputs/brms_May2023/fake-pred-acc_mod.rds")

#### Posterior Distribution ####

# Summary of posterior distribution for fixed and random effects:
# Estimate: mean of posterior distribution
# Est. Error: error associated with mean (standard error)
# CI intervals: if CI intervals encompass 0 => can't be sure effect isn't 0
summary(mod)

plot(conditional_effects(mod), ask = FALSE) #fitted parameters and their CI


#### Prior distribution ####

# What about priors? brms will try to pick reasonable defaults, which you can see:
prior_summary(mod) #can define priors in brm(priors = ...)


# Do priors overwhelm likelihood?
ps <- powerscale_sensitivity(mod) #look at 'diagnosis' column to see if prior isn't appropriate
unique(ps$sensitivity$diagnosis)


#### Model Fit ####

# sample size (Bulk_ESS & Tail_ESS) should > 1000 & rhat < 1.1
summary(mod)

plot(mod) #model convergence (L: does distribution mean match estimate? R: did all values get explored?)

# Posterior predictive check: Does data match model? (could be wrong distribution, not all effects modeled)
pp_check(mod, ndraws = 100) #posterior predictive checks - are predicted values similar to posterior distribution?

# Pairs plots to diagnose sampling problems (should show Gaussian blobs)
pairs(mod)

# Skewness: observed (black line) and simulated (grey distribution) SKEW metric (1000 simulated datasets)
ppc_stat(y = dat1perc$fake.pred.acc,
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
ppc_loo_pit_overlay(dat1perc$fake.pred.acc,
                    posterior_predict(mod),
                    lw = w)

# # CONCLUSION: model fits very well on all counts except for skew. 
# for new simulated pp checks aren't great even when modelling with beta distribution
# however model does converge and ESS and Rhat all look ok, so gonna wait to streamline this model until know better how real data is distributed


#### Models on each taxon-taxon pair ####

# split data into each taxon-taxon pair for separate analysis
# list where each element is a taxa.taxa pair, 29 elements
split.data <- split(dat, dat$taxa.taxa)

# check to make sure variation in CLIMATE, REALM across taxa.taxa in order to fit model
lapply(split.data, class)

# Randomly select 1% of data of each taxon.taxon pair for trial run
split.data_small <- lapply(split.data, sample_frac, .01)


# fit model on each element of list
FAM <- gaussian(link = 'identity')
MODFORM <- bf(fake.pred.acc ~ CLIMATE1 + REALM1 + dist)  #don't need random slopes because included in Q1 models

fitting.model <- function(data) {
  taxon.pair <- unique(data$taxa.taxa)
  FAM <- gaussian(link = 'identity')
  MODFORM <- bf(fake.pred.acc ~ CLIMATE1 + REALM1 + dist)  #don't need random slopes because included in Q1 models
  
  mod <-brm(MODFORM, data = data, family = FAM, seed = 042023, #set seed
            
            chains = 3, iter = 5000, warmup = 1000, cores = 4 #fitting information
            #file = paste0("outputs/brms_May2023/fake-pred-acc_mod_TaxonPair_", taxon.pair, ".rds")
                          )
  }

# throwing error because not enough variation in climate and realm to use as explanatory variable on single taxa.taxa models
lapply(split.data, fitting.model)

#mod <- readRDS(file="outputs/brms_May2023/fake-pred-acc_mod.rds")


#### NEXT STEPS ####

# - include interaction type and treatment once these variables are ready

# - generate fake.pred.acc with different distribution to more realistically simulate what these 
# values are likely to be (done), updated model to use beta distribution because simulated prop values between 0 and 1

# - split dataframe into unique taxon-taxon pairs and run model for each 
# (i.e. importance of each of these factors could be different dependent on taxon-taxon pair)
# done--wrote code to do this but stopped because not enough variation in climate/realm to fit model
# note also that 0 distances and Type:within are perfectly correlated (ie no non 0 distances for within site comparisons)

# - could also include taxon-taxon pair as fixed or random effect, depending on what information we 
# want to get out of it (fixed = want to know the difference; random = want to simply account for 
# different responses), but this adds a lot of degrees of freedom

# when know rough output distribution from Q1 analyses simulate data to reflect this
# and update model specs to reflect this, make sure statistical distribution etc is accurate for distribution of Y data


