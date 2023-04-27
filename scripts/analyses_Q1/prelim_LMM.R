# Aims:
# 1. Define random effect structure
# 2. Create temporal autocorrelation error term -- needed since using temporal change? TAC removed?
# 3. Fit GLM in brms with dummy DF

# Author: Nathalie Chardon
# Date created: 26 April 2023
# Date updated: 26 April 2023 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)


rm(list=ls()) 


# # INPUT FILES # #
dat <- readRDS('data/dummy.dataset.RDS')

# # OUTPUT FILES # #




####################################################################################################

# # DEFINE RANDOM EFFECTS STRUCTURE # # 

####################################################################################################

# Data
dat <- readRDS('data/dummy.dataset.RDS')
  
FAM <- gaussian(link = 'identity')

MODFORM <- brms::bf(Log.prop.change.abun.Gn1 ~ 1 + Log.prop.change.bio.Gn2 + #intercept + fixed effect
                      
                      + (Log.prop.change.bio.Gn2 | SERIES.length) + #rand effect of time series length
                      
                      + (Log.prop.change.bio.Gn2 | PairID/UNIQUE.PAIR.ID)) 
                      #rand effect of genus[i]-genus[j] nested within studyID[i]-studyID[j]

# Fit full model

mod <- brms::brm(MODFORM, data = dat, family = FAM,
                 
                    chains = 3, iter = 5000, warmup = 1000, cores = 4, # fitting information
                 
                    file = 'outputs/brms_April2023/dummy_mod.rds', file_refit = 'on_change') #save


coefs <- coef(mod)[2]

hist(coefs$`PairID:UNIQUE.PAIR.ID`, breaks = 100)


### NOT RUN; STOP 26.4.2023 NC

# Check posterior
summary(mod)
plot(mod)
plot(conditional_effects(mod))

# Check prior
prior_summary(mod) 
ps <- powerscale_sensitivity(mod)
unique(ps$sensitivity$diagnosis)

# Model fit
pp_check(mod, ndraws = 100)
ppc_stat(y = dat$rel_rec, 
         yrep = posterior_predict(mod, ndraws = 1000), stat = "skew")
model_loo <- loo(mod, save_psis = TRUE, cores=4) 
w <- weights(model_loo$psis_object)
ppc_loo_pit_overlay(dat$rel_rec, 
                    posterior_predict(mod), lw = w)



