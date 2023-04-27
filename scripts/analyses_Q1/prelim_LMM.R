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
library(MCMCvis)

rm(list=ls()) 


# # INPUT FILES # #
#dat <- readRDS('data/dummy.dataset.RDS')

# # OUTPUT FILES # #


####################################################################################################

# # DEFINE RANDOM EFFECTS STRUCTURE # # 

####################################################################################################

# Data
#dat <- readRDS('data/dummy.dataset.RDS') #dummy data 

dat <- readRDS('data/preprocessing/log.prop.change.with.meta.RDS')#full cleaned 
names(dat)

#subset to terrestrial and marine realms 
dat_terr<-subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")

dat_marine<-subset(x = dat, subset = REALM1!="Terrestrial" & REALM2!="Terrestrial")

FAM <- gaussian(link = 'identity')

MODFORM <- bf(Log.prop.change.abun.Gn1 ~ Log.prop.change.abun.Gn2 + #intercept + fixed effect
                      
            (Log.prop.change.abun.Gn2 | SERIES.l) + #rand slopes for time series length
                      
            (Log.prop.change.abun.Gn2 | PairID)+    
            
            (Log.prop.change.abun.Gn2 | UNIQUE.PAIR.ID))   

#rand slopes for genus[i]-genus[j] nested within studyID[i]-studyID[j]

# Fit full model
mod<-brm(MODFORM, data = dat_terr, family = FAM, seed = 042023, #set seed
         control = list(adapt_delta=0.99, max_treedepth = 12),    
                   chains = 3, iter = 5000, warmup = 1000, cores = 4) #fitting information
                 
save(mod, file = 'outputs/brms_April2023/full_mod_q1.rds') #save

summary(mod)#look at model outputs
