# Aims:
# 1. Define random effect structure
# 2. Create temporal autocorrelation error term -- needed since using temporal change? TAC removed?
# 3. Fit GLM in brms with dummy DF

# Author: Nathalie Chardon
# Date created: 26 April 2023
# Date updated: 26 April 2023 (NC)

# # LIBRARIES # #
library(dplyr)
library(brms)     
library(MCMCvis)
library(cmdstanr)
library(posterior)
library(bayesplot)

# # INPUT FILES # #
#dat <- readRDS('data/dummy.dataset.RDS')

# # OUTPUT FILES # #

####################################################################################################

# # DEFINE RANDOM EFFECTS STRUCTURE # # 

####################################################################################################

# Data
#dat <- readRDS('data/dummy.dataset.RDS') #dummy data 

dat <- readRDS("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/preprocessing/log.prop.change.full.data.RDS") #full cleaned 
names(dat)

#subset to terrestrial and marine realms 
dat_terr<-subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")

dat_aqua<-subset(x = dat, subset = REALM1!="Terrestrial" & REALM2!="Terrestrial")

FAM <- gaussian(link = 'identity')

MODFORM <- bf(Log.prop.change.abun.Gn1 ~ Log.prop.change.abun.Gn2 + #intercept + fixed effect
                
                (Log.prop.change.abun.Gn2 | SERIES.l) + #rand slopes for time series length
                
                (Log.prop.change.abun.Gn2 | PairID)+    
                
                (Log.prop.change.abun.Gn2 | UNIQUE.PAIR.ID))   

#rand slopes for genus[i]-genus[j] nested within studyID[i]-studyID[j]

# Fit full model
mod.aqua.G1<-brm(MODFORM, data = dat_aqua, family = FAM, seed = 27042023, #set seed
                 control = list(adapt_delta=0.99, max_treedepth = 12),    
                 chains = 2, iter = 50000, warmup = 10000, cores = 40, 
                 backend="cmdstanr", threads = threading(20)) #fitting information


tools::psnice(pid = Sys.getpid(), value = 19)

saveRDS(mod.aqua.G1, file = "/home/shared/isaac.e/R/BioTIME/Data/Outputs/mod.aqua.G1.50K.10K.RDS") #save
mod.terr=readRDS("/home/shared/isaac.e/R/BioTIME/Data/Outputs/mod.terr.RDS")

# extract model output
cc <- coef(mod.aqua.G1)$UNIQUE.PAIR.ID
cc.slopes <- cc[1:dim(cc)[1], 1, dim(cc)[3]] #random slopes for genera-genera ID

hist(cc.slopes, breaks = 100)
abline(v=mean(cc.slopes),col="firebrick",lwd=3)

summary(mod.aqua.G1)#look at model outputs

dim(cc.slopes)

plot(mod.terr.G1)
