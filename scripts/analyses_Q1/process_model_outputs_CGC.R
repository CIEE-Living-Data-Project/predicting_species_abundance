
# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(MCMCvis)
library(tidybayes)
library(ggdist)

#full model----
#w/ corr term
load(file = "outputs/Aug2023/mod_q1.terrestrial_withinstudies.Rdata") 
#w/o corr term 
load(file = 'outputs/Aug2023/mod_q1.terrestrial_withinstudiesv2.Rdata')

#assess convergence issues 
summary(mod)
#convergence issue only on 'cor' term - re-run without 
MCMCtrace(mod) #this takes a while- but mostly looks good!

#Unique study ID
slopes<-coef(object = mod)
slopes<-as.data.frame(slopes$UNIQUE.PAIR.ID)
slopes<-select(slopes, contains("Gn2"))
slopes$UniquePairID<-row.names(slopes)

save(slopes, file = "outputs/Aug2023/randomslopes_q1modelv2.Rdata")

pp_check(mod) #captures mean and variance ok but misses magnitude of mean
#loox<-loo(mod) #too big won't work

#plot results 

