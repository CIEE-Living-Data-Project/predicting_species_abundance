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
#dat <- readRDS('data/preprocessing/dummy.dataset.RDS') #dummy data 
dat <- readRDS('data/data_processing/log.prop.change.full.data.UPDATED.RDS')#full cleaned 6/6/23
names(dat)

#subset to terrestrial and marine realms 
dat_terr<-subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")

dat_marine<-subset(x = dat, subset = REALM1!="Terrestrial" & REALM2!="Terrestrial")

#meta-model----
#run intercept only models for each genus x pairID x unique pairID combination
intmods_terr1=dat_terr%>%group_by(
  PairID, UNIQUE.PAIR.ID)%>%
  do(mod = try(lm(Prop.Change.Gn1 ~ 1,
                       data = .,
  )%>%broom::tidy(.))
  )

intmods_terr2=dat_terr%>%group_by(
  PairID, UNIQUE.PAIR.ID)%>%
  do(mod = try(lm(Prop.Change.Gn2 ~ 1,
                  data = .,
  )%>%broom::tidy(.))
  )

#unnest model columns & combine Gn1 and Gn2 
intmods_terr1=intmods_terr1%>%
  unnest(., cols=mod)%>%
  mutate(genus="Gn1")

intmods_terr2=intmods_terr2%>%
  unnest(., cols=mod)%>%
  mutate(genus="Gn2")

intmods_terr<-rbind(intmods_terr1, intmods_terr2) #82584

#filter anything with se=0, means only 1 value or all same values measured (i.e. no variation)
intmods_terr<-subset(intmods_terr, std.error>0) #81216

#pull wide on unique PairID
intmods_terrw<-group_by(intmods_terr, PairID, UNIQUE.PAIR.ID)%>% select(-statistic, -p.value, -term)%>%
  pivot_wider(names_from = genus, values_from = c(estimate, std.error))


#combine back with random effects 
dat_terrx<-select(dat_terr, -Gn1, -Gn2, -Prop.Change.Gn1, -Prop.Change.Gn2, -YEAR.T, -YEAR.T1) %>%
  distinct(.)

intmods_terrw2<-left_join(intmods_terrw, dat_terrx)%>%distinct(.) #not sure why this has more rows

#going to remove the cross metrics for now as unique pair IDs seem duplicated on this 
intmods_terrw2<-subset(intmods_terrw2, Metric!="CROSS") %>%
  distinct(.)#cuts data in half

#only 10% of data are between studies- filter out 
intmods_terrw2<-subset(intmods_terrw2, Type!="Between") %>%
  distinct(.)


#run hierarchical model with mean, sd from intercept models as joint response & predictors
MODDAT<-  intmods_terrw2
FAM <- gaussian(link = 'identity')

MODFORM<-bf(estimate_Gn1|resp_se(std.error_Gn1, sigma = TRUE)~ me(estimate_Gn2,std.error_Gn2) + 
          (estimate_Gn2 | PairID) +    
          (estimate_Gn2 | SERIES.l))  + set_mecor(FALSE) 
library(cmdstanr)
          
mod<-brm(MODFORM, MODDAT, FAM, #seed = 042023, #set seed
                         control = list(adapt_delta=0.99, max_treedepth = 12),    
                         chains = 3, iter = 2000, warmup = 500, cores = 4) 
         
save(mod, file = 'outputs/brms_June2023/meta_mod_q1.terrestrial_withinstudies.rds') #save          


#model outputs 
summary(mod)
ranef<-ranef(mod)

ranef_terr=as.data.frame(ranef$PairID)
ranef_terr$PairID<-row.names(ranef_terr)

ranef_terr2=as.data.frame(ranef$SERIES.l)
ranef_terr2$SERIES.l<-row.names(ranef_terr2)

#leave one out cross validation 
loo(mod)

loo(mod, moment_match = T, )

#ppchecks 
pp_check(mod, ndraws = 100) #this doesn't look great 


#full model----
#going to remove the cross metrics for now as unique pair IDs seem duplicated on this 
dat_terrx<-subset(dat_terr, Metric!="CROSS") %>%
  distinct(.)#

#only 10% of data are between studies- filter out 
dat_terry<-subset(dat_terrx, Type!="Between") %>%
  distinct(.)

MODDAT<-  dat_terry
FAM <- gaussian(link = 'identity')

library(cmdstanr)
set_cmdstan_path("C:/Users/court/Documents/.cmdstan/cmdstan-2.32.2")
MODFORM<-bf(Prop.Change.Gn1~ Prop.Change.Gn2 + 
              (Prop.Change.Gn2 | PairID) +    
              (Prop.Change.Gn2 | SERIES.l)+
              (Prop.Change.Gn2 | UNIQUE.PAIR.ID)) 


mod<-brm(MODFORM, MODDAT, FAM, #seed = 042023, #set seed
         control = list(adapt_delta=0.99, max_treedepth = 12),    
         chains = 3, iter = 2000, warmup = 500, cores = 4, 
         backend = "cmdstanr") 




