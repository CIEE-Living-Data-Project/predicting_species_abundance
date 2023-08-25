
# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(MCMCvis)
library(tidybayes)
library(ggdist)

rm(list=ls()) 


# Read in data
#dat <- readRDS('data/preprocessing/dummy.dataset.RDS') #dummy data 
#dat <- readRDS('data/data_processing/log.prop.change.full.data.UPDATED.RDS')#full cleaned 6/6/23
#names(dat)
#data with interaction info 
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23

#subset to terrestrial and marine realms 
dat_terr<-subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")

dat_marine<-subset(x = dat, subset = REALM1!="Terrestrial" & REALM2!="Terrestrial")


#full model----
#remove the cross metrics for now as unique pair IDs seem duplicated on this 
dat_terrx<-subset(dat_terr, Metric!="CROSS") %>%
  distinct(.)#

#only 10% of data are between studies- filter out 
dat_terry<-subset(dat_terrx, Type!="Between") %>%
  distinct(.)

MODDAT<-  dat_terry
FAM <- gaussian(link = 'identity')

MODFORM<-bf(Prop.Change.Gn1~ Prop.Change.Gn2 + 
              (Prop.Change.Gn2 | UNIQUE.PAIR.ID)) 

mod<-brm(MODFORM, MODDAT, FAM, #seed = 042023, #set seed
         control = list(adapt_delta=0.99, max_treedepth = 12),    
         chains = 3, iter = 5000, warmup = 500, cores = 4) 
        # backend = "cmdstanr") 

save(mod, file = 'outputs/Aug2023/mod_q1.terrestrial_withinstudies.Rdata') #save (but too big to push to Git)         

#re-run without correlation term-bad rhats
MODFORM2<-bf(Prop.Change.Gn1~ Prop.Change.Gn2 + 
               (Prop.Change.Gn2 || UNIQUE.PAIR.ID)) 

mod<-brm(MODFORM2, MODDAT, FAM, #seed = 042023, #set seed
         control = list(adapt_delta=0.99, max_treedepth = 12),    
         chains = 3, iter = 10000, warmup = 500, cores = 4) 
# backend = "cmdstanr") 

summary(mod)
save(mod, file = 'outputs/Aug2023/mod_q1.terrestrial_withinstudies.Rdata')

#try in lme4----
library(lme4)
library(optimx)

dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')
dat_terr<-subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
dat_terr<-subset(dat_terr, Metric!="CROSS" &  Type!="Between") %>%
  distinct(.)

#filter out any unique Pair IDs where se=0, means only 1 value or all same values measured (i.e. no variation at group level) 
load(file = 'outputs/Aug2023/intercept_only_models.Rdata')
filterpairs<-subset(intmods_terr, std.error==0)%>%select(UNIQUE.PAIR.ID)

dat_terrx<-anti_join(dat_terr, filterpairs)

lmermod<-lmer(Prop.Change.Gn1~ Prop.Change.Gn2 + 
                (Prop.Change.Gn2 | UNIQUE.PAIR.ID), dat_terrx, control=lmerControl(optimizer="optimx", optCtrl=list(method="nlminbwrap"))) 

#singular fit
#trying all optimizers 
glmm_all = allFit(lmermod)


#filter out any unique Pair IDs where n<10 obs

intmods_terr_ct=dat_terr%>%group_by(
  PairID, UNIQUE.PAIR.ID)%>%summarize(count=n())
intmods_terr_ct_filt<-subset(intmods_terr_ct, count<10)%>%select(UNIQUE.PAIR.ID)
intmods_terr_ct_filt$PairID<-NULL

dat_terry<-anti_join(dat_terrx, intmods_terr_ct_filt)

#now re-try
lmermod<-lmer(Prop.Change.Gn1~ Prop.Change.Gn2 + 
                (Prop.Change.Gn2 | UNIQUE.PAIR.ID), dat_terry, control=lmerControl(optimizer="optimx", optCtrl=list(method="nlminb"))) 
rslopes<-ranef(lmermod_optx)
rslopesdf<-as.data.frame(rslopes$UNIQUE.PAIR.ID$Prop.Change.Gn2)


intmods_terr_ct=dat_terr%>%group_by(
  PairID, UNIQUE.PAIR.ID)%>%summarize(count=n())
intmods_terr_ct_filt<-subset(intmods_terr_ct, count<20)%>%select(UNIQUE.PAIR.ID)
intmods_terr_ct_filt$PairID<-NULL

dat_terrz<-anti_join(dat_terry, intmods_terr_ct_filt)


#now re-try
lmermod<-lmer(Prop.Change.Gn1~ Prop.Change.Gn2 + 
                (Prop.Change.Gn2 | UNIQUE.PAIR.ID), dat_terrz, control=lmerControl(optimizer="optimx", optCtrl=list(method="nlminb"))) 
#this had singular fit but not convergence issues...could work 


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

#save these so don't have to re-run 
save(intmods_terr, file = 'outputs/Aug2023/intercept_only_models.Rdata')


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
              (estimate_Gn2 | SERIES.l)) + set_mecor(FALSE) 

mod<-brm(MODFORM, MODDAT, FAM, #seed = 042023, #set seed
         control = list(adapt_delta=0.99, max_treedepth = 12),    
         chains = 3, iter = 5000, warmup = 500, cores = 4) 

save(mod, file = 'outputs/brms_July2023/meta_mod_q1.terrestrial_withinstudies.Rdata') #save          

#run models with interaction, climate info 
#intmods_terrw2<-left_join(intmods_terrw, dat_terrx)%>%distinct(.) #not sure why this has more rows

#going to remove the cross metrics for now as unique pair IDs seem duplicated on this 
#intmods_terrw2<-subset(intmods_terrw2, Metric!="CROSS") %>%
#  distinct(.)#cuts data in half

#interaction info only for between studies- subset 7/20/23 - Emily re-running 
#intmods_terrw3<-subset(intmods_terrw2, Type=="Between") %>%
#  distinct(.)

#run hierarchical model with mean, sd from intercept models as joint response & predictors
MODDAT<-  intmods_terrw2
FAM <- gaussian(link = 'identity')

MODFORM<-bf(estimate_Gn1|resp_se(std.error_Gn1, sigma = TRUE)~ me(estimate_Gn2,std.error_Gn2) + 
              (estimate_Gn2 | SERIES.l) +
              (estimate_Gn2 | CLIMATE1) +
              (estimate_Gn2 | TAXA1) +
              (estimate_Gn2 | interaction_present)+
              (estimate_Gn2 | interaction_benefit)+
              (estimate_Gn2 | interaction_type)) + set_mecor(FALSE) 

mod<-brm(MODFORM, MODDAT, FAM, #seed = 042023, #set seed
         control = list(adapt_delta=0.99, max_treedepth = 12),    
         chains = 3, iter = 5000, warmup = 500, cores = 4) 

save(mod, file = 'outputs/August2023/meta_mod_q1.terrestrial_betweenstudies_interactions.Rdata') #save          
