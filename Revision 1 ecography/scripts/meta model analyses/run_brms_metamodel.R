library(brms)
library(tidyverse)

# modelling timmmmme #####
## set priors ####
priors <-c(prior(normal(0,1), class = b), # set between -3 and 3 for slopes (fixed effects) 
           prior(normal(0,0.33), class = Intercept), # set between -1 and 1 for z score
           prior(cauchy(0,1), class = sd)) # set values b/w (0,1) for SE & random effects 

## define models ####
# model with z scores and total indivs SE as joint response
#load and try with small df first
load(file="Revision 1 ecography/output/prep_data/model_data_small.Rdata")

#if works then try full df
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")


MODFORM.indivs <- bf(z|resp_se(SE.total.indivs, sigma = FALSE) ~ 
                       scale.SERIES.l + treatment_yn_clean + 
                       scale.abs.lat +
                       interaction_present.factor +
                       scale.elev +
                       (1|STUDY_ID) + (1|resolved_taxa_pair))

# model with z scores and total indivs SE as joint response
MODFORM.sp <- bf(z|resp_se(SE.total.sp, sigma = FALSE) ~ 
                   scale.SERIES.l + treatment_yn_clean + 
                   scale.abs.lat +
                   interaction_present.factor +
                   scale.elev +
                   (1|STUDY_ID) + (1|resolved_taxa_pair))

## run models ####

# model 1
start_time <- Sys.time() 
metamod.indivs <- brm(MODFORM.indivs, moddat_small, cores=3, chains=3, prior=priors, 
                      # control = list(adapt_delta=0.9, max_treedepth = 12)
                      iter=3000, family=gaussian, file="Revision 1 ecography/output/meta model/SEindivs_small3K.rmd")
end_time <- Sys.time()
end_time - start_time

#save as Rdata file 
#save(metamod.indivs, file="Revision 1 ecography/output/meta model/metamod.individual.BAD.Rdata")
summary(metamod.indivs)

# model 2
start_time <- Sys.time() 
metamod.sp <- brm(MODFORM.sp, moddat_small,cores=3, chains=3, prior = priors, 
                  # control = list(adapt_delta=0.8, max_treedepth = 12), 
                  iter=5000, family=gaussian, 
                  file="Revision 1 ecography/output/meta model/SEspecies_small5K.rmd")
end_time <- Sys.time()
end_time - start_time

## model checks ####
#sample size (Bulk_ESS & Tail_ESS) should > 1000 & rhat < 1.1
summary(metamod)

plot(metamod.indivs) #model convergence (L: does distribution mean match estimate? R: did all values get explored?)

# posterior predictive checks - are predicted values similar to posterior distribution?
# takes a long time to run
pp_check(metamod.indivs) 

# Pairs plots to diagnose sampling problems (should show Gaussian blobs)
pairs(metamod)

plot(conditional_effects(metamod)) #fitted parameters and their CI

# get coefs for summary tables using fixef(mod)
# this reports mean of posterior
fixef(metamod)
ranef(metamod)

#run subset models####

# note, if models run quickly and we have lots of time, 
