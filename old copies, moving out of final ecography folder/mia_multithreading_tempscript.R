#Mia's changes to run_brms_metamodel.R April 3 2024

library(brms)
library(tidyverse)


## set priors ####
priors <- c(prior(normal(0,0.33),class=b),
            prior(normal(0,0.33),class=Intercept),
            prior(cauchy(0,0.5), class = sd))


## Load the full data ####
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")

#look225 <- subset(moddat, STUDY_ID == "225")
#look39 <- subset(moddat, STUDY_ID == "39")
#removing those two studies
moddat <- filter(moddat, STUDY_ID != "225") 
moddat <- filter(moddat, STUDY_ID != "39")



## define model ####
MODFORM.sp <- bf(z|resp_se(SE.total.sp, sigma = FALSE) ~ 
                   scale.SERIES.l + treatment_yn_clean + 
                   scale.abs.lat +
                   interaction_present.factor +
                   # scale.elev + #removing this term
                   (1|STUDY_ID) + (1|resolved_taxa_pair))

## set up cmdstanr to be able to do multithreading ####
#if this is your first time setting this up, un-comment the following lines. I hope it works easily following these steps!
#install.packages("devtools")
library(devtools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos"))) #this was the only way the installation worked for me

library(cmdstanr)
#install_cmdstan()
# set_cmdstan_path(path = NULL) #I think this isn't needed typically, but I read that this would help solve path problems, but it ended up not helping me
# if you run into issues here, I would recommend restarting R and reinstalling cmdstanr. After doing those two things (and a lot of other fussing, that I haven't transferred to this script), everything just worked
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

#Run the model:
start_time <- Sys.time() 
metamod.sp <- brm(MODFORM.sp, moddat,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
                  # control = list(adapt_delta=0.8, max_treedepth = 12), 
                  iter=10000, warmup = 5000, thin = 5, family=gaussian,
                  file="Revision 1 ecography/output/meta model/SEspecies_full5k_miathreads2.rmd") #please change this to a unique file name so you don't overwrite mine
end_time <- Sys.time() #run this line right after running the model, so that it runs as soon as the model finishes. That way we'll know exactly how long it took.
end_time - start_time
summary(metamod.sp)