#Mia's changes to run_brms_metamodel.R April 3 2024
#I'm only putting the code I altered/added here. All the other code in that script I ran as needed before each of these changes.

## set priors ####
priors <- c(prior(normal(0,0.33),class=b),
            prior(normal(0,0.33),class=Intercept),
            prior(cauchy(0,0.5), class = sd))


#Mia trying full data:
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")

#look225 <- subset(moddat, STUDY_ID == "225")
#look39 <- subset(moddat, STUDY_ID == "39")
#removing those two studies
moddat <- filter(moddat, STUDY_ID != "225") 
moddat <- filter(moddat, STUDY_ID != "39")
unique(moddat$STUDY_ID)

MODFORM.sp <- bf(z|resp_se(SE.total.sp, sigma = FALSE) ~ 
                   scale.SERIES.l + treatment_yn_clean + 
                   scale.abs.lat +
                   interaction_present.factor +
                   # scale.elev + #removing this term
                   (1|STUDY_ID) + (1|resolved_taxa_pair))

#set up cmdstanr 
#install.packages("devtools")
library(devtools)

#install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos"))) #this was the only way the instalation worked for me

library(cmdstanr)
#install_cmdstan()
# set_cmdstan_path(path = NULL) #I think this isn't needed typically, but I read that this would help solve path problems, but it ended up not helping me
# if you run into issues here, I would recommend restarting R and reinstalling cmdstanr. After doing those two things (and a lot of other fussing, that I haven't transferred to this script), everything just worked
cmdstanr::check_cmdstan_toolchain(fix = TRUE)

#mia made changes to this model, including changiing to moddat below.
metamod.sp <- brm(MODFORM.sp, moddat,cores=3, chains=3, backend = "cmdstanr", threads = threading(4), prior = priors, 
                  # control = list(adapt_delta=0.8, max_treedepth = 12), 
                  iter=10000, warmup = 5000, thin = 5, family=gaussian,
                  file="Revision 1 ecography/output/meta model/SEspecies_full5k_miathreads2.rmd")
end_time <- Sys.time()
end_time - start_time
summary(metamod.sp)