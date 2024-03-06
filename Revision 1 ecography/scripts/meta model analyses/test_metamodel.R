library(tidyverse)
library(brms)

#### reading in data and cleaning ####
# current full dataset, update w/ final version
interactions.dat <- readRDS("Revision 1 ecography/output/prep_data/within.study.updated.interactions.020724ENB.RDS") 
abun.dat <- read.csv("Revision 1 ecography/output/prep_data/results.abundance_23Feb24.csv")

# clean up taxa columns
unique(interactions.dat$TAXA) # terrestrial plants, birds, terrestrial inverts, mammals, reptiles and all. All cateogry is only herpetofauna, mainly snakes, but one frog
# lump reptiles and all into herpetofauna category
interactions.dat$TAXA <- gsub("Reptiles", "Herpetofauna", interactions.dat$TAXA)
interactions.dat$TAXA <- gsub("All", "Herpetofauna", interactions.dat$TAXA)
unique(interactions.dat$TAXA)

# clean up organisms column
interactions.dat$ORGANISMS <- gsub("birds|Bird|Birds", "Birds", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub('breeding Birds|breeding bird pairs|Breeding Birds', "Breeding Birds", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub('Plants|plants', "Plants", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub('Small mammal|Small mammals|small mammals', "Small mammals", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub('Insects specifically Coleoptera and Lepidoptera|small mammals|Beetles|Butterflies', "Beetles, Butterflies, Moths", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub("Lizards", "Herpetofauna", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub("rodents", "Small mammals", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub("waterBirds|ducks", "Water birds", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub("Acrididae (grasshoppers)", "Grasshoppers", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub("Acrididae \\(", "", interactions.dat$ORGANISMS)
interactions.dat$ORGANISMS <- gsub("\\)", "", interactions.dat$ORGANISMS)

sort(unique(interactions.dat$ORGANISMS))

# clean up resolved taxa pair column NEEED A RESOLVED TAXA PAIR COLUMN

# clean up habitat column
interactions.dat$HABITAT <- gsub("Chalk grassland|Grassland", "Grassland", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("Boreal forest|Boreal vegetation / Taiga|Scandinavian taiga", "Boreal", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("Tallgrass prairie|Tallgrass prairie gallery forest and riparian edge|Savanna/ Tallgrass prairie|Woodland", "Tallgrass prairie and woodland", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("Praire", "Prairie", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("upper montane forests above 2700 feet in the North", "Alpine", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("Desert Wildlife Refuge|Desert/ Grassland|semiarid thorn scrub/Desert/ Grassland|creosotebush site and grassland site|semiarid thorn scrub", "Desert, grassland, thorn scrub, creosotebush", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("Urban / Desert|Urban / Desert / Riparian / Agricultural|large pine plantatino in Nacogdeoches County. Texa|forests agricultural fields and meadows", "Human modified", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("Northern mixed prairie|Mixed|Forest and grassland", "Mixed prairie and forest", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("Mixed prairie and forest Conifer-Hardwood Forest", "Mixed Conifer-Hardwood Forest", interactions.dat$HABITAT)
interactions.dat$HABITAT <- gsub("birch forest", "Deciduous forest", interactions.dat$HABITAT)

sort(unique(interactions.dat$HABITAT))

# calculate abs latitude metric
interactions.dat$abs.lat <- abs(interactions.dat$LATITUDE)

# combine data sets
interactions.dat_trim <- distinct(select(interactions.dat, c("TS_ID", "STUDY_PLOT", "Gn1", "Gn2", "STUDY_ID", "abs.lat",
                                                    "REALM", "CLIMATE", "HABITAT", "TAXA", "ORGANISMS", "LATITUDE", "LONGITUDE","interaction_present")))
alldat <- left_join(abun.dat, interactions.dat_trim, by = c("TS_ID", "STUDY_PLOT", "Gn1", "Gn2"))

# check to make sure that every pair has interaction info
sum(is.na(alldat$interaction_present))

# add in worldclim data
worldclim <- read.csv("Revision 1 ecography/output/prep_data/worldclim.csv")[,-1]
alldat <- left_join(alldat, worldclim)

# calculate pearson's cor for all time series 
alldat <- alldat %>% 
  group_by(TS_ID) %>%
  mutate(cor=cor(Log.prop.change.Gn1, Log.prop.change.Gn2)) 

# warnings caused by time series length where 0 changes through entire time series

# subset dataset with cleaned disturbance information

#### prep data for model ####
# subset to only relatively strong cors 
alldat_trim <- subset(alldat, cor<0.8 & cor>-0.8) 

hist(alldat_trim$cor)
max(alldat_trim$cor)
min(alldat_trim$cor)

# calculate other metrics for SE
alldat_trim <- alldat_trim %>% 
  group_by(TS_ID) %>% 
  mutate(total.indivs = sum(Abun.T.Gn1)+sum(Abun.T.Gn2),
         total.sp = sum(Rich.T.Gn1)+sum(Rich.T.Gn2),
         abs.total.indivsGn1mGn2 = abs(sum(Abun.T.Gn1) - sum(Abun.T.Gn2)),
         abs.total.spGn1mGn2 = abs(sum(Rich.T.Gn1) - sum(Rich.T.Gn2)))
  
# calculate z scores and SE w/ sample sizes 
alldat_trim <- alldat_trim %>% 
  mutate(z=0.5*log((1+cor)/(1-cor))) #eq 3.11 in meta-analysis book 
  #mutate(CIlow = (-1.96 /sqrt(SERIES.l-3))+ z)%>%
  #mutate(CIhigh=  (1.96 /sqrt(SERIES.l-3))+ z)%>%

# note throws warnings for all SE where sample size is less than 3, but that is expected and okay
# time series length is only one of these for which all n > 3
# takes 30 min
alldat_trim <- alldat_trim %>% 
  mutate(SE.timeseries = (1/sqrt(SERIES.l-3)), # time series length
         SE.total.indivs = (1/sqrt(total.indivs-3)), # total number of individuals
         SE.total.sp = (1/sqrt(total.sp-3)), # total number of species
         SE.abs.total.indivsGn1mGn2 = (1/sqrt(abs.total.indivsGn1mGn2-3)), # abundance difference between genus 1 and 2
         SE.abs.total.spGn1mGn2 = (1/sqrt(abs.total.spGn1mGn2-3)) # # richness difference between genus 1 and 2
         ) #eq 3.12 in meta-analysis book #can update with different sample size definitions here 


# filter df to remove repeat values 
# ie keeping unique cases of z-scores
moddat <- alldat_trim %>% 
  select(-Log.prop.change.Gn1, -Log.prop.change.Gn2, -Abs.change.Gn1, -Abun.T1.Gn1, -Abun.T.Gn2, -Abun.T1.Gn2, 
         -Abs.change.Gn2, -Abun.T.Gn1, -Length.Unique.Values.Gn1, -Length.Unique.Values.Gn2,
         -Rich.T.Gn1, -Rich.T1.Gn1, -Rich.T.Gn2, -Rich.T1.Gn2, -SS.T.Gn1, -SS.T1.Gn1, -SS.T.Gn2, -SS.T1.Gn2,
         -X, -X.1, -YEAR.T, -YEAR.T1) %>%
  distinct(.)
names(moddat)

# look at some to check   
select(moddat, cor, z, SE.timeseries)
hist(moddat$SE.timeseries)

# explore number of rows with various SE metric
# drop the time series length one bc doesn't make sense to also use it as a question we are testing
# trim all data to most restrictive of these SE calculations
# compare models with goodness of fit stats
length(na.omit(moddat$SE.timeseries)) # 411,472
length(na.omit(moddat$SE.total.indivs)) # 411,423
length(na.omit(moddat$SE.total.sp)) #  411,472
length(na.omit(moddat$abs.total.spGn1mGn2)) # 165,082
length(na.omit(moddat$abs.total.indivsGn1mGn2)) # 401,804

# explore how related
# update column numbers when re-run mutate
pairs(moddat[,c(27,28,30,31,32)])


# NEED RESOLVED TAXA PAIRS AND DISTURBANCE INFO
# NOTE THAT DISTURBANCE INFO WILL need to get expanded to >10km for full dataset 
# probably need to go back to the full metadata table methods for each study 
# or just run on clean data subset 
# load("data/prep_biotime/meta_pairs_10km.RData")
# moddat<-left_join(moddat, meta.pairs[, c(1,20,46)])
# names(moddat)

length(unique(moddat$TS_ID))

# make sure interaction_present is encoded as a factor
# scale inputs here bc doesn't work in model formula (maybe why it's so slow???)
moddat$scale.abs.lat <- scale(moddat$abs.lat)
moddat$scale.SERIES.l <- scale(moddat$SERIES.l)  
moddat$interaction_present.factor <- as.factor(moddat$interaction_present)
moddat$scale.elev <- scale(moddat$Elevation)
  
# set some priors 
priors <-c(prior(normal(0,0.33), class = Intercept), #set between -1 and 1 for z score
          prior(cauchy(0,5), class = sd)) #set lower bound 0 for SE, values b/w (0,1)

# model with z scores and SE as joint response
# don't use SE.timeseries
# add interaction between disturbance and time series length
# elevation? Tseas? Pseas?
MODFORM <- bf(z|resp_se(SE.timeseries) ~ 
                scale.SERIES.l + 
                #Resolved.taxa.pair + NEED THIS IN DATA
                scale.abs.lat +
                #treatment_yn + #is there some sort of disturbance yes/no (fertilizer, fire, grazing etc)
                interaction_present.factor + #does GLOBI record these genera as potentially interacting
                scale.elev +
                (1|STUDY_ID))

start_time <- Sys.time() 
metamod <- brm(MODFORM, moddat,
             control = list(adapt_delta=0.8, max_treedepth = 12), cores=3, chains=3, 
             iter=5000, family=gaussian, prior = priors, 
             file="Revision 1 ecography/output/meta model/SEtimeseries length test.rmd")
end_time <- Sys.time()
end_time - start_time
# 1.73074 days
# 5 transitions after warmup that exceeded max treedepth


#sample size (Bulk_ESS & Tail_ESS) should > 1000 & rhat < 1.1
summary(metamod)

plot(metamod) #model convergence (L: does distribution mean match estimate? R: did all values get explored?)

# posterior predictive checks - are predicted values similar to posterior distribution?
# takesa a long time to run
pp_check(metamod, ndraws = 2) 

# Pairs plots to diagnose sampling problems (should show Gaussian blobs)
pairs(metamod)


plot(conditional_effects(metamod)) #fitted parameters and their CI
# need to make sure id is reversible 

# get coefs for summary tables using fixef(mod)
# this reports mean of posterior
fixef(metamod)

##### plot random slopes ####
library(tidybayes)
mod <- metamod

# extract the draws corresponding to posterior distributions of the overall mean and standard deviation of observations
mod %>%
  spread_draws(b_Intercept, sigma) %>%
  head(10)

# median and 95% quantile interval of the variables, we can apply median_qi()
mod %>%
  #spread_draws(b_Intercept, sigma) %>% #wide
  gather_draws(b_Intercept, b_scaleSERIES.l, b_scaleabs.lat, 
               b_scaleElevation, b_interaction_present,
               sigma) %>% #long
  mean_qi()

# summary for all regions [only when have random effect]
mod %>%
  spread_draws(r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR, ]) %>%
  median_qi()

# convergence diagnostics
mod %>%
  #spread_draws(r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR, ]) %>%
  spread_draws(b_Intercept, b_scaleSERIES.l, b_scaleabs.lat, 
               b_scaleElevation, b_interaction_present,
               sigma) %>% 
  summarise_draws()

# maybe we want to know what the overall effect is, ie the global mean as expressed in each region
moddat %>%
  #spread_draws(b_Intercept, r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR,]) %>%
  #median_qi(condition_mean = b_Intercept + r_RESOLVED.TAXA.PAIR) # can change the width of the credible interval
  spread_draws(b_Intercept, b_scaleSERIES.l, b_scaleabs.lat, 
               b_scaleElevation, b_interaction_present) %>%
  median_qi(#tropical_mean = b_Intercept + b_CLIMATE1Tropical,
            #treatmentYES_mean = b_Intercept + b_CLIMATE1Tropical,
            interactionYES_mean = b_Intercept + b_interaction_present) # can change the width of the credible interval

# and plot conditional effects
mod %>%
  #spread_draws(b_Intercept, r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR,]) %>%
  spread_draws(b_Intercept, b_CLIMATE1Tropical) %>%
  median_qi(condition_mean = b_Intercept + b_CLIMATE1Tropical, .width = c(.95, .66)) %>%
  ggplot(aes(y = reorder(b_CLIMATE1, condition_mean), x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, color = "#839496", size = 1) 

# and with density distribution
mod %>%
  spread_draws(b_Intercept, r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR,]) %>%
  mutate(condition_mean = b_Intercept + r_RESOLVED.TAXA.PAIR) %>%
  ggplot(aes(y = reorder(RESOLVED.TAXA.PAIR, condition_mean), x = condition_mean)) +
  geom_vline(xintercept = 0, color = "#839496", size = 1) +
  stat_halfeye(.width = .5, size = 2/3, fill = "#859900")+
  labs(x = expression("Estimate"),
       y = "Associations between taxa pairs") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        text = element_text(family = "Ubuntu")) 

