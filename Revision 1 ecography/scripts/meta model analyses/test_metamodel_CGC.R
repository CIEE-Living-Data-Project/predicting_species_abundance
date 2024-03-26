# script written by CC and GLL March 2024
# calculates correlations and z-scores for genera pairs
# fits metamodel that tests how 4 ecological hypotheses drive relationships between genera-pair
library(tidyverse)
library(brms)

# requires data to be saved locally
# 1. log prop changes
# 2. interactions and taxonomic data
# 3. BioTime metadata: BioTIMEMetadata_24_06_2021.csv

# reading in data and cleaning ####
# interactions and resolved taxa pair column
# update to include corrected genera names from Emily
interactions.dat <- readRDS("Revision 1 ecography/output/prep_data/results_abundance_interactions_taxa_032024ENB.RDS") 

# full dataset log prop change dataset
# update to include file with corrected genera names from Emily
abun.dat <- read.csv("Revision 1 ecography/output/prep_data/results.abundance.csv")
abun.dat$X.1<-NULL
#str(abun.dat)
#names(abun.dat)
#names(interactions.dat)

#pull out only new information not to duplicate on abun.dat
abun.dat_trim<-abun.dat[, c(2,3, 23:30)]%>%distinct(.)
#str(abun.dat_trim)
abun.dat_trim<-distinct_at(abun.dat_trim, "TS_ID", .keep_all = T)
rm(abun.dat)

# add a STUDY_ID column so can join to metadata
#abun.dat$STUDY_ID <- as.character(sapply(strsplit(abun.dat$STUDY_PLOT, "~"), "[", 1))
interactions.dat$STUDY_ID <- as.character(sapply(strsplit(interactions.dat$STUDY_PLOT, "~"), "[", 1))

# meta data
all_meta <- read.csv("data/prep_biotime/BioTIMEMetadata_24_06_2021.csv")
all_meta <- select(all_meta, STUDY_ID, REALM, CLIMATE, HABITAT, ORGANISMS, CENT_LAT, CENT_LONG, NUMBER_OF_SPECIES)
all_meta$STUDY_ID <- as.character(all_meta$STUDY_ID)
str(all_meta)

# clean up organisms column
all_meta$ORGANISMS <- gsub("birds|Bird|Birds", "Birds", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub('breeding Birds|breeding bird pairs|Breeding Birds', "Breeding Birds", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub('Plants|plants', "Plants", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub('Small mammal|Small mammals|small mammals', "Small mammals", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub('Insects specifically Coleoptera and Lepidoptera|small mammals|Beetles|Butterflies', "Beetles, Butterflies, Moths", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub("Lizards", "Herpetofauna", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub("rodents", "Small mammals", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub("waterBirds|ducks", "Water birds", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub("Acrididae (grasshoppers)", "Grasshoppers", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub("Acrididae \\(", "", all_meta$ORGANISMS)
all_meta$ORGANISMS <- gsub("\\)", "", all_meta$ORGANISMS)

# clean up habitat column
all_meta$HABITAT <- gsub("Chalk grassland|Grassland", "Grassland", all_meta$HABITAT)
all_meta$HABITAT <- gsub("Boreal forest|Boreal vegetation / Taiga|Scandinavian taiga", "Boreal", all_meta$HABITAT)
all_meta$HABITAT <- gsub("Tallgrass prairie|Tallgrass prairie gallery forest and riparian edge|Savanna/ Tallgrass prairie|Woodland", "Tallgrass prairie and woodland", all_meta$HABITAT)
all_meta$HABITAT <- gsub("Praire", "Prairie", all_meta$HABITAT)
all_meta$HABITAT <- gsub("upper montane forests above 2700 feet in the North", "Alpine", all_meta$HABITAT)
all_meta$HABITAT <- gsub("Desert Wildlife Refuge|Desert/ Grassland|semiarid thorn scrub/Desert/ Grassland|creosotebush site and grassland site|semiarid thorn scrub", "Desert, grassland, thorn scrub, creosotebush", all_meta$HABITAT)
all_meta$HABITAT <- gsub("Urban / Desert|Urban / Desert / Riparian / Agricultural|large pine plantatino in Nacogdeoches County. Texa|forests agricultural fields and meadows", "Human modified", all_meta$HABITAT)
all_meta$HABITAT <- gsub("Northern mixed prairie|Mixed|Forest and grassland", "Mixed prairie and forest", all_meta$HABITAT)
all_meta$HABITAT <- gsub("Mixed prairie and forest Conifer-Hardwood Forest", "Mixed Conifer-Hardwood Forest", all_meta$HABITAT)
all_meta$HABITAT <- gsub("birch forest", "Deciduous forest", all_meta$HABITAT)


## combine datasets ####
#interactions.dat_trim <- distinct(select(interactions.dat, c("TS_ID", "STUDY_PLOT", "Gn1", "Gn2", "SERIES.l", 
#                                                             "interaction_present", "RESOLVED.TAXA1", "RESOLVED.TAXA2", "resolved_taxa_pair")))
#alldat <- left_join(abun.dat_trim, interactions.dat)
alldat <- left_join(interactions.dat, all_meta)
names(alldat) #look at all column names 

# check to make sure that every pair has interaction info
sum(is.na(alldat$interaction_present))
#intcheck<-subset(alldat, is.na(interaction_present))#1984 obs missing info
#checkIDs<-unique(intcheck$TS_ID) #108 IDs 
#which(is.na(alldat$interaction_present)==TRUE)
#alldat[7010172, ] #??

# calculate abs latitude metric
alldat$abs.lat <- abs(alldat$CENT_LAT)

### add in worldclim data ####
worldclim <- read.csv("Revision 1 ecography/output/prep_data/worldclim.csv")[,-c(1,4,5)]
alldat <- left_join(alldat, worldclim, by=c("CENT_LAT"="LATITUDE", 
                                            "CENT_LONG"="LONGITUDE"))

###  integrate disturbance data ####

# section that creates csv file to populate with cleaned disturbance data,
# ignore during modelling work
# read in metadata from entire study
#all_meta <- read.csv("Revision 1 ecography/raw biotime data/BioTIMEMetadata_24_06_2021.csv")

#meta.info <- distinct(left_join(select(alldat, STUDY_ID, STUDY_PLOT, LATITUDE, LONGITUDE), all_meta))

# add in column for cleaned disturbance from Sam
#load("data/prep_biotime/meta_pairs_10km.RData")
#test <- left_join(meta.info, meta.pairs)
# write.csv(test, "meta.info.csv")

# download final disturbance data from google drive when cleaned
# and bind back to alldat
dist<-read.csv("Revision 1 ecography/output/prep_data/disturbance cleaning.csv")
dist$STUDY_ID<-as.character(dist$STUDY_ID)
alldat<-left_join(alldat, dist)

# prep data for model ####
# calculate pearson's cor for all time series 
alldat <- alldat %>% 
  group_by(TS_ID) %>%
  mutate(cor=cor(Log.prop.change.Gn1, Log.prop.change.Gn2)) 
# warnings caused by time series length where 0 changes through entire time series


# subset to only relatively strong cors - any higher is unrealistic
alldat_trim <- subset(alldat, cor<0.8 & cor>-0.8) 

hist(alldat_trim$cor)
max(alldat_trim$cor)
min(alldat_trim$cor)

#add in richness info from abundance data csv
alldat_trim<-left_join(alldat_trim, abun.dat_trim)#7436718
str(alldat_trim)

#fix studies missing richness data  
#alldat_trimxx<-subset(alldat_trim, is.na(Rich.T.Gn1))#843 with missing info 
alldat_trim<-subset(alldat_trim, !is.na(Rich.T.Gn1))#remove them 

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

# note throws warnings for all SE where sample size is less than 3, but that is expected and okay
# time series length is only one of these for which all n > 3
# note only use the estimates of total sp and individuals to estimate sample size
# differences don't make sense as a sample size
# takes a few min
# eq 3.12 in meta-analysis book 
alldat_trim <- alldat_trim %>% 
  mutate(SE.timeseries = (1/sqrt(SERIES.l-3)), # time series length
         SE.total.indivs = (1/sqrt(total.indivs-3)), # total number of individuals
         SE.total.sp = (1/sqrt(total.sp-3)), # total number of species
         var.total.indivs = (1/(total.indivs-3)), # total number of individuals
         var.total.sp = (1/(total.sp-3)), # total number of species
        ) 

# explore how many studies have very different numbers of sp or individuals between time series
# include plot in supplement
# things with higher error are more weakly correlated 
# takes 10 minutes to render
#par(mfrow=c(1,2), oma=c(3,3,2,0), mar=c(2,4,2,4))
#plot(alldat_trim$cor, alldat_trim$abs.total.spGn1mGn2,
#     ylab="", xlab="", 
#     cex.axis=1.5, cex.lab=1.5, 
#     col=ifelse(alldat_trim$cor<.2 & alldat_trim$cor>-.2 & alldat_trim$abs.total.spGn1mGn2>=175 , "dodgerblue2", 
#                ifelse(alldat_trim$cor<.2 & alldat_trim$cor>-.2 & alldat_trim$abs.total.spGn1mGn2<175, "darkblue",
#                       "grey40")))
#mtext("Difference in species richness \nbetween genera", 2, line=3, cex=1.5)
#mtext("a.", 3, line=1, cex=1.5, at=-1)

#plot(alldat_trim$cor, alldat_trim$abs.total.indivsGn1mGn2, 
#     ylab="", xlab="", cex.axis=1.5,
#     col=ifelse(alldat_trim$cor<.2 & alldat_trim$cor>-.2 & alldat_trim$abs.total.indivsGn1mGn2>=10000 , "dodgerblue2", 
#                ifelse(alldat_trim$cor<.2 & alldat_trim$cor >-.2 & alldat_trim$abs.total.indivsGn1mGn2<10000, "darkblue",
#                       "grey40")))
#mtext("Difference in abundance \nbetween genera", 2, line=3, cex=1.5)
#mtext("b.", 3, line=1, cex=1.5, at=-1)
#mtext("Correlations between genera", 1, line=1, cex=1.5, outer=T)

# filter df to remove repeat values 
# ie keeping unique cases of z-scores

moddat <- alldat_trim %>% 
  select(-Log.prop.change.Gn1, -Log.prop.change.Gn2, -Abs.change.Gn1, -Abun.T1.Gn1, -Abun.T.Gn2, -Abun.T1.Gn2, 
         -Abs.change.Gn2, -Abun.T.Gn1, -Length.Unique.Values.Gn1, -Length.Unique.Values.Gn2,
         -Rich.T.Gn1, -Rich.T1.Gn1, -Rich.T.Gn2, -Rich.T1.Gn2, -SS.T.Gn1, -SS.T1.Gn1, -SS.T.Gn2, -SS.T1.Gn2,
         -X, -YEAR.T, -YEAR.T1) %>%
         #-abs.total.indivsGn1mGn2, -abs.total.spGn1mGn2, -SE.abs.total.indivsGn1mGn2, -SE.abs.total.spGn1mGn2
  distinct(.)

#save so don't need to re-run 
#save(moddat, file="Revision 1 ecography/output/prep_data/all_model_data.Rdata")

#load model data ####
#load(file="Revision 1 ecography/output/prep_data/all_model_data.Rdata")
# look at some to check   
select(moddat, cor, z, SE.timeseries)
hist(moddat$cor)

# explore number of rows with various SE metric
# drop the time series length one bc doesn't make sense to also use it as a question we are testing
# trim all data to most restrictive of these SE calculations
# compare models with goodness of fit stats
length(na.omit(moddat$SE.total.indivs)) # 411,351
length(na.omit(moddat$SE.total.sp)) #  411,400

# if want to run all SE calculations on max restrictive data set then:
# remove any NAs
#length(unique(na.omit(moddat)$STUDY_ID)) # 62 studies
#length(unique(na.omit(moddat)$TS_ID)) #  374557 time series
#moddat <- na.omit(moddat)
moddat<-subset(moddat, !is.na(SE.total.indivs))

# make sure interaction_present is encoded as a factor
# and scale inputs here bc doesn't work in model formula 
moddat$scale.abs.lat <- scale(moddat$abs.lat)
moddat$scale.SERIES.l <- scale(moddat$SERIES.l)  
moddat$interaction_present.factor <- as.factor(moddat$interaction_present)
moddat$scale.elev <- scale(moddat$Elevation)

unique(moddat$treatment_yn_clean)  

moddat$treatment_yn_clean[moddat$treatment_yn_clean=='No'|moddat$treatment_yn_clean=='no']="no"
moddat$treatment_yn_clean[moddat$treatment_yn_clean=='Yes'|moddat$treatment_yn_clean=='yes']="yes"

# run models on subset of data that doesn't have big differences between indivs or sp
# between time series to see if our results are robust
# note that we expect (and can see in the data) that big differences to weaken correlations and add noise, 
# thus including them is a conservative approach
# subset data below like 1000 indivs diff and to less than 30 spp diffs, breaks based on diff x corr plots above

hist(moddat$abs.total.indivsGn1mGn2)
hist(moddat$abs.total.spGn1mGn2)

moddat<-subset(moddat, abs.total.indivsGn1mGn2< 1000)
moddat<-subset(moddat, abs.total.spGn1mGn2 < 30)
unique(moddat$STUDY_ID)

#save AGAIN 
save(moddat, file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")

# modelling timmmmme #####
## set priors ####
#priors <-c(prior(normal(0,0.33), class = Intercept), # set between -1 and 1 for z score
#          prior(normal(0,0.33), class = sd, lb=0)) # set lower bound 0 for SE, values b/w (0,1)

load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")

## define models ####
# model with z scores and total indivs SE as joint response
MODFORM.indivs <- bf(z|resp_se(SE.total.indivs, sigma = FALSE) ~ 
                scale.SERIES.l + treatment_yn_clean + 
                resolved_taxa_pair + 
                scale.abs.lat +
                interaction_present.factor +
                scale.elev +
                (1|STUDY_ID))
# model with z scores and total indivs SE as joint response
MODFORM.sp <- bf(z|resp_se(SE.total.sp, sigma = FALSE) ~ 
                       scale.SERIES.l + treatment_yn_clean + 
                       resolved_taxa_pair + 
                       scale.abs.lat +
                       interaction_present.factor +
                       scale.elev +
                       (1|STUDY_ID))

## run models ####

# model 1
start_time <- Sys.time() 
metamod.indivs <- brm(MODFORM.indivs, moddat, cores=3, chains=3, 
             iter=5000, family=gaussian, file="Revision 1 ecography/output/meta model/SEindivs2.rmd")
end_time <- Sys.time()
end_time - start_time
#save as Rdata file 
#save(metamod.indivs, file="Revision 1 ecography/output/meta model/metamod.individual.BAD.Rdata")
summary(metamod.indivs)
# model 2
start_time <- Sys.time() 
metamod.sp <- brm(MODFORM.sp, moddat,
                      control = list(adapt_delta=0.8, max_treedepth = 12), cores=3, chains=3, 
                      iter=5000, family=gaussian, #prior = priors, 
                      file="Revision 1 ecography/output/meta model/SEspecies.rmd")
end_time <- Sys.time()
end_time - start_time



## model checks ####
#sample size (Bulk_ESS & Tail_ESS) should > 1000 & rhat < 1.1
summary(metamod)

plot(metamod.indivs) #model convergence (L: does distribution mean match estimate? R: did all values get explored?)

# posterior predictive checks - are predicted values similar to posterior distribution?
# takes a long time to run
pp_check(metamod, ndraws = 50) 

# Pairs plots to diagnose sampling problems (should show Gaussian blobs)
pairs(metamod)

plot(conditional_effects(metamod)) #fitted parameters and their CI

# get coefs for summary tables using fixef(mod)
# this reports mean of posterior
fixef(metamod)

#run subset models####

# note, if models run quickly and we have lots of time, 
