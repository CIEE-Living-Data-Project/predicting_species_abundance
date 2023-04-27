##Main Data Cleaning Script
##created and updated by: by GLL, AF, EB, MW, AD, IE
##Script compiled and retested by MW and FK



#Filtering BIOtime (collated pairs)
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping
library(progress)

#load data
load("data/tidy/collated_pairs.RData") #full dataset with species/abundances
load('data/prep_biotime/bio_pairs_10km.RData') #"key" with overlapping studies, and reference
#to species in collated.pairs
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use
#subseted<-collated.pairs[1:10000,]

#load functions
source("scripts/data_processing/Functions_data.processing.R")


# adds taxa ID so that easy to ID which genera belong to what type of animal/plant/fungi
collated.pairs$ID <- as.character(collated.pairs$ID)
meta.pairs$STUDY_ID <- as.character(meta.pairs$STUDY_ID)
collated.pairs <- left_join(collated.pairs, meta.pairs[, c(1, 3,4, 12, 13)], by=c("ID"="STUDY_ID"))


#First, look for species with less than three observation across all studies and years
low_sampling <- collated.pairs %>%
  group_by(GENUS_SPECIES) %>% #group by species
  summarize(n=n()) %>% #count the number of each species
  arrange(n)%>% #arrange in ascending order
  filter(n<3) #only get species with less than three observations
#across all studies and years 

low_sampled_species_list <- low_sampling$GENUS_SPECIES #get list of low sampled species
filtered.collated.pairs <- collated.pairs %>% #get original collated pairs dataset
  filter(!GENUS_SPECIES %in% low_sampled_species_list) #filter out observations with rare species

#######
#Standardizing abundances by sampling effort
#######

#Visualizing variation in sampling effort across years within each study:

temp <- filtered.collated.pairs %>% 
  group_by(ID, YEAR) %>% #group by study, year
  summarize(EFFORT.YEAR = n_distinct(PLOT)) #count number of distinct plots

plot <- ggplot(data = temp, aes(x = YEAR, y = EFFORT.YEAR, group = ID, color = ID)) +
  theme_classic() +
  geom_line();plot
#this plots the number of plots per study (colors) across years


#zooming in on studies other than those 4 studies
plot2 <- temp %>% 
  filter(!(ID %in% c(54, 295, 296, 355))) %>% 
  ggplot(aes(x = YEAR, y = EFFORT.YEAR, group = ID, color = ID)) + 
  theme_classic() +
  geom_line();plot2
#lots of variation in sampling effort (plot number) among years here too

# 1. back at the species level, standardize abundances by sum.allrawdata.ABUNDANCE or BIOMASS / EFFORT.YEAR

collated.pairs_standardized <- filtered.collated.pairs %>% 
  group_by(ID, YEAR) %>% #group by study, year
  summarize(EFFORT.YEAR = n_distinct(PLOT)) %>% #count number of distinct plots 
  left_join(filtered.collated.pairs, ., by = c('ID', 'YEAR')) %>% 
  mutate(ST.ABUN = sum.allrawdata.ABUNDANCE/EFFORT.YEAR) %>% 
  mutate(ST.BIO = sum.allrawdata.BIOMASS/EFFORT.YEAR)

# 2. Calculate mean, median, min, max, sd, and CoV for standardized abundances and biomass

## Only genera and study ID kept (not Lat/Long), which standardizes Lat/Long within each study
## Aggregated species to the genus level


collated.pairs_standardized_summary <- collated.pairs_standardized %>% #read in filtered collated pairs
  group_by(ID,YEAR,GENUS,TAXA, ORGANISMS, CLIMATE, REALM) %>% #group by study, location, genus
  summarize(mean_abun_st=mean(ST.ABUN,na.rm=T), #mean abundance
            median_abun_st=median(ST.ABUN,na.rm=T), #median abundance
            min_abun_st=min(ST.ABUN), #min abundance
            max_abun_st=max(ST.ABUN), #max abundance
            sd_abun_st=sd(ST.ABUN,na.rm=T), #sd abundance (when more than one species per genera)
            CoV_abun_st = sd_abun_st/mean_abun_st, # coef of variation
            mean_bio_st=mean(ST.BIO,na.rm=T), #mean biomass
            median_bio_st=median(ST.BIO,na.rm=T), #median biomass
            min_bio_st=min(ST.BIO), #min biomass
            max_bio_st=max(ST.BIO), #max biomass
            sd_bio_st=sd(ST.BIO,na.rm=T), #sd biomass
            CoV_bio_st = sd_bio_st/mean_bio_st) #coef of variation

#Save file as needed using function below
#save(collated.pairs_standardized_summary, file = "data/cleaned_collated_standardized_MSF.Rdata")

######
#Find overlapping consecutive years for all pairs
######

bio.pairs <- read_csv("data/prep_biotime/bio_pairs_10km.csv")


#Run for BETWEEN studies
pb<-set.prog.bar(nrow(bio.pairs)) #sets progress bar
between.studies.overlap=data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", "Max.Overlap","Type","PAIR.ID")))) #makes an empty dataframe
for (h in 1:nrow(bio.pairs)){
  pb$tick()
  data.pair<-calc.overlap.between(h,collated.pairs_standardized_summary)
  between.studies.overlap<-rbind(between.studies.overlap,data.pair)
  
} #populates dataframe with pairs

#Run for WITHIN studies
pb<-set.prog.bar(length(unique(c(bio.pairs$ID.1,bio.pairs$ID.2)))) #sets progress bar
within.studies.overlap=data.frame(matrix(ncol=5,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", "Max.Overlap","Type","PAIR.ID")))) #makes an empty dataframe
for (h in 1:length(unique(c(bio.pairs$ID.1,bio.pairs$ID.2)))){
  #pb$tick()
  id=unique(c(bio.pairs$ID.1,bio.pairs$ID.2))[h]
  data.pair<-calc.overlap.within(id,collated.pairs_standardized_summary)
  within.studies.overlap<-rbind(within.studies.overlap,data.pair)
  
} #populates dataframe with pairs

#Combine
overlap.all.pairs<-rbind(between.studies.overlap,within.studies.overlap)

#Save as needed
#saveRDS(between.studies.overlap,"data/between.studies.overlap.RDS")


######
# Calculate portion change for species pairs
######




pairs.keep=between.studies.overlap[which(between.studies.overlap$Max.Overlap>9),] #this is generated in the Overlap script



results=data.frame(matrix(ncol=14,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", "Log.prop.change.abun.Gn1","Log.prop.change.abun.Gn2","Log.prop.change.bio.Gn1","Log.prop.change.bio.Gn2","PairID","Type","SERIES.n","SERIES.start","SERIES.end","SERIES.l","YEAR.T","YEAR.T1")))) #makes an empty dataframe
pb<-set.prog.bar(nrow(pairs.keep))
for (i in 1:nrow(pairs.keep)){
  pb$tick()
  wk<-get.log.prop.change(i,collated.pairs_standardized_summary,pairs.keep)
  results<-rbind(results,wk)
  
} #get prop change for every genera pair

#add metadata

meta.data<-make.meta(results,collated.pairs_standardized_summary)
log.prop.change.with.meta<-left_join(results,meta.data) #join meta data with results df

#add unique genera ID col
log.prop.change.with.meta$UNIQUE.PAIR.ID=paste(log.prop.change.with.meta$Gn1,log.prop.change.with.meta$Gn2,log.prop.change.with.meta$PairID,sep="_")

#save
saveRDS(log.prop.change.with.meta,"data/log.prop.change.with.meta.RDS")










