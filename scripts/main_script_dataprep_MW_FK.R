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

# 2. re-calculate mean, median, min, max, sd, and CoV for standardized abundances and biomass

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
#Script to find overlapping consecutive years for all pairs
######

bio.pairs <- read_csv("data/prep_biotime/bio_pairs_10km.csv")

#Load some functions
calc.overlap.between<-function(p,data){
  require(progress)
  require(dplyr)
  
  pair=bio.pairs[p,] #isolate which study pair we are focusing on
  
  
  #get timeseries
  timeseries_1 <- data[which(data$ID==pair$ID.1),] #get all timeseries for ID1
  timeseries_2 <- data[which(data$ID==pair$ID.2),] #get all timeseries for ID2
  
  #get genera lists
  genera1 <- timeseries_1$GENUS %>% unique()
  genera1 = genera1[!is.na(genera1)]
  
  genera2 <- timeseries_2$GENUS %>% unique()
  genera2 = genera2[!is.na(genera2)]
  
  if (length(genera1)==0){return(NULL)} #if no species found exit the function
  if (length(genera2)==0){return(NULL)} #if no species found exit the function
  
  #calc overlap for each genera-genera pairing
  mat<-data.frame("Gn1"=NA,"Gn2"=NA,"Max.Overlap"=NA) #make blank df
  mat<-mat[-1,] #remove 1st row
  for (i in 1:length(genera1)){
    gen1=genera1[i]
    #pb$tick()
    
    for (x in 1:length(genera2)){
      gen2=genera2[x]
      
      #get overlapping years
      years1=timeseries_1$YEAR[which(timeseries_1$GENUS==gen1)]
      years2=timeseries_2$YEAR[which(timeseries_2$GENUS==gen2)]
      
      cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
      
      dat<-data.frame("Period"=NA,"value"=NA)
      for (y in 1:length(cont.overlapping)){
        dat[y,1]=y
        dat[y,2]=length(cont.overlapping[[y]])}
      
      #take max
      mat<-rbind(mat,data.frame("Gn1"=gen1,"Gn2"=gen2,"Max.Overlap"=max(dat$value)))
      
    }
    
  } #for loop to look through genera
  mat$Type="Between" #add some stuff
  mat$PAIR.ID=paste0(pair$ID.1,"_",pair$ID.2) #add pair IDs
  
  return(mat) #exit the function
} #function to calc overlap between pairs BETWEEN two studies
calc.overlap.within<-function(ID,data){
  
  #get timeseries
  timeseries_1 <- data[which(data$ID==ID),]
  
  #get genera lists
  genera1 <- timeseries_1$GENUS %>% unique()
  genera1 = genera1[!is.na(genera1)]
  
  if (length(genera1)==0){return(NULL)}
  
  #then do ID1
  mat.id1<-data.frame("Gn1"=NA,"Gn2"=NA,"Max.Overlap"=NA)
  mat.id1<-mat.id1[-1,]
  for (i in 1:length(genera1)){
    #pb$tick()
    
    gen1=genera1[i]
    
    for (x in 1:length(genera1)){
      gen2=genera1[x]
      
      #get overlapping years
      years1=timeseries_1$YEAR[which(timeseries_1$GENUS==gen1)]
      years2=timeseries_1$YEAR[which(timeseries_1$GENUS==gen2)]
      
      cont.overlapping<-split(years1[which(years1%in%years2)], cumsum(c(1, diff(years1[which(years1%in%years2)]) != 1)))
      
      dat<-data.frame("Period"=NA,"value"=NA)
      for (y in 1:length(cont.overlapping)){
        dat[y,1]=y
        dat[y,2]=length(cont.overlapping[[y]])}
      
      #take max
      mat.id1<-rbind(mat.id1,data.frame("Gn1"=gen1,"Gn2"=gen2,"Max.Overlap"=max(dat$value)))
      
      
    }
    
    
    #print(paste("DONE:",i,"out of",length(genera1)))
    
    
  }
  mat.id1$Type="Within"
  mat.id1$PAIR.ID=ID
  mat.id1<-mat.id1[-which(mat.id1$Gn1==mat.id1$Gn2),]
  
  return(mat.id1)
} #function to calc overlap between pairs WITHIN a single study
set.prog.bar<-function(n_iter){
  #make progress bar
  progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                   total = n_iter,
                   complete = "=",
                   incomplete = "-",
                   current = ">",
                   clear = FALSE,
                   width = 100)
  
} #function to make a progress bar

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




