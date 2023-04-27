#Title: Script to find overlapping consecutive years for all pairs
#Author: Isaac Eckert
#Date: April 26th, 2023

#Libraries
library(dplyr) 
library(tibble) 
library(progress)
library(stringr)

#Load data
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/preprocessing/cleaned_collated_standardized_MSF.Rdata")
bio.pairs <- read.csv("data/prep_biotime/bio_pairs_10km.csv")

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
  mat.id1$PAIR.ID=paste0(ID,"_",ID)
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

#Save
saveRDS(between.studies.overlap,"data/preprocessing/between.studies.overlap.RDS")
saveRDS(within.studies.overlap,"data/preprocessing/within.studies.overlap.RDS")

