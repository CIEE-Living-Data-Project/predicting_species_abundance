#Filtering BIOtime (collated pairs)
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping
library("progress")

#load data
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/tidy/collated_pairs.RData")#collated pairs of overlapping studies
load("~/Library/CloudStorage/OneDrive-McGillUniversity/R Scripts/GitHub/predicting_species_abundance/data/prep_biotime/bio_pairs_10km.RData") #metadata of overlapping studies

cleaned_collated_pairs_EBMW <- read_csv("data/cleaned_collated_pairs_EBMW.csv")

## calc overlap for each pair in each ID
calc.overlap.pair<-function(p,data){

  require(progress)
  require(dplyr)
pair=bio.pairs[p,]

#get timeseries
timeseries_1 <- data[which(data$ID==pair$ID.1)]
timeseries_2 <- data[which(data$ID==pair$ID.2)]

#get genera lists
genera1 <- timeseries_1$GENUS %>% unique()
genera1 = genera1[!is.na(genera1)]

genera2 <- timeseries_2$GENUS %>% unique()
genera2 = genera2[!is.na(genera2)]

if (length(genera1)==0){return(NULL)}
if (length(genera2)==0){return(NULL)}


#print(paste("Comparing study",pair$ID.1,"and study",pair$ID.2))
#print(paste("NUMBER OF GENERA in",pair$ID.1,"=",length(genera1)))
#print(paste("NUMBER OF GENERA in",pair$ID.2,"=",length(genera2)))

#make progress bar
#n_iter <- length(genera1)
#pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",total = n_iter,complete = "=",   # Completion bar character
#                       incomplete = "-", # Incomplete bar character
#                       current = ">",    # Current bar character
#                       clear = FALSE,    # If TRUE, clears the bar when finish
#                       width = 100)      # Width of the progress bar


#start with between studies
mat<-data.frame("Gn1"=NA,"Gn2"=NA,"Max.Overlap"=NA)
mat<-mat[-1,]
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
  
  

  
}
mat$Type="Between"
mat$PAIR.ID=paste0(pair$ID.1,"_",pair$ID.2)

return(mat)
}
calc.overlap.single<-function(ID,data){

  #get timeseries
  timeseries_1 <- data[which(data$ID==ID),]

  #get genera lists
  genera1 <- timeseries_1$GENUS %>% unique()
  genera1 = genera1[!is.na(genera1)]
  
  #print(paste("NUMBER OF GENERA in",ID,"=",length(genera1)))
  
  if (length(genera1)==0){return(NULL)}

  #make progress bar
  #n_iter <- length(genera1)
  #pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
   #                      total = n_iter,
    #                     complete = "=",   # Completion bar character
     #                    incomplete = "-", # Incomplete bar character
      #                   current = ">",    # Current bar character
       #                  clear = FALSE,    # If TRUE, clears the bar when finish
        #                 width = 100)      # Width of the progress bar
  
  
  
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
  mat.id1$PAIR.ID=paste0(pair$ID.1)
  mat.id1<-mat.id1[-which(mat.id1$Gn1==mat.id1$Gn2),]
  
  return(mat.id1)
}
set.prog.bar<-function(n_iter){
  #make progress bar
  n_iter <- n_iter
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = n_iter,
                         complete = "=",
                         incomplete = "-",
                         current = ">",
                         clear = FALSE,
                         width = 100)
  
}

set.prog.bar(nrow(bio.pairs))
results=data.frame("Gn1"=NA,"Gn2"=NA,"Max.Overlap"=NA,"Type"=NA,"PAIR.ID"=NA)
results<-results[-1,]
for (h in 1:nrow(bio.pairs)){
  pb$tick()
  data.pair<-calc.overlap.pair(h,cleaned_collated_pairs_EBMW)
  results<-rbind(results,data.pair)
  
}

#now do within study overlapps
set.prog.bar(length(unique(c(bio.pairs$ID.1,bio.pairs$ID.2))))
results.single=data.frame("Gn1"=NA,"Gn2"=NA,"Max.Overlap"=NA,"Type"=NA,"PAIR.ID"=NA)
results.single<-results.single[-1,]
for (h in 1:length(unique(c(bio.pairs$ID.1,bio.pairs$ID.2)))){
  #pb$tick()
  h=1
  id=unique(c(bio.pairs$ID.1,bio.pairs$ID.2))[h]
  data.pair<-calc.overlap.single(id,cleaned_collated_pairs_EBMW)
  results.single<-rbind(results.single,data.pair)
  
}




length(which(results$Max.Overlap>5))

nrow(results.single)




full<-rbind(results,results.single)
saveRDS(full,"/Users/isaaceckert/Desktop/full.RDS")

head(full)


