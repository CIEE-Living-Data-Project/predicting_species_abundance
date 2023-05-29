#Investigating spatial trends in the bio pairs dataset
# Date created: 26 May 2023


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_


# # LIBRARIES # #
library(tidyverse)
library("dplyr")
library("tibble")
library("readr") 
library("ggplot2")
library("magrittr")
library('rglobi')
library('progress')
library('sf')
library('sp')
library('geosphere')


rm(list=ls()) 


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Step 1. 
#Read in and prepare bio.pairs and collated.pairs


bio.pairs <- read.csv("data/prep_biotime/bio_pairs_10km.csv")
#read in and prepare collated pairs
load("data/tidy/collated_pairs.RData") #collated pairs of overlapping studies

#Set collated pairs as a spatial points data frame for investigation
collated.pairs.spatial <- collated.pairs
xy <- collated.pairs.spatial%>%
  select(LONGITUDE, LATITUDE)
coordinates(collated.pairs.spatial) <- xy
proj4string(collated.pairs.spatial) <- CRS("+proj=longlat +datum=WGS84")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#Step 2. Calculate distances between observations within each study

dummy_data <- collated.pairs %>%
  dplyr::filter(ID=='122')
xy <- dummy_data %>%
  dplyr::select(LONGITUDE, LATITUDE)
coordinates(dummy_data) <- xy
proj4string(dummy_data) <- CRS("+proj=longlat +datum=WGS84")



# Create an empty data frame to store the results
distances <- data.frame()

#Make Isaac's progress bar 
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

pb<-set.prog.bar(length(unique(collated.pairs.spatial$ID))) #sets progress bar
# Iterate over unique study IDs
for (study_id in unique(collated.pairs.spatial$ID)) {
  
  pb$tick()

  # Subset the data frame for the current study ID
  study_data <- collated.pairs.spatial[collated.pairs.spatial$ID == study_id, ]
 study_sf <- st_as_sf(study_data)
  
  # Calculate the average latitude and longitude for the study
  avg_latitude <- mean(study_data$LATITUDE)
  avg_longitude <- mean(study_data$LONGITUDE)
  avg_point <- data.frame(avg_longitude, avg_latitude)
  coordinates(avg_point) <- avg_point
  proj4string(avg_point) <- CRS("+proj=longlat +datum=WGS84")
  avg_sf <- st_as_sf(avg_point)
  
  

  # Calculate the distance to the average point
  distance_to_avg <- st_distance(study_sf, avg_sf) 
  units(distance_to_avg) <- c('km')
  distance_df <- as.data.frame(distance_to_avg)
  distance_df$distance_to_avg <- as.numeric(distance_to_avg)
  
  # Add the study ID as a column in the distance data frame
  distance_df$ID <- study_id
  
  # Append the distances to the 'distances' data frame
  distances <- rbind(distances, distance_df)
}

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Step 3. Interpret the results

#Get some summary stats
mean(distances$distance_to_avg)
#592km
max(distances$distance_to_avg)
#15408 km!
#Total number of studies: 
length(unique(distances$ID))

#(min excluding zeroes)
no_zeroes <- distances %>%
  dplyr::filter(!distance_to_avg==0)
min(no_zeroes$distance_to_avg)
#0.01km, or 10m
#Number of studies with non-zero distances: 
length(unique(no_zeroes$ID))


#Average distances from average across studies
avg_dist_within_study <- distances %>%
  group_by(ID) %>%
  summarize(avg_dist = mean(distance_to_avg), 
            sd = sd(distance_to_avg))


#Pull in realm from meta.pairs to investigate in box plot 
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use
meta.pairs.realm <- meta.pairs %>%
  dplyr::select(STUDY_ID, REALM) 
colnames(meta.pairs.realm) <- c("ID", "REALM")
merged_df <- merge(no_zeroes, meta.pairs.realm, by = "ID")



#Plot the data - not including ones with zero variance 
box_plot_dist <-  merged_df %>%
  ggplot(aes(x=as.character(ID), y=distance_to_avg, group=as.character(ID)))+
  geom_boxplot(aes(colour=as.character(REALM))) +
  xlab("Study ID") +
  ylab("Distance to Average Lat/Lon of Study") +
  labs(colour="Realm")+
  theme_classic()
box_plot_dist

