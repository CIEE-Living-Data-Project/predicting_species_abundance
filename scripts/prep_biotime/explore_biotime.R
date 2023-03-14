# Aims:
# 1. Import raw BioTIME data
# 2. Identify proximal pairs
# 3. Map pairs

# Author: Nathalie Chardon
# Date created: 28 Nov 2022
# Date updated: 15 Dec 2022 (NC)


# # LIBRARIES # #
library(tidyverse)
library(dplyr)


rm(list=ls()) 


# # WORKING DIRECTORIES # #
#raw_dat <- '~/Desktop/Code/predicting_species_abundance_offline/' #WD for NC (biotime query 1.2 GB)
#figs <- '~/Desktop/Code/predicting_species_abundance/figures/explore_biotime_Dec2022/'
#tidy_dat <- '~/Desktop/Code/predicting_species_abundance/data/prep_biotime/'

# # INPUT FILES # #
#setwd(raw_dat)
# biotime.raw <- read.csv('BioTIMEQuery_24_06_2021.csv')
biotimeMeta <- read.csv('data/prep_biotime/BioTIMEMetadata_24_06_2021.csv') #saved in data folder of github repo on local machine 
  

# # OUTPUT FILES # #
setwd(tidy_dat)
#load('data/prep_biotime/bio_pairs_10km.RData') #unique 10 km geographic-overlap years-taxa pairs (explore_biotime.R)
#load('data/prep_biotime/bio_pairs_1km.RData') #unique 1 km geographic-overlap years-taxa pairs (explore_biotime.R)
#load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use in lit review (explore_biotime.R)




####################################################################################################

# # EXPLORE BIOTIME METADATA # # 

####################################################################################################

library(ggplot2)

# Download date: 28 Nov 2022
# Download URL: https://biotime.st-andrews.ac.uk/getFullDownload.php
# R walkthrough: https://biotime.st-andrews.ac.uk/downloads/interactBioTIME.html#using-the-csv-file-from-the-query


# Meta data
dim(biotimeMeta) #381 studies


# Plot data

# draw a basic world map, add "y" or "n" for display of tropics and polar latitudes

drawWorld<-function(lats) {
  world_map<-map_data("world")
  
  g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
  g1<-g1+geom_polygon(data=world_map, aes(x=long, y=lat, group=group), colour="gray60", fill="gray60")
  g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
               panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
               legend.position="none",axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
  
  if(lats=="y") {
    g1<-g1+geom_hline(yintercept=23.5, colour="red")+geom_hline(yintercept =-23.5, colour="red")
    g1<-g1+geom_hline(yintercept=66.5, colour="darkblue")+geom_hline(yintercept =-66.5, colour="darkblue")
  }
  else { return(g1) }
  return(g1)
}

# colors per taxa
taxaCol<-c('#ffffff','#ffffbf','#5e4fa2','#f46d43','#3288bd','#abdda4','#a8c614','#d53e4f','#66c2a5','#e6f598','#fee08b','#9e0142','#fdae61', '#fdae62')

# create plot
points<-drawWorld("y")+geom_point(data=biotimeMeta, aes(x=CENT_LONG, y=CENT_LAT, colour=TAXA, size=TOTAL), alpha=I(0.7))
points<-points+scale_colour_manual(values=taxaCol)+scale_size(range=c(3, 10))
points




####################################################################################################

# # IDENTIFY GEOGRAPHICALLY PROXIMAL PAIRS (1 & 10 km) # # 
# source: https://stackoverflow.com/questions/20982635/identify-points-within-specified-distance-in-r

####################################################################################################

library(fossil)                                     # for earth.dist(...)
library(data.table)


# Set up lat, long data
df <- data.frame(long = biotimeMeta$CENT_LONG, lat = biotimeMeta$CENT_LAT) #only lat, long
df.id <- data.frame(long = biotimeMeta$CENT_LONG, lat = biotimeMeta$CENT_LA, #with ID information
                    row.name = rownames(biotimeMeta), study.ID = biotimeMeta$STUDY_ID)

# Set separate distance (RUN THIS SECTION TWICE, ONCE FOR 1KM AND ONCE FOR 10 KM)
# sep.km   <- 1                                      # critical separation (1 km)
sep.km   <- 10                                      # critical separation (10 km)
m        <- as.matrix(earth.dist(df))               # distance matrix in km

# Identify pairs
coloc    <- data.table(which(m < sep.km, arr.ind=T))  # pairs of observations within critical separation distance

# Edit dataframe
setnames(coloc,c("row","col"),c("Obs.1","Obs.2"))     # rename columns to reflect observation row
coloc    <- coloc[Obs.1<Obs.2,]                       # want only lower triangular part
coloc[,dist:=m[Obs.1,Obs.2],by="Obs.1,Obs.2"] # append distances in km

# Set up data table
observations <- data.table(id=as.integer(rownames(df)),df)
setkey(observations,id)
setkey(coloc,Obs.1)
coloc[observations,c("long.1","lat.1"):=list(long,lat)]
setkey(coloc,Obs.2)
coloc[observations,c("long.2","lat.2"):=list(long,lat)]

dim(coloc) #127 study pairs are within 1 km of each other
dim(coloc) #367 study pairs are within 10 km of each other

# Plot pairs
points<-drawWorld("y")+geom_point(data=coloc, aes(x=long.1, y=lat.1), alpha=I(0.7))
points<-points+scale_colour_manual(values=taxaCol)+scale_size(range=c(3, 10))
points


# Match row numbers to observation ID
coloc$ID.1 <- NA
coloc$ID.2 <- NA

for (i in 1:nrow(coloc)) {
  
  foo <- df.id[df.id$row.name == coloc$Obs.1[i], ] #identify Obs.1 in df.id by rowname
  coloc$ID.1[i] <- foo$study.ID #add study ID to coloc
  
  boo <- df.id[df.id$row.name == coloc$Obs.2[i], ] #identify Obs.2 in df.id by rowname
  coloc$ID.2[i] <- boo$study.ID #add study ID to coloc
}

head(coloc)




####################################################################################################

# # ORGANIZE DATAFRAME # # 

####################################################################################################

# Keep only relevant columns in biotimeMeta
bio <- biotimeMeta %>% select(c('HABITAT', 'TAXA', 'ORGANISMS', 'START_YEAR', 'END_YEAR', 
                                'NUMBER_OF_SPECIES', 'STUDY_ID'))

bio <- bio %>% rename(ID.1 = STUDY_ID) #column names need to match for join
head(bio)

# Match data from ID.1 to coloc
bio.pairs <- left_join(coloc, bio, by = 'ID.1') #join ID.1 info

bio.pairs <- bio.pairs %>% rename('habitat.1' = 'HABITAT', 'taxa.1' = 'TAXA', 'organisms.1' = 'ORGANISMS',
                      'start_year.1' = 'START_YEAR', 'end_year.1' = 'END_YEAR', 
                      'no_species.1' = 'NUMBER_OF_SPECIES') #identify to ID.1
head(bio.pairs)

# Match data from ID.2 to coloc
bio <- bio %>% rename(ID.2 = ID.1) #column names need to match for join
bio.pairs <- left_join(bio.pairs, bio, by = 'ID.2')  #join ID.2 info

bio.pairs <- bio.pairs %>% rename('habitat.2' = 'HABITAT', 'taxa.2' = 'TAXA', 'organisms.2' = 'ORGANISMS',
                                  'start_year.2' = 'START_YEAR', 'end_year.2' = 'END_YEAR', 
                                  'no_species.2' = 'NUMBER_OF_SPECIES') #identify to ID.2

head(bio.pairs)

# Delete irrelevant columns and arrange remaining
bio.pairs <- bio.pairs %>% select(dist, ID.1, ID.2, long.1, long.2, lat.1, lat.2, habitat.1, habitat.2,
                                  taxa.1, taxa.2, organisms.1, organisms.2, no_species.1, no_species.2,
                                  start_year.1, start_year.2, end_year.1, end_year.2)

head(bio.pairs)




####################################################################################################

# # IDENTIFY TEMPORALLY PROXIMAL PAIRS # # 
# source: https://search.r-project.org/CRAN/refmans/DescTools/html/overlaps.html

####################################################################################################

library(DescTools)


# Dataframes with only year info
df1 <- bio.pairs %>% select(start_year.1, end_year.1)
df2 <- bio.pairs %>% select(start_year.2, end_year.2)

head(df1)
head(df2)

# Identify overlap amount
vv.over <- Overlap(df1, df2) #amount of overlap
head(vv.over)

# Add overlap vector to bio.pairs
bio.pairs$overlap.years <- vv.over 

# Keep only overlapping years
bio.pairs <- bio.pairs %>% filter(overlap.years > 0)
dim(bio.pairs)
summary(bio.pairs$overlap.years) #1-51 years overlap, median 11 years, mean 14 years




####################################################################################################

# # ORGANIZE & VISUALIZE DATA # # 

####################################################################################################

# Keep only unique taxa pairs
bio.pairs$unique.taxa <- bio.pairs$taxa.1 == bio.pairs$taxa.2 #identify if taxa match
bio.pairs <- bio.pairs %>% filter(unique.taxa == FALSE) #keep only unique pairs
bio.pairs <- select(bio.pairs, -c(unique.taxa)) #remove unique column

dim(bio.pairs) #68 pairs (1 km); 180 unique pairs (10 km)
head(bio.pairs)
tail(bio.pairs)

# Create unique taxa pair column
bio.pairs$taxa.pairs <- paste(bio.pairs$taxa.1, bio.pairs$taxa.2, sep = '-')


# Plot matching pairs (colored by taxa.pairs and sized by overlap.years)
setwd(figs)
points<-drawWorld("y") + 
  geom_point(data=bio.pairs, aes(x=long.1, y=lat.1, colour = taxa.pairs, size = overlap.years), alpha=I(0.7))
points
# ggsave('overlap_1km.pdf', points)
ggsave('overlap_10km.pdf', points)

# Barplot of taxa pairs
# pdf('barplot_overlap_1km.pdf')
pdf('barplot_overlap_10km.pdf')
par(mai = c(.5,3.5,.1,.1)) #create space under x-axis
barplot(table(bio.pairs$taxa.pairs), las = 2, horiz = T)
dev.off()


# Save dataframe
setwd(tidy_dat)
# save(bio.pairs, file = 'bio_pairs_1km.RData')
# write.csv(bio.pairs %>% arrange(overlap.years), file = 'bio_pairs_1km.csv', row.names = F)
save(bio.pairs, file = 'bio_pairs_10km.RData')
write.csv(bio.pairs %>% arrange(overlap.years), file = 'bio_pairs_10km.csv', row.names = F)




####################################################################################################

# # GENERATE METADATA FILE FOR LIT REVIEW # # 

####################################################################################################

# Data

setwd(tidy_dat)
load('bio_pairs_10km.RData') #unique 10 km geographic-overlap years-taxa pairs (explore_biotime.R)

setwd(raw_dat)
biotimeMeta <- read.csv('BioTIMEMetadata_24_06_2021.csv')


# Select matched studies from metadata

ii <- unique(c(bio.pairs$ID.1, bio.pairs$ID.2)) #unique study IDs for 10 km pairs

meta.pairs <- biotimeMeta %>% 
  filter(STUDY_ID %in% ii) #filter our study IDs from identified pairs


# Add additional columns for lit review

meta.pairs <-  meta.pairs %>% 
  mutate(reviewed_by = NA, .after = STUDY_ID)

meta.pairs <-  meta.pairs %>% 
  mutate(eco_footprint = NA, .after = SUMMARY_METHODS)

meta.pairs <-  meta.pairs %>% 
  mutate(detectability = NA, .after = eco_footprint)


# Check license on datasets

unique(meta.pairs$LICENSE)

# https://biotime.st-andrews.ac.uk/usageGuidelines.php => all open access
# CC-by - attribution and credit
# ODbl - attribution and share alike
# ODC-by - attribution and outline where changes were made
# PDDL - no restrictions


# Save
setwd(tidy_dat)
save(meta.pairs, file = 'meta_pairs_10km.RData')
write.csv(meta.pairs, file = 'meta_pairs_10km.csv', row.names = F)




####################################################################################################

# # EXPLORE BIOTIME SPECIES DATA # # 

####################################################################################################

# # Nathalie double checking studies----
# Exclude 492 because data is only biomass

# Data
biotime.raw <- read.csv('BioTIMEQuery_24_06_2021.csv')
head(biotime.raw)

# Look at one pair from bio.pairs
head(bio.pairs)

sp1 <- filter(biotime.raw, STUDY_ID == 492)
sp2 <- filter(biotime.raw, STUDY_ID == 58)

summary(sp1$sum.allrawdata.BIOMASS)
summary(sp1$sum.allrawdata.ABUNDANCE)

table(sp1$GENUS_SPECIES) #different species

summary(sp2)
table(sp2$GENUS_SPECIES)


# Fungi data = from sequencing studies?
fung.pairs <- filter(bio.pairs, taxa.pairs )
fungi <- filter(biotime.raw, STUDY_ID == 461)


# # Courtney double checking studies----
biotime.raw <- read.csv('data/prep_biotime/BioTIMEQuery_24_06_2021.csv')
#studies reviewed by CC
#334 is all zeroes (biomass only) so removed 
IDs=c(295, 296, 300, 301, 305, 307, 311, 313, 332, 333, 458, 459, 460, 461,462,463,464,465) 
check<-filter(biotime.raw, STUDY_ID %in% IDs)

#add in taxa info 
biotimeMeta <- read.csv('data/prep_biotime/BioTIMEMetadata_24_06_2021.csv') 
names(biotimeMeta)
meta<-select(biotimeMeta, STUDY_ID, REALM, CLIMATE, TAXA, HABITAT)

check<-left_join(check, meta)

#plot abundances @ taxa level over time for each study  
ggplot(data=check, aes(y=sum.allrawdata.ABUNDANCE, x=as.factor(YEAR)))+
   geom_point(aes(color=TAXA)) + facet_wrap(~as.factor(STUDY_ID), scales='free_y')+ theme_bw() +
  ggtitle("ABUNDANCE") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#look at set pairs 
load("data/prep_biotime/bio_pairs_10km.RData")
checkpairs<-filter(bio.pairs, ID.1 %in% IDs)%>%filter(ID.2 %in% IDs)

checkpairs$ID.1
checkpairs$ID.2


