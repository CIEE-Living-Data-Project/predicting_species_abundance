
library(tidyverse)
library(dplyr)
library(ggplot2)

# # INPUT FILES # #
biotime.raw <- read.csv('data/prep_biotime/BioTIMEQuery_24_06_2021.csv') #in gitignore - must have saved locally
biotimeMeta <- read.csv('data/prep_biotime/BioTIMEMetadata_24_06_2021.csv') 
load('data/prep_biotime/bio_pairs_10km.RData') #unique 10 km geographic-overlap years-taxa pairs (explore_biotime.R)

### check which studies we actually use ######
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23
dat_terr <- subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
dat_terr <- subset(dat_terr, Metric!="CROSS" &  Type!="Between") %>%
  select(-c("Prop.Change.Gn1", "Prop.Change.Gn2", "YEAR.T","YEAR.T1", "SERIES.start", "SERIES.end")) %>% 
  distinct(.)

study_id_used <- do.call(rbind, strsplit(unique(dat_terr$PairID), split="_"))[,1] # 26 studies

### Plot full data on map ####

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

### filter BioTime raw to only types of data that want in study ####
##### ORGANIZE DATAFRAME #####

# Keep only relevant columns in biotimeMeta
bio <- biotimeMeta %>% select(c("REALM", "CLIMATE",'HABITAT', 'TAXA', 'ORGANISMS', 'START_YEAR', 'END_YEAR', 
                                'NUMBER_OF_SPECIES', 'STUDY_ID', "CENT_LAT","CENT_LONG", "TREATMENT", "TREAT_COMMENTS"))

#bio <- bio %>% rename(ID.1 = STUDY_ID) #column names need to match for join
head(bio)

# Match data from ID.1 to bio.raw
bio.single <- left_join(biotime.raw, bio, by = 'STUDY_ID') #join ID.1 info

head(bio.single)



#### IDENTIFY TEMPORALLY PROXIMAL PAIRS ####
# source: https://search.r-project.org/CRAN/refmans/DescTools/html/overlaps.html

####################################################################################################

# Dataframes with only year info
df1 <- bio.single %>% select(START_YEAR, END_YEAR)

head(df1)

# Identify overlap amount
vv.over <- df1[,2]-df1[,1] #amount of overlap
hist(vv.over)

# Add overlap vector to bio.pairs
bio.single$overlap.years <- vv.over 

# Keep only 10 overlapping years
bio.single.10 <- bio.single %>% 
  filter(overlap.years > 9)
dim(bio.single.10)
summary(bio.single.10$overlap.years) # 10-129 years overlap, median 29 years, mean 32.3 years

# only abundance
# Exclude 492 because data is only biomass
# 334 is all zeroes (biomass only) so removed 
bio.single.10 <- bio.single.10[-which(bio.single.10$STUDY_ID==492), ]
bio.single.10 <- bio.single.10[-which(bio.single.10$STUDY_ID==334), ]

# terrestrial
bio.single.10_terr <- subset(bio.single.10, REALM=="Terrestrial")
length(unique(bio.single.10_terr$STUDY_ID)) # now we have 127 studies not 26...

# taxa category "All" is all herps so rename to amphibs
bio.single.10_terr$TAXA <- ifelse(bio.single.10_terr$TAXA=="All", "Amphibians", bio.single.10_terr$TAXA)


##### plot cleaned data on map #####

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
taxaCol<-c('#ffffbf','#5e4fa2','#fdae61','#a8c614','#66c2a5','#9e0142')

# create plot
points<-drawWorld("y")+geom_point(data=distinct(select(bio.single.10_terr, CENT_LONG, CENT_LAT, TAXA)), 
                                              aes(x=CENT_LONG, y=CENT_LAT, colour=TAXA), size=2,alpha=I(0.5))
points<-points+scale_colour_manual(values=taxaCol) #+scale_size(range=c(3, 10))
points


unique(biotimeMeta$TAXA)
unique(bio.single.10_terr$TAXA)

