# Aims:
# 1. Visualize BioTIME data in space and complete the following tasks sprinkled throughout code below
# 2. TASK 1: Calculate correlations and visualize time series of species pairs
# 3. TASK 2: Identify study pairs that have no data in their overlapping time series
# 4. TASK 3: Explore the effect of treatments on time series data
# 5. Ensure that R and GitHub are running smoothly through these exercises

# Authors: Ryan Langendorf, Courtney Collins, Nathalie Chardon
# Date created: 13 Mar 2023
# Date updated: 23 Mar 2023 (NC)

# # LIBRARIES # #
library(tidyverse)
library("dplyr")
library("tibble")
library("readr") 
library("ggplot2")
library("magrittr")

rm(list=ls()) 


# # INPUT FILES # #
load("data/tidy/collated_pairs.RData") #collated pairs of overlapping studies
load('data/prep_biotime/bio_pairs_10km.RData') #metadata of overlapping studies




####################################################################################################

# # UPDATE R AND RSTUDIO # # 

####################################################################################################

## Please make sure you have R version 4.2.3 installed. 

# Check your version: 
R.version

# Update your version: https://cran.rstudio.com

# Note that for macOS, the newest version can either be loaded for silicon (M1 or newer) or Intel chips. 
# Check what kind you have: Apple symbol in menu bar --> About this Mac --> look at "Processor" info


## Please make sure you have RStudio version 1.4.1106 installed. 

# Check your version: RStudio in menu bar --> About RStudio --> look at Version info (2023.03.0+386 will show on Mac)

# Update your version: https://support--rstudio-com.netlify.app/products/rstudio/download/

## Here is a useful blog about the installation if you have questions: 
## https://www.r-bloggers.com/2022/01/how-to-install-and-update-r-and-rstudio/




####################################################################################################

# # VISUALIZE DATA # # 

####################################################################################################

## Load Data
#these are the datasets in bioTime that overlap both in time (>=1 year) and space (within a 10km distance) 
load("data/tidy/collated_pairs.RData")
load('data/prep_biotime/bio_pairs_10km.RData') #metadata of overlapping studies


#let's look at the data to get a sense of their structure 
head(bio.pairs)

#each row reflects a pair of 2 time series that overlap (1D.1 and 1D.2) 
#Each time series in a pair is defined by a unique taxon in a unique location 
#For example the first pair is of Marine plants & Marine invertebrates (ID 120 & 122) 
#Each taxon can have multiple species in it, e.g. there are 260 species of Marine plants in the timeseries 1D=120
#if we look closer at the organisms.1 columns we see these are all tropical algae
#You will also notice that some IDs are in more than one pair, e.g. ID 459 (Birds) is in 4 of the 6 pairs shown here 


# Let's look at how these studies are distributed across the globe

# First draw a basic world map, add "y" or "n" for display of tropics and polar latitudes

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

# Now let's use the above function to plot these studies across the globe
(gplot <- drawWorld("y") + 
    geom_point(data=bio.pairs, aes(x=long.1, y=lat.1, colour = taxa.pairs, size = overlap.years), 
               alpha=I(0.7)))

# Let's also look at the distribution of taxa pairs in these data
par(mai = c(.5,3.5,.1,.1)) #create space under x-axis
barplot(table(bio.pairs$taxa.pairs), las = 2, horiz = T)




####################################################################################################

# # TASK 1: SPECIES PAIRS # # 

####################################################################################################

###PLOT SPP PAIRS----
## So let's explore a particular pair of time series with lots of overlapping data
pair=164

pair_1_ID <- bio.pairs$ID.1[pair]
pair_2_ID <- bio.pairs$ID.2[pair]

timeseries_1 <- collated.pairs %>% dplyr::filter(ID == pair_1_ID)
timeseries_2 <- collated.pairs %>% dplyr::filter(ID == pair_2_ID)


years_overlap <- unique(timeseries_1$YEAR)[unique(timeseries_1$YEAR) %in% unique(timeseries_2$YEAR)] %>% sort()

years_overlap #so these timeseries have overlapping data for 30 years!


## Let's make a list of the species in these two timeseries
timeseries_1_species <- timeseries_1$SPECIES %>% unique()
timeseries_1_species = timeseries_1_species[!is.na(timeseries_1_species)]
timeseries_1_species_length <- sapply(timeseries_1_species, function(x){
  timeseries_1 %>%
    dplyr::filter(SPECIES == x) %>%
    dplyr::select(YEAR) %>%
    unique() %>%
    unlist() %>%
    length()
})

timeseries_2_species <- timeseries_2$SPECIES %>% unique()
timeseries_2_species = timeseries_2_species[!is.na(timeseries_2_species)]
timeseries_2_species_length <- sapply(timeseries_2_species, function(x){
  timeseries_2 %>%
    dplyr::filter(SPECIES == x) %>%
    dplyr::select(YEAR) %>%
    unique() %>%
    unlist() %>%
    length()
})

## Now we will randomly select a species from each timeseries, with at least 10 observations
sp1 = timeseries_1_species[sample(x = which(timeseries_1_species_length > 10), size = 1)]
sp2 = timeseries_2_species[sample(x = which(timeseries_2_species_length > 10), size = 1)]


#let's look at the data for one of the species 
sp1_data <- timeseries_1 %>%
  dplyr::filter(
    SPECIES == sp1,
    YEAR %in% years_overlap
  ) 
head(sp1_data)
#notice that there are multiple observations of this species within the same year with different abundance values
#these are due to sampling occurring at different spatial locations within the study area (for example plots)

#Because we want to look at how the abundances of species pairs co-vary over time, for now, let's take the average 
#of all spatial replicates within each year for each species. 

#TASK 1: There are many reasons why this might not be the most accurate approach, try to think of some different ideas for 
#how we can best approach this issue and bring them to the working group

sp1_data <- timeseries_1 %>%
  dplyr::filter(
    SPECIES == sp1,
    YEAR %in% years_overlap
  ) %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarise(
    Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
      mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
      mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(YEAR)%>%
  distinct(.)

sp2_data <- timeseries_2 %>%
  dplyr::filter(
    SPECIES == sp2,
    YEAR %in% years_overlap
  ) %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarize(
    Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
      mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
      mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(YEAR)

## So you will remember, the overlap we calculated above is at the taxon level, so for all species in the timeseries 
#but each pair of species may not overlap across all the years 
#so we need to filter the data by only the years that these species overlap in sampling 
sp1_data %<>% dplyr::filter(YEAR %in% sp2_data$YEAR)
sp2_data %<>% dplyr::filter(YEAR %in% sp1_data$YEAR)

## Calculate correlation between two species
cor_pair <- cor(sp1_data$Abundance, sp2_data$Abundance)

## Plot our timeseries
fig_data <- dplyr::bind_rows(
  sp1_data %>% dplyr::mutate(Species = sp1),
  sp2_data %>% dplyr::mutate(Species = sp2)
)

#plot correlation over time 
fig_data %>%
  ggplot() +
  theme_bw() +
  geom_line(
    aes(
      x = YEAR,
      y = Abundance,
      color = Species
    ), 
    linewidth = 3
  ) +
  labs(
    title = paste0("Correlation = ", round(cor_pair, 2))
  ) + xlab("Year")

#How does it look? Are your species highly correlated over time? positively or negatively? 
#let's look back at the metadata and see what kind of organisms we are looking at 

bio.pairs%>%subset(ID.1==pair_1_ID&ID.2==pair_2_ID)

#so you can see we are looking at small mammals and plants in tallgrass prairie 
#seems like biologically we might expect these groups to co-vary? does your species pair reflect this?
#if you want to dig deeper, let's look at the full species names 
timeseries_1%>%filter(grepl(sp1, SPECIES))%>%select(GENUS_SPECIES)%>%distinct(.)
timeseries_2%>%filter(grepl(sp2, SPECIES))%>%select(GENUS_SPECIES)%>%distinct(.)

#feel free to research more about these individual species to inform why they may or may not correlate over time



####################################################################################################

# # TASK 2: CHECK FOR ZEROS IN DATA FOR TAXA PAIRS IN OVERLAPPING YEARS # # 

####################################################################################################

# Now let's make the same graph as above, but at the taxa level.

# Abundance data for taxa 1 from the pair identified at the start of the last section
taxon1_data <- timeseries_1 %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarise(
    Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
      mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
      mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(YEAR)%>%
  distinct(.)

# Abundance data for taxa 2 from the pair identified at the start of the last section
taxon2_data <- timeseries_2 %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarise(
    Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
      mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
      mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(YEAR)%>%
  distinct(.)

# Calculate correlation between two taxa
cor_pair <- cor(taxon1_data$Abundance, taxon2_data$Abundance)

# Plot our timeseries
fig_data <- dplyr::bind_rows(
  taxon1_data %>% dplyr::mutate(Taxon = bio.pairs$taxa.1[pair]),
  taxon2_data %>% dplyr::mutate(Taxon = bio.pairs$taxa.2[pair])
)

# Plot correlation over time 
fig_data %>%
  ggplot() +
  theme_bw() +
  geom_line(
    aes(
      x = YEAR,
      y = Abundance,
      color = Taxon
    ), 
    linewidth = 3
  ) +
  labs(
    title = paste0("Correlation = ", round(cor_pair, 2))
  ) + xlab("Year")

# How does this graph compare to the species-level graph you made? Is the correlation higher or lower?

# This example had abundance data for both taxa groups in all years. What about taxa pairs that 
# have many 0s in their overlapping years? Let's look for an example of such a pair. 

zeros <- data.frame(study_id = NA, years = NA) #initialize empty dataframe

for (i in 1:nrow(bio.pairs)) {
  
  # create dataframe for each timeseries in a pair
  timeseries_1 <- collated.pairs %>% dplyr::filter(ID == bio.pairs$ID.1[i])
  timeseries_2 <- collated.pairs %>% dplyr::filter(ID == bio.pairs$ID.2[i])
  
  # summarize mean abundance per year for taxon 1
  taxon1_data <- timeseries_1 %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      Abundance1 = sum( ## This assumes there will only ever be either abundance or biomass data
        mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
        mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
        na.rm = TRUE
      )
    )
  
  # summarize mean abundance per year for taxon 2
  taxon2_data <- timeseries_2 %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarise(
      Abundance2 = sum( ## This assumes there will only ever be either abundance or biomass data
        mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
        mean(sum.allrawdata.BIOMASS, na.rm = TRUE),
        na.rm = TRUE
      )
    )
  
  # keep only overlapping years
  tax12_data <- inner_join(taxon1_data, taxon2_data, by = 'YEAR')
  
  # identify and record any years of study ID that have 0s in abundance or biomass
  if (0 %in% tax12_data$Abundance1) {
    
    zeros <- rbind(zeros, 
                   c(bio.pairs$ID.1[i], paste(min(tax12_data$YEAR), max(tax12_data$YEAR), sep = '_')))
  }
}

# Look at zeros
zeros

# Turns out there are no overlapping years with zeros in abundance or biomass! That's great news, 
# and will make our analyses easier. 

# TASK 2: Think about some reasons why having 0s in taxa 1 and taxa 2 of overlapping pairs of years
# could present some difficulties in our analyses of asking whether changes in abundance in taxa in 
# one year can predict changes in abundance of another taxa.




####################################################################################################

# # TASK 3: TREATMENTS # # 

####################################################################################################

###TREATMENTS----
#Now let's consider the broader stud(ies) and experiments in which these data were collected 
#We can read in a dataframe with additional metadata for each study from BioTime 
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use

#filter for our unique IDs   
pair=164

pair_1_ID <- bio.pairs$ID.1[pair]
pair_2_ID <- bio.pairs$ID.2[pair]

meta.pairs<-filter(meta.pairs, STUDY_ID==pair_1_ID|STUDY_ID==pair_2_ID)

#let's look at the columns that describe the methods of these studies 
head(meta.pairs$METHODS)

head(meta.pairs$GENERAL_TREAT)

#We can see for ID 311 that there were seven fire and grazing treatments with 2 trap-lines per treatment
#however we have not accounted for these different treatments in any way when plotting the species correlations above

##This could potentially cause issues in our interpretations of the different time series
#TASK 3: Try to think of some different ideas for how we might account for this and bring them to the working group

#now let's see where these studies were located 
head(meta.pairs$CENT_LAT)
head(meta.pairs$CENT_LONG)

#So it appears these studies were in the exact same location, and upon further investigation 
head(meta.pairs$DATA_SOURCE)
#we see that they are both from the Konza Prairie LTER, but they don't both report the same treatments. 
#So we might need to do some more digging to sort out these inconsistencies... 

#See if you can find any other inconsistencies in the meta-data for these 2 studies, or look through the entire 
#meta.pairs dataframe to see if anything else looks concerning! 
