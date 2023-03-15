## Load Packages
library("dplyr")
library("tibble")
library("readr")
library("ggplot2")
library("magrittr")

## Load Data
#these are the datasets in bioTime that overlap both in time (>=1 year) and space (within a 10km distance) 
load("data/tidy/collated_pairs.RData")
bio_pairs_10km <- read.csv("data/prep_biotime/bio_pairs_10km.csv") 

#let's look at the data to get a sense of their structure 
head(bio_pairs_10km)

#each row reflects a pair of 2 time series that overlap (1D.1 and 1D.2) 
#Each time series in a pair is defined by a unique taxon in a unique location 
#For example the first pair is of Marine plants & Marine invertebrates (ID 120 & 122) 
#Each taxon can have multiple species in it, e.g. there are 260 species of Marine plants in the timeseries 1D=120
#if we look closer at the organisms.1 columns we see these are all tropical algae
#You will also notice that some IDs are in more than one pair, e.g. ID 459 (Birds) is in 4 of the 6 pairs shown here 

###PLOT SPP PAIRS----
## So let's explore a particular pair of time series with lots of overlapping data
pair=164

pair_1_ID <- bio_pairs_10km$ID.1[pair]
pair_2_ID <- bio_pairs_10km$ID.2[pair]

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

#There are many reasons why this might not be the most accurate approach, try to think of some different ideas for 
#how we can best approach this issue and bring them to the working group (TASK 1)

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

bio_pairs_10km%>%subset(ID.1==pair_1_ID&ID.2==pair_2_ID)

#so you can see we are looking at small mammals and plants in tallgrass prairie 
#seems like biologically we might expect these groups to co-vary? does your species pair reflect this?
#if you want to dig deeper, let's look at the full species names 
timeseries_1%>%filter(grepl(sp1, SPECIES))%>%select(GENUS_SPECIES)%>%distinct(.)
timeseries_2%>%filter(grepl(sp2, SPECIES))%>%select(GENUS_SPECIES)%>%distinct(.)

#feel free to research more about these individual species to inform why they may or may correlate over time



###TREATMENTS----
#Now let's consider the broader stud(ies) and experiments in which these data were collected 
#We can read in a dataframe with additional metadata for each study from BioTime 
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use

#filter for our unique IDs   
meta.pairs<-filter(meta.pairs, STUDY_ID==pair_1_ID|STUDY_ID==pair_2_ID)

#let's look at the columns that describe the methods of these studies 
head(meta.pairs$METHODS)

head(meta.pairs$GENERAL_TREAT)

#We can see for ID 311 that there were seven fire and grazing treatments with 2 trap-lines per treatment
#however we have not accounted for these different treatments in any way when plotting the species correlations above

##This could potentially cause issues in our interpretations of the different time series
#Try to think of some different ideas for how we might account for this and bring them to the working group (TASK 3)

#now let's see where these studies were located 
head(meta.pairs$CENT_LAT)
head(meta.pairs$CENT_LONG)

#So it appears these studies were in the exact same location, and upon further investigation 
head(meta.pairs$DATA_SOURCE)
#we see that they are both from the Konza Prairie LTER, but they don't both report the same treatments. 
#So we might need to do some more digging to sort out these inconsistencies... 

#See if you can find any other inconsistencies in the meta-data for these 2 studies, or look through the entire 
#meta.pairs dataframe to see if anything else looks concerning! 