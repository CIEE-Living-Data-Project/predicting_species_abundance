#Filtering BIOtime (collated pairs)
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping

#load data
load("data/tidy/collated_pairs.RData") #full dataset with species/abundances
load('data/prep_biotime/bio_pairs_10km.RData') #"key" with overlapping studies, and reference
                                              #to species in collated.pairs
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use
#subseted<-collated.pairs[1:10000,]

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


#Wrangle the collated.pairs into genus organization
working<-filtered.collated.pairs %>% #read in filtered collated pairs
  group_by(ID,YEAR,LATITUDE,LONGITUDE,GENUS) %>% #group by study, location, genus
  summarize(mean_abun=mean(sum.allrawdata.ABUNDANCE,na.rm=T), #mean abundance
            median_abun=median(sum.allrawdata.ABUNDANCE,na.rm=T), #median abundance
            min_abun=min(sum.allrawdata.ABUNDANCE), #min abundance
            max_abun=max(sum.allrawdata.ABUNDANCE), #max abundance
            sd_abun=sd(sum.allrawdata.ABUNDANCE,na.rm=T), #sd abundance (when more than one species per genera)
            mean_bio=mean(sum.allrawdata.BIOMASS,na.rm=T), #mean biomass
            median_bio=median(sum.allrawdata.BIOMASS,na.rm=T), #median biomass
            min_bio=min(sum.allrawdata.BIOMASS), #min biomass
            max_bio=max(sum.allrawdata.BIOMASS), #max biomass
            sd_bio=sd(sum.allrawdata.BIOMASS,na.rm=T)) #sd biomass 

#calculate sampling effort
working <- collated.pairs %>% 
  group_by(ID, YEAR) %>% #group by study, year
  summarize(EFFORT.YEAR = n_distinct(PLOT)) %>% #count number of distinct plots 
  left_join(working, ., by = c('ID', 'YEAR')) #join with new dataset above by year and ID

write.csv(working, "data/cleaned_collated_pairs_EBMW.csv")


####################################################################################################

## FINISHED DATA CLEANING - VISUALIZATION AND EXPLORATION BELOW

####################################################################################################

#Double checking genus low sampling - did removing rare species remove rare genera? 
low_sampling_cleaned <- working %>%
  group_by(GENUS) %>%
  summarize(n=n()) %>%
  arrange(n)%>%
  filter(n<3) %>%
  nrow() %>%
  print()

low_sampling_collated <- collated.pairs %>%
  group_by(GENUS) %>%
  summarize(n=n()) %>%
  arrange(n)%>%
  filter(n<3)%>%
  nrow() %>%
  print()

#removed almost half the rare genera! 



#histogram of mean abundance across all studies and years
number_of_zeroes <- working %>%
  ggplot(aes(x=mean_abun)) + 
  geom_histogram()+
  #xlim(0,10)+ #uncomment to see more detail in plot
  #ylim(0, 500)+
  labs(x='Mean abundance of all species', y='Count')+
  theme_classic();number_of_zeroes
#looks like some datasets have very large abundances, while some have very small 

#histogram of mean abundance across all studies and years
number_of_zeroes <- working %>%
  ggplot(aes(x=mean_bio)) + 
  geom_histogram()+
  #xlim(0,10)+ #uncomment to see more detail in plot
  #ylim(0, 500)+
  labs(x='Mean biomass of all species', y='Count')+
  theme_classic();number_of_zeroes
#looks like some datasets have very large biomass, while some have very small 

#see if there are any zero mean abundance values
num_zeroes_abun_mean <- working %>%
  filter(mean_abun==0) # %>%
  #nrow()
#none
#see if there are any zero median abundance values 
num_zeroes_abun_med <- working %>%
  filter(median_abun==0) %>%
  nrow()
#none


#Create a plot showing a sample of 10 studies, and their mean abundance

#Select 10 random studies
selected_studies <- sample(unique(working$ID), 10)
ggplot(working[working$ID %in% selected_studies, ], aes(x = mean_abun)) + 
  geom_histogram(binwidth = 0.5, color = "black", fill = "lightblue") + 
  facet_wrap(~ID, nrow = 2) +
  ylim(0,500)+
  xlim(0, 500)+
  labs(title = "Mean Abundance by Study ID", x = "Mean Abundance", y = "Frequency")+
theme_classic()

#What are these studies? 
bio_sets <- bio.pairs %>%
  filter(ID.1 %in% selected_studies)


#Sampling effort  of all studies as histogram
ggplot(working, aes(x = EFFORT.YEAR)) + 
  geom_histogram(binwidth = 10, color = "black", fill = "lightblue") + 
  labs(x = "Sampling effort", y = "Frequency")+
  theme_classic()

#What had sampling efforts greater than 500?
high_sampling_effort <- working %>%
  filter(EFFORT.YEAR > 400) 
high_studies <- unique(high_sampling_effort$ID)
#These studies are: 
bio_sets <- bio.pairs %>%
  filter(ID.1 %in% high_studies)
collated_sets <- working %>%
  filter(ID %in% high_studies, 
         EFFORT.YEAR>400)
collated_sets <- working %>%
  filter(ID %in% high_studies) %>%
  nrow()
  




#Correlation over time map - sample to see if data cleaned properly

#Sample pair - daphnia and their algal food source (found to be predator by 
                                                    #Globi package)

pair=142

pair_1_ID <- bio.pairs$ID.1[pair] #Pair 1
pair_2_ID <- bio.pairs$ID.2[pair] #Pair 2

timeseries_1 <- working %>% dplyr::filter(ID == pair_1_ID) #get data for genus
timeseries_2 <- working %>% dplyr::filter(ID == pair_2_ID) #get data for genus

#How many years overlap? 
years_overlap <- unique(timeseries_1$YEAR)[unique(timeseries_1$YEAR) %in% 
                                             unique(timeseries_2$YEAR)] %>% sort()


#Summarize daphnia data in each year (already done in cleaned dataset)
sp1_data <- timeseries_1 %>%
  dplyr::filter(
    GENUS == 'Daphnia',
    YEAR %in% years_overlap
  ) %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarise(
    Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
      mean(mean_abun, na.rm = TRUE),
      mean(mean_bio, na.rm = TRUE),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(YEAR)%>%
  distinct(.)


#Summarize algae data in each year (already done in cleaned dataset)
sp2_data <- timeseries_2 %>%
  dplyr::filter(
    GENUS == 'Ankistrodesmus',
    YEAR %in% years_overlap
  ) %>%
  dplyr::group_by(YEAR) %>%
  dplyr::summarize(
    Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
      mean(mean_abun, na.rm = TRUE),
      mean(mean_bio, na.rm = TRUE),
      na.rm = TRUE
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(YEAR)


## So you will remember, the overlap we calculated above is at the taxon level, so for all species in the timeseries 
#but each pair of species may not overlap across all the years 
#so we need to filter the data by only the years that these species overlap in sampling 
sp1_data %<>% dplyr::filter(YEAR %in% sp2_data$YEAR)
sp1_data <- sp1_data %>%
  mutate(scaled_abundance = scale(Abundance))
sp2_data %<>% dplyr::filter(YEAR %in% sp1_data$YEAR)
sp2_data <- sp2_data %>%
  mutate(scaled_abundance = scale(Abundance))

## Calculate correlation between two species
cor_pair <- cor(sp1_data$scaled_abundance[,1], sp2_data$scaled_abundance[,1])

## Plot our timeseries
fig_data <- dplyr::bind_rows(
  sp1_data %>% dplyr::mutate(Species = 'Daphnia'),
  sp2_data %>% dplyr::mutate(Species = 'Anistrodesmus')
)

#plot correlation over time 
fig_data %>%
  ggplot() +
  theme_bw() +
  geom_line(
    aes(
      x = YEAR,
      y = scaled_abundance[,1],
      color = Species
    ), 
    linewidth = 3
  ) +
  labs(
    title = paste0("Correlation = ", round(cor_pair, 2))
  ) + xlab("Year")






