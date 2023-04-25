#Filtering BIOtime (collated pairs)
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping

#load data
#these are the datasets in bioTime that overlap both in time (>=1 year) and space (within a 10km distance) 
load("data/tidy/collated_pairs.RData")
load('data/prep_biotime/bio_pairs_10km.RData') #metadata of overlapping studies
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use
#subseted<-collated.pairs[1:10000,]

working<-collated.pairs %>%
  group_by(ID,YEAR,LATITUDE,LONGITUDE,GENUS) %>%
  summarize(mean_abun=mean(sum.allrawdata.ABUNDANCE,na.rm=T),
            median_abun=median(sum.allrawdata.ABUNDANCE,na.rm=T),
            min_abun=min(sum.allrawdata.ABUNDANCE),
            max_abun=max(sum.allrawdata.ABUNDANCE),
            sd_abun=sd(sum.allrawdata.ABUNDANCE,na.rm=T),
            mean_bio=mean(sum.allrawdata.BIOMASS,na.rm=T),
            median_bio=median(sum.allrawdata.BIOMASS,na.rm=T),
            min_bio=min(sum.allrawdata.BIOMASS),
            max_bio=max(sum.allrawdata.BIOMASS),
            sd_bio=sd(sum.allrawdata.BIOMASS,na.rm=T))

#sampling effort
working <- collated.pairs %>%
  group_by(ID, YEAR) %>%
  summarize(EFFORT.YEAR = n_distinct(PLOT)) %>%
  left_join(working, ., by = c('ID', 'YEAR'))


#Zeroes
#Example studies
#Mean, median, min/max? 
#Genus distribution
#Outliers in sampling effort - why?
  #Snail one

number_of_zeroes <- working %>%
  ggplot(aes(x=mean_abun)) + 
  geom_histogram()+
  #xlim(0,10)+
  ylim(0, 500)+
  labs(x='Mean abundance of all species', y='Count')+
  theme_classic();number_of_zeroes

num_zeroes_abun_mean <- working %>%
  filter(!is.na(mean_abun)) 
num_zeroes_abun_med <- working %>%
  filter(median_abun==0) %>%
  nrow()

#Create a plot showing a sample of 10 studies, and their mean abundance

#Select 10 studies
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

#What's going on with 274?
set_274 <- working %>%
  filter(ID==274) 
original_274 <- collated.pairs %>%
  filter(ID==274)

#Sampling effort 
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
  






#Correlation over time map

pair=142

pair_1_ID <- bio.pairs$ID.1[pair]
pair_2_ID <- bio.pairs$ID.2[pair]

timeseries_1 <- working %>% dplyr::filter(ID == pair_1_ID)
timeseries_2 <- working %>% dplyr::filter(ID == pair_2_ID)

years_overlap <- unique(timeseries_1$YEAR)[unique(timeseries_1$YEAR) %in% unique(timeseries_2$YEAR)] %>% sort()



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

ggplot(sp2_data, aes(x=PLOT, y=sum.allrawdata.BIOMASS ))+
  geom_point()

#Ok, so if there is plot data could we do a mixed effects model? Incorporate differences among plots? 

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


low_sampling <- working %>%
  group_by(GENUS, ID) %>%
  summarize(n=n()) %>%
  arrange(n)%>%
  filter(n<3)

low_sampling_genus <- working %>%
  group_by(GENUS) %>%
  summarize(n=n()) %>%
  arrange(n) %>%
  filter(n<3)




