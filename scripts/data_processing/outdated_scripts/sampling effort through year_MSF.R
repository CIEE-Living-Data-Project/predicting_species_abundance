#### script that calculates sampling effort through year 
#### authors: Gavia Lertzman-Lepofsky, Haley Branch
#### 25 April 2023
#### further work done by Mia, Sophia, and Faraz 26 April 2023, indicated below


# libraries
library(dplyr)
library(tidyr)

rm(list=ls()) 



#### DATA #####
load("./data/tidy/collated_pairs.RData") # collated pairs of overlapping studies
load('./data/prep_biotime/bio_pairs_10km.RData') # metadata of overlapping studies
load("./data/prep_biotime/meta_pairs_10km.RData") # metadata on study methods

# piped statement that calculates summary statistics per genus per year
collated.pairs_genus <- collated.pairs %>%
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
#This swallows the plot-specific data into a single genus row per study


# calculates sampling effort per year by summing number of plots 
# but this doesn't account for how sampling effort might be spread through year
collated.pairs_genus1 <- collated.pairs %>%
  group_by(ID, YEAR) %>%
  summarize(EFFORT.YEAR = n_distinct(PLOT)) %>%
  left_join(collated.pairs_genus, ., by = c("ID", "YEAR"))
#number of plots sampled per study per year is propagated across all study/year genera because of left_join

### ----- New code from Mia, Sophia, Faraz  -----------------

#JUST a dataset with study ID, year, and number of plots sampled:
collated.pairs_genus2 <- collated.pairs %>%
  group_by(ID, YEAR) %>%
  summarize(EFFORT.YEAR = n_distinct(PLOT))
  
plot <- ggplot(data = collated.pairs_genus2, aes(x = YEAR, y = EFFORT.YEAR, group = ID, color = ID)) +
  theme_classic() +
  geom_line();plot
#this plots the number of plots per study (colors) across years

collated.pairs_genus2 %>%
  filter(EFFORT.YEAR >= 200) %>% 
  summarize(big_studies = n_distinct(ID))
#gives the IDs of each study that had greater than 200 plots in a year (4 studies)

#what were those 4 studies? Use STUDY_ID to get a feel for what habitats etc
df <- meta.pairs %>%
  filter(STUDY_ID %in% c(54, 295, 296, 355)) %>% 
  head()
#snail study, among others. HAS_PLOTS is confusing - what is S?

#zooming in on studies other than those 4 studies
plot2 <- collated.pairs_genus2 %>% 
  filter(!(ID %in% c(54, 295, 296, 355))) %>% 
  ggplot(aes(x = YEAR, y = EFFORT.YEAR, group = ID, color = ID)) + 
  theme_classic() +
  geom_line();plot2
#lots of variation in sampling effort (plot number) among years here too


#Next steps: 
#Standardize the abundance of each genus by the number of plots sampled PER year in each study
#What does "HAS_PLOT" mean in meta.data?



### -------------------------------




# 1. number of months that sampling effort occurs over
# counts number of months of sampling effort, not organized by season
collated.pairs_month <- collated.pairs %>% 
  filter(is.na(MONTH)!=TRUE) %>%   # not study has month, so this is only for sites that have month info
  group_by(ID, YEAR) %>%
  summarize(no.months=n_distinct(MONTH))
# tested:
# length(unique(subset(collated.pairs, ID==54 & YEAR=="1991")$MONTH)), returns same value as dplyr pipe

# visualize: number of months that have at least one plot surveyed
hist(collated.pairs_month$no.months)


# 2. now calculate sampling effort through year but grouped into season: winter, spring, summer, fall
# each a 3 month chunk of time, which months are winter/spring/summer/fall is determined by lattitude

# assign hemisphere
# tested, correctly assigned head(subset(collated.pairs, LATITUDE>=-26 & LATITUDE<=26))
collated.pairs$HEMISPHERE <- ifelse(collated.pairs$LATITUDE>=23.43626, "north", 
                                    ifelse(collated.pairs$LATITUDE<=-23.43626, "south", "equatorial"))

# assign seasonality based on hemisphere classification
# north
collated.pairs_north <- subset(collated.pairs, HEMISPHERE=="north")
collated.pairs_north$SEASON <- ifelse(as.character(collated.pairs_north$MONTH)=="6" | as.character(collated.pairs_north$MONTH)=="7" | as.character(collated.pairs_north$MONTH)=="8",  "summer", 
                                  ifelse( as.character(collated.pairs_north$MONTH)=="9" | as.character(collated.pairs_north$MONTH)=="10" | as.character(collated.pairs_north$MONTH)=="11", "fall",
                                       ifelse(as.character(collated.pairs_north$MONTH)=="12" | as.character(collated.pairs_north$MONTH)=="1" | as.character(collated.pairs_north$MONTH)=="2", "winter", 
                                              ifelse(as.character(collated.pairs_north$MONTH)=="3" | as.character(collated.pairs_north$MONTH)=="4" | as.character(collated.pairs_north$MONTH)=="5","spring", NA))))

# south
collated.pairs_south <- subset(collated.pairs, HEMISPHERE=="south")
collated.pairs_south$SEASON <- ifelse(as.character(collated.pairs_south$MONTH)=="12" | as.character(collated.pairs_south$MONTH)=="1" | as.character(collated.pairs_south$MONTH)=="2",  "summer", 
                                      ifelse( as.character(collated.pairs_south$MONTH)=="3" | as.character(collated.pairs_south$MONTH)=="4" | as.character(collated.pairs_south$MONTH)=="5", "fall",
                                              ifelse(as.character(collated.pairs_south$MONTH)=="6" | as.character(collated.pairs_south$MONTH)=="7" | as.character(collated.pairs_south$MONTH)=="8", "winter", 
                                                     ifelse(as.character(collated.pairs_south$MONTH)=="9" | as.character(collated.pairs_south$MONTH)=="10" | as.character(collated.pairs_south$MONTH)=="11","spring", NA))))

# equatorial
# no true seasons, so just filled with equator
collated.pairs_equatorial <- subset(collated.pairs, HEMISPHERE=="equatorial")
collated.pairs_equatorial$SEASON <- "equator"

# combine three hemisphere data sets back together
collated.pairs_seasonal <- rbind(collated.pairs_south, collated.pairs_north, collated.pairs_equatorial)

# count number of months that have a plot surveyed in it
# can visualize then which seasons have have sampling efforts
collated.pairs_seasonal_summary <- collated.pairs_seasonal %>% 
  filter(is.na(MONTH)!=TRUE) %>% 
  group_by(ID, YEAR, SEASON) %>%
  summarize(no.months.sampled = n_distinct(MONTH))

# visualize sampling efforts per season
ggplot(collated.pairs_seasonal_summary, aes(no.months.sampled)) + 
  geom_histogram() +
  facet_wrap(~SEASON, scales="free")

# convert this to wide format where SEASON each has own column, 
# this would allow seasonality of sampling to be merged back into sampling effort
# note that because no true seasons in equator, yearly sampling effort is not binned to season
collated.pairs_seasonal_summary_wide <- spread(collated.pairs_seasonal_summary, SEASON, no.months.sampled)

# bind back to full data if wanted 
#collated.pairs_seasonal_genus <- left_join(collated.pairs_genus, collated.pairs_seasonal_summary)

