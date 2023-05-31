##Main Data Cleaning Script
##created and updated by: by GLL, AF, EB, MW, AD, IE
##Script compiled and retested by MW and FK



#Filtering BIOtime (collated pairs)
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping
library(progress)
library('rglobi')



#load data
load("data/prep_biotime/collated_pairs.RData") #full dataset with species/abundances
load('data/prep_biotime/bio_pairs_10km.RData') #"key" with overlapping studies, and reference
#to species in collated.pairs
load('data/prep_biotime/meta_pairs_10km.RData') #biotime metadata for 10 km pairs to use
#subseted<-collated.pairs[1:10000,]

#load functions
source("scripts/data_processing/Functions_data.processing.R")


# adds taxa ID so that easy to ID which genera belong to what type of animal/plant/fungi
collated.pairs$ID <- as.character(collated.pairs$ID)
meta.pairs$STUDY_ID <- as.character(meta.pairs$STUDY_ID)
collated.pairs <- left_join(collated.pairs, meta.pairs[, c(1, 3,4, 12, 13)], by=c("ID"="STUDY_ID"))


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

#######
#Standardizing abundances by sampling effort
#######

#Visualizing variation in sampling effort across years within each study:

temp <- filtered.collated.pairs %>% 
  group_by(ID, YEAR) %>% #group by study, year
  summarize(EFFORT.YEAR = n_distinct(PLOT)) #count number of distinct plots

plot <- ggplot(data = temp, aes(x = YEAR, y = EFFORT.YEAR, group = ID, color = ID)) +
  theme_classic() +
  geom_line();plot
#this plots the number of plots per study (colors) across years


#zooming in on studies other than those 4 studies
plot2 <- temp %>% 
  filter(!(ID %in% c(54, 295, 296, 355))) %>% 
  ggplot(aes(x = YEAR, y = EFFORT.YEAR, group = ID, color = ID)) + 
  theme_classic() +
  geom_line();plot2
#lots of variation in sampling effort (plot number) among years here too

# 1. back at the species level, standardize abundances by sum.allrawdata.ABUNDANCE or BIOMASS / EFFORT.YEAR

collated.pairs_standardized <- filtered.collated.pairs %>% 
  group_by(ID, YEAR) %>% #group by study, year
  summarize(EFFORT.YEAR = n_distinct(PLOT)) %>% #count number of distinct plots 
  left_join(filtered.collated.pairs, ., by = c('ID', 'YEAR')) %>% 
  mutate(ST.ABUN = sum.allrawdata.ABUNDANCE/EFFORT.YEAR) %>% 
  mutate(ST.BIO = sum.allrawdata.BIOMASS/EFFORT.YEAR)

# 2. Calculate mean, median, min, max, sd, and CoV for standardized abundances and biomass

## Only genera and study ID kept (not Lat/Long), which standardizes Lat/Long within each study
## Aggregated species to the genus level
## tested to make sure quantile calculation is correct
collated.pairs_standardized_summary <- collated.pairs_standardized %>% #read in filtered collated pairs
  group_by(ID,YEAR,GENUS,TAXA, ORGANISMS, CLIMATE, REALM) %>% #group by study, location, genus
  summarize(mean_abun_st=mean(ST.ABUN,na.rm=T), #mean abundance
            median_abun_st=median(ST.ABUN,na.rm=T), #median abundance
            min_abun_st=min(ST.ABUN), # min abundance
            avg_abun_min_5perc = mean(ST.ABUN[ST.ABUN <= quantile(ST.ABUN, .05, na.rm = T)], 
                                     na.rm=T), # calculate specified quantile (in this case lowest 5% of data), and then calc mean of these values
            max_abun_st=max(ST.ABUN), # max abundance
            avg_abun_max_5perc = mean(ST.ABUN[ST.ABUN <= quantile(ST.ABUN, .95, na.rm = T)], 
                                     na.rm=T), # mean of top 5% values, abundance
            sd_abun_st=sd(ST.ABUN,na.rm=T), #sd abundance (when more than one species per genera)
            CoV_abun_st = sd_abun_st/mean_abun_st, # coef of variation
            mean_bio_st=mean(ST.BIO,na.rm=T), #mean biomass
            median_bio_st=median(ST.BIO,na.rm=T), #median biomass
            min_bio_st=min(ST.BIO), #min biomass
            avg_bio_min_5perc = mean(ST.BIO[ST.BIO <= quantile(ST.BIO, .05, na.rm = T)], 
                                         na.rm=T), # mean lowest 5% biomass
            max_bio_st=max(ST.BIO), #max biomass
            avg_bio_max_5perc = mean(ST.ABUN[ST.BIO <= quantile(ST.BIO, .95, na.rm = T)], 
                                          na.rm=T), # mean of top 5% values, biomass
            sd_bio_st=sd(ST.BIO,na.rm=T), #sd biomass
            CoV_bio_st = sd_bio_st/mean_bio_st) #coef of variation

#Save file as needed using function below
#save(collated.pairs_standardized_summary, file = "data/prep_biotime/collated.pairs_standardized_summary_GLL.Rdata")

######
#Find overlapping consecutive years for all pairs
######

bio.pairs <- read_csv("data/prep_biotime/bio_pairs_10km.csv")


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

#Save as needed
#saveRDS(between.studies.overlap,"data/data_processing/between.studies.overlap.RDS")

#Load previously saved version
readRDS(./data/preprocessing/between.studies.overlap)


######
# Calculate portion change for species pairs
######




pairs.keep=between.studies.overlap[which(between.studies.overlap$Max.Overlap>9),] #this is generated in the Overlap script



results=data.frame(matrix(ncol=14,nrow=0, dimnames=list(NULL, c("Gn1", "Gn2", "Log.prop.change.abun.Gn1","Log.prop.change.abun.Gn2","Log.prop.change.bio.Gn1","Log.prop.change.bio.Gn2","PairID","Type","SERIES.n","SERIES.start","SERIES.end","SERIES.l","YEAR.T","YEAR.T1")))) #makes an empty dataframe
pb<-set.prog.bar(nrow(pairs.keep))
for (i in 1:nrow(pairs.keep)){
  pb$tick()
  wk<-get.log.prop.change(i,collated.pairs_standardized_summary,pairs.keep)
  results<-rbind(results,wk)
  
} #get prop change for every genera pair

#add metadata

meta.data<-make.meta(results,collated.pairs_standardized_summary)
log.prop.change.with.meta<-left_join(results,meta.data) #join meta data with results df

#add unique genera ID col
log.prop.change.with.meta$UNIQUE.PAIR.ID=paste(log.prop.change.with.meta$Gn1,log.prop.change.with.meta$Gn2,log.prop.change.with.meta$PairID,sep="_")

#save
saveRDS(log.prop.change.with.meta,"data/data_processing/log.prop.change.with.meta.RDS")



######
# Assign interaction types between pairs
######

#Read in the updated modelling file
log.prop.change.with.meta <- readRDS("data/preprocessing/log.prop.change.full.data.RDS")
log_change <- log.prop.change.with.meta
head(log_change)

#get list of genus pairs from log.prop.change 
genus_pairs <- apply(log_change[, c("Gn1", "Gn2")], 1, 
                     function(x) paste(sort(x), collapse = "-")) #merge genus columns
unique_pairs <- unique(genus_pairs) #get unique genera pairs
split_pairs <- data.frame(Gn1 = sapply(strsplit(unique_pairs, "-"), "[", 1),
                          Gn2 = sapply(strsplit(unique_pairs, "-"), "[", 2)) #split pairs back up

#The existing data is too large to run Globi on all together
# Split into 5 groups of roughly equal size...
#To limit the amount of data being run through Globi

#Split the data into groups
groups <- cut(seq_along(1:nrow(split_pairs)), breaks = 5, labels = FALSE)

# Split groups into separate data frames
df_split <- split(split_pairs, groups)
names(df_split) <- paste0("split_", names(df_split))

#Read list into environment
list2env(df_split, envir = .GlobalEnv)

#Read progress bar function written by Isaac
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


df_list <- list()


#Split 1: 

#Write progress bar function for split 1: 
pb<-set.prog.bar(nrow(split_1)) #sets progress bar
#Run the Globi function on split 1: 
for (i in 1:nrow(split_1)) {
  pb$tick()
  gn1 <- split_1[i, "Gn1"]
  gn2 <- split_1[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_1 <- do.call(rbind, df_list)
write.csv(pair_interactions_1, "data/data_processing/pair_interactions_1.csv")




#Write progress bar function for split 2: 
pb<-set.prog.bar(nrow(split_2)) #sets progress bar
#Run the Globi function on split 2: 
for (i in 1:nrow(split_2)) {
  pb$tick()
  gn1 <- split_2[i, "Gn1"]
  gn2 <- split_2[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}
pair_interactions_2 <- do.call(rbind, df_list)
write.csv(pair_interactions_2, "data/data_processing/pair_interactions_2.csv")



#Write progress bar function for split 3: 
pb<-set.prog.bar(nrow(split_3)) #sets progress bar
#Run the Globi function on split 3: 
for (i in 1:nrow(split_3)) {
  pb$tick()
  gn1 <- split_3[i, "Gn1"]
  gn2 <- split_3[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_3 <- do.call(rbind, df_list)
write.csv(pair_interactions_3, "data/data_processing/pair_interactions_3.csv")


#Write progress bar function for split 4:
pb<-set.prog.bar(nrow(split_4)) #sets progress bar
#Run the Globi function on split 4: 
for (i in 1:nrow(split_4)) {
  pb$tick()
  gn1 <- split_4[i, "Gn1"]
  gn2 <- split_4[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_4 <- do.call(rbind, df_list)
write.csv(pair_interactions_4, "data/data_processing/pair_interactions_4.csv")




#Write progress bar function for split 5:
pb<-set.prog.bar(nrow(split_5)) #sets progress bar
#Run the Globi function on split 5:
for (i in 1:nrow(split_5)) {
  pb$tick()
  gn1 <- split_5[i, "Gn1"]
  gn2 <- split_5[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_5 <- do.call(rbind, df_list)
write.csv(pair_interactions_5, "data/data_processing/pair_interactions_5.csv")


#Merge all the interaction splits back into a single dataframe
all_interactions <- rbind(pair_interactions_1, pair_interactions_2, pair_interactions_3, 
                          pair_interactions_4, pair_interactions_5)
colnames(all_interactions) <- c("Gn1", "Gn2", "interaction")
write.csv(all_interactions, "data/data_processing/genus_interaction_list.csv")


#Read in all_interactions_2 
all_interactions <- read.csv("data/data_processing/genus_interaction_list.csv")

#Remove interaction pairs that are mistakes (pathogen, host, hemiparasite)
#These were identified with manual inspection
mistake_interactions <- c("hasHost", "pathogenOf", "hemiparasiteOf", 'hostOf', 'pathogenOf', 
                          'hasVector', 'parasiteOf')
all_interactions_2<- all_interactions %>%
  filter(!interaction %in% mistake_interactions)
# write.csv(all_interactions_2, "data/data_processing/genus_interaction_list.csv")

#Get number of times each pair shows up in the data
summary_interactions <- all_interactions_2 %>%
  group_by(Gn1, Gn2) %>%
  summarize(n=n()) %>%
  arrange(desc(n))

#Get list of unique pairs 
unique_interactions_1 <-summary_interactions %>%
  select(Gn1, Gn2) %>%
  distinct()


#Next steps: start collapsing interaction types into larger groups
#Key: 
#predator-prey (including herbivory)
#eatenBy, eats, preysOn, preyedUponBy

#mutualism
#mutualistOf

#dispersal 
#hasDispersalVector, dispersalVectorOf

#visits
#flowersVisitedBy, #visits, #visitsFlowersOf, #visitedBy

#uncategorized_interaction
#interactsWith
#All others were removed 


# Use case_when to classify the interactions into broader types 
cased_interactions<- all_interactions_2%>%
  mutate(
    interaction_type = case_when(
      interaction %in%  c("preysOn", "preyedUponBy", "eatenBy", "eats") ~ "predator_prey",
      interaction %in% c("mutualistOf") ~ "mutualism",
      interaction %in% c("hasDispersalVector", "dispersalVectorOf") ~ "dispersal",
      interaction %in% c("interactsWith", "flowersVisitedBy", "visits", "visitsFlowersOf", "visitedBy") ~ "uncategorized_interaction",
      TRUE ~ NA_character_
    ),
    .keep = "unused"
  )

# Remove pairs whose characters don't match anything
cased_interactions <- cased_interactions %>%
  filter(!is.na(interaction_type))
cased_interaction_unique <- cased_interactions %>% distinct()

#Check how many pairs still have multiple interactions assigned
summary_interactions_cased<- cased_interaction_unique %>%
  group_by(Gn1, Gn2) %>%
  summarize(n=n()) %>%
  arrange() %>%
  left_join(cased_interaction_unique, by=c('Gn1', "Gn2"))


cased_interactions_filtered_2<- cased_interactions %>%
  group_by(Gn1, Gn2) %>%
  mutate(num_interactions=n())



#Identify whether there is an interaction in pairs for log_change

#Get pairs 
distinct_pairs <- cased_interactions_filtered_2 %>%
  select(Gn1, Gn2) %>%
  distinct()

#Read in interactions 
log_change <- readRDS("data/data_processing/log.prop.change.with.meta.RDS")
head(log_change)

# Create new column called interaction_present and set all values to 0
log_change$interaction_present <- 0

# Loop through each row in log_change and check if it matches any pairs with interactions
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

pb<-set.prog.bar(nrow(log_change)) #sets progress bar
for (i in 1:nrow(log_change)) {
  pb$tick()
  if (paste(log_change$Gn1[i], log_change$Gn2[i], sep = "-") %in% 
      paste(distinct_pairs$Gn1, distinct_pairs$Gn2, sep = "-")
      | paste(log_change$Gn2[i], log_change$Gn1[i], sep = "-") %in% 
      paste(distinct_pairs$Gn1, distinct_pairs$Gn2, sep = "-")) {
    log_change$interaction_present[i] <- 1
  }
}

colnames(log_change)[26] <- "interaction_found"

saveRDS(log_change, "data/data_processing/log.prop.change.interactions.RDS")




#Do positive negative neutral 
log_change_interaction <- readRDS("data/data_processing/log.prop.change.interactions.RDS")


pos_neg_interactions<- cased_interactions_filtered_2%>%
  mutate(
    interaction_benefit = case_when(
      interaction_type %in%  c("predator_prey") ~ "negative",
      interaction_type %in% c("mutualism", "dispersal") ~ "positive",
      interaction_type %in% c("uncategorized_interaction") ~ "neutral",
      TRUE ~ NA_character_
    ),
    .keep = "unused"
  ) %>%
  distinct()

#Assign positive/negative interactions 
dummy_data <- log_change_interaction %>%
  sample_n(10000)

#Create vectors to store
positive_interaction <- vector("integer", nrow(log_change_interaction))
negative_interaction <- vector("integer", nrow(log_change_interaction))
neutral_interaction <- vector("integer", nrow(log_change_interaction))


# Iterate over each row in the log_prop_change_interaction dataset

pb<-set.prog.bar(nrow(log_change_interaction)) #sets progress bar
#take interaction types, and hyphenate if there are multiple
for (i in 1:nrow(log_change_interaction)) {
  pb$tick()
  
  # Get the genus pair from the current row
  gn1 <- log_change_interaction$Gn1[i]
  gn2 <- log_change_interaction$Gn2[i]
  
  # Find the corresponding row(s) in the pos_neg_interactions dataset
  match_rows <- (pos_neg_interactions$Gn1 == gn1 & pos_neg_interactions$Gn2 == gn2) |
    (pos_neg_interactions$Gn1 == gn2 & pos_neg_interactions$Gn2 == gn1)
  
  if (any(match_rows)) {
    # If matching row(s) are found, retrieve interaction types
    interaction_benefit <- pos_neg_interactions$interaction_benefit[match_rows]
    
    # Hyphenate multiple interaction types or assign "NA" if none
    if (length(interaction_benefit) > 0) {
      hyphenated_interaction <- paste(interaction_benefit, collapse = "-")
    } else {
      hyphenated_interaction <- "NA"
    }
    
    # Assign the hyphenated interaction to the corresponding column
    log_change_interaction$interaction_benefit[i] <- hyphenated_interaction
  } else {
    # If no matching row is found, assign "NA" to the interaction type
    log_change_interaction$interaction_benefit[i] <- "NA"
  }
}


#saveRDS(log_change_interaction, "data/preprocessing/log.prop.change.interactions.hyphenated.RDS")


#interaction types detailed - hyphenated

#Read in hyphenated data from positive/negative code above (if starting from here)
#log_change_interaction <- readRDS("data/preprocessing/log.prop.change.interactions.hyphenated.RDS")
#Create a dummy dataset to test function below
dummy_data <- sample_n(log_change_interaction, 10000)


pb <- set.prog.bar(nrow(log_change_interaction)) # Sets progress bar
# Create an empty vector to store the interaction types for each row
interaction_types_detailed <- vector("character", length = nrow(log_change_interaction))
# Define a for loop to process each row, and assign interaction types to each pair
for (i in 1:nrow(log_change_interaction)) {
  pb$tick()
  
  # Get the genus pair from the current row
  gn1 <- log_change_interaction$Gn1[i]
  gn2 <- log_change_interaction$Gn2[i]
  
  # Find the corresponding row(s) in the cased_interactions_filtered_2 dataset (with both possible orderings)
  match_rows <- (cased_interactions_filtered_2$Gn1 == gn1 & cased_interactions_filtered_2$Gn2 == gn2) |
    (cased_interactions_filtered_2$Gn1 == gn2 & cased_interactions_filtered_2$Gn2 == gn1)
  
  if (any(match_rows)) {
    # If matching row(s) are found, hyphenate the interaction types into a single string
    interaction_types <- paste(cased_interactions_filtered_2$interaction_type[match_rows], collapse = "-")
  } else {
    # If no matching row is found, assign NA
    interaction_types <- NA
  }
  
  # Assign the interaction types to the corresponding index in the vector
  interaction_types_detailed[i] <- interaction_types
}

# Assign the interaction types to the log_change_interaction dataset as a new column
log_change_interaction$interaction_types_detailed <- interaction_types_detailed


