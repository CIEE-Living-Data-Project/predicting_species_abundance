# Date created: 26 Apr 2023
# Date updated: 26 Apr 2023 (NC)

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
library(progress)


rm(list=ls()) 


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

log.prop.change.with.meta <- readRDS("data/log.prop.change.with.meta.RDS")
log_change <- log.prop.change.with.meta

head(log.prop.change.with.meta)

#list of bio pairs with at least 10 overlapping years
genus_pairs <- apply(log_change[, c("Gn1", "Gn2")], 1, function(x) paste(sort(x), collapse = "-"))
unique_pairs <- unique(genus_pairs)
split_pairs <- data.frame(Gn1 = sapply(strsplit(unique_pairs, "-"), "[", 1),
                          Gn2 = sapply(strsplit(unique_pairs, "-"), "[", 2))

#short sample dataset: 
# Split into 5 groups of roughly equal size
groups <- cut(seq_along(1:12623), breaks = 5, labels = FALSE)

# Split data into separate data frames
df_split <- split(split_pairs, groups)
names(df_split) <- paste0("split_", names(df_split))


list2env(df_split, envir = .GlobalEnv)


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

#Run for BETWEEN studies

df_list <- list()


pb<-set.prog.bar(nrow(split_1)) #sets progress bar
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
write.csv(pair_interactions_1, "data/pair_interactions_1.csv")





pb<-set.prog.bar(nrow(split_2)) #sets progress bar
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
write.csv(pair_interactions_2, "data/pair_interactions_2.csv")




pb<-set.prog.bar(nrow(split_3)) #sets progress bar
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
write.csv(pair_interactions_3, "data/pair_interactions_3.csv")



pb<-set.prog.bar(nrow(split_4)) #sets progress bar
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
write.csv(pair_interactions_4, "data/pair_interactions_4.csv")





pb<-set.prog.bar(nrow(split_5)) #sets progress bar
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
write.csv(pair_interactions_5, "data/pair_interactions_5.csv")

all_interactions <- rbind(pair_interactions_1, pair_interactions_2, pair_interactions_3, 
                          pair_interactions_4, pair_interactions_5)
colnames(all_interactions) <- c("Gn1", "Gn2", "interaction")
write.csv(all_interactions, "data/genus_interaction_list.csv")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_


##START HERE IF YOU JUST WANT TO WORK WITH THE DATA##
##Note: all removal of mistake interactions below has already been done on the uploaded
#dataset

#Read in the completed dataset from above
all_interactions <- read.csv("data/preprocessing/genus_interaction_list.csv")
all_interactions <- all_interactions %>%
  select(!X)

parasite <- all_interactions %>%
  filter(interaction=="pathogenOf")

#Remove interaction pairs that are mistakes (pathogen, host, hemiparasite
#- see original all interactions above)

mistake_interactions <- c("hasHost", "pathogenOf", "hemiparasiteOf", 'hostOf', 'pathogenOf', 
                          'hasVector', 'parasiteOf')
all_interactions_2<- all_interactions %>%
  filter(!interaction %in% mistake_interactions)
# write.csv(all_interactions_2, "data/genus_interaction_list.csv")

summary_interactions <- all_interactions_2 %>%
  group_by(Gn1, Gn2) %>%
  summarize(n=n()) %>%
  arrange(desc(n))

#unique rows summary_interactions: 
unique_interactions_1 <-summary_interactions %>%
  select(Gn1, Gn2) %>%
  distinct()


#Next steps: start collapsing interaction types into larger groups, 
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


# Use case_when to classify the columns
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

#Let's check how many have multiple interactions still assigned

summary_interactions_cased<- cased_interaction_unique %>%
  group_by(Gn1, Gn2) %>%
  summarize(n=n()) %>%
  arrange() %>%
  left_join(cased_interaction_unique, by=c('Gn1', "Gn2"))


# 
# # Remove rows with "unclassified_interaction" that are associated with any other interaction type for each genus pair
# cased_interactions_filtered <- summary_interactions_cased %>%
# filter(!(n>1 & interaction_type =="uncategorized_interaction"))
# cased_interactions_filtered <- summary_interactions_cased %>%
#   filter(!(n>1 & interaction_type =="visits"))

cased_interactions_filtered_2<- cased_interaction_unique %>%
  group_by(Gn1, Gn2) %>%
  mutate(num_interactions=n())

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Adding interactions into the log.prop.change dataset
#Yes/No
#Positive/negative/neutral
#Type of interaction

#Yes/no interaction
#Get pairs 
distinct_pairs <- cased_interactions_filtered_2 %>%
  select(Gn1, Gn2) %>%
  distinct()

#Read in interactions 
log_change <- readRDS("data/preprocessing/log.prop.change.with.meta.RDS")
head(log_change)

# Create a new column in df2 called "interaction" and initialize all values to 0
log_change$interaction_present <- 0

# Loop through each row in df2 and check if it matches any pairs in df1
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

saveRDS(log_change, "data/preprocessing/log.prop.change.interactions.RDS")
saveRDS(cased_interactions_filtered_2, "data/preprocessing/all.interactions.genus.pairs.RDS")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_


#Do positive negative neutral 
log_change_interaction <- readRDS("data/preprocessing/log.prop.change.interactions.RDS")


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
for (i in 1:nrow(log_change_interaction)) {
  pb$tick()
  
  # Get the genus pair from the current row
  gn1 <- log_change_interaction$Gn1[i]
  gn2 <- log_change_interaction$Gn2[i]
  
  # Find the corresponding row(s) in the pos_neg_interactions dataset
  match_rows <- (pos_neg_interactions$Gn1 == gn1 & pos_neg_interactions$Gn2 == gn2) |
    (pos_neg_interactions$Gn1 == gn2 & pos_neg_interactions$Gn2 == gn1)
  
  if (any(match_rows)) {
    # If matching row(s) are found, check for each interaction type and assign 1 or 0 accordingly
    positive_interaction[i] <- as.integer("positive" %in% pos_neg_interactions$interaction_benefit[match_rows])
    negative_interaction[i] <- as.integer("negative" %in% pos_neg_interactions$interaction_benefit[match_rows])
    neutral_interaction[i] <- as.integer("neutral" %in% pos_neg_interactions$interaction_benefit[match_rows])
  } else {
    # If no matching row is found, assign 0 to all interaction types
    positive_interaction[i] <- 0
    negative_interaction[i] <- 0
    neutral_interaction[i] <- 0
  }
}

log_change_interaction$positive_interaction <- positive_interaction
log_change_interaction$negative_interaction <- negative_interaction
log_change_interaction$neutral_interaction <- neutral_interaction


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#interaction type


# Create empty vectors to store the interaction types
uncategorized_interaction <- vector("integer", nrow(log_change_interaction))
predator_prey_interaction <- vector("integer", nrow(log_change_interaction))
mutualism_interaction <- vector("integer", nrow(log_change_interaction))
dispersal_interaction <- vector("integer", nrow(log_change_interaction))

# Iterate over each row in the log_change_interaction dataset
for (i in 1:nrow(log_change_interaction)) {
  # Get the genus pair from the current row
  gn1 <- log_change_interaction$Gn1[i]
  gn2 <- log_change_interaction$Gn2[i]
  
  # Find the corresponding row(s) in the cased_interactions_filtered_2 dataset (with both possible orderings)
  match_rows <- (cased_interactions_filtered_2$Gn1 == gn1 & cased_interactions_filtered_2$Gn2 == gn2) |
    (cased_interactions_filtered_2$Gn1 == gn2 & cased_interactions_filtered_2$Gn2 == gn1)
  
  if (any(match_rows)) {
    # If matching row(s) are found, check for each interaction type and assign 1 or 0 accordingly
    uncategorized_interaction[i] <- as.integer("uncategorized_interaction" %in% cased_interactions_filtered_2$interaction_type[match_rows])
    predator_prey_interaction[i] <- as.integer("predator_prey" %in% cased_interactions_filtered_2$interaction_type[match_rows])
    mutualism_interaction[i] <- as.integer("mutualism" %in% cased_interactions_filtered_2$interaction_type[match_rows])
    dispersal_interaction[i] <- as.integer("dispersal" %in% cased_interactions_filtered_2$interaction_type[match_rows])
  } else {
    # If no matching row is found, assign 0 to all interaction types
    uncategorized_interaction[i] <- 0
    predator_prey_interaction[i] <- 0
    mutualism_interaction[i] <- 0
    dispersal_interaction[i] <- 0
  }
}

# Assign the interaction types to the log_change_interaction dataset as new columns
log_change_interaction$uncategorized_interaction <- uncategorized_interaction
log_change_interaction$predator_prey_interaction <- predator_prey_interaction
log_change_interaction$mutualism_interaction <- mutualism_interaction
log_change_interaction$dispersal_interaction <- dispersal_interaction


#Write new log change interaction: 
saveRDS(log_change_interaction, "data/preprocessing/log.prop.change.interactions.RDS")
log_change_interaction <- readRDS("data/preprocessing/log.prop.change.interactions.RDS")


interesting_interactions <- log_change_interaction %>%
  filter(predator_prey_interaction==1 & dispersal_interaction==0) %>%
  filter(Gn1=='Achillea', Gn2=="Melanoplus") %>%
  select(Gn1, Gn2, YEAR.T, Log.prop.change.bio.Gn1, Log.prop.change.abun.Gn2)

# Reshape the data into a longer format
data_long <- interesting_interactions %>%
  pivot_longer(cols=starts_with('log.prop'), 
               names_prefix='Log.prop.change.',
               names_to = "Variable", 
               values_to = "Abundance") 


data_long %>%
  ggplot(aes(x=YEAR.T, y=Abundance))+
  geom_point(aes(colour=Variable))+
  theme_classic()+
  geom_smooth(aes(colour=Variable),method='loess', span=0.1)
