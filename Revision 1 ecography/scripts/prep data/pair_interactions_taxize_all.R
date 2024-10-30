# Date created: 26 Apr 2023
# Date updated: 26 Apr 2023 (NC)
# Date updated: 18 Sept 2024 (ENB)

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

#Step 1. 
    #Assign interactions to the list of unique genus pairs from the modelling dataset


#read in the full dataset used for Q1 and Q2 analysis
log.prop.change.with.meta <-read.csv("Revision 1 ecography/output/prep_data/results.abundance.csv")
log_change <- log.prop.change.with.meta

head(log.prop.change.with.meta)

# get list of unique bio pair genera 
genus_pairs <- apply(log_change[, c("Gn1", "Gn2")], 1, function(x) paste(sort(x), collapse = "-"))
unique_pairs <- unique(genus_pairs)
split_pairs <- data.frame(Gn1 = sapply(strsplit(unique_pairs, "-"), "[", 1),
                          Gn2 = sapply(strsplit(unique_pairs, "-"), "[", 2))

#Since the genera pairs are very large, they will crash Globi as is
#So, split them into twentyfour pairs of roughly equal size, so
#we can split the globi processing up 

# Split into 24 groups of roughly equal size
groups <- cut(seq_along(1:nrow(split_pairs)), breaks = 24, labels = FALSE)

# Turn the splits into data frames
df_split <- split(split_pairs, groups)
names(df_split) <- paste0("split_", names(df_split))

#Bring the splits into the environment
list2env(df_split, envir = .GlobalEnv)

#Make progress bar 
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

#Create a blank list for the for loop below 
df_list <- list()

#Loop for assigning genera interactions: loop 1
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



#Loop for assigning genera interactions: loop 2
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


#Loop for assigning genera interactions: loop 3
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

#Loop for assigning genera interactions: loop 4
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



#Loop for assigning genera interactions: loop 5
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

#Loop for assigning genera interactions: loop 6
pb<-set.prog.bar(nrow(split_6)) #sets progress bar
for (i in 1:nrow(split_6)) {
  pb$tick()
  gn1 <- split_6[i, "Gn1"]
  gn2 <- split_6[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_6 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 7
pb<-set.prog.bar(nrow(split_7)) #sets progress bar
for (i in 1:nrow(split_7)) {
  pb$tick()
  gn1 <- split_7[i, "Gn1"]
  gn2 <- split_7[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_7 <- do.call(rbind, df_list)


#Loop for assigning genera interactions: loop 8
pb<-set.prog.bar(nrow(split_8)) #sets progress bar
for (i in 1:nrow(split_8)) {
  pb$tick()
  gn1 <- split_8[i, "Gn1"]
  gn2 <- split_8[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_8 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 9
pb<-set.prog.bar(nrow(split_9)) #sets progress bar
for (i in 1:nrow(split_9)) {
  pb$tick()
  gn1 <- split_9[i, "Gn1"]
  gn2 <- split_9[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_9 <- do.call(rbind, df_list)


#Loop for assigning genera interactions: loop 10
pb<-set.prog.bar(nrow(split_10)) #sets progress bar
for (i in 1:nrow(split_10)) {
  pb$tick()
  gn1 <- split_10[i, "Gn1"]
  gn2 <- split_10[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_10 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 11
pb<-set.prog.bar(nrow(split_11)) #sets progress bar
for (i in 1:nrow(split_11)) {
  pb$tick()
  gn1 <- split_11[i, "Gn1"]
  gn2 <- split_11[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_11 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop interactions_12
pb<-set.prog.bar(nrow(split_12)) #sets progress bar
for (i in 1:nrow(split_12)) {
  pb$tick()
  gn1 <- split_12[i, "Gn1"]
  gn2 <- split_12[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_12 <- do.call(rbind, df_list)


#Loop for assigning genera interactions: loop 13
pb<-set.prog.bar(nrow(split_13)) #sets progress bar
for (i in 1:nrow(split_13)) {
  pb$tick()
  gn1 <- split_13[i, "Gn1"]
  gn2 <- split_13[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_13 <- do.call(rbind, df_list)


#Loop for assigning genera interactions: loop interactions_14
pb<-set.prog.bar(nrow(split_14)) #sets progress bar
for (i in 1:nrow(split_14)) {
  pb$tick()
  gn1 <- split_14[i, "Gn1"]
  gn2 <- split_14[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_14 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 15
pb<-set.prog.bar(nrow(split_15)) #sets progress bar
for (i in 1:nrow(split_15)) {
  pb$tick()
  gn1 <- split_15[i, "Gn1"]
  gn2 <- split_15[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_15 <- do.call(rbind, df_list)


#Loop for assigning genera interactions: loop interactions_16
pb<-set.prog.bar(nrow(split_16)) #sets progress bar
for (i in 1:nrow(split_16)) {
  pb$tick()
  gn1 <- split_16[i, "Gn1"]
  gn2 <- split_16[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_16 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 17
pb<-set.prog.bar(nrow(split_17)) #sets progress bar
for (i in 1:nrow(split_17)) {
  pb$tick()
  gn1 <- split_17[i, "Gn1"]
  gn2 <- split_17[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_17 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop interactions_18
pb<-set.prog.bar(nrow(split_18)) #sets progress bar
for (i in 1:nrow(split_18)) {
  pb$tick()
  gn1 <- split_18[i, "Gn1"]
  gn2 <- split_18[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_18 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 19
pb<-set.prog.bar(nrow(split_19)) #sets progress bar
for (i in 1:nrow(split_19)) {
  pb$tick()
  gn1 <- split_19[i, "Gn1"]
  gn2 <- split_19[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_19 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop interactions_20
pb<-set.prog.bar(nrow(split_20)) #sets progress bar
for (i in 1:nrow(split_20)) {
  pb$tick()
  gn1 <- split_20[i, "Gn1"]
  gn2 <- split_20[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_20 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 21
pb<-set.prog.bar(nrow(split_21)) #sets progress bar
for (i in 1:nrow(split_21)) {
  pb$tick()
  gn1 <- split_21[i, "Gn1"]
  gn2 <- split_21[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_21 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop interactions_22
pb<-set.prog.bar(nrow(split_22)) #sets progress bar
for (i in 1:nrow(split_22)) {
  pb$tick()
  gn1 <- split_22[i, "Gn1"]
  gn2 <- split_22[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_22 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop 23
pb<-set.prog.bar(nrow(split_23)) #sets progress bar
for (i in 1:nrow(split_23)) {
  pb$tick()
  gn1 <- split_23[i, "Gn1"]
  gn2 <- split_23[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_23 <- do.call(rbind, df_list)

#Loop for assigning genera interactions: loop interactions_24
pb<-set.prog.bar(nrow(split_24)) #sets progress bar
for (i in 1:nrow(split_24)) {
  pb$tick()
  gn1 <- split_24[i, "Gn1"]
  gn2 <- split_24[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_24 <- do.call(rbind, df_list)



#Bind the results of all loops together 
all_interactions <- rbind(pair_interactions_1, pair_interactions_2, pair_interactions_3, 
                          pair_interactions_4, pair_interactions_5 ,pair_interactions_6, 
                          pair_interactions_7, pair_interactions_8, pair_interactions_9, pair_interactions_10, 
                          pair_interactions_11, pair_interactions_12, pair_interactions_13, 
                          pair_interactions_14, pair_interactions_15, pair_interactions_16, 
                          pair_interactions_17, pair_interactions_18, pair_interactions_19, 
                          pair_interactions_20, pair_interactions_21, pair_interactions_22, 
                          pair_interactions_23, pair_interactions_24)
colnames(all_interactions) <- c("Gn1", "Gn2", "interaction")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#manually inspect the types of interactions
unique(all_interactions$interaction)
#figure out how many distinct pairs have an interaction
all_interactions_pairs <- all_interactions %>%
  select(-interaction)
all_interactions_pairs <- distinct(all_interactions_pairs)

#Upon manual inspection, some interactions are clearly mistakes due to
#ambiguity in the genera
#will manually inspect
#Unfortunately, pollination also seems to be a mistake - it's always plant to plant
mistake_interactions <- c("hasHost", "pathogenOf", "hemiparasiteOf", 'hostOf', 'pathogenOf', 
                          'hasVector', 'parasiteOf','ectoparasiteOf', 'hasParasite', 'hasPathogen', "endoparasiteOf", 'pollinates')
mistake_rows <- all_interactions %>%
  filter(interaction %in% mistake_interactions)

#The brown headed cowbird, Molothrus, has some parasitic interactions
#so we will keep those interactions and remove all others 
all_interactions_2<- all_interactions %>%
  filter(!interaction %in% mistake_interactions)
#put molothrus rows back in
molothrus <- mistake_rows %>%
  filter(Gn1 == "Molothrus" | Gn2 == "Molothrus")
all_interactions_2 <- rbind(all_interactions_2, molothrus)

#Get number of interactions assigned to each genera
summary_interactions <- all_interactions_2 %>%
  group_by(Gn1, Gn2) %>%
  summarize(n=n()) %>%
  arrange(desc(n))

#unique interactions for each row:  
unique_interactions_1 <-summary_interactions %>%
  select(Gn1, Gn2) %>%
  distinct()


#Next steps: collapsing interaction types into larger groups, 
#Key: 
#predator-prey (including herbivory)
  #eatenBy, eats, preysOn, preyedUponBy
#mutualism
  #mutualistOf
#dispersal 
  #hasDispersalVector, dispersalVectorOf
#pollination
  #pollinates
#parasitism
  #hostOf, "hhemiparasiteOf, #hasParasite, #parasiteOf
#uncategorized_interaction
  #interactsWith, "ecologicallyRelatedTo",   flowersVisitedBy, #visits, #visitsFlowersOf, #visitedBy

# Use case_when to classify the columns
cased_interactions<- all_interactions_2%>%
  mutate(
    interaction_type = case_when(
      interaction %in%  c("preysOn", "preyedUponBy", "eatenBy", "eats", "killedBy") ~ "predator_prey",
      interaction %in% c("mutualistOf") ~ "mutualism",
      interaction %in% c("hasDispersalVector", "dispersalVectorOf") ~ "dispersal",
      interaction %in% c("pollinates") ~ "pollination", 
      interaction %in% c("hostOf", "hemiparasiteOf", "hasParasite", "parasiteOf", "hasHost") ~ "parasitism", 
      interaction %in% c("interactsWith", "flowersVisitedBy", "visits", "visitsFlowersOf", "visitedBy", "adjacentTo", 
                         "ecologicallyRelatedTo", "coOccursWith") ~ "uncategorized_interaction",
      TRUE ~ interaction
    ),
    .keep = "unused"
  )
unique(cased_interactions$interaction_type)

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

cased_interactions_filtered_2<- cased_interaction_unique %>%
  group_by(Gn1, Gn2) %>%
  mutate(num_interactions=n())


#These need to be cut down even further
#If the interactions contain "dispersal" and "predator_prey", this is a dispersal
#relationship 
#And, if it contains "uncategorized interaction" and "predator_prey", clearly this is predator prey
#So, we can remove the extra rows 


cased_interactions_filtered_3 <- cased_interactions_filtered_2 %>%
  group_by(Gn1, Gn2) %>%
  filter(
    !("dispersal" %in% interaction_type & n_distinct(interaction_type) > 1) |
      interaction_type == "dispersal"
  ) %>%
  ungroup()

cased_interactions_filtered_3 <- cased_interactions_filtered_3 %>%
  group_by(Gn1, Gn2) %>%
  filter(
    !("predator_prey" %in% interaction_type & "uncategorized_interaction" %in% interaction_type) |
      interaction_type == "predator_prey"
  ) %>%
  ungroup()

cased_interactions_filtered_3 <- cased_interactions_filtered_3 %>%
  group_by(Gn1, Gn2) %>%
  filter(
    !("parasitism" %in% interaction_type & "uncategorized_interaction" %in% interaction_type) |
      interaction_type == "parasitism"
  ) %>%
  ungroup()

cased_interactions_filtered_3 <- cased_interactions_filtered_3 %>%
  group_by(Gn1, Gn2) %>%
  filter(
    !("mutualism" %in% interaction_type & "uncategorized_interaction" %in% interaction_type) |
      interaction_type == "mutualism"
  ) %>%
  ungroup()

#Check for new number of interactions for each pair
cased_interactions_filtered_3<- cased_interactions_filtered_3 %>%
  group_by(Gn1, Gn2) %>%
  mutate(new_num_interactions=n())
table(cased_interactions_filtered_3$new_num_interactions)

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Add interaction categories: 
  #Interaction Yes (1)/No(0)



#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#Yes/no interaction
#Get pairs 
distinct_pairs <- cased_interactions_filtered_3 %>%
  select(Gn1, Gn2) %>%
  distinct()


# Create a new column in log_change called "interaction" and initialize all values to 0
log_change$interaction_present <- 0

# Loop through each row in log_change and check if it matches any pairs in our interaction dataset
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

#All pairs assigned interaction Y/N

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Assign Classes to each pair, using the taxize library

#Load the taxize package
library(taxize)

#Get the class names for each group 
taxa1<-tax_name(query = unique(log_change$Gn1), get = "class", db = "itis", messages = FALSE, callopts=list(http_version = 0))
taxa2<-tax_name(query = unique(log_change$Gn2)[-which(unique(log_change$Gn2)%in%unique(log_change$Gn1))], get = "class", db = "itis", 
                messages = FALSE, callopts=list(http_version = 0))
taxa<-rbind(taxa1,taxa2)

#add new taxa col
taxa$Gn1=taxa$query
taxa$Gn2=taxa$query
names(taxa)[which(names(taxa)=="class")]<-"RESOLVED.TAXA1"
taxa$RESOLVED.TAXA2=taxa$RESOLVED.TAXA1

log_change.taxa<-left_join(log_change,taxa[,c("Gn1","RESOLVED.TAXA1")])
log_change.taxa<-left_join(log_change.taxa,taxa[,c("Gn2","RESOLVED.TAXA2")])


#Generate a resolved_taxa_pair column to match other datasets
sorted_words <- apply(log_change.taxa[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
log_change.taxa$resolved_taxa_pair <- sorted_words
unique(log_change.taxa$resolved_taxa_pair)



#When we do this, we see some erroneous NA names 
#E.g. no marine datasets but see chondrichytes
#We can correct these since they are labelled in the BioTIME metadata

#Investigate Acentrosomata
acentro <- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Acentrosomata")
unique(acentro$Gn1)
unique(acentro$STUDY_PLOT)
#Replace all Acentrosomata with Insecta (from metadata)
log_change.taxa$RESOLVED.TAXA1 <- gsub("Acentrosomata", "Insecta", log_change.taxa$RESOLVED.TAXA1)
log_change.taxa$RESOLVED.TAXA2 <- gsub("Acentrosomata", "Insecta", log_change.taxa$RESOLVED.TAXA2)


#Investigate euchelicerata
euchel <- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Euchelicerata")
unique(euchel$Gn1)
unique(euchel$STUDY_PLOT)
#Looks like there are some erronerous names: WIWA, 
#And a mislabelled moth genus - Dioryctria 
#Fix this 
#WIWA should be aves based on metadata 
log_change.taxa$RESOLVED.TAXA1 <- ifelse(log_change.taxa$Gn1 == "WIWA", 
       log_change.taxa$RESOLVED.TAXA1 <- gsub("Euchelicerata", "Aves", log_change.taxa$RESOLVED.TAXA1),
       log_change.taxa$RESOLVED.TAXA1)
log_change.taxa$RESOLVED.TAXA2 <- ifelse(log_change.taxa$Gn2 == "WIWA", 
                                         log_change.taxa$RESOLVED.TAXA2 <- gsub("Euchelicerata", "Aves", log_change.taxa$RESOLVED.TAXA2),
                                         log_change.taxa$RESOLVED.TAXA2)

#These should be insecta, as per the metadata
log_change.taxa$RESOLVED.TAXA1 <- ifelse(log_change.taxa$Gn1 == "Dioryctria", 
                                         log_change.taxa$RESOLVED.TAXA1 <- gsub("Euchelicerata", "Insecta", log_change.taxa$RESOLVED.TAXA1),
                                         log_change.taxa$RESOLVED.TAXA1)
log_change.taxa$RESOLVED.TAXA2 <- ifelse(log_change.taxa$Gn2 == "Dioryctria", 
                                         log_change.taxa$RESOLVED.TAXA2 <- gsub("Euchelicerata", "Insecta", log_change.taxa$RESOLVED.TAXA2),
                                         log_change.taxa$RESOLVED.TAXA2)



#Investigate chondrichyes 
chondro <- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Chondrichthyes")
unique(chondro$Gn1)
unique(chondro$STUDY_PLOT)
#Looks like one mislabelled genus - Prays 
#What else is in this study? 
study_249 <- log_change.taxa[grepl("249", log_change.taxa$STUDY_PLOT), ]
unique(study_249$Gn1)
#Prays - is a lepidopteran according to metadata
#Relabel as Insecta
log_change.taxa$RESOLVED.TAXA1 <- gsub("Chondrichthyes", "Insecta", log_change.taxa$RESOLVED.TAXA1)
log_change.taxa$RESOLVED.TAXA2 <- gsub("Chondrichthyes", "Insecta", log_change.taxa$RESOLVED.TAXA2)

#Next, double check Reptilia
repts<- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Reptilia")
unique(repts$Gn1)
unique(repts$STUDY_PLOT)
#319 and 316 are fine - 225 needs to be changed to Aves according to metadata
log_change.taxa$RESOLVED.TAXA2 [grepl("225", log_change.taxa$STUDY_PLOT )] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("225", log_change.taxa$STUDY_PLOT)] <- "Aves"



#Check Teleostei 
teleost<- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Teleostei")
unique(teleost$Gn1)
unique(teleost$STUDY_PLOT)
#Study 225 should be all birds according to metadata
log_change.taxa$RESOLVED.TAXA2 [grepl("225", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("225", log_change.taxa$STUDY_PLOT)] <- "Aves"


#Check gastropoda 
gastropod<- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Gastropoda")
unique(gastropod$Gn1)
unique(gastropod$STUDY_PLOT)
#Study 382 should be all mammals according to metadata
log_change.taxa$RESOLVED.TAXA2 [grepl("382", log_change.taxa$STUDY_PLOT )] <- "Mammalia"
log_change.taxa$RESOLVED.TAXA1 [grepl("382", log_change.taxa$STUDY_PLOT)] <- "Mammalia"

#Check Gnetopsida 
gnetops<- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Gnetopsida")
unique(gnetops$Gn1)
unique(gnetops$STUDY_PLOT)
#Looks ok! 

#Check bryopsida 
Bryops<- log_change.taxa %>%
  filter(RESOLVED.TAXA1 =="Bryopsida")
unique(Bryops$Gn1)
unique(Bryops$STUDY_PLOT)
#This should be birds according to metadata
#Study 333
log_change.taxa$RESOLVED.TAXA2 [grepl("333", log_change.taxa$STUDY_PLOT )] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("333", log_change.taxa$STUDY_PLOT)] <- "Aves"

#Create new resolved taxa list 
sorted_words <- apply(log_change.taxa[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
log_change.taxa$resolved_taxa_pair <- sorted_words
table(log_change.taxa$resolved_taxa_pair)



#Check rows with NA for Class 
na_classes <-  log_change.taxa[grepl("NA", log_change.taxa$resolved_taxa_pair), ]
#show for which study it is NA
unique(na_classes$STUDY_PLOT)
#Seeing study 195, 221, 240, 248, 249, 300, 308, 
#340, 366, 375, 380, 39, 413, 414, 416, 420, 46, 471, 
#54

#Study 195 - should all be birds according to metadata
log_change.taxa$RESOLVED.TAXA2 [grepl("195", log_change.taxa$STUDY_PLOT )] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("195", log_change.taxa$STUDY_PLOT)] <- "Aves"

#Study 221 - boreal vegetation, so we're gonna have to dive deeper
#NOTE: study 221 was ultimately corrected in results.abundance code
#As genus/species names were incorrect 
#And then subsequently removed from analyses due to inconsistencies
study_221_na <- log_change %>%
  filter(STUDY_PLOT == "221~Was_NA") 
unique(study_221_na$Gn1)


#Study 240
#Terrestrial plants, need to investigate 
study_240_na <- log_change.taxa %>%
  filter(STUDY_PLOT == "240~Was_NA") %>%
  filter(is.na(RESOLVED.TAXA1))
unique(study_240_na$Gn1)
#Should be Magnoliopsida and Monocot
#Will manually correct
log_change.taxa$RESOLVED.TAXA1 [grepl("Chamaesyce", log_change.taxa$Gn1)] <- "Magnoliopsida"
log_change.taxa$RESOLVED.TAXA2 [grepl("Chamaesyce", log_change.taxa$Gn2)] <- "Magnoliopsida"

log_change.taxa$RESOLVED.TAXA1 [grepl("Lesquerella", log_change.taxa$Gn1)] <- "Magnoliopsida"
log_change.taxa$RESOLVED.TAXA2 [grepl("Lesquerella", log_change.taxa$Gn2)] <- "Magnoliopsida"

log_change.taxa$RESOLVED.TAXA1 [grepl("Pleuraphis", log_change.taxa$Gn1)] <- "Monocots"
log_change.taxa$RESOLVED.TAXA2 [grepl("Pleuraphis", log_change.taxa$Gn2)] <- "Monocots"

#Study 248 
#Plants again - need to investigate further
study_248_na <- log_change.taxa %>%
  filter(grepl("248", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
unique(study_248_na$Gn1)
#This is a eudicot 
log_change.taxa$RESOLVED.TAXA1 [grepl("Stenosiphon", log_change.taxa$Gn1)] <- "Eudicots"
log_change.taxa$RESOLVED.TAXA2 [grepl("Stenosiphon", log_change.taxa$Gn2)] <- "Eudicots"

#Study 249
study_249_na <- log_change.taxa %>%
  filter(grepl("249", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
unique(study_249_na$Gn1)
#Lepidopterans - update to Insecta
log_change.taxa$RESOLVED.TAXA1[grepl("249", log_change.taxa$STUDY_PLOT) & log_change.taxa$RESOLVED.TAXA1 == "NA"] <- "Insecta"
log_change.taxa$RESOLVED.TAXA2[grepl("249", log_change.taxa$STUDY_PLOT) & log_change.taxa$RESOLVED.TAXA2 == "NA"] <- "Insecta"


#What is monocots insecta class group? 
monocots_insecta <- log_change.taxa %>%
  filter(resolved_taxa_pair=="Monocots.Insecta")
unique(monocots_insecta$STUDY_PLOT)
#This also needs fixing - these should all be Monocots
unique(monocots_insecta$Gn2)
log_change.taxa$RESOLVED.TAXA2[grepl("Scleropogon", log_change.taxa$Gn2)] <- "Monocots"
log_change.taxa$RESOLVED.TAXA1[grepl("Scleropogon", log_change.taxa$Gn1)] <- "Monocots"


#Study 300
study_300_na <- log_change.taxa %>%
  filter(grepl("300", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
unique(study_300_na$Gn1)
#These are lepidopterans - should be updated to Insecta
log_change.taxa$RESOLVED.TAXA1[grepl("Harmonia", log_change.taxa$Gn1)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA2[grepl("Harmonia", log_change.taxa$Gn2)] <- "Insecta"

#Study 308
study_308_na <- log_change.taxa %>%
  filter(grepl("308", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#Should all be mammalia 
log_change.taxa$RESOLVED.TAXA1[grepl("308", log_change.taxa$STUDY_PLOT)] <- "Mammalia"
log_change.taxa$RESOLVED.TAXA2[grepl("308", log_change.taxa$STUDY_PLOT)] <- "Mammalia"

#Study 340
study_340_na <- log_change.taxa %>%
  filter(grepl("340", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#Plants 
unique(study_340_na$Gn1)
#Looks like we're good here

#Study 366
study_366_na <- log_change.taxa %>%
  filter(grepl("366", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#Sould be updated to mammals
log_change.taxa$RESOLVED.TAXA1[grepl("366", log_change.taxa$STUDY_PLOT)] <- "Mammalia"
log_change.taxa$RESOLVED.TAXA2[grepl("366", log_change.taxa$STUDY_PLOT)] <- "Mammalia"

#Study 375 
study_375_na <- log_change.taxa %>%
  filter(grepl("375", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#Probably beetles (Insecta), but let's check 
unique(study_375_na$Gn1)
#All beetles - update to Insecta based on metadata
log_change.taxa$RESOLVED.TAXA1[grepl("Phelotrupes", log_change.taxa$Gn1)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA2[grepl("Phelotrupes", log_change.taxa$Gn2)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA1[grepl("Platydracus", log_change.taxa$Gn1)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA2[grepl("Platydracus", log_change.taxa$Gn2)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA1[grepl("Staphylinus", log_change.taxa$Gn1)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA2[grepl("Staphylinus", log_change.taxa$Gn2)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA1[grepl("Eusilpha", log_change.taxa$Gn1)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA2[grepl("Eusilpha", log_change.taxa$Gn2)] <- "Insecta"

#Studies remaining: 
#380, 39, 413, 414, 416, 420, 46, 471, 
#54

#Study 380
study_380_na <- log_change.taxa %>%
  filter(grepl("380", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#Probably moths and butterflies, will double check 
butterflies <- unique(study_380_na$Gn1)
butterflies <- paste(butterflies, collapse = "|")
# Use grepl with the combined pattern to update to Insecta based on metadata
log_change.taxa$RESOLVED.TAXA2[grepl(butterflies, log_change.taxa$Gn2)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA1[grepl(butterflies, log_change.taxa$Gn1)] <- "Insecta"

#Study 39
study_39_na <- log_change.taxa %>%
  filter(grepl("\\b39\\b", STUDY_PLOT) & is.na(RESOLVED.TAXA1))
log_change.taxa$RESOLVED.TAXA1[grepl("\\b39\\b", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2[grepl("\\b39\\b", log_change.taxa$STUDY_PLOT)] <- "Aves"
#NOTE: this study was corrected in new results.abundance code
#And then subsequently removed in analyses due to inconsistencies in genus names 

#Study 413
study_413_na <- log_change.taxa %>%
  filter(grepl("413", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#should all be birds from metadata- update to Aves
log_change.taxa$RESOLVED.TAXA1[grepl("413", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2[grepl("413", log_change.taxa$STUDY_PLOT)] <- "Aves"

#414
study_414_na <- log_change.taxa %>%
  filter(grepl("414", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#Also all birds from metadata- update to Aves
log_change.taxa$RESOLVED.TAXA1[grepl("414", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2[grepl("414", log_change.taxa$STUDY_PLOT)] <- "Aves"

#416
study_416_na <- log_change.taxa %>%
  filter(grepl("416", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#More birds from metadata- update to Aves
log_change.taxa$RESOLVED.TAXA1[grepl("416", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2[grepl("416", log_change.taxa$STUDY_PLOT)] <- "Aves"

#420
study_420_na <- log_change.taxa %>%
  filter(grepl("420", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#More birds from metadata- update to Aves
log_change.taxa$RESOLVED.TAXA1[grepl("420", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2[grepl("420", log_change.taxa$STUDY_PLOT)] <- "Aves"

#46
study_46_na <- log_change.taxa %>%
  filter(grepl("\\b46\\b", STUDY_PLOT) & is.na(RESOLVED.TAXA1))
#More birds from metadata- update to Aves
log_change.taxa$RESOLVED.TAXA1[grepl("\\b46\\b", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2[grepl("\\b46\\b", log_change.taxa$STUDY_PLOT)] <- "Aves"

#471
study_471_na <- log_change.taxa %>%
  filter(grepl("471", STUDY_PLOT)) %>%
  filter(is.na(RESOLVED.TAXA1))
#Grasses from metadata - need to double check all is resolved
unique(study_471_na$Log.prop.change.Gn2)
#Looks like it all got sorted

#54
#invertebrates and snails
study_54_na <- log_change.taxa %>%
  filter(grepl("\\b54\\b", STUDY_PLOT) & is.na(RESOLVED.TAXA1))
unique(study_54_na$Gn1)
#Should be Gastropod from metadata- update
log_change.taxa$RESOLVED.TAXA1[grepl("Gaeotis", log_change.taxa$Gn1)] <- "Gastropoda"
log_change.taxa$RESOLVED.TAXA2[grepl("Gaeotis", log_change.taxa$Gn2)] <- "Gastropoda"




#Create a new resolved taxa pair column with corrected names 
sorted_words <- apply(log_change.taxa[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
log_change.taxa$resolved_taxa_pair <- sorted_words
table(log_change.taxa$resolved_taxa_pair)

#looks like we still have a few NAs, where are they? 
na_classes <-  log_change.taxa[grepl("NA", log_change.taxa$resolved_taxa_pair), ]
unique(na_classes$STUDY_PLOT)
#249, 54, 375
study_54_na <- log_change.taxa %>%
  filter(grepl("\\b54\\b", STUDY_PLOT))
unique(study_54_na$Gn2)
#Need to update to Gastropoda from metadata
log_change.taxa$RESOLVED.TAXA1[grepl("Nenia", log_change.taxa$Gn1)] <- "Gastropoda"
log_change.taxa$RESOLVED.TAXA2[grepl("Nenia", log_change.taxa$Gn2)] <- "Gastropoda"

#Study 249
study_249_na <- log_change.taxa %>%
  filter(grepl("249", STUDY_PLOT))
unique(study_249_na$Gn2)
#More moths - update to Insecta from metadata
butterflies <- unique(study_249_na$Gn2)
butterflies <- paste(butterflies, collapse = "|")
# Use grepl with the combined pattern
log_change.taxa$RESOLVED.TAXA2[grepl(butterflies, log_change.taxa$Gn2)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA1[grepl(butterflies, log_change.taxa$Gn1)] <- "Insecta"

#Study 375
study_375_na <- log_change.taxa %>%
  filter(grepl("375", STUDY_PLOT))
unique(study_375_na$Gn2)
#This one is tricky, because it has birds and beetles
#So let's do the beetles first - update to Insecta
beetles <-  c("Pterostichus", "Silpha", "Phelotrupes", "Staphylinus", "Synuchus", 
                            "Platydracus", "Coleoptera", "Eusilpha", "Drusilla", "Panelus")
beetles <- paste(beetles, collapse = "|")
# Use grepl with the combined pattern
log_change.taxa$RESOLVED.TAXA2[grepl(beetles, log_change.taxa$Gn2)] <- "Insecta"
log_change.taxa$RESOLVED.TAXA1[grepl(beetles, log_change.taxa$Gn1)] <- "Insecta"
#Birds - update to Aves
birds <- c("Columba", "Corvus", "Cyanocitta", "Eremophila", "Geothlypis", 
                          "Hirundo", "Melanerpes", "Melospiza", "Molothrus", "Passer", 
                          "Passerina", "Phasianus", "Quiscalus", "Spiza", "Spizella", 
                          "Sturnella", "Sturnus", "Toxostoma", "Turdus", "Tyrannus", 
                          "Zenaida")
birds <- paste(birds, collapse = "|")
log_change.taxa$RESOLVED.TAXA2[grepl(birds, log_change.taxa$Gn2)] <- "Aves"
log_change.taxa$RESOLVED.TAXA1[grepl(birds, log_change.taxa$Gn1)] <- "Aves"

#Double checking aves mammalia: 
aves_mammalia <- log_change.taxa %>%
  filter(resolved_taxa_pair=="Aves.Mammalia")
#Study 195 should all be birds according to metadata
log_change.taxa$RESOLVED.TAXA2 [grepl("195", log_change.taxa$STUDY_PLOT )] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("195", log_change.taxa$STUDY_PLOT)] <- "Aves"

#Double check - why is aves.magnoliopsida so overrepresented? 
aves_aves <- log_change.taxa %>%
  filter(resolved_taxa_pair=="Aves.Magnoliopsida")
table(aves_aves$STUDY_PLOT)
#Fix 363, 339, 475
#All birds 
log_change.taxa$RESOLVED.TAXA2 [grepl("363", log_change.taxa$STUDY_PLOT )] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("363", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2 [grepl("339", log_change.taxa$STUDY_PLOT )] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("339", log_change.taxa$STUDY_PLOT)] <- "Aves"
log_change.taxa$RESOLVED.TAXA2 [grepl("475", log_change.taxa$STUDY_PLOT )] <- "Aves"
log_change.taxa$RESOLVED.TAXA1 [grepl("475", log_change.taxa$STUDY_PLOT)] <- "Aves"
#It seems we just have many large bird studies in our dataset

#Correcting some erroneous Insecta.Magnoliopsida
insecta_magno <- log_change.taxa %>%
  filter(resolved_taxa_pair=="Insecta.Magnoliopsida")
table(insecta_magno$STUDY_PLOT)
#336 should be all plants, what's going on 
unique(insecta_magno$Gn1)
#Boerhaavia should be Eudicot 
log_change.taxa$RESOLVED.TAXA1[grepl("Boerhaavia", log_change.taxa$Gn1)] <- "Eudicots"

#Double checking Aves.Amphibia
amphibia_aves <- log_change.taxa %>%
  filter(grepl("173", STUDY_PLOT))
table(amphibia_aves$STUDY_PLOT)
unique(amphibia_aves$Gn1)
#looks to be correct 

#Fixing some errors from study 173
#Some should be Reptilia
repts <- c("Coluber", "Diadophis", "Eumeces", "Tantilla", "Cnemidophorus", "Uta", "Thamnophis")
repts <- paste(repts, collapse = "|")
log_change.taxa$RESOLVED.TAXA2[grepl(repts, log_change.taxa$Gn2)] <- "Reptilia"
log_change.taxa$RESOLVED.TAXA1[grepl(repts, log_change.taxa$Gn1)] <- "Reptilia"


#Check that it's all resolved: 
sorted_words <- apply(log_change.taxa[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
log_change.taxa$resolved_taxa_pair <- sorted_words
table(log_change.taxa$resolved_taxa_pair)
#All fixed! (Minus study 221, which was removed from analysis)

#Export 
saveRDS(log_change.taxa, "Revision 1 ecography/output/results_abundance_interactions_taxa_032024ENB.RDS")



