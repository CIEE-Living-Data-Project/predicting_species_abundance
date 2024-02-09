# Date created: 26 Apr 2023
# Date updated: 26 Apr 2023 (NC)
# Date updated: 29 Jan 2024 (ENB)

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
log.prop.change.with.meta <-read.csv("~/Documents/Work and Career/LDP/Working Group/within.study.updated.data.050224.csv")
log_change <- log.prop.change.with.meta

head(log.prop.change.with.meta)

# get list of unique bio pair genera 
genus_pairs <- apply(log_change[, c("Gn1", "Gn2")], 1, function(x) paste(sort(x), collapse = "-"))
unique_pairs <- unique(genus_pairs)
split_pairs <- data.frame(Gn1 = sapply(strsplit(unique_pairs, "-"), "[", 1),
                          Gn2 = sapply(strsplit(unique_pairs, "-"), "[", 2))

#Since the genera pairs are very large, they will crash Globi as is
#So, I split them into twentyfour pairs of roughly equal size, so
#we can split the globi processing up 

# Split into 24 groups of roughly equal size
groups <- cut(seq_along(1:nrow(split_pairs)), breaks = 24, labels = FALSE)

# Turn the splits into data frames
df_split <- split(split_pairs, groups)
names(df_split) <- paste0("split_", names(df_split))

#Bring the splits into the environment
list2env(df_split, envir = .GlobalEnv)

#Make Isaac's progress bar 
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
write.csv(pair_interactions_1, "data/preprocessing/pair_interactions_1.csv")




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
write.csv(pair_interactions_2, "data/preprocessing/pair_interactions_2.csv")



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
write.csv(pair_interactions_3, "data/preprocessing/pair_interactions_3.csv")


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
write.csv(pair_interactions_4, "data/preprocessing/pair_interactions_4.csv")




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
write.csv(pair_interactions_5, "data/preprocessing/pair_interactions_5.csv")


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
write.csv(pair_interactions_6, "data/preprocessing/pair_interactions_6.csv")

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
write.csv(pair_interactions_7, "data/preprocessing/pair_interactions_7.csv")



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
write.csv(pair_interactions_8, "data/preprocessing/pair_interactions_8.csv")


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
write.csv(pair_interactions_9, "data/preprocessing/pair_interactions_9.csv")



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
write.csv(pair_interactions_10, "data/preprocessing/pair_interactions_10.csv")

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
write.csv(pair_interactions_11, "data/preprocessing/pair_interactions_11.csv")


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
write.csv(pair_interactions_12, "data/preprocessing/pair_interactions_12.csv")


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
write.csv(pair_interactions_13, "data/preprocessing/pair_interactions_13.csv")



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
write.csv(pair_interactions_14, "data/preprocessing/pair_interactions_14.csv")


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
write.csv(pair_interactions_15, "data/preprocessing/pair_interactions_15.csv")



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
write.csv(pair_interactions_16, "data/preprocessing/pair_interactions_16.csv")


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
write.csv(pair_interactions_17, "data/preprocessing/pair_interactions_17.csv")


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
write.csv(pair_interactions_18, "data/preprocessing/pair_interactions_18.csv")


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
write.csv(pair_interactions_19, "data/preprocessing/pair_interactions_19.csv")


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
write.csv(pair_interactions_20, "data/preprocessing/pair_interactions_20.csv")


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
write.csv(pair_interactions_21, "data/preprocessing/pair_interactions_21.csv")


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
write.csv(pair_interactions_22, "data/preprocessing/pair_interactions_22.csv")


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
write.csv(pair_interactions_23, "data/preprocessing/pair_interactions_23.csv")


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
write.csv(pair_interactions_24, "data/preprocessing/pair_interactions_24.csv")




#Bind the results of all five loops together 
all_interactions <- rbind(pair_interactions_1, pair_interactions_2, pair_interactions_3, 
                          pair_interactions_4, pair_interactions_5 ,pair_interactions_6, 
                          pair_interactions_7, pair_interactions_8, pair_interactions_9, pair_interactions_10, 
                          pair_interactions_11, pair_interactions_12, pair_interactions_13, 
                          pair_interactions_14, pair_interactions_15, pair_interactions_16, 
                          pair_interactions_17, pair_interactions_18, pair_interactions_19, 
                          pair_interactions_20, pair_interactions_21, pair_interactions_22, 
                          pair_interactions_23, pair_interactions_24)
colnames(all_interactions) <- c("Gn1", "Gn2", "interaction")
write.csv(all_interactions, "data/preprocessing/genus_interaction_list_ENB_012924.csv")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_


##START HERE IF YOU JUST WANT TO WORK WITH THE DATA##
##Note: all removal of mistake interactions below has already been done on the uploaded
#dataset

#Read in the completed dataset from above
all_interactions <- read.csv("data/preprocessing/genus_interaction_list_ENB_012924.csv")
all_interactions <- all_interactions %>%
  select(!X)

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
#In general, these appear to be implausible
#EXCEPT in the case of the brown headed cowbird, Molothrus
#so we will keep those interactions and remove all others 
all_interactions_2<- all_interactions %>%
  filter(!interaction %in% mistake_interactions)
#put molothrus rows back in
molothrus <- mistake_rows %>%
  filter(Gn1 == "Molothrus" | Gn2 == "Molothrus")
all_interactions_2 <- rbind(all_interactions_2, molothrus)

#Get number o
summary_interactions <- all_interactions_2 %>%
  group_by(Gn1, Gn2) %>%
  summarize(n=n()) %>%
  arrange(desc(n))

#unique rows summary_interactions: 
unique_interactions_1 <-summary_interactions %>%
  select(Gn1, Gn2) %>%
  distinct()


#Next steps: start collapsing interaction types into larger groups, 
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
#If the interactions contain "mutualism" and "predator_prey", clearly this is a dispersal
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

#Did not get down to 1, but given we only want presence/absence it's fine for now 

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Adding interactions into the log.prop.change dataset
#interaction info being added: 
  #Interaction Yes (1)/No(0)
  #Positive/negative/neutral
  #Type of interaction


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#Yes/no interaction
#Get pairs 
distinct_pairs <- cased_interactions_filtered_3 %>%
  select(Gn1, Gn2) %>%
  distinct()

#Read in interactions 
log_change <- read.csv("~/Documents/Work and Career/LDP/Working Group/within.study.updated.data.050224.csv")
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

saveRDS(log_change, "data/data_processing/log.prop.change.interactions.020724ENB.RDS")
write.csv(log_change, "data/data_processing/within.study.updated.interactions.020724ENB.csv")

saveRDS(cased_interactions_filtered_3, "data/preprocessing/all.interactions.genus.pairs.012924.ENB.RDS")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#Positive negative neutral (hyphenated)

log_change_interaction <- readRDS("data/data_processing/log.prop.change.interactions.RDS")


pos_neg_interactions<- cased_interactions_filtered_3%>%
  mutate(
    interaction_benefit = case_when(
      interaction_type %in%  c("predator_prey", "parasitism") ~ "negative",
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


saveRDS(log_change_interaction, "data/data_processing/log.prop.change.interactions.RDS")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Write new log change interaction: 
log_change_interaction <- readRDS("data/data_processing/log.prop.change.interactions.RDS")


#do interaction type detailed
pb<-set.prog.bar(nrow(log_change_interaction)) #sets progress bar
#take interaction types, and hyphenate if there are multiple
for (i in 1:nrow(log_change_interaction)) {
  pb$tick()
  
  # Get the genus pair from the current row
  gn1 <- log_change_interaction$Gn1[i]
  gn2 <- log_change_interaction$Gn2[i]
  
  # Find the corresponding row(s) in the cased_interactions_filtered_3 dataset
  match_rows <- (cased_interactions_filtered_3$Gn1 == gn1 & cased_interactions_filtered_3$Gn2 == gn2) |
    (cased_interactions_filtered_3$Gn1 == gn2 & cased_interactions_filtered_3$Gn2 == gn1)
  
  if (any(match_rows)) {
    # If matching row(s) are found, retrieve interaction types
    interaction_type <- cased_interactions_filtered_3$interaction_type[match_rows]
    
    # Hyphenate multiple interaction types or assign "NA" if none
    if (length(interaction_type) > 0) {
      hyphenated_interaction <- paste(interaction_type, collapse = "-")
    } else {
      hyphenated_interaction <- "NA"
    }
    
    # Assign the hyphenated interaction to the corresponding column
    log_change_interaction$interaction_type[i] <- hyphenated_interaction
  } else {
    # If no matching row is found, assign "NA" to the interaction type
    log_change_interaction$interaction_type[i] <- "NA"
  }
}

#checking everything worked
unique(log_change_interaction$interaction_type)

#remove the "detailed" column
log_change_interaction$interaction_benefit <- log_change_interaction_hyphenated$interaction_benefit

saveRDS(log_change_interaction, "data/data_processing/log.prop.change.interactions.RDS")


#non-na rows with interactions
summarize_rows <- log_change_interaction %>%
  filter(interaction_present==1)
nrow(summarize_rows)/nrow(log_change_interaction)
unique(log_change_interaction$Gn1, log_change_interaction$Gn2)