# Date created: 26 Apr 2023
# Date updated: 26 Apr 2023 (NC)
# Date updated: 18 Mar 2024 (ENB)

#ONLY FOR STUDIES 221 AND 39

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
log.prop.change.with.meta <-read.csv("~/Documents/Work and Career/LDP/Working Group/results.abundance_221_39.csv")
log_change <- log.prop.change.with.meta

head(log.prop.change.with.meta)

# get list of unique bio pair genera 
genus_pairs <- apply(log_change[, c("Gn1", "Gn2")], 1, function(x) paste(sort(x), collapse = "-"))
unique_pairs <- unique(genus_pairs)
split_pairs <- data.frame(Gn1 = sapply(strsplit(unique_pairs, "-"), "[", 1),
                          Gn2 = sapply(strsplit(unique_pairs, "-"), "[", 2))


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
pb<-set.prog.bar(nrow(split_pairs)) #sets progress bar
for (i in 1:nrow(split_pairs)) {
  pb$tick()
  gn1 <- split_pairs[i, "Gn1"]
  gn2 <- split_pairs[i, "Gn2"]
  interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                                           targettaxon=gn2)
  unique <- unique(interactions$interaction_type)
  new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                       genus_2 = rep(gn2, each = length(unique)))
  pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
  df_list[[i]] <- pairs
}

pair_interactions_1 <- do.call(rbind, df_list)
#write.csv(pair_interactions_1, "data/preprocessing/pair_interactions_1.csv")



colnames(pair_interactions_1) <- c("Gn1", "Gn2", "interaction")
all_interactions <- pair_interactions_1
#write.csv(all_interactions, "data/preprocessing/genus_interaction_list_ENB_012924.csv")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_


##START HERE IF YOU JUST WANT TO WORK WITH THE DATA##
##Note: all removal of mistake interactions below has already been done on the uploaded
#dataset

#Read in the completed dataset from above

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
log_change <- log_change

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

#saveRDS(log_change, "Revision 1 ecography/output/results_abundance_interactions_022824ENB.RDS")
##write.csv(log_change, "data/data_processing/within.study.updated.interactions.020724ENB.csv")

##saveRDS(cased_interactions_filtered_3, "data/preprocessing/all.interactions.genus.pairs.012924.ENB.RDS")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#Add resolved taxa pairs - to class level, and corrected 
#Study 221 - boreal vegetation - all Eudicots 
study_221 <- log_change %>%
  filter(STUDY_ID == "221")

study_221$RESOLVED.TAXA1 <- c("Eudicots")
study_221$RESOLVED.TAXA2 <- c("Eudicots")
sorted_words <- apply(study_221[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
study_221$resolved_taxa_pair <- sorted_words
table(study_221$resolved_taxa_pair)


#Plot study 221 just for fun 
#Create new column with genus pair for labelling
sorted_words <- apply(study_221[, c('Gn1', 'Gn2')], 1, function(x) paste(x, collapse = "."))
study_221$genus_pair <- sorted_words

#Pivot 221 longer
study_221_longer <- study_221 %>%
  pivot_longer(cols = c("Log.prop.change.Gn1", "Log.prop.change.Gn2"), 
               names_to = "genera", 
               values_to = "log_prop_change")


study_221_longer %>%
  ggplot(aes(x = YEAR.T, y = log_prop_change, group=interaction(genus_pair, genera), 
             color = genus_pair, shape = genera)) + 
  geom_point() + 
  geom_jitter()+
  geom_smooth(se=FALSE, span = 0.5, aes(linetype = genera)) +
  facet_wrap(~genus_pair)+
  labs(x="Year", y = "Log(Proportional change in abundance)", 
       shape = "Genera 1 or 2", linetype = "Genera 1 or 2", 
       color = "Genus grouping")+
  theme_classic() 


#Add study 221 back into log_change.taxa
#Remove initial study 221 from log_change
log_change <- log_change %>%
  filter(!STUDY_ID=="221")
log_change <- bind_rows(study_221, log_change) 
#Arrange in order of increasing study 
log_change <- log_change %>%
  arrange(X)

#Same for study 39
study_39 <- log_change %>%
  filter(STUDY_ID == "39")

study_39$RESOLVED.TAXA1 <- c("Aves")
study_39$RESOLVED.TAXA2 <- c("Aves")
sorted_words <- apply(study_39[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
study_39$resolved_taxa_pair <- sorted_words
table(study_39$resolved_taxa_pair)


#Plot study 39 just for fun 
#Create new column with genus pair for labelling
sorted_words <- apply(study_39[, c('Gn1', 'Gn2')], 1, function(x) paste(x, collapse = "."))
study_39$genus_pair <- sorted_words

#Pivot 39 longer
study_39_longer <- study_39 %>%
  pivot_longer(cols = c("Log.prop.change.Gn1", "Log.prop.change.Gn2"), 
               names_to = "genera", 
               values_to = "log_prop_change")


study_39_longer %>%
  ggplot(aes(x = YEAR.T, y = log_prop_change, group=interaction(genus_pair, genera), 
             color = genus_pair, shape = genera)) + 
  geom_point() + 
  geom_jitter()+
  geom_smooth(se=FALSE, span = 0.5, aes(linetype = genera)) +
  facet_wrap(~genus_pair)+
  labs(x="Year", y = "Log(Proportional change in abundance)", 
       shape = "Genera 1 or 2", linetype = "Genera 1 or 2", 
       color = "Genus grouping")+
  theme_classic() 


#Add study 39 back into log_change.taxa
#Remove initial study 39 from log_change
log_change <- log_change %>%
  filter(!STUDY_ID=="39")
log_change <- bind_rows(study_39, log_change) 
#Arrange in order of increasing study 
log_change <- log_change %>%
  arrange(X)

#Add back into the full dataframe 
full_log_change <- readRDS("Revision 1 ecography/output/results_abundance_interactions_taxa_030724ENB_2.RDS")
full_log_change <- full_log_change %>%
 filter(!grepl("221", STUDY_PLOT)) %>%
  filter(!grepl("39", STUDY_PLOT))
#Keep only intersecting columns from full log change
common_cols <- intersect(names(full_log_change), names(log_change))
log_change <- log_change[, common_cols]
full_log_change <- bind_rows(full_log_change, log_change)

#Check that worked 
study_221_rbind <- full_log_change %>%
  filter(grepl("221", STUDY_PLOT))
study_39_rbind <- full_log_change %>%
  filter(grepl("39", STUDY_PLOT))
#Save final dataset
saveRDS(full_log_change, "Revision 1 ecography/output/results_abundance_interactions_taxa_031824ENB.RDS")



#Export 
#saveRDS(log_change.taxa, "Revision 1 ecography/output/results_abundance_interactions_taxa_030724ENB_2.RDS")


#Reassign the genera 

