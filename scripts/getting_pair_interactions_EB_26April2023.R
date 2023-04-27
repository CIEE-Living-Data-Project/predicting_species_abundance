# Date created: 26 Apr 2023
# Date updated: 26 Apr 2023 (NC)

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




log.prop.change.with.meta <- readRDS("data/log.prop.change.with.meta.RDS")
log_change <- log.prop.change.with.meta

head(log.prop.change.with.meta)

#list of bio pairs with at least 10 overlapping years
genus_pairs <- apply(log_change[, c("Gn1", "Gn2")], 1, function(x) paste(sort(x), collapse = "-"))
unique_pairs <- unique(genus_pairs)
split_pairs <- data.frame(Gn1 = sapply(strsplit(unique_pairs, "-"), "[", 1),
                          Gn2 = sapply(strsplit(unique_pairs, "-"), "[", 2))

#short sample dataset: 
short_pairs <- head(split_pairs, 100)

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


pb<-set.prog.bar(nrow(split_pairs)) #sets progress bar
for (i in 1:nrow(split_pairs)) {
  pb$tick()
  gn1 <- short_pairs[i, "Gn1"]
  gn2 <- short_pairs[i, "Gn2"]
interactions <- get_interactions_by_taxa(sourcetaxon=gn1, 
                         targettaxon=gn2)
unique <- unique(interactions$interaction_type)
new_df <- data.frame(genus_1 = rep(gn1, each = length(unique)), 
                  genus_2 = rep(gn2, each = length(unique)))
pairs <- cbind(new_df, char = rep(unique, length.out = nrow(new_df)))
df_list[[i]] <- pairs
}

pair_interactions <- do.call(rbind, df_list)

