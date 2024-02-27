# Date created: 22 Feb 2024 (ENB)
# Date updated:  

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
#Load in the new abundance data 
results_abundance <-read.csv("~/Documents/Work and Career/LDP/Working Group/results.abundance.csv", 
                             )


#Check everything loaded in properly and is named properly: 
head(results_abundance)
colnames(results_abundance) #need to remove a column called X, which is just row numbers
summary(results_abundance)
results_abundance <- results_abundance %>%
  select(-X)
colnames(results_abundance)

#Great! Now let's do some variable checks - are the values what we expect? 
#Gn1 - 
head(unique(results_abundance$Gn1))
unique(head(results_abundance$Gn2, 100))
unique(head(results_abundance$YEAR.T1), 100)
unique(head(results_abundance$YEAR.T), 100)

#Do histograms and QQplots of the numerical variables we're interested in
#Get a subsample of the data for processing time 
set.seed(4477)
sampled_results_abundance <- slice_sample(results_abundance, n=500000)

hist(sampled_results_abundance$Log.prop.change.Gn1, breaks = 100)
qqnorm(sampled_results_abundance$Log.prop.change.Gn1, pch = 1, frame = FALSE)
qqline(sampled_results_abundance$Log.prop.change.Gn1, col = "steelblue", lwd = 2)
#We have a quite heavy tailed QQplot 

hist(sampled_results_abundance$Log.prop.change.Gn2, 
     breaks = 100)
qqnorm(sampled_results_abundance$Log.prop.change.Gn2, pch = 1, frame = FALSE)
qqline(sampled_results_abundance$Log.prop.change.Gn2, col = "steelblue", lwd = 2)
#Quite a heavy tailed distribution here too 



hist(results_abundance$Abs.change.Gn1, breaks = 100)
min(results_abundance$Abs.change.Gn1)
max(results_abundance$Abs.change.Gn1)
mean(results_abundance$Abs.change.Gn1)


#Looks like we have a strange outlier here - an Abs change of -2161. Where is that, 
#and what group is it in? 
#This is in our bird study 
qqnorm(sampled_results_abundance$Abs.change.Gn1, pch = 1, frame = FALSE)
qqline(sampled_results_abundance$Abs.change.Gn1, col = "steelblue", lwd = 2)
#As expected, quite a heavy tailed distribution here too

#Same here as above
hist(sampled_results_abundance$Abs.change.Gn2, breaks = 100)
qqnorm(sampled_results_abundance$Abs.change.Gn2, pch = 1, frame = FALSE)
qqline(sampled_results_abundance$Abs.change.Gn2, col = "steelblue", lwd = 2)


#Ok, now abundances (this might be where things get weird)
hist(sampled_results_abundance$Abun.T.Gn1, breaks = 100)
hist(sampled_results_abundance$Abun.T1.Gn1, breaks = 100)
#So in general, we see absolute abundances that are low. Can we zoom in on this a bit? 
hist(sampled_results_abundance$Abun.T.Gn1, breaks = 1000, 
     xlim = c(0, 100))
#So we do have some low sample sizes here - let's check the QQ
qqnorm(sampled_results_abundance$Abun.T.Gn1, pch = 1, frame = FALSE)
qqline(sampled_results_abundance$Abun.T.Gn1, col = "steelblue", lwd = 2)
#Looks like what I would expect 
#Gn2 should also just be the reverse of this 

########################################################
#ok, so now checking out those flip-flopping abundances 
#Get a function that pulls out all the patterns in the dataframe 
#(Adapted from stack overflow - https://stackoverflow.com/questions/4617407/finding-repeated-patterns-in-r)
head_results_abundance <- head(results_abundance, 50000)
column_matrix <- t(matrix(head_results_abundance$Log.prop.change.Gn1, ncol = 5000))
pasteSort <- function( x ) do.call(paste, as.list( sort( x) ) )
#Test on a smaller subet
pairs <- c( apply(column_matrix, 1, function(row) apply( combn(row, 4), 2, pasteSort)) )
pairFreqs <- as.data.frame(table(pairs))
pairFreqs <- pairFreqs %>%
  arrange(desc(Freq))
pairFreqs[ pairFreqs > 1 ]
#Split pairFreqs into four columns 
pairFreqs_split <- pairFreqs %>%
  separate(col = pairs, into = c("num_1", "num_2", "num_3", "num_4"), sep = " ", convert = TRUE) 

sorted_pairs <- pairFreqs_split %>%
  mutate(across(starts_with("num"), sort))

summary_pairs <- sorted_pairs %>%
  group_by(num_1, num_2, num_3, num_4) %>%
  summarise(total_freq = sum(Freq))

#Let's search for one pattern in particular:
pairFreqs$pairs <- as.character(pairFreqs$pairs)
pairFreqs <- pairFreqs %>%
  arrange(desc(Freq))
pattern_test <- pairFreqs$pairs[1]
pattern_test <- gsub(" ", ", ", pattern_test)
pattern_vector <- as.numeric(unlist(strsplit(pattern_test, ", ")))

# pattern_vector_test <- c(0.6, 0.7, 0.8, 0.9)

matching_indices <- c()

# Iterate over the dataframe to find matching segments
for (i in 1:(length(head_results_abundance$Log.prop.change.Gn1) - length(pattern_vector) + 1)) {
  segment <- head_results_abundance$Log.prop.change.Gn1[i:(i + length(pattern_vector) - 1)]
  if(all(segment == pattern_vector)) {
    matching_indices <- c(matching_indices, i:(i + length(pattern_vector) - 1))
  }
}

# Filter the dataframe using the matching indices
strange_values <- head_results_abundance[matching_indices, ]


