log_change_interaction <- readRDS("data/data_processing/log.prop.change.interactions.RDS")


log_change_interaction_meta <- readRDS("data/data_processing/log.prop.change.with.meta.RDS")

load('data/prep_biotime/meta_pairs_10km.RData')



library(stringr)

# Sample data frame
df <- data.frame(id = 1:5,
                 text = c("This is a sample text",
                          "Another example of text",
                          "Text with a word",
                          "No match here",
                          "Some more sample text"))

# Word to search for
word <- "treatment"

# Query the text column for the word
result <- which(str_detect(meta.pairs$METHODS, "treatment"), )

# Output the result
print(result)


## #1 had a couple clearcuts but not really stated as a treatment
meta.pairs$METHODS[9]
length(which(!is.na(meta.pairs$TREATMENT))) #3
length(which((meta.pairs$GENERAL_TREAT==TRUE))) #2
