log_change_interaction <- readRDS("data/data_processing/log.prop.change.interactions.RDS")


log_change_interaction_meta <- readRDS("data/data_processing/log.prop.change.with.meta.RDS")

load('data/prep_biotime/meta_pairs_10km.RData')



library(stringr)

# Words to search for
word <- paste(c("treatment", "manipulat*", "clearcut", "clear-cut", "burn*", "logged"), collapse = '|')


# Query the text column for the word
result <- which(str_detect(meta.pairs$METHODS, word), )
result_negate <- which(str_detect(meta.pairs$METHODS, word, negate = TRUE), )

# Output the result
print(result)
length(result)

print(result_negate)

## #1 had a couple clearcuts but not really stated as a treatment
meta.pairs$METHODS[13]
length(which(!is.na(meta.pairs$TREATMENT))) #3 (8, 12, 33)
length(which((meta.pairs$GENERAL_TREAT==TRUE))) #2 (8,12)


meta.pairs$METHODS[52]
no.treat <- c(2,3,4,5,6,7,9,10,11,13,14,15,16,17,18,19,20,21,22,23,26,29,30,34,35,36,37,38,39,49,51,52)
length(no.treat) #32

yes.treat <- c(1,8,12,24,25,27,28,31,32,33,40,41,42,43,44,45,46,47,48,50,53,54,55,56,57,58,59,
               60,61,62,63,64,65,66,67,68,69,70,71,72,73)
length(yes.treat) #41

meta.pairs$treatment_yn <- rep("", times = length(meta.pairs$STUDY_ID))
meta.pairs$treatment_yn[no.treat] <- "no"
meta.pairs$treatment_yn[yes.treat] <- "yes"

## 9,10 were a laboratory treatment, probably want to exclude
## 29 -- old agricultural field abandonment -- is this maybe a succession treatment? 
    # different fiels were abandoned at various times, no control that i can tell
    # for now classified as no.treat

## treatment notes, will probably need to manually enter
# 1 - clearcut (logging)
# 8 - prescribed burning
# 12 - logging, not captured by str_detect...
# 24 - mentions treatment transects, but not what the treatment is (something with succession)
# 25 - prescribed burning
# 27 - tillage and crop systems
# 28 - fire-grazing treatments
# 31 - mentions treatment transects, but not what the treatment is
# 32 - fertilization
# 33 - small mammal exclosures
# 40 - logging
# 41 - logging
# 42 - logging
# 43 - logging
# 44 - logging
# 45 - logging
# 46 - logging
# 47 - logging
# 48 - prescribed burning
# 50 - warming
# 53 - prescribed burning
# 54 - prescribed burning
# 55 - prescribed burning
# 56 - prescribed burning
# 57 - prescribed burning
# 58 - prescribed burning
# 59 - prescribed burning
# 60 - prescribed burning and grazing
# 61 - prescribed burning
# 62 - prescribed burning and grazing
# 63 - prescribed burning
# 64 - prescribed burning and grazing
# 65 - prescribed burning
# 66 - logging
# 67 - logging
# 68 - logging
# 69 - logging
# 70 - logging
# 71 - logging
# 72 - logging
# 73 - logging


meta.pairs$treatment_descr <- rep("", times = length(meta.pairs$STUDY_ID))

meta.pairs$treatment_descr[1] <-  "clearcut"

meta.pairs$treatment_descr[c(24,31)] <- c("unclear", "unclear")

meta.pairs$treatment_descr[27] <-  "tillage and crop systems"
meta.pairs$treatment_descr[28] <-  "fire-grazing treatments"
meta.pairs$treatment_descr[32] <-  "fertilization"
meta.pairs$treatment_descr[33] <-  "small mammal exclosures"
meta.pairs$treatment_descr[50] <-  "warming"

meta.pairs$treatment_descr[c(12,40,41,42,43,44,45,46,47,
                             66,67,68,69,70,71,72,73)] <- c(rep("logging", times = 17))

meta.pairs$treatment_descr[c(8,25,48,53,54,55,56,57,58,59,61,63,65)] <- c(rep("prescribed burning", times = 13))

meta.pairs$treatment_descr[c(60,62,64)] <- c(rep("prescribed burning and grazing", times = 3))

meta.pairs$treatment_descr[which(meta.pairs$treatment_descr == "")] <- NA


meta.pairs$treatment_simplified <- rep("", times = length(meta.pairs$STUDY_ID))

meta.pairs$treatment_simplified[c(1,12,40,41,42,43,44,45,46,47,66,67,68,69,70,71,72,73)] <- c(rep("logging", times = 9))
meta.pairs$treatment_simplified[c(8,25,48,53,54,55,56,57,58,59,61,63,65, 
                                  60,62,64)] <- c(rep("prescribed burning", times = 16))
meta.pairs$treatment_simplified[c(24,31, 27,28,32,33,50)] <- c(rep("other", times = 7))

meta.pairs$treatment_simplified[which(meta.pairs$treatment_simplified == "")] <- NA
