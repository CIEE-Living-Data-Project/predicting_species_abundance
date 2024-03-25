library(metafor)
library(tidyverse)

#load data for metafor model 
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")
memory.limit(999999999)

metamodel.indivs<-rma(z ~ scale.SERIES.l + treatment_yn_clean + 
                        resolved_taxa_pair + 
                        scale.abs.lat +
                        interaction_present.factor +
                        scale.elev, # z is the effect size from each time series pair (correlation)
                        #random= ~1|STUDY_ID, 
                        sei=SE.total.indivs, # measure of standard error for each z value 
                        method = "REML", # Using  a REML estimator which is common for random effect meta-analyses
                        data = moddat) 

