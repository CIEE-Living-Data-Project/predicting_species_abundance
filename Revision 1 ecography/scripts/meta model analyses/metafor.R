library(metafor)
library(tidyverse)

#load data for metafor model 
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")
memory.limit(999999999)

#clean up moddat
moddat<-select(moddat, z, scale.SERIES.l, treatment_yn_clean, 
                 resolved_taxa_pair, 
                 scale.abs.lat, 
                 interaction_present.factor, 
                 scale.elev, SE.total.indivs, SE.total.sp, var.total.indivs, var.total.sp, STUDY_ID)
names(moddat)

metamodel.indivs<-rma.mv(z ~ scale.SERIES.l + # z is the effect size from each time series pair (correlation)
                        treatment_yn_clean + 
                        resolved_taxa_pair + 
                        scale.abs.lat +
                        interaction_present.factor +
                        scale.elev, 
                        V=var.total.indivs, # measure of variance for each z value w/ indiv sample sizes 
                        random = ~1|STUDY_ID, 
                        method = "REML", # REML estimator- common for random effect meta-analyses
                        data = moddat) 
metamodel.indivs

metamodel.sp<-rma.mv(z ~ scale.SERIES.l + # z is the effect size from each time series pair (correlation)
                          treatment_yn_clean + 
                          resolved_taxa_pair + 
                          scale.abs.lat +
                          interaction_present.factor +
                          scale.elev, 
                          V=var.total.sp, # measure of variance for each z value w/ spp sample sizes
                          random = ~1|STUDY_ID, 
                          method = "REML", # REML estimator- common for random effect meta-analyses
                          data = moddat) 
metamodel.sp