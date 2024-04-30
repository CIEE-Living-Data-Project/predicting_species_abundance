library(tidyverse)
library(lmerTest)

#load data---- 
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")
#remove studies not using 
moddat <- filter(moddat, STUDY_ID != "225") 
moddat <- filter(moddat, STUDY_ID != "39")

#run hierarchical model----
lmermod<-lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + (1|STUDY_ID) +
                 (1| resolved_taxa_pair),  dat=moddat)
save(lmermod, file="Revision 1 ecography/output/lmer_model_final.Rdata")

#model summary and R2---- 
summary(lmermod)
MuMIn::r.squaredGLMM(lmermod)
performance::r2(lmermod)

#model k fold cross validation---- 
load(file="Revision 1 ecography/output/lmer_model_final.Rdata")
library(cv)
#fold 10 times with n/10 dataset each time 
kfoldcv<-cv(lmermod, k=10) 
save(kfoldcv, file="Revision 1 ecography/output/lmer_model_kfoldstudy.Rdata")

#fold n times (63 study IDs) leaving out 1 study each time 
kfoldcv_study<-cv(lmermod, clusterVariables="STUDY_ID") 
save(kfoldcv_study, file="Revision 1 ecography/output/lmer_model_kfoldstudy.Rdata")

#fold n times (16 taxa categories) leaving out 1 category each time 
kfoldcv_tax<-cv(lmermod, clusterVariables="resolved_taxa_pair")
save(kfoldcv_tax, file="Revision 1 ecography/output/lmer_model_kfoldtax.Rdata")

#posterior predictive check plots- Fig 5a, b
library(performance)
library(see)
#w random effects- Fig 5a
performance::check_predictions(lmermod, re_formula = NULL)
#w/o random effects-Fig 5b
performance::check_predictions(lmermod, re_formula = NA)


