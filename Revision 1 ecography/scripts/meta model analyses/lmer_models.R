
library(tidyverse)
library(lmerTest)
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")

lmermod<-lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + (1|STUDY_ID) +
                 (1| resolved_taxa_pair),  dat=moddat)
summary(lmermod)
MuMIn::r.squaredGLMM(lmermod)
save(lmermod, file="Revision 1 ecography/output/lmer_model_final.Rdata")

#add total species # back in 

load(file="Revision 1 ecography/output/prep_data/all_model_data.Rdata")
moddat2<-select(moddat, TS_ID, total.sp, total.indivs)
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")

moddatx<-left_join(moddat, moddat2)

lmermod2<-lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor +
                          scale(total.sp) + (1|STUDY_ID) + (1| resolved_taxa_pair),  dat=moddatx)
summary(lmermod2) ##spp highly correlated with time serires length 
MuMIn::r.squaredGLMM(lmermod)

lmermod3<-lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor +
                           scale(SE.total.sp) + (1|STUDY_ID) + (1| resolved_taxa_pair),  dat=moddatx)
summary(lmermod3) #SE highly correlated with time series length 
MuMIn::r.squaredGLMM(lmermod3)

moddatx$zw<-moddatx$z*moddatx$SE.total.sp
hist(moddatx$zw)

lmermod4<-lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor +
                           (1|STUDY_ID) + (1| resolved_taxa_pair),  dat=moddatx)
summary(lmermod4)  
MuMIn::r.squaredGLMM(lmermod3)

#try with weighted z values?? NOT finished 

# 1. Calculate weights
#load back in raw data with sample sizes 
load(file="Revision 1 ecography/output/prep_data/all_model_data.Rdata")

moddat$weights<-moddat$total.sp-3

# 2. Weighted mean of z for each study ID
moddat <- group_by(moddat, STUDY_ID)%>%
          mutate(zbar= sum(z*weights)/sum(weights))%>%ungroup(.)

unique(moddat$zbar)#62- same as # studies

#3. Calculate Q- the weighted heterogeneity among the observed Z values
moddat<-rowwise(moddat)%>%mutate(Q= sum(weights * (z- zbar)^2))

#4. Calculate Tau squared- the within study sampling variance
moddat<-group_by(moddat, STUDY_ID)%>%#rowwise(moddat)%>%
  mutate(tau.sq =(Q - (total.sp - 1))/(sum(weights) - sum(weights^2)/sum(weights)))

#5. Calculate random weights
moddat<-mutate(ranweights <- 1/(variance + tau.sq)

#6. Calculte random Weighted mean of z for each study ID
moddat <- group_by(moddat, STUDY_ID)%>%
  mutate(zbar_ran= sum(z*ranweights)/sum(ranweights))%>%ungroup(.)


moddat_small<-group_by(moddat, resolved_taxa_pair)%>%slice_sample(n=10)%>%ungroup(.)

moddat_small<-select(moddat_small, TS_ID, STUDY_ID, STUDY_PLOT,  Gn1, Gn2, z, zbar, NUMBER_OF_SPECIES, total.sp)
