library(tidyverse)
library(lmerTest)

# load data ####
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")
#remove studies not using 
moddat <- filter(moddat, STUDY_ID != "225") 
moddat <- filter(moddat, STUDY_ID != "39")

# run hierarchical model ####
lmermod <- lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + 
                          (1|STUDY_ID) +
                 (1| resolved_taxa_pair),  dat=moddat)
save(lmermod, file="Revision 1 ecography/output/lmer_model_final.Rdata")


# run model with only one random intercept to test differences in magnitude of conditional R2 driven by each term
lmermod.studyOnly <- lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + 
                          (1|STUDY_ID),  dat=moddat)
# results match with full model

lmermod.taxaOnly <- lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + 
                                       (1| resolved_taxa_pair),  dat=moddat)
# interesting, when no study ID then latitude is significant


# model summary and R2 ####
## full model ####
summary(lmermod)
MuMIn::r.squaredGLMM(lmermod)
performance::r2(lmermod)
# Conditional R2: 0.298
# Marginal R2: 0.003

## compare R2 between variations of random effects models ####
# very similar conditional between two models, 2% higher in study ID model
# marginal .07% higher in study ID model
MuMIn::r.squaredGLMM(lmermod.studyOnly)
MuMIn::r.squaredGLMM(lmermod.taxaOnly)

# model k fold cross validation ####
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


# null model shuffle to get p-values ####
# to account for non-independence of pairwise comparison 
# (1000 shuffles for faster run time, should be 10K for pubs)
func.pvalShuffle.acrossPlots <- function(resp.dat, times) {
  mod.real<-lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + (1|STUDY_ID) +
                             (1| resolved_taxa_pair), dat=moddat)
  
  core.function<-function(x){
    n.ord<-sample(1:length(resp.dat))
    mod.rand<-lmerTest::lmer(resp.dat[n.ord]~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + (1|STUDY_ID) +
                               (1| resolved_taxa_pair), dat=moddat)
    drop1(mod.rand, test="F")[ ,5]
  }
  null <- sapply(1:times, core.function) # output each row is a term in model, each column is iteration
  
  Fobs <- drop1(mod.real, test='F')[,5]
  p <- apply(null>Fobs, 1, sum)/times # count how many times null F is greater than obs and divide by iterations, columns are terms in model
  data.frame(Frand=apply(null, 1, mean), Fobs, p)
  
}
# all singular fits due to lack of stdev across taxonomic pair
# 2.485368 hrs
start.time <- Sys.time()
null.out.acrossPlots <- func.pvalShuffle.acrossPlots(moddat$z, 1000)
end.time <- Sys.time()
end.time - start.time
# 
# # shuffle response within plot, largely consistent with shuffle across plots
# func.pvalShuffle.withinPlots <- function(resp.dat, times) {
#   mod.real<-lmerTest::lmer(z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + (1|STUDY_ID) +
#                              (1| resolved_taxa_pair), dat=moddat)
#   
#   core.function<-function(x){
#     shuffle.z <- moddat %>% 
#       group_by(STUDY_PLOT) %>% 
#       mutate(n.ord = sample(1:length(z))) %>%  # this works, checked
#       select(STUDY_PLOT, z, n.ord) %>% 
#       group_by(STUDY_PLOT) %>% 
#       arrange(-desc(n.ord), .by_group = TRUE)
#     
#     mod.rand<-lmerTest::lmer(shuffle.z$z~scale.SERIES.l+ treatment_yn_clean + scale.abs.lat + interaction_present.factor + (1|STUDY_ID) +
#                                (1| resolved_taxa_pair), dat=moddat)
#     drop1(mod.rand, test="F")[ ,5]
#   }
#   null <- sapply(1:times, core.function) # output each row is a term in model, each column is iteration
#   
#   Fobs <- drop1(mod.real, test='F')[,5]
#   p <- apply(null>Fobs, 1, sum)/times # count how many times null F is greater than obs and divide by iterations, columns are terms in model
#   data.frame(Frand=apply(null, 1, mean), Fobs, p)
#   
# }
# 
# # ~40 warnings about singular fits
# null.out.withinPlots <- func.pvalShuffle.withinPlots(moddat, 1000)

