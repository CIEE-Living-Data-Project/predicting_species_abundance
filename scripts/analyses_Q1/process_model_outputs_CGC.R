
# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(MCMCvis)
library(tidybayes)
library(ggdist)

load(file = "outputs/brms_July2023/meta_mod_q1.terrestrial_withinstudies.Rdata") 


#meta model outputs---- 
summary(mod)
posterior_summary(mod)

get_variables(mod)
draws<-spread_draws(mod,b_Intercept, bsp_meestimate_Gn2std.error_Gn2, sdme_meestimate_Gn2, sigma) 

#random slopes by study ID and time series length 
#study ID
slopes<-coef(object = mod)
slopes<-as.data.frame(slopes$PairID)
slopes<-select(slopes, contains("Gn2"))
slopes$PairID<-row.names(slopes)

#merge back with metadata  
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23
dat_terr<-subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
clim<-select(dat_terr,PairID, CLIMATE1, CLIMATE2)  %>%
  distinct(.)
slopes<-left_join(slopes, clim)

#plot
#actual values- captures direction of relationship- temperate slightly more likely to positively co-occur over time (competition??) 
ggplot(data = slopes, aes(y=Estimate.estimate_Gn2, x=CLIMATE1))+
  geom_boxplot()+ theme_bw()
#absolute values-captures strength of correlation -no real difference by climate  
ggplot(data = slopes, aes(y=abs(Estimate.estimate_Gn2), x=CLIMATE1))+
  geom_boxplot()+ theme_bw()

#Series length 
slopes2<-coef(object = mod)
slopes2<-as.data.frame(slopes2$SERIES.l)
slopes2<-select(slopes2, contains("Gn2"))
slopes2$SERIES.l<-row.names(slopes2)

#plot
#actual values- captures direction of relationship -more likely to negatively co-occur over time (competition??) 
ggplot(data = slopes2, aes(y=Estimate.estimate_Gn2, x=as.numeric(SERIES.l)))+
  geom_point()+ theme_bw()+
  geom_smooth()+
  geom_smooth(method = "lm")
#absolute values-captures strength of correlation -no real change in strength over time 
ggplot(data = slopes2, aes(y=abs(Estimate.estimate_Gn2), x=as.numeric(SERIES.l)))+
  geom_point()+ theme_bw()+
  geom_smooth()+
  geom_smooth(method = "lm")

ranef_terr=as.data.frame(ranef$PairID)
ranef_terr$PairID<-row.names(ranef_terr)

ranef_terr2=as.data.frame(ranef$SERIES.l)
ranef_terr2$SERIES.l<-row.names(ranef_terr2)

#leave one out cross validation 
loo1<-loo(mod)
save(loo1, file = 'outputs/brms_July2023/looCV_withinstudies_meta.Rdata') #save          

#loo2<-loo(mod, moment_match = T )

#ppchecks 
ppcheck<-pp_check(mod, ndraws = 100) #this doesn't look great 

ppcheck<-pp_check(mod,type = "error_hist", ndraws = 100) #this doesn't look great 

MCMCvis::MCMCtrace(mod)


#predictive accuracy 
load("outputs/brms_July2023/meta_mod_q1.terrestrial_withinstudies.Rdata")

# Extract the posterior samples of the random slopes
random_slopes <- posterior_samples(mod, pars = "r_PairID")

# Make predictions using the posterior samples
predictions <- posterior_predict(model, newdata = your_data)

# Calculate the predictive accuracy for each observation
accuracy <- abs(predictions - your_data$y)



