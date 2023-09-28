# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(MCMCvis)
library(tidybayes)
library(ggdist)
library(marginaleffects)

#full model----
#w/ corr term
load(file = "outputs/Aug2023/mod_q1.terrestrial_withinstudies.Rdata") 
#w/o corr term 
load(file = 'outputs/Aug2023/mod_q1.terrestrial_withinstudiesv2.Rdata')

moddat<-mod$data

#plot model predline  
#use range of values from original dataset    
nd<-expand_grid(Prop.Change.Gn1  = NA, 
                Prop.Change.Gn2  = c(-7, -5, -3, -1,0, 1, 3, 5))

#or predict Gn1 at AVERAGE Gn2 
#nd <- with(moddat, expand.grid(Prop.Change.Gn2=mean(Prop.Change.Gn2), 
#                               Prop.Change.Gn=NA)) 

#or predict Gn1 with a subset of Gn2 spanning the data range 
nd<-expand_grid(Prop.Change.Gn1  = NA, 
                Prop.Change.Gn2  = seq(min(moddat$Prop.Change.Gn2), 
                                       max(moddat$Prop.Change.Gn2), 
                                       length.out = 20))

#moddat<-mod$data
#ndx<- expand_grid(Prop.Change.Gn1  = unique(moddat$Prop.Change.Gn1),
#                 Prop.Change.Gn2  = unique(moddat$Prop.Change.Gn2))%>%
#  slice(which(row_number() %% 20000 == 1))#make much smaller n~2000

#conditional effects 
#w/measurement error but not random error 
pred<-predictions(mod, nd, re_formula = NA) |> #setting RE=NA 
  posterior_draws()
gc()

predr<-predictions(mod, nd, re_formula = NULL) |> #setting RE=NULL (need UNIQUE ID in new data)
  posterior_draws()

plot(pred$Prop.Change.Gn2, pred$estimate) #mean of all draws per x 

plot(pred$Prop.Change.Gn2, pred$draw) #all draws

#w/o measurement error 
fit<-fitted(mod, re_formula = NA, newdata = nd,
        ndraws=100,
        summary = FALSE)
fit<-as.data.frame(fit)


#epred<-posterior_epred(mod, nd, re_formula = NA)|> #setting RE=NA here because otherwise error ribbon is weird
#  posterior_draws()

#preds<-predictions(mod, re_formula = NA) #marginal effects version using full dataset 

#plot model predline over raw data 
#I'm not sure if I should plot the draws or the estimates here??
ggplot(pred, aes(x = Prop.Change.Gn2, y = estimate)) + 
  geom_line(aes(color="red"))+
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.5)+
  #stat_lineribbon()+ #scale_fill_brewer() +
  #geom_smooth(method='gam', formula= y ~ s(x, bs = "cs", fx = TRUE, k = 2), se = T)+
  #geom_ribbon(aes(ymin = (conf.low*1.11)+2,
  #                ymax = (conf.high*1.11)+2), alpha=0.2, outline.type = "both")+
  #geom_point(data=flowdat, aes(x = (doy*14)+172, y=(value*1.11)+2), alpha=0.2)+
  # plot raw data
geom_point(data=moddat, aes(x = Prop.Change.Gn2, y=Prop.Change.Gn1,alpha=0.2))+
  labs(x = "Prop yearly abundance change genus x ",
       y = "Prop yearly abundance change genus y") + theme_bw()

#try with marginal effects plot_predictions 
#https://marginaleffects.com/dev/articles/predictions.html#bayesian-models
plot_predictions(mod, condition = "Prop.Change.Gn2")


#combine with raw data for pred vs obs 
predall<-select(pred, -Prop.Change.Gn1)%>%left_join(., select(moddat, -UNIQUE.PAIR.ID))

data=moddat

#predict over new data 
#use range of original dataset
hist(mod$data$Prop.Change.Gn1)
hist(mod$data$Prop.Change.Gn2)
mean(mod$data$Prop.Change.Gn1)
sd(mod$data$Prop.Change.Gn1)

mean(mod$data$Prop.Change.Gn2)
sd(mod$data$Prop.Change.Gn2)

#simulate with normal dist 
nd2 <- expand_grid(Prop.Change.Gn1  = rnorm(35, mean=0, sd=0.5),
                   Prop.Change.Gn2 = (35, mean=0, sd=0.5))



#OLD CODE----

# # FITTING ISSUES

# NA's -> these are biomass values

# divergent transitions -> NEED TO DEAL WITH THIS (http://mc-stan.org/misc/warnings.html#divergent-transitions-after-warmup)
# low ESS -> NEED TO DEAL WITH THIS

coefs <- coef(mod)$`PairID:UNIQUE.PAIR.ID`[1:52, 1, 2] #random int & slopes
coef(mod)
hist(coefs, breaks = 100)

loo<-loo(mod)
loo
### NOT RUN; STOP 26.4.2023 NC

# Check posteriors
summary(mod)
plot(mod)
plot(conditional_effects(mod))
MCMCvis::MCMCtrace(mod) #plots all fixed & radom slope posteriors in a separate pdf 

# Check prior
prior_summary(mod) 
ps <- powerscale_sensitivity(mod)
unique(ps$sensitivity$diagnosis)

# Model fit
pp_check(mod, ndraws = 100)
ppc_stat(y = dat$rel_rec, 
         yrep = posterior_predict(mod, ndraws = 1000), stat = "skew")
model_loo <- loo(mod, save_psis = TRUE, cores=4) 
w <- weights(model_loo$psis_object)
ppc_loo_pit_overlay(dat$rel_rec, 
                    posterior_predict(mod), lw = w)



#k fold cross validation 
## LOO Script
## add meta info
library(loo)

#load dummy model 
load('outputs/brms_April2023/dummy_mod2.rds')


kfold_group<-kfold_split_grouped(K = 5, x =  mod$data$UNIQUE.PAIR.ID)


kfold_out <- rstanarm::kfold(mod, K=10, folds = kfold_group)

datx<-left_join(dat, kfold_strat)

kfold_strat
head(dat$UNIQUE.PAIR.ID)
summary(mod)
coef <- as.data.frame(coef(mod)$'PairID:UNIQUE.PAIR.ID')
coef2 <- as.data.frame(coef(mod)$'UNIQUE.PAIR.ID')
coef3 <- as.data.frame(coef(mod)$'PairID')


