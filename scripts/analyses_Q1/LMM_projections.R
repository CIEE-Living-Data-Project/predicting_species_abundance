
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


