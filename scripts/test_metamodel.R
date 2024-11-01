library(tidyverse)
library(brms)

# current full dataset, update w/ final version
alldat <- readRDS("data/data_processing/within.study.updated.interactions.020724ENB.RDS") 

names(alldat)
alldat <- alldat %>% 
  group_by(TS_ID) %>%
  mutate(cor=cor(Log.prop.change.Gn1, Log.prop.change.Gn2)) #calculate pearson's cor for all time series 

hist(alldat$cor)

#play with test dataset
alldat_test <- subset(alldat, cor<0.8 & cor>-0.8) #subset to only relatively strong cors 

hist(alldat_test$cor)

max(alldat_test$cor)
min(alldat_test$cor)

#calculate z scores and SE w/ sample sizes 
# here we use the length of the time series (per plot per study per genus pair)
# other options include:
# total number of individuals recorded across time series for both genera
# abs difference in total abundance of indivs between the genera
# total number of species recorded across both genera
# abs diff in species between the genera 
alldat_test<-mutate(alldat_test, z=0.5*log((1+cor)/(1-cor)))%>% #eq 3.11 in meta-analysis book 
  #mutate(CIlow = (-1.96 /sqrt(SERIES.l-3))+ z)%>%
  #mutate(CIhigh=  (1.96 /sqrt(SERIES.l-3))+ z)%>%
  mutate(SE=(1/sqrt(SERIES.l-3))) #eq 3.12 in meta-analysis book #can update with different sample size defintions here 


#look at some to check   
select(alldat_test, cor, z, SE)
hist(alldat_test$SE)
#model with z scores and SE as joint response

# filter df to remove repeat values 
# ie keeping unique cases of z-scores
head(alldat_test)
moddat<- select(alldat_test,-Log.prop.change.Gn1, -Log.prop.change.Gn2, -X,  -YEAR.T, -YEAR.T1)%>%distinct(.)
names(moddat)


# add in metadata info
# this will need to get expanded to >10km for full dataset 
# probably need to go back to the full metadata table methods for each study 
load("data/prep_biotime/meta_pairs_10km.RData")
moddat<-left_join(moddat, meta.pairs[, c(1,20,46)])
names(moddat)

#set some priors 
priors<-c(prior(normal(0,0.33), class = Intercept), #set between -1 and 1 for z score
          prior(normal(0,0.33), lb=0, class = sd)) #set lower bound 0 for SE, values b/w (0,1)

MODFORM <- bf(z|resp_se(SE) ~ 
                scale(SERIES.l) + # if two disjoint time series, adding the together to get total length
                #TAXA1 +
                scale(LATITUDE)+
                #CLIMATE1 + #*interaction_present  + #identical to column CLIMATE2 ie all within study comparisons are within the same climate
                treatment_yn + #is there some sort of disturbance yes/no (fertilizer, fire, grazing etc)
                #interaction_benefit +
                #does GLOBI record these genera as potentially interacting
                interaction_present + 
                #Centrality_Betweenness_Edge +
                (1|STUDY_ID))#(1|RESOLVED.TAXA.PAIR) # want this to be study and plot ID??

metamod <- brm(MODFORM, moddat,
             control = list(adapt_delta=0.8, max_treedepth = 11), cores=3, chains=3, 
             iter=5000, family=gaussian, prior = priors)
