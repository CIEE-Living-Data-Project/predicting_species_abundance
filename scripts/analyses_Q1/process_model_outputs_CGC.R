
# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(MCMCvis)
library(tidybayes)
library(ggdist)

#full model----
#w/ corr term
load(file = "outputs/Aug2023/mod_q1.terrestrial_withinstudies.Rdata") 
#w/o corr term 
load(file = 'outputs/Aug2023/mod_q1.terrestrial_withinstudiesv2.Rdata')

#raw data from model for plotting
load(file = 'outputs/Aug2023/Q1_rawdata_full.Rdata')

#assess convergence issues 
summary(mod)
#convergence issue only on 'cor' term - re-run without 
#MCMCtrace(mod) #this takes a while- but mostly looks good!

#Unique study ID
slopes<-coef(object = mod)
slopes<-as.data.frame(slopes$UNIQUE.PAIR.ID)
slopes<-select(slopes, contains("Gn2"))
slopes$UniquePairID<-row.names(slopes)

#save(slopes, file = "outputs/Aug2023/randomslopes_q1modelv2.Rdata")

pp_check(mod) #captures mean and variance ok but misses magnitude of mean
#loox<-loo(mod) #too big won't work


# Check the normal distribution of random effect
qqnorm(slopes$Estimate.Prop.Change.Gn2, main = "Normal Q-Q plot of random slopes",  bty="n")
qqline(slopes$Estimate.Prop.Change.Gn2, col = 2, lwd = 3, lty = 3) # a bit off on the tails 

#plot results 


#plot random slopes 
bars<-ggplot(data=slopes, aes(y = UniquePairID, x=Estimate.Prop.Change.Gn2)) + 
  geom_pointrange(aes(xmin=Q2.5.Prop.Change.Gn2, xmax=Q97.5.Prop.Change.Gn2), size=0.01, alpha=0.5, color='grey')+ 
  geom_vline(xintercept = 0, color='black', lty=2)+
   theme(axis.text.y=element_blank(),  #remove y axis labels
                axis.ticks.y=element_blank()) + xlab("Group level slope estimates")

#add posteriors on top 
#these are equivalent 
postmod<-as_tibble(fixef(mod, summary=F))
postmod2<-as_tibble(posterior_samples(mod, pars = "b_Prop.Change.Gn2"))

transparent_theme <- theme( axis.title.y = element_blank(), axis.title.x = element_blank(),
                             axis.text.x = element_blank(), 
                             axis.text.y = element_blank(), 
                             axis.ticks.x = element_blank(), 
                             axis.ticks.y = element_blank(), 
                             panel.grid = element_blank(), 
                             axis.line = element_blank(), 
                             panel.background = element_blank(), 
                             plot.background = element_blank(), 
                             panel.border = element_blank(), 
                             panel.grid.major = element_blank(), 
                             panel.grid.minor = element_blank(), 
                             panel.grid.major.x = element_blank(),
                             panel.grid.major.y = element_blank())


hist<-ggplot(postmod2,aes(x=b_Prop.Change.Gn2))+ geom_density()+   transparent_theme
grob_hist<-ggplotGrob(hist)

bars + annotation_custom(grob = grob_hist,  xmin=min(postmod2$b_Prop.Change.Gn2), xmax=max(postmod2$b_Prop.Change.Gn2))


moddat<-mod$data

moddat%>%
add_predicted_draws(mod) 


