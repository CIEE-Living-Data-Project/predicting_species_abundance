# # LIBRARIES # #
library(tidyverse)
library(dplyr)
library(brms)
library(MCMCvis)
library(tidybayes)
library(ggdist)

#QA QC full model----
#w/ corr term
load(file = "outputs/Aug2023/mod_q1.terrestrial_withinstudies.Rdata") 
#w/o corr term 
#load(file = 'outputs/Aug2023/mod_q1.terrestrial_withinstudiesv2.Rdata')
moddat<-mod$data

#raw data from model for plotting
#load(file = 'outputs/Aug2023/Q1_rawdata_full.Rdata')

#assess convergence issues---- 
summary(mod)
#convergence issue only on 'cor' term - re-run without 
#MCMCtrace(mod) #this takes a while- but mostly looks good!

#pull out Unique study ID slopes----
slopes<-coef(object = mod)
slopes<-as.data.frame(slopes$UNIQUE.PAIR.ID)
slopes<-select(slopes, contains("Gn2"))
slopes$UniquePairID<-row.names(slopes)

#save(slopes, file = "outputs/Aug2023/randomslopes_q1modelv2.Rdata")

#run ppchecks ----
ppc<-pp_check(mod) #captures mean and variance ok but misses magnitude of mean
#save(ppc, file='outputs/Aug2023/ppcheck_modq1.Rdata')

load(file = "outputs/Aug2023/ppcheck_modq1.Rdata")

#re-plot for predictive accuracy 
ppcdata<-ppc$data
obsdata<-subset(ppcdata, is_y==TRUE)
preddata<-subset(ppcdata, is_y==FALSE) #10 predictions for every 1 observation 

obsdatax<-select(obsdata, y_id, value)%>%rename(value0=value)
  
alldata<-left_join(obsdatax, preddata)

rm(ppc)
plot(alldata$value~alldata$value0)
save(alldata, file='outputs/Sep2023/Q1_ppc_data.Rdata')

#make plot nicer
pdf(file = "figures/ppc_plot.pdf", width = 8, height = 8)
ggplot(data=alldata, aes(y = value, x=value0)) + 
  geom_point()+ 
  stat_lineribbon()+
  ylab("Predicted") + xlab("Observed") +
  scale_fill_brewer() + theme_bw()
dev.off()


#run correlation test 
ct<-cor.test(alldata$value, alldata$value0)

           
# Check the normal distribution of random effects----
qqnorm(slopes$Estimate.Prop.Change.Gn2, main = "Normal Q-Q plot of random slopes",  bty="n")
qqline(slopes$Estimate.Prop.Change.Gn2, col = 2, lwd = 3, lty = 3) # a bit off on the tails 


#Kfold cross validation----
#https://avehtari.github.io/modelselection/rats_kcv.html#5_K-fold_cross-validation
#use kfold split grouped to preserves group structure within folds but assigns groups randomly to folds 
grps<-loo::kfold_split_grouped(K = 10, x = moddat$UNIQUE.PAIR.ID)
print(grps)
tab<-table(moddat$UNIQUE.PAIR.ID, grps)
colSums(tab)# look pretty even ~49k each 

options(future.globals.maxSize = 8000 * 1024^2)# give a lot of space so kf doesn't crash 


kf<-kfold(mod, folds = grps, save_fits = T ) #running very slowly


#loo-CV subsample- 
#Keep getting ERROR: can't allocate vector this size for both of these
#https://mc-stan.org/loo/articles/loo2-large-data.html
#try this? 
loo_ss_1 <-
  brms::loo_subsample(mod)
print(loo_ss_1)

r2<-loo_R2(mod)


#plot random slopes ----
bars<-ggplot(data=slopes, aes(y = UniquePairID, x=Estimate.Prop.Change.Gn2)) + 
  geom_pointrange(aes(xmin=Q2.5.Prop.Change.Gn2, xmax=Q97.5.Prop.Change.Gn2), size=0.01, alpha=0.5, color='grey')+ 
  geom_vline(xintercept = 0, color='black', lty=2)+
   theme(axis.text.y=element_blank(),  #remove y axis labels
                axis.ticks.y=element_blank()) + xlab("Group level slope estimates")

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



