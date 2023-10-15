#### Script to re-create+format MS main plots ####


# script created 17 Sept 2023 by MTW
# last updated ______________________

#### Exploratory Figures not used in MS ####
# Goals for Sept 17th and following days: 
# 1) Try to decide on a clean visualization of the genus-genus slopes. Decide whether this should be a panel of Fig 2 or 
# 2) Conceptual figure - in Illustrator, not R
# 3) Start working on cleaning up 

# Loading output from the Q1 model - version 1 is the correct one to be using, not v2. 
load(file = "outputs/Aug2023/randomslopes_q1model.Rdata")
View(slopes)
#go back and look- what are the "Q2.5..." columns? 

#Pull in the BioTIME data:
# From Emily's script, "figure_generation_Q2.R":
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23
head(dat)
dat_select <- dat %>%
  select(UNIQUE.PAIR.ID, TAXA1, TAXA2, CLIMATE1, CLIMATE2, REALM1, REALM2, SERIES.l, 
         interaction_present, interaction_benefit, interaction_type)
dat_select <- dat_select %>%
  rename(UniquePairID = UNIQUE.PAIR.ID)
slopes_join <- left_join(slopes, dat_select, by='UniquePairID')
head(slopes_join)

#From Gavia's script, "analysis with Q1 outputs.R":
dat_terr <- subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
dat_terr <- subset(dat_terr, Metric!="CROSS" &  Type!="Between") %>%
  select(-c("Prop.Change.Gn1", "Prop.Change.Gn2", "YEAR.T","YEAR.T1", "SERIES.start", "SERIES.end")) %>% 
  distinct(.)

#An important difference here to move forward with is that Gavia subset for terrestrial and within-study, which we want to do for all figures since that was a decision we made for the analysis.

#now my new code is below:
head(dat_terr)
dat_terr <- dat_terr %>%
  rename(UniquePairID = UNIQUE.PAIR.ID)
slopes_terr <- left_join(slopes, dat_terr, by='UniquePairID')
head(slopes_terr)

# Create a density plot of taxa 2 association strengths
#What is a given taxa2's strength of association with taxa 1, and how are these association strengths distributed across the dataset? Are they more positive than negative, do some pairs have more/less variation in their associations than other pairs?
terr_density <- slopes_terr %>%
  ggplot(aes(x = Estimate.Prop.Change.Gn2, group = RESOLVED.TAXA2, fill = RESOLVED.TAXA2)) +
  geom_density(adjust=1.5, alpha=.4) +
  xlim(-1.1, 1.1)+
  labs(x = "Strength of association with Class 1", y = "Frequency", fill = "Class 2") +
  theme_classic() +
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)) +
  #faceting by Climate (temperate vs tropical) and taxa 1. 
  facet_wrap(~ CLIMATE1 + RESOLVED.TAXA1);terr_density 

#faceted by climate only
terr_density_Climate <- slopes_terr %>%
  ggplot(aes(x = Estimate.Prop.Change.Gn2, group = RESOLVED.TAXA2, fill = RESOLVED.TAXA2)) +
  geom_density(adjust=1.5, alpha=.4) +
  xlim(-1.1, 1.1)+
  labs(x = "Strength of association", y = "Frequency", fill = "Class 2") +
  theme_classic() +
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)) +
  #faceting by Climate (temperate vs tropical) and taxa 1. 
  facet_wrap(~ CLIMATE1);terr_density_Climate 

#Color Brewer - get the hex codes for the "paired" color set.We're considering these because they're colorblind friendly and fairly B&W freindly.
library(RColorBrewer)

brewer.pal(8, "Set2") #gives the hex values for each color.



#### Figure 5 plot ####
#load the Q2 model output, emailed by Gavia on 9/17/2023 and added locally to my WG folder (since it's too big to push/pull)
# load("C:/Users/miatw/Desktop/predicting_species_abundance/Q2.model.wTaxa.fixed.wTreatment.interactionYNxClimate.Rdata") = this is the old model

# Notes from Gavia's email 10/10/23:
#1. Continuous predictors are now centred and scaled so that they are directly comparable, these will need to be back transformed for ease of interpretation
#2. In the predictions model, the response is inverse transformed to meet assumptions of normality, so coefficients will need to be back transformed before interpretation 

#The models that were run are of the form ~ scale(SERIES.l.new) + scale(abs.lat) + treatment_yn + interaction_present + RESOLVED.TAXA.PAIR. After a cursory examination, there are some slight differences between the models, namely if the interaction_present term has credible intervals that cross 0. So this will be interesting to talk about.

#For the figures, proceed with the associations model, as it is less weird than the other one.
load("C:/Users/miatw/Desktop/predicting_species_abundance/Q2.model.wTaxa.fixed.wTreatment.abs.lat.scaled.Rdata")
Q2mod <- Q2mod.assoc
head(Q2mod$data)
Q2mod$data
Q2mod$formula
Q2mod$fit 

# figure 5: association strength x time = Estimate.Prop.Change.Gn2 x SERIES.1.new (? or should I use the scale(Series.1.new)? ) ...and try to color by latitude? Facet by latitude?

lat_values <- unique(round(Q2mod$data$abs.lat, digits = 0)); lat_values 
# 6 different latitudes
## Actually I want to try to plot latitude as a 3rd axis, but first I will get the association x time 2D plot:

#Set up workspace with Emily's code from figure_generation_Q2.R, then:

#Hypothesis: As time series length increases, we should see an increase in association strength.

#Overall:
ggplot(slopes_join_stats_all, 
       aes(x = SERIES.l, y = emmean)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", linewidth = 1, se = FALSE) +
  labs(x = "Time series length (years)", 
       y = "Association strength") +
  #ggtitle("") +
  theme_classic()
#We actually see a decrease now?? Is that right?? With the previous model, I thought we were getting positive? Have I done this incorrectly?

#How does this trend differ across latitudes?
slopes_join_stats_all$abs.lat <- as.factor(slopes_join_stats_all$abs.lat)
str(slopes_join_stats_all)

fig5_time_lat <- ggplot(slopes_join_stats_all, 
       aes(x = SERIES.l, y = emmean, colour = abs.lat)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", linewidth = 1, se = FALSE) +
  labs(x = "Time series length (years)", 
       y = "Association strength") +
  #ggtitle("") +
  theme_classic(); fig5_time_lat #that's incredibly ugly, but some latitudes show a positive effect, others are negative

#by Taxa1
ggplot(slopes_join_stats_all, 
       aes(x = SERIES.l, y = emmean, colour = TAXA1)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", linewidth = 1, se = FALSE) +
  labs(x = "Time series length (years)", 
       y = "Association strength") +
  #ggtitle("") +
  theme_classic() #hmm. These are just the rough taxa groups of genus 1, but I was curious

#across taxa pairs?
ggplot(slopes_join_stats_all, 
       aes(x = SERIES.l, y = emmean, colour = resolved_taxa_pair)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", linewidth = 1, se = FALSE) +
  labs(x = "Time series length (years)", 
       y = "Association strength") +
  #ggtitle("") +
  theme_classic() # this is too messy to use, but interesting to see... within a pair, as the time series increases, there lots of instances of neg effects, and some neutral, a couple are positive.

#But am I even using the right emmeans? I'm pretty sure I am, but I'm going to double check before trying to interpret this too deeply.

fig5_time_pairs <- ggplot(slopes_join_stats_all, 
       aes(x = SERIES.l, y = emmean, colour = resolved_taxa_pair)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", linewidth = 1, se = FALSE) +
  labs(x = "Time series length (years)", 
       y = "Association strength") +
  #ggtitle("") +
  theme_classic() + 
  facet_wrap(~abs.lat) #...Weird! Only pairs at 34 latitude show a change (mostly very negative!) in association strength as time increases. Is there something different about studies at this latitude? (or again, do I need to calculate the emmeans differently than Emily did...)

#Also I want to be using SERIES.1.new... that's used in Emily's emmeans, but needs to be added (from Gavia's latest model output) to the slopes_join_stats_all df, then I could use it. Although I'm not sure how much of a difference there is between that and the original SERIES.1? That's a Q for Gavia.

#here it is with all the raw estimates from the model instead of the marginal means
fig5_time_pairs_raw <- ggplot(slopes_join_stats_all, 
       aes(x = SERIES.l, y = Estimate.Prop.Change.Gn2, colour = resolved_taxa_pair)) + 
  geom_point(size = 2, alpha = 0.7) + 
  geom_smooth(method = "lm", linewidth = 1, se = FALSE) +
  labs(x = "Time series length (years)", 
       y = "Association strength") +
  #ggtitle("") +
  theme_classic() + 
  facet_wrap(~abs.lat) 

fig5_time_pairs
fig5_time_pairs_raw

#AH, I wonder if the solution is just that I need to backtransform the data, and then recalculate the means. I'm going to sleep and will look more into this in the morning. 
