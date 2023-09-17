#### Script to re-create+format MS main plots ####


# script created 17 Sept 2023 by MTW
# last updated ______________________

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


