#Figures with Q1 outputs
# Created by: ENB
# Created: 13 Sept 2023 (ENB)
# Last Modified: 13 Sept 2023 (ENB)


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 0: Read in neccesary packages
rm(list=ls()) 

# # LIBRARIES # #
library(tidyverse)
library("dplyr")
library("tibble")
library("readr") 
library("ggplot2")
library("magrittr")
library(progress)




#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 1: Read in random slopes, get a sense of the data, and set up Q2 analysis 
load("outputs/Aug2023/randomslopes_q1modelv2.Rdata")
head(slopes)
colnames(slopes)
hist(slopes$Estimate.Prop.Change.Gn2)
hist(slopes$Est.Error.Prop.Change.Gn2)

#This histogram is probably something we want to visualize, so let's write that up in 
#ggplot now

# Create a histogram using ggplot2
plain_histogram <- slopes %>%
ggplot(aes(x = Estimate.Prop.Change.Gn2)) +
  geom_histogram(binwidth = 0.05, fill = "grey", color = "black") +
  xlim(-1.1, 1.1)+
  
  # Customize the theme to theme_classic
  theme_classic() +
  
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  
  # Add labels to axes
  labs(
    x = "Random slope of interaction",
    y = "Frequency"
  )
plain_histogram

#So it appears that each pair has been assigned several metrics from Q1
# So we need to read in the pair information to do further analyses 
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23
head(dat)
dat_select <- dat %>%
  select(UNIQUE.PAIR.ID, TAXA1, TAXA2, CLIMATE1, CLIMATE2, REALM1, REALM2, SERIES.l, 
         interaction_present, interaction_benefit, interaction_type)
dat_select <- dat_select %>%
  rename(UniquePairID = UNIQUE.PAIR.ID)
slopes_join <- left_join(slopes, dat_select, by='UniquePairID')
slopes_join <- distinct(slopes_join)
head(slopes_join)

#Since these are all within studies, are there any cases where the 1/2 versions don't match? 
identical(slopes_join$TAXA1, slopes_join$TAXA2)
identical(slopes_join$CLIMATE1, slopes_join$CLIMATE2)
identical(slopes_join$REALM1, slopes_join$REALM2)

#remove duplicate columns
slopes_join <- slopes_join %>%
  select(-TAXA2, -CLIMATE2, -REALM2)
colnames(slopes_join)


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 2: Try out some visualizations to see how we end up 

#First thing we want to visualize: interaction present vs interaction absent
#We could do this using violin plots... 
presence_absence_violin <- slopes_join %>%
  ggplot(aes(x = as.character(interaction_present), y=Estimate.Prop.Change.Gn2)) +
  geom_violin() +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red", position = position_dodge(0.75)) +  # Add mean points
  stat_boxplot(geom = "errorbar", width = 0.25, color = "blue", position = position_dodge(0.75)) +  # Add boxplot
  ylim(-1.1, 1)+
  
  
  # Customize the theme to theme_classic
  theme_classic() +
  
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  ) +
  
  # Add labels to axes
  labs(
    x = "Interaction present",
    y = "Random slope of interaction"
  )
presence_absence_violin

climate_density <- slopes_join %>%
  ggplot(aes(x=Estimate.Prop.Change.Gn2, group=CLIMATE1)) + 
  geom_density(aes(fill=CLIMATE1), alpha=0.7, adjust=0.5)+
  xlim(-1.1, 1)+
  labs(x="Random slope of interaction", y="Density of observations", 
       fill="Climate")+
  theme_classic()+
  
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
climate_density


type_density <- slopes_join %>%
  filter(!interaction_type=="uncategorized_interaction") %>%
  ggplot(aes(x=Estimate.Prop.Change.Gn2, group=interaction_type)) + 
  geom_density(aes(fill=interaction_type), alpha=0.7, adjust=0.5)+
  labs(x="Random slope of interaction", y="Density of observations", 
       fill="Interaction type")+
  xlim(-1.1, 1)+
  theme_classic()+
  
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
type_density



benefit_density <- slopes_join %>%
  filter(!interaction_benefit=="NA") %>%
  ggplot(aes(x=Estimate.Prop.Change.Gn2, group=interaction_benefit)) + 
  geom_density(aes(fill=interaction_benefit), alpha=0.5, adjust=0.5)+
  labs(x="Random slope of interaction", y="Density of observations", 
       fill="Interaction benefit")+
  xlim(-1.1, 1)+
  theme_classic()+
  
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
benefit_density



taxa_density <- slopes_join %>%
  ggplot(aes(x=Estimate.Prop.Change.Gn2, group=TAXA1)) + 
  geom_density(aes(fill=TAXA1), alpha=0.7, adjust=0.5)+
  xlim(-1.1, 1)+
  labs(x="Random slope of interaction", y="Density of observations", 
       fill="Taxa")+
  theme_classic()+
  
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
taxa_density

#Interaction with series length and slope
series_slopes_plot <- slopes_join %>%
  ggplot(aes(x=Estimate.Prop.Change.Gn2, y=as.numeric(SERIES.l))) + 
  geom_point() +
  theme_classic()+
  
  # Increase the size of axis titles and labels
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12)
  )
series_slopes_plot

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 3: Investigate Gavia's linear model
load("~/Documents/Work and Career/LDP/Working Group/Q2.model.wTaxa.fixed.wTreatment.interactionYNxClimate.Rdata")
Q2mod
head(Q2mod$data)
posterior_samples <- posterior_samples(Q2mod, pars = c("b", "sd"))
#calculate the mean and standard error of the mean 
# Calculate the mean of each column
means <- colMeans(posterior_samples)

# Calculate the standard deviation of each column
sds <- apply(posterior_samples, 2, sd)

# Calculate the standard error of each column
ses <- sds / sqrt(nrow(posterior_samples))

summary_df <- data.frame(
  effect = names(means),  # Column names as 'effect'
  mean = means,           # Mean values
  sd = sds,               # Standard deviation values
  se = ses # Standard error values
)
rownames(summary_df) <- NULL

#Get a taxa pair dataset
taxa_pair_betas <- subset(summary_df, grepl("RESOLVED.TAXA.PAIR", effect))
taxa_pair_betas$effect <- gsub(".*RESOLVED.TAXA.PAIR", "", taxa_pair_betas$effect)
taxa_pair_betas <- taxa_pair_betas[!grepl("NA", taxa_pair_betas$effect), ]



#Now plot these 
Q2_plot <- ggplot(taxa_pair_betas, aes(x = reorder(effect, -mean), y = mean)) +
  geom_point() +  # Add points
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +  # Add error bars
  labs(x = "Effect", y = "Mean") +  # Label axes
  theme_classic() +  # Use a minimal theme (you can customize this)
  coord_flip()  # Flip the axes
Q2_plot


#So looks like something we could do with this is model the taxa-level 
