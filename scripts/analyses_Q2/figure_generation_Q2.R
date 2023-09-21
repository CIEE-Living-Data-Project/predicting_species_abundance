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
load("outputs/Aug2023/randomslopes_q1model.Rdata")
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
#read in the disturbance data
load("data/prep_biotime/meta_pairs_10km.RData")
meta.pairs$STUDY_ID <- as.character(meta.pairs$STUDY_ID)
dat <- left_join(dat, meta.pairs[, c(1,46)],
                        by=c("ID1" = "STUDY_ID"))

dat_select <- dat %>%
  select(UNIQUE.PAIR.ID, TAXA1, TAXA2, CLIMATE1, CLIMATE2, REALM1, REALM2, SERIES.l, 
         RESOLVED.TAXA1, RESOLVED.TAXA2, PairID,
         interaction_present, interaction_benefit, interaction_type, treatment_yn)


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

#Get taxa pairs similar to Gavia's 
sorted_words <- apply(slopes_join[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
slopes_join$resolved_taxa_pair <- sorted_words
unique(slopes_join$resolved_taxa_pair)


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
  labs(x="Strength of association of taxa pairs", y="Density of observations", 
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

#Just investigate eudicots eudicots
slopes_BG <- slopes_join %>%
  filter(resolved_taxa_pair=="Bivalvia.Gastropoda") 
unique(slopes_BG$PairID)
#54_54
slopes_GB <- slopes_join %>%
  filter(resolved_taxa_pair=="Gastropoda.Bivalvia") 
unique(slopes_GB$PairID)
#54_54

slopes_MP <- slopes_join %>%
  filter(resolved_taxa_pair=="Magnoliopsida.Pinopsida") 
unique(slopes_MP$PairID)
#355, 240
slopes_PM <- slopes_join %>%
  filter(resolved_taxa_pair=="Pinopsida.Magnoliopsida") 
unique(slopes_PM$PairID)
#355,240


#filter to only temperate, and no disturbance 
slopes_baseline <- slopes_join %>%
  filter(CLIMATE1=="Temperate", treatment_yn=="no")

ggplot(slopes_baseline, aes(x = Estimate.Prop.Change.Gn2)) +
  geom_density(color = "black", aes(fill=resolved_taxa_pair)) +
  facet_wrap(~resolved_taxa_pair , ncol = 2) +
  labs(x = "Value", y = "Frequency") +
  xlim(-1, 1)+
  theme_classic()
  

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 3: Investigate Gavia's linear model
load("~/Documents/Work and Career/LDP/Working Group/Q2.model.wTaxa.fixed.wTreatment.interactionYNxClimate.Rdata")
Q2mod
head(Q2mod$data)
Q2mod$formula
Q2mod$fit
posterior_samples <- posterior_samples(Q2mod, pars = c("b", "sd"))
coef_bays <- coef(Q2mod)

#Using emmeans, extract the marginal effects
climate_means <- Q2mod %>%
  emmeans(~CLIMATE1 )
taxa_means <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR )
#taxa means at different climates
taxa_means_clim <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR + CLIMATE1,
          level=0.95,
          at = list(interaction_present = '0', 
        treatment_yn = 'no')
          )
#plot the taxa means 


#save as dataframe 
taxa_means_clim_df <- as.data.frame(taxa_means_clim)
#rename some variables
taxa_means_clim_df <- taxa_means_clim_df %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR)
#merge with slopes join
slopes_join_stats <- left_join(slopes_join, taxa_means_clim_df, by=c('resolved_taxa_pair', 'CLIMATE1'))

df_sampled <- slopes_join_stats %>%
  group_by(resolved_taxa_pair, CLIMATE1) %>%
  mutate(row_number = row_number()) %>%  # Add a row_number column
  filter(row_number <= min(100, n())) %>%  # Filter to keep up to 100 points per group
  ungroup() %>%  # Ungroup the data frame
  select(-row_number)  

custom_colors <- c("Temperate" = "#2E3C4C", "Tropical" = "#006400")

df_sampled %>%
  filter(interaction_present == '0',
         treatment_yn == 'no') %>%
ggplot(aes(x = reorder(resolved_taxa_pair, -emmean), y=as.numeric(Estimate.Prop.Change.Gn2), group=CLIMATE1)) +
  geom_point(colour = 'grey', alpha=0.2,
             position = position_jitter(width = 0.08, height = 0), stroke=1) +
  #geom_hline(yintercept = 0.236448, linetype = "dotted", color = "black")+
  geom_errorbar(aes(x = resolved_taxa_pair, ymin = lower.HPD, ymax = upper.HPD, color = CLIMATE1),
                position = position_dodge(width = 0.3), width = 1) +
  geom_point(aes(x=resolved_taxa_pair, y=emmean, colour = CLIMATE1), size = 2)+
  labs(x = "Resolved taxa pair", y = "Strength of association", color = 'Climate') +
  scale_color_manual(values = custom_colors)+
  ylim(-1, 1)+
  coord_flip()+
  theme_classic()

#with density
slopes_join_stats %>%
  filter(interaction_present == '0', treatment_yn == 'no') %>%
  ggplot(aes(x = reorder(resolved_taxa_pair, -emmean), y = as.numeric(Estimate.Prop.Change.Gn2), fill = CLIMATE1)) +
  geom_violin(alpha = 0.4, bw=0.05) +  # Violin plot by CLIMATE1
  geom_errorbar(aes(x = resolved_taxa_pair, ymin = lower.HPD, ymax = upper.HPD, color = CLIMATE1),
                position = position_dodge(width = 0.3), width = 1) +
  geom_point(aes(x = resolved_taxa_pair, y = emmean, colour = CLIMATE1), size = 2) +
  labs(x = "Resolved taxa pair", y = "Strength of association", colour = 'Climate') +
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors, guide = guide_legend(override.aes = list(fill = "grey"))) +
  ylim(-1, 1) +
  coord_flip() +
  theme_classic()+
  guides(fill = "none") 

#Do the same thing but with marginal effects for climate, interaction, disturbance, etc
climate_means <- as.data.frame(climate_means)
interaction_means <- Q2mod %>%
  emmeans(~interaction_present )
interaction_means <- as.data.frame(interaction_means)
treatment_means <- Q2mod %>%
  emmeans(~treatment_yn)
treatment_means <- as.data.frame(treatment_means)





#Now plot these 
Q2_plot <- ggplot(taxa_pair_betas, aes(x = reorder(effect, -mean), y = mean)) +
  geom_point() +  # Add points
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +  # Add error bars
  labs(x = "Taxa pair", y = "Strength of Association") +  # Label axes
  theme_classic() +  # Use a minimal theme (you can customize this)
  coord_flip()  # Flip the axes
Q2_plot


#So looks like something we could do with this is model the taxa-level 
