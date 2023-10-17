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
library("emmeans")
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
dat <- left_join(dat, meta.pairs[, c(1,20, 46)],
                        by=c("ID1" = "STUDY_ID"))

dat_select <- dat %>%
  select(UNIQUE.PAIR.ID, TAXA1, TAXA2, CLIMATE1, CLIMATE2, REALM1, REALM2, SERIES.l, 
         RESOLVED.TAXA1, RESOLVED.TAXA2, ORGANISMS1, ORGANISMS2, UNIQUE.PAIR.ID, treatment_yn, CENT_LAT, 
         interaction_present, interaction_benefit, interaction_type, interaction_present)


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

#Resolve taxa names in slopes_join 
slopes_join$RESOLVED.TAXA1 <- ifelse(
  is.na(slopes_join$RESOLVED.TAXA1),
  ifelse(
    slopes_join$ORGANISMS1 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(slopes_join$ORGANISMS1 == "birds", "Aves", 
           ifelse(slopes_join$ORGANISMS1 == "rodents", "Mammalia", NA)
    )
  ),
  slopes_join$RESOLVED.TAXA1
)



slopes_join$RESOLVED.TAXA2 <- ifelse(
  is.na(slopes_join$RESOLVED.TAXA2),
  ifelse(
    slopes_join$ORGANISMS2 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(slopes_join$ORGANISMS2 == "birds", "Aves", 
           ifelse(slopes_join$ORGANISMS2 == "rodents", "Mammalia", NA)
    )
  ),
  slopes_join$RESOLVED.TAXA2
)



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
  filter(CLIMATE1=="Temperate", interaction_present=="0")

ggplot(slopes_baseline, aes(x = Estimate.Prop.Change.Gn2)) +
  geom_density(color = "black", aes(fill=resolved_taxa_pair)) +
  facet_wrap(~resolved_taxa_pair , ncol = 2) +
  labs(x = "Value", y = "Frequency") +
  xlim(-1, 1)+
  theme_classic()
  

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 3: Investigate Gavia's linear model
load("~/Documents/Work and Career/LDP/Working Group/Q2.model.wTaxa.fixed.wTreatment.abs.lat.scaled.Rdata")
Q2mod <- Q2mod.assoc
head(Q2mod$data)
Q2mod$formula
Q2mod$fit
#posterior_samples <- posterior_samples(Q2mod, pars = c("b", "sd"))

#get unique series sample values 

#Using emmeans, extract the marginal effects
lat_means <- Q2mod %>%
  emmeans(~abs.lat )
taxa_means <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR, 
          level=0.95,
          at = list(
                    interaction_present = '0', 
                    treatment_yn = 'no'))
#Get unique latitude values 
lat_values <- unique(round(Q2mod$data$abs.lat, digits = 0))
#taxa means at different climates
taxa_means_clim <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR + abs.lat + SERIES.l.new,
          level=0.95,
          at = list(abs.lat = c(39, 34, 42, 45, 18, 44),
                SERIES.l.new = c(10, 15, 22, 30, 45), 
            interaction_present = '0', 
        treatment_yn = 'no')
          )
#plot the taxa means 

#get the opposite means
taxa_means_clim_opposite <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR + abs.lat + SERIES.l.new,
          level=0.95,
          at = list(abs.lat = c(39, 34, 42, 45, 18, 44), 
                    SERIES.l.new = c(10, 15, 22, 30, 45), 
            interaction_present = '0', 
                    treatment_yn = 'yes')
  )



#save as dataframe 
taxa_means_clim_df <- as.data.frame(taxa_means_clim)
taxa_means_clim_opposite_df <- as.data.frame(taxa_means_clim_opposite)
#rename some variables
taxa_means_clim_df <- taxa_means_clim_df %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR)
taxa_means_clim_df$treatment_yn <- c("no")
#add a rounded latitude column to slopes_join 
slopes_join <- slopes_join %>%
  mutate(abs.lat = round(CENT_LAT, digits = 0))

taxa_means_clim_opposite_df <- taxa_means_clim_opposite_df %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR)
#add a column indicating treatment is yes
taxa_means_clim_opposite_df$treatment_yn <- c("yes")

taxa_means_clim_all <- rbind(taxa_means_clim_df, taxa_means_clim_opposite_df)
#rename the slopes column 
taxa_means_clim_all <- taxa_means_clim_all %>%
  rename(assigned_sl = SERIES.l.new)

#Add a new column in slopes_join saying whether the average value of their 
#series length is closer to the min, mean, or max
average_series_length_slopes <- slopes_join %>%
  group_by(resolved_taxa_pair) %>%
  summarize(mean_sl = median(SERIES.l))

reference_values <- c(10, 15, 22, 30, 45)
average_series_length_slopes$assigned_sl <- sapply(average_series_length_slopes$mean_sl, function(x) {
  closest_value <- reference_values[which.min(abs(x - reference_values))]
  closest_value
})

#Build the values back into slopes_join 
slopes_join <- left_join(slopes_join, average_series_length_slopes, by=c('resolved_taxa_pair'))
slopes_join_stats_all <- left_join(slopes_join, taxa_means_clim_all, by=c('resolved_taxa_pair', 'abs.lat' ,'treatment_yn', 'assigned_sl'))

#Filter the data for a violin plot 
slopes_join_stats_all_filtered <- slopes_join_stats_all %>%
  filter(interaction_present == '0') %>%
  group_by(resolved_taxa_pair, abs.lat, treatment_yn) %>%
  filter(n() >= 2) %>%
  ungroup()

#Add custom theme
my.theme<-theme(axis.text=element_text(size=12),
                axis.title = element_text(size = 14),
                legend.text=element_text(size=10),
                legend.title = element_text(size=12),
                plot.title = element_text(face="bold",size=14,margin=margin(0,0,20,0),hjust = 0.5),
                axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
                axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))

#Add a latitude grouping variable
# Define the hex codes for light green and dark green
green_colors <- colorRampPalette(c("#e9f0ff", "#0d3da8"))(6)

slopes_join_stats_all_filtered$abs.lat <- factor(slopes_join_stats_all_filtered$abs.lat, levels = c(18, 34, 39, 42, 44, 45))

supp.labs <- c("No Disturbance", "Disturbance")
names(supp.labs) <- c("no", "yes")
figure_4 <- slopes_join_stats_all_filtered %>%
  filter(interaction_present == '0') %>%
  ggplot(aes(x = reorder(resolved_taxa_pair, -as.numeric(abs.lat)), y = as.numeric(Estimate.Prop.Change.Gn2),
             fill = abs.lat)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")+
  geom_violin(alpha = 0.95, bw=0.04, trim=FALSE, position = position_dodge(width = 1), width = 1) +  # Violin plot by CLIMATE1
  geom_errorbar(aes(x = resolved_taxa_pair, ymin = lower.HPD, ymax = upper.HPD, group = abs.lat),
                color = 'black', position = position_dodge(width = 1), width = 1) +
  geom_point(aes(x = resolved_taxa_pair, y = emmean, group = abs.lat),  size = 3, 
             position = position_dodge(width = 1), colour = 'black') +
  labs(x = "Taxonomic group", y = "Strength of association", fill = 'Latitude') +
  scale_color_manual(values = green_colors) +
  scale_fill_manual(values = green_colors) +
  ylim(-1, 1) +
  coord_flip() +
  facet_grid(~treatment_yn, 
             labeller = labeller(treatment_yn = supp.labs)) + 
  theme_bw()+
  #guides(fill = "none") +
  theme(
    panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
    panel.grid.minor.x = element_blank()   # Hide minor x-axis grid lines
  )+
  my.theme
figure_4
ggsave("figures/figure_4.png", plot = figure_4, width = 11, height = 7, units = 'in')

#Rearrange the taxa groups to be by plant-plant, plant-animal, and animal-animal
unique(slopes_join_stats_all$resolved_taxa_pair)
interaction_list <- c(
  "Eudicots.Eudicots", "Eudicots.Gnetopsida",
  "Eudicots.Magnoliopsida", "Eudicots.Monocots", "Eudicots.Pinopsida", "Gnetopsida.Eudicots",
  "Gnetopsida.Monocots", "Gnetopsida.Magnoliopsida", "Monocots.Eudicots", "Monocots.Gnetopsida", 
  "Monocots.Monocots", "Monocots.Magnoliopsida", "Monocots.Pinopsida", "Magnoliopsida.Eudicots", 
  "Magnoliopsida.Gnetopsida",  "Magnoliopsida.Magnoliopsida",
  "Magnoliopsida.Monocots", "Magnoliopsida.Pinopsida", "Pinopsida.Magnoliopsida",
  "Pinopsida.Monocots", "Pinopsida.Eudicots", "Bryopsida.Aves", "Aves.Bryopsida", 
  "Aves.Aves", "Bivalvia.Gastropoda",
  "Gastropoda.Bivalvia","Gastropoda.Gastropoda", "Insecta.Insecta", "Mammalia.Mammalia"
)

slopes_join_stats_all$resolved_taxa_pair <- factor(slopes_join_stats_all$resolved_taxa_pair, 
                                                   levels = interaction_list)

slopes_join_stats_all$abs.lat <- factor(slopes_join_stats_all$abs.lat, levels = c(18, 34, 39, 42, 44, 45))
figure_4_noviolin <- slopes_join_stats_all %>%
  filter(interaction_present == '0') %>%
  mutate(resolved_taxa_pair = fct_reorder(resolved_taxa_pair, CENT_LAT, .fun='max')) %>%
  ggplot(aes(x = resolved_taxa_pair, y = as.numeric(Estimate.Prop.Change.Gn2))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")+
  #geom_violin(alpha = 0.95, bw=0.04, trim=FALSE, position = position_dodge(width = 1), width = 1) +  # Violin plot by CLIMATE1
  geom_errorbar(aes(x = resolved_taxa_pair, ymin = lower.HPD, ymax = upper.HPD, group = abs.lat),
                color = 'black', position = position_dodge(width = 1), width = 1) +
  geom_point(aes(x = resolved_taxa_pair, y = emmean, group = abs.lat, color = abs.lat),  size = 3, 
             position = position_dodge(width = 1)) +
  labs(x = "Taxonomic groups", y = "Strength of association", colour = "Latitude") +
  scale_color_manual(values = green_colors) +
  scale_fill_manual(values = green_colors) +
  ylim(-0.3, 0.65) +
  coord_flip() +
  facet_grid(~treatment_yn, 
             labeller = labeller(treatment_yn = supp.labs)) + 
  theme_bw()+
  #guides(fill = "none") +
  theme(
    panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
    panel.grid.minor.x = element_blank()   # Hide minor x-axis grid lines
  )+
  my.theme
figure_4_noviolin
ggsave("figures/figure_4_noviolin.png", plot = figure_4_noviolin, width = 11, height = 7, units = 'in')




custom_colors_2 <- c("yes" = "maroon", "no" = "darkblue")
slopes_join_stats_all %>%
  filter(interaction_present == '0') %>%
  ggplot(aes(x = reorder(resolved_taxa_pair, -emmean), y = as.numeric(Estimate.Prop.Change.Gn2),
             fill = treatment_yn)) +
  geom_violin(alpha = 0.35, bw=0.05, position='dodge') +  # Violin plot by CLIMATE1
  geom_errorbar(aes(x = resolved_taxa_pair, ymin = lower.HPD, ymax = upper.HPD),
                color = 'black', position = position_dodge(width = 0.7), width = 1) +
  geom_point(aes(x = resolved_taxa_pair, y = emmean), color = 'black', size = 2, 
             position = position_dodge(width = 0.7)) +
  labs(x = "Resolved taxa pair", y = "Strength of association", colour = 'Disturbance') +
  scale_color_manual(values = custom_colors_2) +
  scale_fill_manual(values = custom_colors_2, guide = guide_legend(override.aes = list(fill = "grey"))) +
  ylim(-1, 1) +
  coord_flip() +
  theme_classic()+
  guides(fill = "none") 




#Do the same thing but with marginal effects for climate, interaction, disturbance, etc
climate_means <- Q2mod %>%
  emmeans(~CLIMATE1, 
          level=0.95,
          at = list(treatment_yn = 'no', 
                    interaction_present = '0'))
climate_means <- as.data.frame(climate_means)
interaction_means <- Q2mod %>%
  emmeans(~interaction_present, 
          level=0.95,
          at = list(CLIMATE1 = 'Temperate', 
                    treatment_yn = 'no'))
interaction_means <- as.data.frame(interaction_means)
treatment_means <- Q2mod %>%
  emmeans(~treatment_yn, 
          at = list(interaction_present = '0', 
                    CLIMATE1 = 'Temperate'))
treatment_means <- as.data.frame(treatment_means)

treatment_means$treatment_yn <- as.character(treatment_means$treatment_yn)
slopes_join$interaction_present <- as.character(slopes_join$interaction_present)

slopes_join_treatment_stats <- left_join(slopes_join, treatment_means, by=c('treatment_yn'))
df_sampled_treatment <- slopes_join_treatment_stats %>%
  group_by(interaction_present, treatment_yn, resolved_taxa_pair, CLIMATE1) %>%
  mutate(row_number = row_number()) %>%  # Add a row_number column
  filter(row_number <= min(150, n())) %>%  # Filter to keep up to 100 points per group
  ungroup() %>%  # Ungroup the data frame
  select(-row_number)  

#get some example pairs to highlight 
sample_pairs_treatment <-  df_sampled_treatment %>%
  filter(interaction_present == '0', CLIMATE1 == 'Temperate') 

df_sampled_treatment %>%
  filter(interaction_present == '0', CLIMATE1 == 'Temperate') %>%
  ggplot(aes(x = treatment_yn, y = as.numeric(Estimate.Prop.Change.Gn2))) +
  geom_point(alpha = 0.1, position = position_jitter(width = 0.1), colour = 'darkblue')+
  geom_errorbar(aes(x =treatment_yn, ymin = lower.HPD, ymax = upper.HPD),
                position = position_dodge(width = 0.3), width = 0.3, colour='black') +
  geom_point(aes(x = treatment_yn, y = emmean), size = 4, colour='black') +
  labs(x = "Treatment", y = "Strength of association") +
  ylim(-1, 1) +
  coord_flip() +
  theme_classic()+
  guides(fill = "none") 


interaction_means$interaction_present <- as.character(interaction_means$interaction_present)
slopes_join$interaction_present <- as.character(slopes_join$interaction_present)
slopes_join_interaction_stats <- left_join(slopes_join, interaction_means, by=c('interaction_present'))
df_sampled_interaction <- slopes_join_interaction_stats %>%
  group_by(interaction_present, resolved_taxa_pair, CLIMATE1) %>%
  mutate(row_number = row_number()) %>%  # Add a row_number column
  filter(row_number <= min(50, n())) %>%  # Filter to keep up to 100 points per group
  ungroup() %>%  # Ungroup the data frame
  select(-row_number)  

#get some example pairs to highlight 
sample_pairs_interaction <-  df_sampled_interaction %>%
  filter(treatment_yn=="no", CLIMATE1 == 'Temperate') 

df_sampled_interaction %>%
  filter(treatment_yn=="no", CLIMATE1 == 'Temperate') %>%
  ggplot(aes(x = interaction_present, y = as.numeric(Estimate.Prop.Change.Gn2))) +
  geom_point(alpha = 0.1, position = position_jitter(width = 0.1), colour = 'darkblue')+
  geom_errorbar(aes(x = interaction_present, ymin = lower.HPD, ymax = upper.HPD),
                position = position_dodge(width = 0.3), width = 0.3, colour='black') +
  geom_point(aes(x = interaction_present, y = emmean), size = 4, colour='black') +
  labs(x = "Interaction present", y = "Strength of association") +
  ylim(-1, 1) +
  coord_flip() +
  theme_classic()+
  guides(fill = "none")







#Now plot these 
Q2_plot <- ggplot(taxa_pair_betas, aes(x = reorder(effect, -mean), y = mean)) +
  geom_point() +  # Add points
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2) +  # Add error bars
  labs(x = "Taxa pair", y = "Strength of Association") +  # Label axes
  theme_classic() +  # Use a minimxaal theme (you can customize this)
  coord_flip()  # Flip the axes
Q2_plot


#So looks like something we could do with this is model the taxa-level 
