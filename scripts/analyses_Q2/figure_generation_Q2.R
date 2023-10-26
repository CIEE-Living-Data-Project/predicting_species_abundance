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

  

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 3: Investigate Gavia's linear model
load("~/Documents/Work and Career/LDP/Working Group/Q2.model.wTaxa.fixed.wTreatment.abs.lat.scale.Rdata")
Q2mod <- Q2mod
head(Q2mod$data)
Q2mod$formula
Q2mod$fit

#Get the Q2 dat 
Q2_dat <- Q2mod$data

#get unique series sample values 

#Using emmeans, extract the marginal effects
lat_means <- Q2mod %>%
  emmeans(~lat.scale )
taxa_means <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR, 
          level=0.95,
          at = list(
                    interaction_present = '0', 
                    treatment_yn = 'no'))
#Get unique latitude values 
lat_unique <- unique(Q2mod$data$lat.scale)
lat_values <- unique(round(Q2mod$data$lat.scale, digits = 1))
lat_abs_values <- unique(round(slopes_join$CENT_LAT, digits = 1))

hist(Q2mod$data$lat.scale)
hist(slopes_join$CENT_LAT)
#Get scale values
max(Q2mod$data$series.scale)
min(Q2mod$data$series.scale)
scale_values <- round(Q2mod$data$series.scale, digits = 1)
scale_values_list <- c(-1.55, -1, -0.5, 0, 0.5,1, 1.5, 2, 2.5, 3)

#taxa means at different climates
taxa_means_clim <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR + lat.scale + series.scale,
          level=0.95,
          at = list(lat.scale = lat_values,
                series.scale = scale_values_list, 
            interaction_present = '0', 
        treatment_yn = 'no')
          )

#get the opposite means
taxa_means_clim_opposite <- Q2mod %>%
  emmeans(~RESOLVED.TAXA.PAIR + lat.scale + series.scale,
          level=0.95,
          at = list(lat.scale = lat_values,
                    series.scale = scale_values_list, 
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

taxa_means_clim_opposite_df <- taxa_means_clim_opposite_df %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR)
#add a column indicating treatment is yes
taxa_means_clim_opposite_df$treatment_yn <- c("yes")

taxa_means_clim_all <- bind_rows(taxa_means_clim_df, taxa_means_clim_opposite_df)
#rename the slopes column 
taxa_means_clim_all <- taxa_means_clim_all %>%
  rename(scale.lat = lat.scale)

#Add a new column in slopes_join saying whether the average value of their 
#series length is closer to the min, mean, or max
average_series_length_slopes <- Q2_dat %>%
  group_by(RESOLVED.TAXA.PAIR) %>%
  summarize(mean_sl = median(series.scale))

average_series_length_slopes$assigned_sl <- sapply(average_series_length_slopes$mean_sl, function(x) {
  closest_value <- scale_values_list[which.min(abs(x - scale_values_list))]
  closest_value
})

#Build the values back into slopes_join 
Q2_join_join <- left_join(Q2_dat, average_series_length_slopes, by=c('RESOLVED.TAXA.PAIR'))
Q2_join_join <- Q2_join_join %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR )
Q2_join_join <- Q2_join_join %>%
  mutate(lat.scale = round(lat.scale, digits = 1))
taxa_means_clim_all <- taxa_means_clim_all %>%
  rename(assigned_sl = series.scale) %>%
  rename(lat.scale = scale.lat)
slopes_join_stats_all <- left_join(Q2_join_join, taxa_means_clim_all, 
                                   by=c('resolved_taxa_pair', 'lat.scale' ,'treatment_yn', 'assigned_sl'))

#Unscale the latitude 
#Replace slopes_join resolved taxa bivalvia
slopes_join$resolved_taxa_pair <- gsub("Gastropoda.Bivalvia", "Gastropoda.Gastropoda", slopes_join$resolved_taxa_pair)
slopes_join$resolved_taxa_pair <- gsub("Bivalvia.Gastropoda", "Gastropoda.Gastropoda", slopes_join$resolved_taxa_pair)

slopes_join_stats_all <- left_join(slopes_join_stats_all, slopes_join[, c(1,2, 15, 19)], 
                                   by=c("Estimate.Prop.Change.Gn2", 'Est.Error.Prop.Change.Gn2', "resolved_taxa_pair"))




#Add a latitude grouping variable
# Define the hex codes for light green and dark green
green_colors <- colorRampPalette(c("#e6f6ff", "#006199"))(6)

#Make a abs.lat column
slopes_join_stats_all <- slopes_join_stats_all %>%
  filter(!is.na(emmean)) %>%
  mutate(abs.lat = round(CENT_LAT, digits = 0))
slopes_join_stats_all$abs.lat <- factor(slopes_join_stats_all$abs.lat, levels = c(18, 34, 39, 42, 44, 45))
#Get the average latitude of each group
average_lat <- slopes_join_stats_all %>%
  group_by(resolved_taxa_pair) %>%
  summarize(mean_lat = mean(CENT_LAT)) %>%
  arrange(desc(mean_lat))
average_lat


#Rearrange the taxa groups to be by plant-plant, plant-animal, and animal-animal
unique(slopes_join_stats_all$resolved_taxa_pair)
interaction_list <- c(
  "Gnetopsida.Monocots",   "Monocots.Gnetopsida","Magnoliopsida.Gnetopsida", "Gnetopsida.Magnoliopsida",
  "Eudicots.Gnetopsida","Gnetopsida.Eudicots","Eudicots.Eudicots", 
  "Eudicots.Pinopsida", "Pinopsida.Eudicots","Monocots.Eudicots", "Eudicots.Monocots", "Monocots.Pinopsida", "Pinopsida.Monocots",
  "Magnoliopsida.Eudicots", "Eudicots.Magnoliopsida", "Pinopsida.Magnoliopsida", "Magnoliopsida.Pinopsida",
  "Monocots.Magnoliopsida", "Magnoliopsida.Monocots", "Magnoliopsida.Magnoliopsida", "Monocots.Monocots",
  "Bryopsida.Aves", "Aves.Bryopsida", 
  "Gastropoda.Gastropoda","Mammalia.Mammalia", "Aves.Aves", "Insecta.Insecta" 
)

slopes_join_stats_all$resolved_taxa_pair <- factor(slopes_join_stats_all$resolved_taxa_pair, 
                                                   levels = interaction_list)

slopes_join_stats_all$abs.lat <- factor(slopes_join_stats_all$abs.lat, levels = c(18, 34, 39, 42, 44, 45))

#Add custom theme
my.theme<-theme(axis.text=element_text(size=12),
                axis.title = element_text(size = 14),
                legend.text=element_text(size=10),
                legend.title = element_text(size=12),
                plot.title = element_text(face="bold",size=14,margin=margin(0,0,20,0),hjust = 0.5),
                axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
                axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))


figure_4_noviolin <- slopes_join_stats_all %>%
  filter(interaction_present == '0') %>%
  #mutate(resolved_taxa_pair = fct_reorder(resolved_taxa_pair, CENT_LAT, .fun='max')) %>%
  ggplot(aes(x = resolved_taxa_pair, y = as.numeric(Estimate.Prop.Change.Gn2))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")+
  #geom_violin(alpha = 0.95, bw=0.04, trim=FALSE, position = position_dodge(width = 1), width = 1) +  # Violin plot by CLIMATE1
  geom_errorbar(aes(x = resolved_taxa_pair, ymin = lower.HPD, ymax = upper.HPD, group = abs.lat),
                color = 'black', position = position_dodge(width = 1), width = 1) +
  geom_point(aes(x = resolved_taxa_pair, y = emmean, group = abs.lat, color = abs.lat),  size = 3, 
             position = position_dodge(width = 1)) +
  labs(x = "Taxonomic category", y = "Strength of association", colour = "Latitude") +
  scale_color_manual(values = green_colors) +
  scale_fill_manual(values = green_colors) +
  ylim(-0.3, 0.65) +
  coord_flip() +
  facet_grid(~treatment_yn) + 
  theme_bw()+
  #guides(fill = "none") +
  theme(
    panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
    panel.grid.minor.x = element_blank(), 
    strip.text.x = element_blank()# Hide minor x-axis grid lines
  )+
  my.theme
figure_4_noviolin
ggsave("figures/figure_4_noviolin_10262023.png", plot = figure_4_noviolin, width = 11, height = 7, units = 'in')
ggsave("figures/figure_4_noviolin_10262023.pdf", plot = figure_4_noviolin, width = 11, height = 7, units = 'in')



#Do again but with the predictive model 
#Part 3: Investigate Gavia's linear model
load("~/Documents/Work and Career/LDP/Working Group/Q2.predictions.model.invtransform.scale.Rdata")
head(Q2predictions.mod$data)
Q2predictions.mod$formula
Q2predictions.mod$fit
Q2predictions_dat <- Q2predictions.mod$data


#Using emmeans, extract the marginal effects
lat_means_predictions <- Q2predictions.mod %>%
  emmeans(~lat.scale )
taxa_means_predictions <- Q2predictions.mod %>%
  emmeans(~RESOLVED.TAXA.PAIR, 
          level=0.95,
          at = list(
            interaction_present = '0', 
            treatment_yn = 'no'))
#Get unique latitude values 
lat_unique <- unique(Q2predictions.mod$data$lat.scale)
lat_values <- unique(round(Q2predictions.mod$data$lat.scale, digits = 1))
lat_abs_values <- unique(round(slopes_join$CENT_LAT, digits = 1))

hist(Q2predictions.mod$data$lat.scale)


#taxa means at different climates
taxa_means_clim_predictions <- Q2predictions.mod %>%
  emmeans(~RESOLVED.TAXA.PAIR + lat.scale + series.scale,
          level=0.95,
          at = list(lat.scale = lat_values,
                    series.scale = scale_values_list, 
                    interaction_present = '0', 
                    treatment_yn = 'no')
  )
#plot the taxa means 

#get the opposite means
taxa_means_clim_predictions_opposite <- Q2predictions.mod %>%
  emmeans(~RESOLVED.TAXA.PAIR + lat.scale + series.scale,
          level=0.95,
          at = list(lat.scale = lat_values,
                    series.scale = scale_values_list, 
                    interaction_present = '0', 
                    treatment_yn = 'yes')
  )




#save as dataframe 
taxa_means_clim_predictions_df <- as.data.frame(taxa_means_clim_predictions)
taxa_means_clim_predictions_opposite_df <- as.data.frame(taxa_means_clim_predictions_opposite)

#rename some variables
taxa_means_clim_predictions_df <- taxa_means_clim_predictions_df %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR)
taxa_means_clim_predictions_df$treatment_yn <- c("no")

taxa_means_clim_predictions_opposite_df <- taxa_means_clim_predictions_opposite_df %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR)
#add a column indicating treatment is yes
taxa_means_clim_predictions_opposite_df$treatment_yn <- c("yes")

taxa_means_clim_predictions_all <- bind_rows(taxa_means_clim_predictions_df, taxa_means_clim_predictions_opposite_df)
#rename the slopes column 
taxa_means_clim_predictions_all <- taxa_means_clim_predictions_all %>%
  rename(assigned_sl = series.scale) %>%
  rename(scale.lat = lat.scale)

#Add a new column in slopes_join saying whether the average value of their 
#series length is closer to the min, mean, or max
average_series_length_slopes <- Q2predictions_dat %>%
  group_by(RESOLVED.TAXA.PAIR) %>%
  summarize(mean_sl = median(series.scale))

average_series_length_slopes$assigned_sl <- sapply(average_series_length_slopes$mean_sl, function(x) {
  closest_value <- scale_values_list[which.min(abs(x - scale_values_list))]
  closest_value
})

#Build the values back into slopes_join 
Q2predictions_join <- left_join(Q2predictions_dat, average_series_length_slopes, by=c('RESOLVED.TAXA.PAIR'))
Q2predictions_join <- Q2predictions_join %>%
  rename(resolved_taxa_pair = RESOLVED.TAXA.PAIR )
Q2predictions_join <- Q2predictions_join %>%
  mutate(lat.scale = round(lat.scale, digits = 1))
taxa_means_clim_predictions_all <- taxa_means_clim_predictions_all %>%
  #rename(assigned_sl = series.scale) %>%
  rename(lat.scale = scale.lat)
slopes_join_predictions_all <- left_join(Q2predictions_join, taxa_means_clim_predictions_all,
                                         by=c('resolved_taxa_pair', 'lat.scale' ,'treatment_yn', 'assigned_sl'))

#get latitude data from original data
load("outputs/Sep2023/Q1_ppc_data.Rdata")
head(pred_estimates_pairid)
pred_estimates_pairid <- pred_estimates_pairid %>%
  rename(UniquePairID = UNIQUE.PAIR.ID)
pred_estimates_withlat <- left_join(pred_estimates_pairid,slopes_join[, c(1,2, 5, 15, 19)],  
                                    by=c("UniquePairID"))

slopes_join_predictions_all <- left_join(slopes_join_predictions_all, pred_estimates_withlat[, c(5,6, 13)], 
                                   by=c("mean_diff", "sd_diff"))



#Add a latitude grouping variable
# Define the hex codes for light green and dark green
green_colors <- colorRampPalette(c("#e6f6ff", "#006199"))(6)

#Make a abs.lat column
slopes_join_predictions_all <- slopes_join_predictions_all %>%
  filter(!is.na(emmean)) %>%
  mutate(abs.lat = round(CENT_LAT, digits = 0))
slopes_join_predictions_all$abs.lat <- factor(slopes_join_predictions_all$abs.lat, levels = c(18, 34, 39, 42, 44, 45))



#Rearrange the taxa groups to be by plant-plant, plant-animal, and animal-animal

slopes_join_predictions_all$resolved_taxa_pair <- factor(slopes_join_predictions_all$resolved_taxa_pair, 
                                                   levels = interaction_list)

slopes_join_predictions_all$abs.lat <- factor(slopes_join_predictions_all$abs.lat, levels = c(18, 34, 39, 42, 44, 45))


figure_4_noviolin_predictions <- slopes_join_predictions_all %>%
  filter(interaction_present == '0') %>%
  #mutate(resolved_taxa_pair = fct_reorder(resolved_taxa_pair, CENT_LAT, .fun='max')) %>%
  ggplot(aes(x = resolved_taxa_pair, y = as.numeric(mean_diff))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")+
  #geom_violin(alpha = 0.95, bw=0.04, trim=FALSE, position = position_dodge(width = 1), width = 1) +  # Violin plot by CLIMATE1
  geom_errorbar(aes(x = resolved_taxa_pair, ymin = lower.HPD, ymax = upper.HPD, group = abs.lat),
                color = 'black', position = position_dodge(width = 1), width = 1) +
  geom_point(aes(x = resolved_taxa_pair, y = emmean, group = abs.lat, color = abs.lat),  size = 3, 
             position = position_dodge(width = 1)) +
  labs(x = "Taxonomic category", y = "Predictive accuracy (1/MAE)", colour = "Latitude") +
  scale_color_manual(values = green_colors) +
  scale_fill_manual(values = green_colors) +
  #ylim(-0.3, 0.65) +
  coord_flip() +
  facet_grid(~treatment_yn) + 
  theme_bw()+
  #guides(fill = "none") +
  theme(
    panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
    panel.grid.minor.x = element_blank(), 
    strip.text.x = element_blank()# Hide minor x-axis grid lines
  )+
  my.theme
figure_4_noviolin_predictions
ggsave("figures/figure_4_noviolin_predictions_10262023.png", plot = figure_4_noviolin_predictions, width = 11, height = 7, units = 'in')
ggsave("figures/figure_4_noviolin_predictions_10262023.pdf", plot = figure_4_noviolin_predictions, width = 11, height = 7, units = 'in')



#Figure 5
#Pick four representative groups: two plant, one animal, one animal-plant
#Aves/Aves
#Insecta/Insecta
#Eudicots/Magnoliopsida
#Aves/Bryopsida 

#Get the different series lengths for the different taxa groups
table(slopes_join$resolved_taxa_pair, slopes_join$SERIES.l)

#Get our slopes table down to the groups in question
selected_pairs <- c("Aves.Aves", "Insecta.Insecta","Mammalia.Mammalia", 
                    "Aves.Bryopsida", 
                    "Bryopsida.Aves", "Eudicots.Eudicots",  
                    "Monocots.Magnoliopsida", "Magnoliopsida.Monocots", 
                    "Magnoliopsida.Magnoliopsida")
select_groups <- slopes_join %>%
  filter(resolved_taxa_pair %in% selected_pairs)
#Add an abs.lat column
select_groups$abs_lat <- round(select_groups$CENT_LAT, digits =0)
#Plot the relationships vs the series lengths 
select_groups$factor_groups <- paste(select_groups$abs_lat, select_groups$treatment_yn, sep = "-")
unique(select_groups$factor_groups)

#now get four blue colours and four orange colours
green_colors_2 <- colorRampPalette(c("#b3e3ff", "#005180"))(4)
orange_colours <- colorRampPalette(c("#ffcf66", "#996900"))(4)
all_colours <- c(green_colors_2,orange_colours )

#Set up the linetype factoring
factor_list <- c("18-no", "34-no", "39-no", "45-no", 
                 "34-yes", "39-yes","42-yes", "44-yes")
select_groups$factor_groups <- factor(select_groups$factor_groups, 
                                                   levels = factor_list)

#Come up with a factor table 

#Get summary r-squared values for each group
# Calculate the overall R-squared for your linear model

series_length_plot <- select_groups %>%
  ggplot(aes(x=SERIES.l, y=Estimate.Prop.Change.Gn2, color = factor(factor_groups), 
             shape = factor(treatment_yn))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey")+
  geom_point(colour ='darkgrey', alpha = 0.5, size = 2,
             position = position_jitter(width = 0.1, height = 0)) + 
  geom_smooth(method = "lm", se = TRUE, 
              aes(group = factor_groups, fill = factor_groups), size = 2) +
  facet_wrap(~factor(resolved_taxa_pair, levels = selected_pairs))+
  scale_color_manual(values = all_colours) +
  scale_fill_manual(values = all_colours) +
  ylim(-1, 1)+
  labs(x = "Length of study (years)", y ="Strength of Association", shape = "Disturbance", 
       color = "Latitude - Disturbance", fill ="Latitude - Disturbance" )+
  theme_classic()+
  my.theme+ 
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  theme(legend.key.width = unit(2, "cm"))
series_length_plot
ggsave("figures/figure_5_10262023.png", plot = series_length_plot, width = 11, height = 7, units = 'in')
ggsave("figures/figure_5_10262023.pdf", plot = series_length_plot, width = 11, height = 7, units = 'in')



