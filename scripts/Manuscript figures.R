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
library(tidybayes)
library(brms)
#library(brmstools)
library(bayesplot)
library(cowplot)


##### Read in Associations, predictions model and Q1 model ####
load("outputs/Oct2023/Q2.model.wTaxa.fixed.wTreatment.abs.lat.scale.Rdata")
load("outputs/Aug2023/randomslopes_q1model.Rdata")
load("outputs/Oct2023/Q2.predictions.model.invtransform.scale.Rdata")

####### read in data and clean #######
# this is same as pre-processing step for when fit Q2 model
#data with interaction info 
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23
dat_terr <- subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
dat_terr <- subset(dat_terr, Metric!="CROSS" &  Type!="Between") %>%
  select(-c("Prop.Change.Gn1", "Prop.Change.Gn2", "YEAR.T","YEAR.T1", "SERIES.start", "SERIES.end")) %>% 
  distinct(.)


slopes.meta <- left_join(slopes, dat_terr, by=c("UniquePairID"="UNIQUE.PAIR.ID"))


# fixing inconsistency in names
slopes.meta$ORGANISMS1 <- ifelse(slopes.meta$ORGANISMS1=="Plants", "plants", slopes.meta$ORGANISMS1)
slopes.meta$ORGANISMS2 <- ifelse(slopes.meta$ORGANISMS2=="Plants", "plants", slopes.meta$ORGANISMS2)

slopes.meta$ORGANISMS1 <- ifelse(slopes.meta$ORGANISMS1=="Acrididae (grasshoppers)", "grasshoppers", slopes.meta$ORGANISMS1)
slopes.meta$ORGANISMS2 <- ifelse(slopes.meta$ORGANISMS2=="Acrididae (grasshoppers)", "grasshoppers", slopes.meta$ORGANISMS2)

slopes.meta$ORGANISMS1 <- ifelse(slopes.meta$ORGANISMS1=="Grasshoppers", "grasshoppers", slopes.meta$ORGANISMS1)
slopes.meta$ORGANISMS2 <- ifelse(slopes.meta$ORGANISMS2=="Grasshoppers", "grasshoppers", slopes.meta$ORGANISMS2)


# Fix: “Oleacina” was assigned the classification of “Bivalvia”, 
# but it’s a terrestrial snail in the Gastropoda
slopes.meta$RESOLVED.TAXA1 [slopes.meta$Gn1 == "Oleacina" ] <- "Gastropoda"
slopes.meta$RESOLVED.TAXA2 [slopes.meta$Gn2 == "Oleacina" ] <- "Gastropoda"


# only pair id where repeat series
meas.two.n <- subset(slopes.meta, SERIES.n==2)$UniquePairID

slopes.series.2 <- slopes.meta[slopes.meta$UniquePairID %in% meas.two.n, ]

# add series lengths together
slopes.meta2 <- slopes.meta %>% 
  group_by(UniquePairID) %>% 
  mutate(SERIES.l.new = ifelse(n()==2, sum(SERIES.l), SERIES.l)) %>% 
  select(-c(SERIES.l, SERIES.n)) %>%
  distinct()
# and back now to the original number of rows of slopes

# if the genera have no known interactions, include that info in the interaction type column
slopes.meta2$interaction_type <- ifelse(slopes.meta2$interaction_present=="0", "no_interaction", slopes.meta2$interaction_type)

# create a resolved taxa pair id column
slopes.meta2$RESOLVED.TAXA.PAIR <- paste0(slopes.meta2$RESOLVED.TAXA1, ".",slopes.meta2$RESOLVED.TAXA2)


#### read in centrality measures ####
# use edge centralities  since they are immediately comparable to the slopes

edge.centrality <- readRDS("outputs/Aug2023/centrality_edges.RDS")


slopes.meta3 <- left_join(slopes.meta2, edge.centrality,
                          by = c("Estimate.Prop.Change.Gn2", "Est.Error.Prop.Change.Gn2", "Q2.5.Prop.Change.Gn2", "Q97.5.Prop.Change.Gn2", "UniquePairID", "Gn1", "Gn2"))
# note that if no known GLOBI interaction, edge centrality is NA, losing a lot of data


#### read in disturbance data ####
# last 3 columns: treatment_yn, treatment_desc, and treatment_simplified. 
# note that not all treatments are necessarily disturbances,
# e.g., there’s at least one fertilization treatment

load("data/prep_biotime/meta_pairs_10km.RData")

meta.pairs$STUDY_ID <- as.character(meta.pairs$STUDY_ID)

slopes.meta4 <- left_join(slopes.meta3, meta.pairs[, c(1,20,46)],
                          by=c("ID1" = "STUDY_ID"))

slopes.meta4$interaction_present <- as.factor(slopes.meta4$interaction_present)

slopes.meta4$interaction_benefit <- ifelse(slopes.meta4$interaction_benefit=="NA", "no interaction", slopes.meta4$interaction_benefit)

# absolute latitude
slopes.meta4$abs.lat <- abs(slopes.meta4$CENT_LAT)

### fix NA resolved taxa pairs ####
# code chunk by EB

#Check the NA values and their corresponding organism values 
dat_na <- slopes.meta4 %>%
  filter(grepl("\\.NA|NA\\.", RESOLVED.TAXA.PAIR)) 

#double check that the NAs match
unique(dat_na$ORGANISMS1)
unique(dat_na$ORGANISMS2)
#Adjust names below 

# taxa 1
slopes.meta4$RESOLVED.TAXA1 <- ifelse(
  is.na(slopes.meta4$RESOLVED.TAXA1),
  ifelse(
    slopes.meta4$ORGANISMS1 %in% c("insects", "grasshoppers"),
    "Insecta",
    ifelse(slopes.meta4$ORGANISMS1 == "birds", "Aves", 
           ifelse(slopes.meta4$ORGANISMS1 == "rodents", "Mammalia", NA))
  ),
  slopes.meta4$RESOLVED.TAXA1
)


# taxa 2
slopes.meta4$RESOLVED.TAXA2 <- ifelse(
  is.na(slopes.meta4$RESOLVED.TAXA2)==TRUE,
  ifelse(
    slopes.meta4$ORGANISMS2 %in% c("insects", "grasshoppers"),
    "Insecta",
    ifelse(slopes.meta4$ORGANISMS2 == "birds", "Aves", 
           ifelse(slopes.meta4$ORGANISMS2 == "rodents", "Mammalia", NA)
    )
  ),
  slopes.meta4$RESOLVED.TAXA2
)

sorted_words <- apply(slopes.meta4[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
slopes.meta4$RESOLVED.TAXA.PAIR <- sorted_words
unique(slopes.meta4$RESOLVED.TAXA.PAIR)
table(slopes.meta4$RESOLVED.TAXA.PAIR)

# drop taxa groups that only have one observation
slopes.meta5 <- subset(slopes.meta4, RESOLVED.TAXA.PAIR!="Monocots.Gnetopsida" & RESOLVED.TAXA.PAIR!="Gnetopsida.Monocots")

#double check if there are any more NAs
dat_na_check <- slopes.meta5 %>%
  filter(grepl("\\.NA|NA\\.", RESOLVED.TAXA.PAIR)) 


###### Calculations for results section #####
paste0("b_", rownames(fixef(Q2mod))[6:29])
Q2mod %>%
  spread_draws(b_Intercept, b_RESOLVED.TAXA.PAIRAves.Bryopsida,
               b_RESOLVED.TAXA.PAIRBryopsida.Aves,            
               b_RESOLVED.TAXA.PAIREudicots.Eudicots, 
               b_RESOLVED.TAXA.PAIREudicots.Gnetopsida,      
               b_RESOLVED.TAXA.PAIREudicots.Magnoliopsida, 
               b_RESOLVED.TAXA.PAIREudicots.Monocots,         
               b_RESOLVED.TAXA.PAIREudicots.Pinopsida,
               b_RESOLVED.TAXA.PAIRGastropoda.Gastropoda,     
               b_RESOLVED.TAXA.PAIRGnetopsida.Eudicots,
               b_RESOLVED.TAXA.PAIRGnetopsida.Magnoliopsida,
               b_RESOLVED.TAXA.PAIRInsecta.Insecta, 
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Eudicots,
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Gnetopsida, 
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Magnoliopsida,
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Monocots, 
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Pinopsida,  
               b_RESOLVED.TAXA.PAIRMammalia.Mammalia,
               b_RESOLVED.TAXA.PAIRMonocots.Eudicots, 
               b_RESOLVED.TAXA.PAIRMonocots.Magnoliopsida, 
               b_RESOLVED.TAXA.PAIRMonocots.Monocots,        
               b_RESOLVED.TAXA.PAIRMonocots.Pinopsida, 
               b_RESOLVED.TAXA.PAIRPinopsida.Eudicots,        
               b_RESOLVED.TAXA.PAIRPinopsida.Magnoliopsida, 
               b_RESOLVED.TAXA.PAIRPinopsida.Monocots) %>%
  mean_qi(aves.aves = b_Intercept,
          Aves.Bryopsida = b_Intercept + b_RESOLVED.TAXA.PAIRAves.Bryopsida,
          Bryopsida.Aves = b_Intercept + b_RESOLVED.TAXA.PAIRBryopsida.Aves,
          Eudicots.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Eudicots,
          Eudicots.Gnetopsida = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Gnetopsida,
          Eudicots.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Magnoliopsida,
          Eudicots.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Monocots,
          Eudicots.Pinopsida = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Pinopsida,
          Gastropoda.Gastropoda = b_Intercept + b_RESOLVED.TAXA.PAIRGastropoda.Gastropoda,
          Gnetopsida.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIRGnetopsida.Eudicots,
          Gnetopsida.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRGnetopsida.Magnoliopsida,
          Insecta.Insecta = b_Intercept + b_RESOLVED.TAXA.PAIRInsecta.Insecta,
          Magnoliopsida.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Eudicots,
          Magnoliopsida.Gnetopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Gnetopsida,
          Magnoliopsida.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Magnoliopsida,
          Magnoliopsida.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Monocots,
          Magnoliopsida.Pinopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Pinopsida,
          Mammalia.Mammalia = b_Intercept + b_RESOLVED.TAXA.PAIRMammalia.Mammalia,
          Monocots.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMonocots.Magnoliopsida,
          Monocots.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIRMonocots.Monocots,
          Monocots.Pinopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMonocots.Pinopsida,
          Pinopsida.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIRPinopsida.Eudicots,
          Pinopsida.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRPinopsida.Magnoliopsida,
          Pinopsida.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIRPinopsida.Monocots
          )  %>% 
  select(seq(1, 75, by=3)) %>% 
  tidyr::pivot_longer(1:25) %>% 
  summarize(min = min(value), max=max(value))
 
Q2predictions.mod %>%
  spread_draws(b_Intercept, b_RESOLVED.TAXA.PAIRAves.Bryopsida,
               b_RESOLVED.TAXA.PAIRBryopsida.Aves,            
               b_RESOLVED.TAXA.PAIREudicots.Eudicots, 
               b_RESOLVED.TAXA.PAIREudicots.Gnetopsida,      
               b_RESOLVED.TAXA.PAIREudicots.Magnoliopsida, 
               b_RESOLVED.TAXA.PAIREudicots.Monocots,         
               b_RESOLVED.TAXA.PAIREudicots.Pinopsida,
               b_RESOLVED.TAXA.PAIRGastropoda.Gastropoda,     
               b_RESOLVED.TAXA.PAIRGnetopsida.Eudicots,
               b_RESOLVED.TAXA.PAIRGnetopsida.Magnoliopsida,
               b_RESOLVED.TAXA.PAIRInsecta.Insecta, 
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Eudicots,
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Gnetopsida, 
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Magnoliopsida,
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Monocots, 
               b_RESOLVED.TAXA.PAIRMagnoliopsida.Pinopsida,  
               b_RESOLVED.TAXA.PAIRMammalia.Mammalia,
               b_RESOLVED.TAXA.PAIRMonocots.Eudicots, 
               b_RESOLVED.TAXA.PAIRMonocots.Magnoliopsida, 
               b_RESOLVED.TAXA.PAIRMonocots.Monocots,        
               b_RESOLVED.TAXA.PAIRMonocots.Pinopsida, 
               b_RESOLVED.TAXA.PAIRPinopsida.Eudicots,        
               b_RESOLVED.TAXA.PAIRPinopsida.Magnoliopsida, 
               b_RESOLVED.TAXA.PAIRPinopsida.Monocots) %>%
  mean_qi(aves.aves = b_Intercept,
          Aves.Bryopsida = b_Intercept + b_RESOLVED.TAXA.PAIRAves.Bryopsida,
          Bryopsida.Aves = b_Intercept + b_RESOLVED.TAXA.PAIRBryopsida.Aves,
          Eudicots.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Eudicots,
          Eudicots.Gnetopsida = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Gnetopsida,
          Eudicots.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Magnoliopsida,
          Eudicots.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Monocots,
          Eudicots.Pinopsida = b_Intercept + b_RESOLVED.TAXA.PAIREudicots.Pinopsida,
          Gastropoda.Gastropoda = b_Intercept + b_RESOLVED.TAXA.PAIRGastropoda.Gastropoda,
          Gnetopsida.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIRGnetopsida.Eudicots,
          Gnetopsida.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRGnetopsida.Magnoliopsida,
          Insecta.Insecta = b_Intercept + b_RESOLVED.TAXA.PAIRInsecta.Insecta,
          Magnoliopsida.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Eudicots,
          Magnoliopsida.Gnetopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Gnetopsida,
          Magnoliopsida.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Magnoliopsida,
          Magnoliopsida.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Monocots,
          Magnoliopsida.Pinopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMagnoliopsida.Pinopsida,
          Mammalia.Mammalia = b_Intercept + b_RESOLVED.TAXA.PAIRMammalia.Mammalia,
          Monocots.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMonocots.Magnoliopsida,
          Monocots.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIRMonocots.Monocots,
          Monocots.Pinopsida = b_Intercept + b_RESOLVED.TAXA.PAIRMonocots.Pinopsida,
          Pinopsida.Eudicots = b_Intercept + b_RESOLVED.TAXA.PAIRPinopsida.Eudicots,
          Pinopsida.Magnoliopsida = b_Intercept + b_RESOLVED.TAXA.PAIRPinopsida.Magnoliopsida,
          Pinopsida.Monocots = b_Intercept + b_RESOLVED.TAXA.PAIRPinopsida.Monocots
  )  %>% 
  select(seq(1, 75, by=3)) %>% 
  tidyr::pivot_longer(1:25) %>% 
  summarize(min = min(value), max=max(value))


########## Main text figures #########
# Figure generation

#Add custom theme
my.theme<-theme(axis.text=element_text(size=20),
                axis.title = element_text(size = 25),
                legend.text=element_text(size=10),
                legend.title = element_text(size=12),
                plot.title = element_text(face="bold",size=14,margin=margin(0,0,20,0),hjust = 0.5),
                axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
                axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))


###### Figure 2. Q1 results #######
# Q1 predictive accuracy
load(file = "outputs/Sep2023/Q1_ppc_data.Rdata")

a<-ggplot(data = pred_estimates, aes(x=value0, y=mean_pred, colour=diff))+
  #geom_smooth(method='lm')+
  geom_abline(slope=1, intercept=0, color="darkblue", lty=2)+
  geom_point(alpha=1)+ 
  ylab(" Predicted log change in abundance") + xlab("Observed log change in abundance")  +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  scale_colour_gradient(
    low = "#8DA0CB",
    #mid = "blue",
    high = "#ff0000",
    #midpoint = 0.15,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour", name="Residuals")+
  annotate("text", label="R^2 = 0.38
             Pearson's R = 0.62
             95% CI (0.608, 0.639)
             t = 61.869, df = 6023
             p-value < 2.2e-16",

           x=5.2, y=-1.2, size=5, hjust=1) 

# "R^2 == 0.38
# Pearson's R == 0.62
#               95% CI (0.608, 0.639)
#               t = 61.869, df = 6023
#               p-value < 2.2e-16"


b<-ggplot(subset(pred_estimates,diff<0.25), aes(x=value0, y=mean_pred, colour=diff))+
  #geom_smooth(method='lm')+
  geom_abline(slope=1, intercept=0, color="black", lty=2)+
  geom_point(alpha=1)+
  ylab("Predicted log change \nin abundance") + 
  #xlab("Observed log change in abundance")  +
  xlab("")  +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15)) +
  scale_colour_gradient(
    low = "#8DA0CB",
    #mid = "blue",
    high = "#Ff0000",
    #midpoint = 0.15,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
  aesthetics = "colour", name="Residuals")#+
  # theme(axis.title = element_text(size = 10), 
  #       legend.title = element_text(size = 10))

check2.1<-subset(pred_estimates, diff<0.25)#2870/6025 ~47% 
c<-ggplot(check2.1, aes(x=value0)) + 
  geom_histogram(aes(y=..count..), colour="white", fill="#AFA7AE")+
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  geom_density(alpha=.2, fill="#AFA7AE") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  xlab("Observed log change in abundance")+
  ylab("Number of observations")+ xlim(-1,1) #+
  # theme(axis.title = element_text(size = 10))

top_row <- plot_grid(b, align = "h", axis = "l", ncol = 1)

bottom_row <- plot_grid(c, align = "h", axis = "l", ncol = 2, rel_widths = c(1, .19))

bc <- plot_grid(top_row, bottom_row, ncol = 1)

abc <- plot_grid(a, bc, ncol=2, rel_widths = c(3, 2))
#pdf(file = "figures/Figure2.pdf", width = 12, height = 8)

######  Figure 3. Q2 coef plot  #######

# with posterior distribution
# associations model
fig3 <- Q2mod %>%
  spread_draws(b_Intercept, b_series.scale, b_lat.scale, b_treatment_ynyes, b_interaction_present1) %>%
  tidyr::pivot_longer(cols=4:8, names_to = "draw.name") %>% 
  ggplot(aes(y = draw.name, x = value)) +
  #geom_vline(xintercept = 0, color = "#839496", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey")+
  stat_halfeye(aes(fill=draw.name), 
               .width = c(0.66, 0.95),
               point_size=3, 
               point_interval = "mean_qi",
               interval_size_range=c(1,2)) +
  labs(x = expression("Association strength"),
       y = "Model parameters") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        #text = element_text(family = "Ubuntu")
        ) +
  theme_bw() +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  scale_y_discrete(
    labels = c("b_Intercept" = "Intercept",
               "b_series.scale" = "Time series length",
               "b_lat.scale" = "Latitude",
               "b_treatment_ynyes" = "Disturbance:yes",
               "b_interaction_present1" = "Interaction:yes"),
    limits = c("b_series.scale", 
               "b_treatment_ynyes",
               "b_interaction_present1",
               "b_lat.scale", 
               "b_Intercept" )) +
  guides(fill = "none") +
  scale_fill_manual(values=c("#009E73", "grey50","#0072B2" ,"#E69F00" ,"#CC79A7"))
fig3

######  Figure 4. Q2 coef plot  #######


######  Figure 5. Q2 linear predictor plots  #######



### Sup mat figures ####
###### Barplot of genera by taxa in only terrestrial ####
# added by GLL for figure 2 (methods conceptual figure)
# only terrestrial
# only within

dat_terr1 <- subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
dat_terr2 <- subset(dat_terr1, Metric!="CROSS" &  Type!="Between") %>%
  select(-c("Prop.Change.Gn1", "Prop.Change.Gn2", "YEAR.T","YEAR.T1", "SERIES.start", "SERIES.end")) %>% 
  distinct(.)

# Fix: “Oleacina” was assigned the classification of “Bivalvia”, 
# but it’s a terrestrial snail in the Gastropoda
dat_terr2$RESOLVED.TAXA1 [dat_terr2$Gn1 == "Oleacina" ] <- "Gastropoda"
dat_terr2$RESOLVED.TAXA2 [dat_terr2$Gn2 == "Oleacina" ] <- "Gastropoda"



dat.terr3 <- dat_terr2

dat.terr3$RESOLVED.TAXA.PAIR <- paste0(dat.terr3$RESOLVED.TAXA1, ".",dat.terr3$RESOLVED.TAXA2)

# fix resolved taxa names
#Check the NA values and their corresponding organism values 
dat_na2 <- dat.terr3 %>%
  filter(grepl("\\.NA|NA\\.", RESOLVED.TAXA.PAIR)) 

#double check that the NAs match
unique(dat_na2$ORGANISMS1)
unique(dat_na2$ORGANISMS2)
#Adjust names below 

# taxa 1
dat.terr3$RESOLVED.TAXA1 <- ifelse(
  is.na(dat.terr3$RESOLVED.TAXA1),
  ifelse(
    dat.terr3$ORGANISMS1 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat.terr3$ORGANISMS1 == "birds", "Aves", 
           ifelse(dat.terr3$ORGANISMS1 == "rodents", "Mammalia", NA)
    )
  ),
  dat.terr3$RESOLVED.TAXA1
)


# taxa 2
dat.terr3$RESOLVED.TAXA2 <- ifelse(
  is.na(dat.terr3$RESOLVED.TAXA2)==TRUE,
  ifelse(
    dat.terr3$ORGANISMS2 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat.terr3$ORGANISMS2 == "birds", "Aves", 
           ifelse(dat.terr3$ORGANISMS2 == "rodents", "Mammalia", NA)
    )
  ),
  dat.terr3$RESOLVED.TAXA2
)

sorted_words2 <- apply(dat.terr3[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
dat.terr3$RESOLVED.TAXA.PAIR <- sorted_words2
unique(dat.terr3$RESOLVED.TAXA.PAIR)
table(dat.terr3$RESOLVED.TAXA.PAIR)

data.frame(table(dat.terr3$RESOLVED.TAXA.PAIR))


library(forcats)
library(rcartocolor)

options(scipen = 999) #converts to nice numbers on axes
# colors per taxa
taxaCol <- rcartocolor::carto_pal(11, "Safe")

ggplot(data.frame(table(dat.terr3$RESOLVED.TAXA1)), 
       aes(x=fct_reorder(Var1, Freq, .desc=TRUE), 
           y=log(Freq), fill=Var1)) +
  #geom_bar(stat = "identity", fill="#66CC66") +
  #geom_col( fill="#3382BF") +
  geom_col() +
  labs(x = "Taxanomic category", y = "Number of genera") +
  theme_classic(base_size = 25) + 
  theme(legend.position = "none",
        axis.text=element_text(size=25)) +
  theme(aspect.ratio = 0.4) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values = taxaCol) 

###### World map where size of point corresponds to length of time series ####
drawWorld<-function(lats) {
  world_map<-map_data("world")
  
  g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
  g1<-g1+geom_polygon(data=world_map, aes(x=long, y=lat, group=group), colour="gray60", fill="gray60")
  g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
               panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
               legend.position="none",axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
  
  if(lats=="y") {
    g1<-g1+geom_hline(yintercept=23.5, colour="red")+geom_hline(yintercept =-23.5, colour="red")
    g1<-g1+geom_hline(yintercept=66.5, colour="darkblue")+geom_hline(yintercept =-66.5, colour="darkblue")
  }
  else { return(g1) }
  return(g1)
}

# Now let's use the above function to plot these studies across the globe
(gplot <- drawWorld("y") + 
    geom_point(data=left_join(dat.terr3, meta.pairs, by=c("ID1" = "STUDY_ID")), 
               aes(x=CENT_LONG, y=CENT_LAT, colour = RESOLVED.TAXA1, size = SERIES.l), 
               alpha=0.1)) +
  scale_colour_manual(values=taxaCol)


###### Predictions model coef plot, main hypotheses #####
Q2predictions.mod %>%
  spread_draws(b_Intercept, b_series.scale, b_lat.scale, b_treatment_ynyes, b_interaction_present1) %>%
  tidyr::pivot_longer(cols=4:8, names_to = "draw.name") %>% 
  ggplot(aes(y = draw.name, x = value)) +
  #geom_vline(xintercept = 0, color = "#839496", size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey")+
  stat_halfeye(aes(fill=draw.name), 
               .width = c(0.66, 0.95),
               point_size=3, 
               point_interval = "mean_qi",
               interval_size_range=c(1,2)) +
  labs(x = expression("Predictive accuracy (1/MAE)"),
       y = "Model parameters") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  scale_y_discrete(
    labels = c("b_Intercept" = "Intercept",
               "b_series.scale" = "Time series length",
               "b_lat.scale" = "Latitude",
               "b_treatment_ynyes" = "Disturbance:yes",
               "b_interaction_present1" = "Interaction:yes"),
    limits = c("b_series.scale", 
               "b_treatment_ynyes",
               "b_interaction_present1",
               "b_lat.scale", 
               "b_Intercept" )) +
  guides(fill = "none") +
  scale_fill_manual(values=c("#009E73", "grey50","#0072B2" ,"#E69F00" ,"#CC79A7"))


###### Barplot genera pair by latitude ######
ggplot(data=slopes.meta5, aes(x=abs.lat)) +
  #geom_vline(xintercept = 43, linetype = "dashed", color = "darkgrey")+
  geom_histogram(binwidth = 1.2, bins = 40, fill="grey50", colour="white") +
  labs(x = expression("Latitude"),
       y = "Frequency") +
  scale_y_log10() +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  guides(fill = "none") #+
  #scale_fill_manual(values=c("#66C2A5", "grey50","#8DA0CB" ,"#FC8D62" ,"#E78AC3"))




#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
###############################################################################
#Emily's Figure 4, Sup Fig 4, and Fig 5 Code


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


