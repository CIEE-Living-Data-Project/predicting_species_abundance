#### script to run Q2 analyses #######
#### what drives variation in the associations between log proportional abundance changes betweeen genera[i] and genera[j]
#### where associations are the random slopes from Q1 analysis

# script created 25 August 2023 by GLL
# last updated ______________________

library(brms)
library(dplyr)
library(ggplot2)
library(tidybayes)

#### read in model output from Q1 ####
# only terrestrial, only within study comparisons
load("~/Documents/14 U of T/CIEE/predicting_species_abundance/outputs/Aug2023/randomslopes_q1model.Rdata")
#load("~/Documents/14 U of T/CIEE/predicting_species_abundance/outputs/Aug2023/intercept_only_models.Rdata")

head(slopes)
hist(slopes$Estimate.Prop.Change.Gn2)
hist(slopes$Est.Error.Prop.Change.Gn2)


#### read in unique pair id metadata ####
#data with interaction info 
dat <- readRDS('data/data_processing/log.prop.change.interactions.RDS')#full cleaned 6/6/23
dat_terr <- subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
dat_terr <- subset(dat_terr, Metric!="CROSS" &  Type!="Between") %>%
  select(-c("Prop.Change.Gn1", "Prop.Change.Gn2", "YEAR.T","YEAR.T1", "SERIES.start", "SERIES.end")) %>% 
  distinct(.)


slopes.meta <- left_join(slopes, dat_terr, by=c("UniquePairID"="UNIQUE.PAIR.ID"))
# some unique pair ids duplicated bc have 2 series, series 1 and series 2
# so about 1000 more rows than original slopes dataframe 
# estimated slopes are identical per unique pair id across both series 1 and series 2
# for those unique pair IDs that have two instances of a series
# therefore add series length together for both series

# fixing inconsistency in names
slopes.meta$ORGANISMS1 <- ifelse(slopes.meta$ORGANISMS1=="Plants", "plants", slopes.meta$ORGANISMS1)
slopes.meta$ORGANISMS2 <- ifelse(slopes.meta$ORGANISMS2=="Plants", "plants", slopes.meta$ORGANISMS2)

slopes.meta$ORGANISMS1 <- ifelse(slopes.meta$ORGANISMS1=="Acrididae (grasshoppers)", "grasshoppers", slopes.meta$ORGANISMS1)
slopes.meta$ORGANISMS2 <- ifelse(slopes.meta$ORGANISMS2=="Acrididae (grasshoppers)", "grasshoppers", slopes.meta$ORGANISMS2)

slopes.meta$ORGANISMS1 <- ifelse(slopes.meta$ORGANISMS1=="Grasshoppers", "grasshoppers", slopes.meta$ORGANISMS1)
slopes.meta$ORGANISMS2 <- ifelse(slopes.meta$ORGANISMS2=="Grasshoppers", "grasshoppers", slopes.meta$ORGANISMS2)


##### exploratory plots ####
# series info
ggplot(slopes.meta, aes(SERIES.l, Estimate.Prop.Change.Gn2, 
                        colour=as.factor(interaction_present))) + geom_point()

ggplot(slopes.meta, aes(as.factor(SERIES.n), Estimate.Prop.Change.Gn2, 
                        colour=as.factor(interaction_present))) + geom_point()

# only pair id where repeat series
meas.two.n <- subset(slopes.meta, SERIES.n==2)$UniquePairID

slopes.series.2 <- slopes.meta[slopes.meta$UniquePairID %in% meas.two.n, ]

ggplot(slopes.series.2, aes(as.factor(SERIES.n), Estimate.Prop.Change.Gn2)) + geom_point()
# okay slopes are identical for series that have two instances 

# add series lengths together
slopes.meta2 <- slopes.meta %>% 
  group_by(UniquePairID) %>% 
  mutate(SERIES.l.new = ifelse(n()==2, sum(SERIES.l), SERIES.l)) %>% 
  select(-c(SERIES.l, SERIES.n)) %>%
  distinct()
# and back now to the original number of rows of slopes
  

# dist
ggplot(slopes.meta2, aes(dist, Estimate.Prop.Change.Gn2)) + geom_point() # all 0 because its a within studies model

# climate
ggplot(slopes.meta2, aes(CLIMATE1, Estimate.Prop.Change.Gn2)) + geom_boxplot() 
ggplot(slopes.meta2, aes(CLIMATE2, Estimate.Prop.Change.Gn2)) + geom_boxplot() # identical bc no within study comparison across climates

# birds vs mammals vs inverts vs plants
ggplot(slopes.meta2, aes(TAXA1, Estimate.Prop.Change.Gn2)) + geom_boxplot() 
ggplot(slopes.meta2, aes(TAXA2, Estimate.Prop.Change.Gn2)) + geom_boxplot() # identical

# more specific common name eg birds, woodland vegetation, plants
# seems like a silly way of splitting the data
# lots of inconsistencies eg grasshoppers vs insects, rodents vs small mammals
ggplot(slopes.meta2, aes(ORGANISMS1, Estimate.Prop.Change.Gn2)) + geom_boxplot() 
ggplot(slopes.meta2, aes(ORGANISMS2, Estimate.Prop.Change.Gn2)) + geom_boxplot() # identical

# scientific classification eg aves
ggplot(slopes.meta2, aes(RESOLVED.TAXA1, Estimate.Prop.Change.Gn2)) + geom_boxplot() 
ggplot(slopes.meta2, aes(RESOLVED.TAXA2, Estimate.Prop.Change.Gn2)) + geom_boxplot() 

# interactions
ggplot(slopes.meta2, aes(as.factor(interaction_present), Estimate.Prop.Change.Gn2)) + geom_boxplot() 
ggplot(slopes.meta2, aes(as.factor(interaction_benefit), Estimate.Prop.Change.Gn2)) + geom_boxplot() 
ggplot(slopes.meta2, aes(as.factor(interaction_type), Estimate.Prop.Change.Gn2)) + geom_boxplot() 

# if the genera have no known interactions, include that info in the interaction type column
slopes.meta2$interaction_type <- ifelse(slopes.meta2$interaction_present=="0", "no_interaction", slopes.meta2$interaction_type)

# create a resolved taxa pair id column
slopes.meta2$RESOLVED.TAXA.PAIR <- paste0(slopes.meta2$RESOLVED.TAXA1, ".",slopes.meta2$RESOLVED.TAXA2)

# think about if this makes sense to make pairs reciprocal
# I don't think so? random slopes are vastly different depending on direction
# plot random slopes first
# NA.Insecta
# NA.Aves
# NA.Mammalia
# Eudicots.Magnoliopsida
# Monocots.Magnoliopsida
# Gastropoda.Bivalvia
# Pinopsida.Magnoliopsida

#### read in centrality measures ####
# use edge centralities  since they are immediately comparable to the slopes

edge.centrality <- readRDS("/Users/Gavia/Documents/14 U of T/CIEE/predicting_species_abundance/outputs/Aug2023/centrality_edges.RDS")


slopes.meta3 <- left_join(slopes.meta2, edge.centrality,
                          by = c("Estimate.Prop.Change.Gn2", "Est.Error.Prop.Change.Gn2", "Q2.5.Prop.Change.Gn2", "Q97.5.Prop.Change.Gn2", "UniquePairID", "Gn1", "Gn2"))
# note that if no known GLOBI interaction, edge centrality is NA, losing a lot of data


#### read in disturbance data ####
# last 3 columns: treatment_yn, treatment_desc, and treatment_simplified. 
# note that not all treatments are necessarily disturbances,
# e.g., there’s at least one fertilization treatment

load("data/prep_biotime/meta_pairs_10km.RData")

meta.pairs$STUDY_ID <- as.character(meta.pairs$STUDY_ID)

slopes.meta4 <- left_join(slopes.meta3, meta.pairs[, c(1,46)],
                          by=c("ID1" = "STUDY_ID"))

slopes.meta4$interaction_present <- as.factor(slopes.meta4$interaction_present)

#### fit model ####
FAM <- gaussian(link = 'identity')

# centrality measures
# interaction as yes/no
MODFORM <- bf(Estimate.Prop.Change.Gn2|resp_se(Est.Error.Prop.Change.Gn2, sigma = TRUE) ~ 
                SERIES.l.new + # if two disjoint time series, adding the together to get total length
                #TAXA1 +
                CLIMATE1 + #identical to column CLIMATE2 ie all within study comparisons are within the same climate
                treatment_yn + #is there some sort of disturbance yes/no (fertilizer, fire, grazing etc)
                interaction_present + #does GLOBI record these genera as potentially interacting
                #interaction_type #no effect here, so removing bc adds unneseccary complexity. also small sample sizes (eg 1) within some categories reduce power
                #Centrality_Betweenness_Edge +
                RESOLVED.TAXA.PAIR 
                #(1|RESOLVED.TAXA.PAIR) 
                ) 

Q2mod <- brm(MODFORM, slopes.meta4, FAM, seed = 042023, 
         control = list(adapt_delta=0.99, max_treedepth = 12),    
         chains = 4, iter = 5000, warmup = 1000, cores = 4) 

save(Q2mod, file = 'outputs/Aug2023/Q2.model.wTaxa.fixed.wTreatment.interactionYN.Rdata')      


# output models
# Q2.model.wTaxa: ~ SERIES.l.new + TAXA1 + CLIMATE1 + interaction_type + (1|RESOLVED.TAXA.PAIR) 
# Q2.model.wTaxa.random.wTreatment ~ SERIES.l.new + CLIMATE1 + treatment_yn + interaction_type + (1|RESOLVED.TAXA.PAIR) 
# Q2.model.wTaxa.fixed.wTreatment ~ SERIES.l.new  + CLIMATE1 +treatment_yn + interaction_type + RESOLVED.TAXA.PAIR
# Q2.model.wTaxa.fixed.wTreatment.interactionYN ~ SERIES.l.new  + CLIMATE1 + treatment_yn + interaction_present + RESOLVED.TAXA.PAIR


# notes on models
# all converge, diagnostics look good
# seems like resolved taxa matters for prediction, so using as fixed effects
# interaction type does not, so simplifying things just to include interaction yes/no
# series length, climate, treatment all have an effect on strength of genera associations,
# regardless of what else include in model
# intercept CI cross zero except in models with resolved taxa as fixed effect
# incorporating centrality will require subsetting data, so waiting to do this
# prediction vs inference? were we gonna compare results? 

load("/Users/Gavia/Documents/14 U of T/CIEE/predicting_species_abundance/outputs/Aug2023/Q2.model.wTaxa.fixed.wTreatment.interactionYN.Rdata")

#sample size (Bulk_ESS & Tail_ESS) should > 1000 & rhat < 1.1
summary(Q2mod)

plot(Q2mod) #model convergence (L: does distribution mean match estimate? R: did all values get explored?)


# Posterior predictive check: Does data match model? (could be wrong distribution, not all effects modeled)
pp_check(Q2mod, ndraws = 100) #posterior predictive checks - are predicted values similar to posterior distribution?

# Pairs plots to diagnose sampling problems (should show Gaussian blobs)
pairs(mod)


plot(conditional_effects(Q2mod)) #fitted parameters and their CI
# need to make sure id is reversible ie monocot.dicot is same as dicot.monocot


##### plot random slopes ####
mod <- Q2mod

# extract the draws corresponding to posterior distributions of the overall mean and standard deviation of observations
mod %>%
  spread_draws(b_Intercept, sigma) %>%
  head(10)

# median and 95% quantile interval of the variables, we can apply median_qi()
mod %>%
  #spread_draws(b_Intercept, sigma) %>% #wide
  gather_draws(b_Intercept, b_SERIES.l.new, b_CLIMATE1Tropical, 
               b_treatment_ynyes, b_interaction_present1,
               sigma) %>% #long
  median_qi()

# summary for all regions [only when have random effect]
mod %>%
  spread_draws(r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR, ]) %>%
  median_qi()

# convergence diagnostics
mod %>%
  #spread_draws(r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR, ]) %>%
  spread_draws(b_Intercept, b_SERIES.l.new, b_CLIMATE1Tropical, 
               b_treatment_ynyes, b_interaction_present1, sigma) %>% 
  summarise_draws()

# maybe we want to know what the overall effect is, ie the global mean as expressed in each region
mod %>%
  #spread_draws(b_Intercept, r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR,]) %>%
  #median_qi(condition_mean = b_Intercept + r_RESOLVED.TAXA.PAIR) # can change the width of the credible interval
  spread_draws(b_Intercept, b_CLIMATE1Tropical, 
               b_treatment_ynyes, b_interaction_present1) %>%
  median_qi(tropical_mean = b_Intercept + b_CLIMATE1Tropical,
            treatmentYES_mean = b_Intercept + b_CLIMATE1Tropical,
            interactionYES_mean = b_Intercept + b_CLIMATE1Tropical) # can change the width of the credible interval
  
# and plot conditional effects
mod %>%
  #spread_draws(b_Intercept, r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR,]) %>%
  spread_draws(b_Intercept, b_CLIMATE1Tropical) %>%
  median_qi(condition_mean = b_Intercept + b_CLIMATE1Tropical, .width = c(.95, .66)) %>%
  ggplot(aes(y = reorder(b_CLIMATE1, condition_mean), x = condition_mean, xmin = .lower, xmax = .upper)) +
  geom_pointinterval() +
  geom_vline(xintercept = 0, color = "#839496", size = 1) 

# and with density distribution
mod %>%
  spread_draws(b_Intercept, r_RESOLVED.TAXA.PAIR[RESOLVED.TAXA.PAIR,]) %>%
  mutate(condition_mean = b_Intercept + r_RESOLVED.TAXA.PAIR) %>%
  ggplot(aes(y = reorder(RESOLVED.TAXA.PAIR, condition_mean), x = condition_mean)) +
  geom_vline(xintercept = 0, color = "#839496", size = 1) +
  stat_halfeye(.width = .5, size = 2/3, fill = "#859900")+
  labs(x = expression("Estimate"),
       y = "Associations between taxa pairs") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        text = element_text(family = "Ubuntu")) 

##### plot conditional effects ####
# series
# taxa
# climate
# interaction type
mod.DF.series.l <- conditional_effects(mod, re_formula = NULL)[[1]]
mod.DF.taxa <- conditional_effects(mod, re_formula = NULL)[[2]]


ggplot() +
  # error
  # geom_ribbon(data=mod.DF.series.l, 
  #             aes(SERIES.l.new, Estimate.Prop.Change.Gn2, ymin = lower__, ymax = upper__), alpha=.2) + 
  # # fit lines
  # geom_line(data=mod.DF.series.l, 
  #           aes(effect1__, estimate__), 
  #           size=1.2) +
  # # raw data
  geom_point(data=slopes.meta4, 
             aes(SERIES.l.new, Estimate.Prop.Change.Gn2, colour=treatment_yn, group=treatment_yn), 
             size=4, alpha=.8) +
  facet_wrap(~interaction_type) +
  theme_classic(base_size = 25) +
  labs(y="Length of time series", 
       x="Estimate of association",
       col="", group="") +
  theme(plot.background = element_blank(),
        strip.text.x = element_text(size=25),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), #r, t, l, b
        legend.position="none") 

