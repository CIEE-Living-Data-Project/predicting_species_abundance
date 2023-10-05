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

slopes.meta4$interaction_benefit <- ifelse(slopes.meta4$interaction_benefit=="NA", "no interaction", slopes.meta4$interaction_benefit)


### fix NA resolved taxa pairs ####
# code chunk by EB

# taxa 1
slopes.meta4$RESOLVED.TAXA1 <- ifelse(
  is.na(slopes.meta4$RESOLVED.TAXA1),
  ifelse(
    slopes.meta4$ORGANISMS1 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
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
    slopes.meta4$ORGANISMS2 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(slopes.meta4$ORGANISMS2 == "birds", "Aves", 
           ifelse(slopes.meta4$ORGANISMS2 == "rodents", "Mammalia", NA)
    )
  ),
  slopes.meta4$RESOLVED.TAXA2
)

sorted_words <- apply(slopes.meta4[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
slopes.meta4$resolved_taxa_pair <- sorted_words
unique(slopes.meta4$resolved_taxa_pair)
table(slopes.meta4$resolved_taxa_pair)


#### fit model ####
FAM <- gaussian(link = 'identity')

# centrality measures
# interaction as yes/no
MODFORM <- bf(Estimate.Prop.Change.Gn2|resp_se(Est.Error.Prop.Change.Gn2, sigma = TRUE) ~ 
                SERIES.l.new + # if two disjoint time series, adding the together to get total length
                #TAXA1 +
                CLIMATE1 + #*interaction_present  + #identical to column CLIMATE2 ie all within study comparisons are within the same climate
                treatment_yn + #is there some sort of disturbance yes/no (fertilizer, fire, grazing etc)
                interaction_benefit +
                #does GLOBI record these genera as potentially interacting
                #interaction_type #no effect here, so removing bc adds unneseccary complexity. also small sample sizes (eg 1) within some categories reduce power
                #Centrality_Betweenness_Edge +
                RESOLVED.TAXA.PAIR 
                #(1|RESOLVED.TAXA.PAIR) 
                ) 

Q2mod <- brm(MODFORM, slopes.meta4, FAM, seed = 042023, 
         control = list(adapt_delta=0.99, max_treedepth = 12),    
         chains = 4, iter = 10000, warmup = 500, cores = 4) 

save(Q2mod, file = 'outputs/Aug2023/Q2.model.wTaxa.fixed.wTreatment.interactionBenefit.Rdata')      

# re-run with interaction between interaction_present*CLIMATE1 or type of interaction

# output models
# Q2.model.wTaxa: ~ SERIES.l.new + TAXA1 + CLIMATE1 + interaction_type + (1|RESOLVED.TAXA.PAIR) 
# Q2.model.wTaxa.random.wTreatment ~ SERIES.l.new + CLIMATE1 + treatment_yn + interaction_type + (1|RESOLVED.TAXA.PAIR) 
# Q2.model.wTaxa.fixed.wTreatment ~ SERIES.l.new  + CLIMATE1 +treatment_yn + interaction_type + RESOLVED.TAXA.PAIR
# Q2.model.wTaxa.fixed.wTreatment.interactionYN ~ SERIES.l.new  + CLIMATE1 + treatment_yn + interaction_present + RESOLVED.TAXA.PAIR
# Q2.model.wTaxa.fixed.wTreatment.interactionYNxClimate ~ SERIES.l.new  + CLIMATE1*interaction_present + treatment_yn  + RESOLVED.TAXA.PAIR


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

##### plot conditional effects with random effects incorporated ####
# series
# taxa
# climate
# interaction type
mod.DF.series.l <- conditional_effects(mod, re_formula = NULL)[[1]]
mod.DF.taxa <- conditional_effects(mod, re_formula = NULL)[[2]]


ggplot() +
  # error
  geom_ribbon(data=mod.DF.series.l,
              aes(SERIES.l.new, Estimate.Prop.Change.Gn2, ymin = lower__, ymax = upper__), alpha=.2) +
  # fit lines
  geom_line(data=mod.DF.series.l,
            aes(effect1__, estimate__),
            size=1.2) +
  # raw data
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


#### plot conditional effects when no random effects structure ####

##### coefficient plots ####
color_scheme_set("blue")
mcmc_intervals(as.matrix(mod),
           pars = c("b_Intercept", "b_SERIES.l.new","b_CLIMATE1Tropical", "b_treatment_ynyes", "b_interaction_present1"),
           prob = 0.8, 
           prob_outer = .9,
           point_est = "mean", 
           point_size = 3) +
  scale_y_discrete(
    labels = c("b_Intercept" = "Intercept",
               "b_SERIES.l.new" = "Time series length",
               "b_CLIMATE1Tropical" = "Climate:tropical",
               "b_treatment_ynyes" = "Treatment:yes",
               "b_interaction_present1" = "Interaction:yes"),
    limits = c("b_interaction_present1",
               "b_treatment_ynyes",
               "b_CLIMATE1Tropical", 
               "b_SERIES.l.new", 
               "b_Intercept"
               ))

# just time series length, hard to see on previous
mcmc_areas(as.matrix(mod),
           pars = c( "b_SERIES.l.new")) +
  scale_y_discrete(
    labels = c("b_SERIES.l.new" = "Time series length"))


# just taxa
# baseline is Aves.Aves
color_scheme_set("blue")
mod.taxa.coefs <- as.matrix(mod)[,6:40]
colnames(mod.taxa.coefs) <- gsub(".*PAIR","", colnames(mod.taxa.coefs)) # rename
colnames(mod.taxa.coefs) <- gsub("NA","Unclassified", colnames(mod.taxa.coefs)) 

mcmc_intervals(mod.taxa.coefs,
               prob = 0.8, 
               prob_outer = .9,
               point_est = "mean")

##### plots with fitted lines ####
# cool, holds everything at intercept value
# but not what we want in this case. want to be able to specify value combinations that
# actually exist in the data
series.lengthDF <- conditional_effects(mod, 
                                     re_formula = NULL)[[1]]

ggplot() +
  # error
  geom_ribbon(data=series.lengthDF, 
              aes(SERIES.l.new, Estimate.Prop.Change.Gn2, ymin = lower__, ymax = upper__,
                  group = CLIMATE1, fill = CLIMATE1), alpha=.2) + 
  # fit lines
  geom_line(data=series.lengthDF, 
            aes(effect1__, estimate__), 
            size=1.2) +
  # raw data
  #geom_point(data=site.alpha.intraDF, 
             #aes(elev, convex.vol, pch=habitat, colour=habitat, group=habitat), 
             #size=4, alpha=.8) +
  # manual colour scheme
  #scale_color_manual(values=c( "#877800", "#004A39")) + 
  #scale_fill_manual(values=c( "#877800", "#004A39"))  +
  #scale_shape_manual(values=c(17, 19))+
  theme_classic(base_size = 25) +
  labs(y="Association", 
       x="Series length"#,
       #col="", group=""
       ) +
  theme(plot.background = element_blank(),
        strip.text.x = element_text(size=0),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), #r, t, l, b
        legend.position="none") 


# want to be able to specify data
predict.Q2 <- predict(mod) # same as posterior_predict (?) equiv to random draws from posterior normal, not what we want
predict.Q2B <- posterior_epred(mod) # expectation of posterior
predict.Q2C <- posterior_linpred(mod) # expectation of linear predictor. identical to expectation of posterior bc gaussian model

library(emmeans)
rg <- ref_grid(mod)
# marginal effects calculated for us
em <- emmeans(rg, specs=c("CLIMATE1", "treatment_yn", "interaction_present", "RESOLVED.TAXA.PAIR"))
summary(em, point.est = mean)
plot(em) # gives coefficient plot, note that not all these conditions actually exist in the data


epred <- emmeans(mod, 
                 specs=c("CLIMATE1", "treatment_yn", "interaction_present", "RESOLVED.TAXA.PAIR"), 
                 epred = TRUE)
summary(epred, point.est = mean)


new_data <- expand.grid(elev = seq(min(subset(site.alpha.intraDF, island=="Hispaniola")$elev), 
                                   max(subset(site.alpha.intraDF, island=="Hispaniola")$elev), 
                                   by = 5),
                        habitat = c("Forest", "Anthro"))



ggplot(slopes.meta5) +
  # error
  geom_ribbon(aes(SERIES.l.new, Estimate.Prop.Change.Gn2, ymin = Q2.5, ymax = Q97.5,
                  group = RESOLVED.TAXA.PAIR, fill = RESOLVED.TAXA.PAIR), alpha=.2) + 
  # fit lines
  geom_line(aes(SERIES.l.new, Estimate), 
            size=1.2) +
  facet_grid(treatment_yn~CLIMATE1)
  # raw data
  #geom_point(data=site.alpha.intraDF, 
  #aes(elev, convex.vol, pch=habitat, colour=habitat, group=habitat), 
  #size=4, alpha=.8) +
  # manual colour scheme
  #scale_color_manual(values=c( "#877800", "#004A39")) + 
  #scale_fill_manual(values=c( "#877800", "#004A39"))  +
  #scale_shape_manual(values=c(17, 19))+
  theme_classic(base_size = 25) +
  labs(y="Association", 
       x="Series length"#,
       #col="", group=""
  ) +
  theme(plot.background = element_blank(),
        strip.text.x = element_text(size=0),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), #r, t, l, b
        legend.position="none")

  
# not quite doing what I want, but closer
# want each marginal slope plotted for each of the taxa.taxa pairs
# try by just pulling out specifc taxa.pairs--this fixes it!
library(rphylopic)
  
# insecta.insecta
img <- pick_phylopic(name = "libellula pulchella", n = 1)
A <- slopes.meta4 %>% 
    subset(RESOLVED.TAXA.PAIR=="Insecta.Insecta" & CLIMATE1=="Temperate" & treatment_yn=="yes") %>%   
  #data_grid(flipper_length_mm = seq_range(flipper_length_mm, n = 100)) |> 
    group_by(RESOLVED.TAXA.PAIR) %>%
    add_linpred_draws(mod, ndraws = 100) %>% 
    ggplot(aes(x = SERIES.l.new), colour=RESOLVED.TAXA.PAIR) +
    stat_lineribbon(aes(y = .linpred), .width = 0.95,
                    alpha = 0.5
                    ) +
    facet_wrap(treatment_yn~CLIMATE1, scales = "free") +
    geom_point(data = subset(slopes.meta4, RESOLVED.TAXA.PAIR=="Insecta.Insecta" & CLIMATE1=="Temperate" & treatment_yn=="yes"),
               aes(y = Estimate.Prop.Change.Gn2), size = 1, alpha = 0.7) +
  theme_classic(base_size = 25) +
  labs(y="Association", 
       x="Series length") +
  theme(plot.background = element_blank(),
        #strip.text.x = element_text(size=25),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), #r, t, l, b
        legend.position="none") +
  add_phylopic(x = 12, y = 1.3, img = img, alpha = .8, ysize = .2) +
  add_phylopic(x = 16, y = 1.3, img = img, alpha = .8, ysize = .2) 

# mamallia.mamallia
#img <- pick_phylopic(name = "Giraffa camelopardalis", n = 2)
uuid <- get_uuid(name = "Giraffa camelopardalis", n = 1)
img <- get_phylopic(uuid = uuid)
B <- slopes.meta4 %>% 
  subset(RESOLVED.TAXA.PAIR=="Mammalia.Mammalia" & CLIMATE1=="Temperate") %>%   
  #data_grid(flipper_length_mm = seq_range(flipper_length_mm, n = 100)) |> 
  group_by(RESOLVED.TAXA.PAIR) %>%
  add_linpred_draws(mod, ndraws = 200) %>% 
  ggplot(aes(x = SERIES.l.new), colour=RESOLVED.TAXA.PAIR) +
  stat_lineribbon(aes(y = .linpred), .width = 0.95,
                  alpha = 0.5
  ) +
  facet_wrap(treatment_yn~CLIMATE1) +
  geom_point(data = subset(slopes.meta4, RESOLVED.TAXA.PAIR=="Mammalia.Mammalia" & CLIMATE1=="Temperate"),
             aes(y = Estimate.Prop.Change.Gn2), size = 1, alpha = 0.7) +
  theme_classic(base_size = 25) +
  labs(y="Association", 
       x="Series length") +
  theme(plot.background = element_blank(),
        #strip.text.x = element_text(size=25),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), #r, t, l, b
        legend.position="none") +
  add_phylopic(x = 12, y = .8, img = img, alpha = .8, ysize = .1) +
  add_phylopic(x = 16, y = .8, img = img, alpha = .8, ysize = .1) 

# Pinopsida.Monocots
#img <- pick_phylopic(name = "Giraffa camelopardalis", n = 2)
uuid <- get_uuid(name = "Welwitschia mirabilis", n = 1)
img <- get_phylopic(uuid = uuid)
uuid2 <- get_uuid(name = "Veratrum", n = 1)
img2 <- get_phylopic(uuid = uuid2)
C <- slopes.meta4 %>% 
  subset(RESOLVED.TAXA.PAIR=="Pinopsida.Monocots" & CLIMATE1=="Temperate" & treatment_yn=="no") %>%   
  #data_grid(flipper_length_mm = seq_range(flipper_length_mm, n = 100)) |> 
  group_by(RESOLVED.TAXA.PAIR) %>%
  add_linpred_draws(mod, ndraws = 200) %>% 
  ggplot(aes(x = SERIES.l.new), colour=RESOLVED.TAXA.PAIR) +
  stat_lineribbon(aes(y = .linpred), .width = 0.95,
                  alpha = 0.5
  ) +
  facet_wrap(treatment_yn~CLIMATE1) +
  geom_point(data = subset(slopes.meta4, RESOLVED.TAXA.PAIR=="Pinopsida.Monocots" & CLIMATE1=="Temperate" & treatment_yn=="no"),
             aes(y = Estimate.Prop.Change.Gn2), size = 1, alpha = 0.7) +
  theme_classic(base_size = 25) +
  labs(y="Association", 
       x="Series length") +
  theme(plot.background = element_blank(),
        #strip.text.x = element_text(size=25),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), #r, t, l, b
        legend.position="none") +
  add_phylopic(x = 12, y = .8, img = img2, alpha = .8, ysize = .1) +
  add_phylopic(x = 12.5, y = .8, img = img, alpha = .8, ysize = .1) 


# Pinopsida.Monocots
#img <- pick_phylopic(name = "Giraffa camelopardalis", n = 2)
uuid <- get_uuid(name = "Quercus robur", n = 1)
img <- get_phylopic(uuid = uuid)
uuid2 <- get_uuid(name = "Gnetum gnemon", n = 1)
img2 <- get_phylopic(uuid = uuid2)
D <- slopes.meta4 %>% 
  subset(RESOLVED.TAXA.PAIR=="Magnoliopsida.Gnetopsida" & CLIMATE1=="Temperate" ) %>%   
  #data_grid(flipper_length_mm = seq_range(flipper_length_mm, n = 100)) |> 
  group_by(RESOLVED.TAXA.PAIR) %>%
  add_linpred_draws(mod, ndraws = 200) %>% 
  ggplot(aes(x = SERIES.l.new), colour=RESOLVED.TAXA.PAIR) +
  stat_lineribbon(aes(y = .linpred), .width = 0.95,
                  alpha = 0.5
  ) +
  facet_wrap(treatment_yn~CLIMATE1) +
  geom_point(data = subset(slopes.meta4, RESOLVED.TAXA.PAIR=="Magnoliopsida.Gnetopsida" & CLIMATE1=="Temperate" ),
             aes(y = Estimate.Prop.Change.Gn2), size = 1, alpha = 0.7) +
  theme_classic(base_size = 25) +
  labs(y="Association", 
       x="Series length") +
  theme(plot.background = element_blank(),
        #strip.text.x = element_text(size=25),
        plot.margin = unit(c(1, 1, 1, 1), "cm"), #r, t, l, b
        legend.position="none") +
  add_phylopic(x = 12, y = .8, img = img, alpha = .8, ysize = .1) +
  add_phylopic(x = 12.5, y = .8, img = img2, alpha = .8, ysize = .1) 

library(cowplot)
plot_grid(A,B,C,D, labels=c("A", "B", "C", "D"))

slopes.meta4 %>%
    add_predicted_draws(Q2mod, ndraws = 10) %>%  # adding the posterior distribution
    ggplot(aes(x = SERIES.l.new, y = Estimate.Prop.Change.Gn2)) +  
    stat_lineribbon(aes(y = .prediction), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    geom_point(data = slopes.meta4, aes(SERIES.l.new, y = Estimate.Prop.Change.Gn2),
               colour = "darkseagreen4", size = 3) +   # raw data
    scale_fill_brewer(palette = "Greys") +
    ylab("Association") + 
    xlab("Series length") +
    theme_bw() +
    theme(legend.title = element_blank())


slopes.meta4 %>%
    group_by(RESOLVED.TAXA.PAIR) %>%
    add_predicted_draws(Q2mod, ndraws = 100) %>%
  ggplot(aes(x = SERIES.l.new, y = Estimate.Prop.Change.Gn2, 
             color = ordered(RESOLVED.TAXA.PAIR), fill = ordered(RESOLVED.TAXA.PAIR))) +  
  # geom_point(data = slopes.meta4, aes(SERIES.l.new, y = Estimate.Prop.Change.Gn2),
  #            colour = "darkseagreen4", size = 3, alpha=.3) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95, .80),  # regression line and CI
                  alpha = 0.5, colour="black") +
  facet_wrap(CLIMATE1~treatment_yn) +   # raw data
  #scale_fill_brewer(palette = "Set2") +
  #scale_color_brewer(palette = "Dark2") +
  ylab("Association") + 
  xlab("Series length") +
  theme_classic() +
  theme(legend.title = element_blank(),
        )


