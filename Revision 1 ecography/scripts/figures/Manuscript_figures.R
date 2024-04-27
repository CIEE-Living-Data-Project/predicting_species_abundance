# packages
rm(list=ls()) 

# # LIBRARIES # #
library(tidyverse)
library("ggplot2")
library("magrittr")
library(cowplot)


# Read in data and model ####
load(file="Revision 1 ecography/output/prep_data/model_data_final.Rdata")
moddat <- subset(moddat, STUDY_ID!=39 & STUDY_ID!=221 & STUDY_ID!=225)
load(file="Revision 1 ecography/output/lmer_model_final.Rdata")
model <- lmermod


# Main text figures ####
## set custom theme ####
my.theme <- theme(axis.text=element_text(size=20),
                axis.title = element_text(size = 25),
                legend.text=element_text(size=20),
                legend.title = element_text(size=25),
                plot.title = element_text(face="bold",size=14,margin=margin(0,0,20,0),hjust = 0.5),
                axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0)),
                axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)))


## Figure 2 inset ####
#### Plot z-scores ####
# plot as histogram of raw data
moddat %>% 
  ggplot(aes(x=z)) + 
  geom_histogram(colour="black", bins=29, fill="#8DA0CB") +
  theme_bw()+
  xlab("Fisher's z-scores")+
  ylab("Number of genus pairs")+
  geom_vline(aes(xintercept=0),
             color="black", linetype="dashed", size=1)+
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme 


### Barplot of genera by taxonomic category ####
library(forcats)
library(rcartocolor)

options(scipen = 999) #converts to nice numbers on axes
# colors per taxa
taxaCol <- rcartocolor::carto_pal(11, "Safe")

a <- data.frame(table(moddat$RESOLVED.TAXA2))
b <- data.frame(table(moddat$RESOLVED.TAXA1))
c <- full_join(a,b, by="Var1")
c[10,2] <- 0
c$sum <- c$Freq.x+c$Freq.y

ggplot(c, 
       aes(x=fct_reorder(Var1, sum, .desc=TRUE), 
           y=log(sum), fill=Var1)) +
  #geom_bar(stat = "identity", fill="#66CC66") +
  #geom_col( fill="#3382BF") +
  geom_col() +
  labs(x = "Taxanomic category", y = "log Number of genera") +
  theme_classic(base_size = 25) + 
  theme(legend.position = "none",
        axis.text=element_text(size=25)) +
  theme(aspect.ratio = 0.4) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values = taxaCol) 

### World map where col of point corresponds to length of time series ####
drawWorld<-function(lats) {
  world_map<-map_data("world")
  
  g1<-ggplot()+coord_fixed()+xlab("")+ylab("")
  g1<-g1+geom_polygon(data=world_map, aes(x=long, y=lat, group=group), colour="gray60", fill="gray60")
  g1<-g1+theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
               panel.background=element_rect(fill="white", colour="white"), axis.line=element_line(colour="white"),
               legend.position="none",axis.ticks=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())
  
  if(lats=="y") {
    g1<-g1+geom_hline(yintercept=23.5, colour="darkred")+geom_hline(yintercept =-23.5, colour="darkred")
    g1<-g1+geom_hline(yintercept=66.5, colour="darkblue")+geom_hline(yintercept =-66.5, colour="darkblue")
  }
  else { return(g1) }
  return(g1)
}

# Now let's use the above function to plot these studies across the globe
(gplot <- drawWorld("y") + 
    geom_point(data=distinct(select(moddat, LONGITUDE, LATITUDE, SERIES.l, STUDY_ID)), 
               aes(x=LONGITUDE, y=LATITUDE, 
                   #colour = RESOLVED.TAXA1, 
                   colour = SERIES.l, 
                   size = 1
               ), 
               alpha=0.4)) 


## Figure 3. model accurary ####
### raw correlations ####
# histogram of pearson cor
hist <- moddat %>% 
  mutate(group=ifelse(cor>=.3, "positive", 
                      ifelse(cor<=-.3, "negative", "neutral"))) %>% 
  ggplot(aes(x=cor, fill=group)) + 
  geom_histogram(colour="black", bins=29) +
  #geom_histogram(aes(y=..count..., fill=group), colour="black") + #, fill="#8DA0CB"
  #geom_histogram(aes(y=..density..), colour="black", fill="white")+
  #geom_density(alpha=.2, fill="#8DA0CB") +
  theme_bw()+
  xlab("Correlation")+
  ylab("Number of genus pairs")+
  geom_vline(aes(xintercept=0),
             color="black", linetype="dashed", size=1)+
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0),
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  scale_fill_manual(values = c("#66C2A5", "#8DA0CB","#FC8D62"))

ggsave("Revision 1 ecography/output/figures/raw data histogram.pdf", hist, 
       width=8, height=5, units="in")


##  Figure 4. coef plot  #######

# main figure
# time series and interactions error so small can't see, make insets
# better to use confidence intervals?
total <- summary(model)$coefficients %>% 
  data.frame() %>% 
  mutate(term = rownames(.)) %>% 
  ggplot(aes(Estimate, term, xmin=Estimate-`Std..Error`, xmax=Estimate+`Std..Error`, colour=term)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", size=2) +
  geom_pointrange(linewidth=2, size=1.8) +
  labs(x = expression("Estimate"),
       y = "Model parameters") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  scale_y_discrete(
    labels = c("(Intercept)" = "Intercept",
      "scale.SERIES.l" = "Time series length",
      "scale.abs.lat" = "Latitude",
      "interaction_present.factor1" = "GLoBI interactions",
      "treatment_yn_cleanyes" = "Disturbance"),
    limits = c("scale.SERIES.l", 
               "treatment_yn_cleanyes",
               "interaction_present.factor1",
               "scale.abs.lat",
               "(Intercept)" 
    )) +
  guides(col = "none") +
  scale_color_manual(values=c("grey30", "#009E73","#0072B2" ,"#CC79A7","#E69F00" ))


# lat inset

time.inset <- summary(model)$coefficients %>% 
  data.frame() %>% 
  mutate(term = rownames(.)) %>% 
  filter(term=="scale.SERIES.l") %>% 
  ggplot(aes(Estimate, term, xmin=Estimate-`Std..Error`, xmax=Estimate+`Std..Error`, colour=term)) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", size=2) +
  geom_pointrange(linewidth=2, size=1.2) +
  labs(x = expression(""),
       y = "") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  scale_y_discrete(labels=c("scale.SERIES.l" = ""),
    limits = c("scale.SERIES.l")) +
  scale_x_continuous(breaks=c(.0142, 0.0148, .0153)) +
  guides(col = "none") +
  scale_color_manual(values=c("#CC79A7" ))

# GLoBI inset
int.inset <- summary(model)$coefficients %>% 
  data.frame() %>% 
  mutate(term = rownames(.)) %>% 
  filter(term=="interaction_present.factor1") %>% 
  ggplot(aes(Estimate, term, xmin=Estimate-`Std..Error`, xmax=Estimate+`Std..Error`, colour=term)) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", size=2) +
  geom_pointrange(linewidth=2, size=1.2) +
  labs(x = expression(""),
       y = "") +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  scale_y_discrete(labels=c("interaction_present.factor1" = ""),
                   limits = c("interaction_present.factor1")) +
  guides(col = "none") +
  scale_color_manual(values=c("#009E73" )) +
  scale_x_continuous(breaks=c(.0070, .0083, .0095)) 


fig4 <- total + 
  annotation_custom(
    ggplotGrob(time.inset), 
    #xmin = .10, xmax = .2, ymin = .6, ymax = 1.6) +
    xmin = .12, xmax = .22, ymin = .6, ymax = 1.6) +
  annotation_custom(
    ggplotGrob(int.inset), 
    #xmin = .05, xmax = .15, ymin = 2.6, ymax = 3.6)
    xmin = .12, xmax = .22, ymin = 2.6, ymax = 3.6) 

ggsave("Revision 1 ecography/output/figures/figure4.pdf", 
       fig4, width=15, height=8, units="in")

##  Figure 5. ranef plots  #######

str(dd <- as.data.frame(ranef(model)))
ranef_study <- dd %>% 
  filter(grpvar=="STUDY_ID") %>% 
  ggplot(aes(condval, grp,
             xmin=condval -1.96*condsd, 
             xmax=condval + 1.96*condsd)) +
  geom_pointrange(linewidth=2, size=1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = 0.1666366, linetype = "dashed", color = "#8DA0CB", size=2)+
  geom_pointinterval() +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  labs(x = "Condition mean",
       y = "Study ID")

ranef_taxa <- dd %>% 
  filter(grpvar=="resolved_taxa_pair") %>% 
  ggplot(aes(condval, grp,
             xmin=condval -1.96*condsd, 
             xmax=condval +1.96*condsd)) +
  geom_pointrange(linewidth=2, size=1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = 0.1666366, linetype = "dashed", color = "#8DA0CB", size=2)+
  geom_pointinterval() +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)
        #text = element_text(family = "Ubuntu")
  ) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  labs(x = "Condition mean",
       y = "Taxonomic categories") +
  scale_y_discrete(
    limits = c("Insecta.Insecta", 
               "Reptilia.Reptilia",
               "Gastropoda.Gastropoda",
               "Mammalia.Mammalia",
               "Aves.Aves",
               "Amphibia.Reptilia",
               "Monocots.Monocots",     
               "Monocots.Magnoliopsida",
               "Magnoliopsida.Monocots",    
               "Magnoliopsida.Magnoliopsida",
               "Magnoliopsida.Pinopsida"    ,
               "Pinopsida.Pinopsida",
               "Pinopsida.Magnoliopsida",
               "Gnetopsida.Magnoliopsida",
               "Magnoliopsida.Gnetopsida",
               "Gnetopsida.Monocots")) 

fig5 <- plot_grid(ranef_study, ranef_taxa)
ggsave("Revision 1 ecography/output/figures/figure5.pdf", fig5, width=20, height=20, units="in")

##  Figure 6. linear predictor plots  #######
# 4x4 panel plot highlighting random effect from 
# Monocots.Monocots
# Aves.Aves
# Insecta.Insecta
# Pinopsida.Pinopsida
# for an average site

# predict new data across time series length
# disturbance=0,1, interaction=0, latitude=0

new.data <- expand.grid(scale.SERIES.l = seq(min(moddat$scale.SERIES.l), 
                                                     max(moddat$scale.SERIES.l), by = .1),
                        scale.abs.lat = 0,
                        interaction_present.factor=c("0","1"),
                        treatment_yn_clean= c("no", "yes"), 
                        resolved_taxa_pair=c("Aves.Aves", "Monocots.Monocots", "Insecta.Insecta", "Pinopsida.Pinopsida"), 
                        STUDY_ID=c("54") # chosen bc doesn't exert strong effect, close to zero
                        ) 
new.predict <- cbind(predict(model, newdata=new.data), new.data)

# unscale time.series length
unscale <-  function (x, x.avg, s) { # x is transformed variable, x.avg is mean before scale and s is sd before scale
  x*s + x.avg
}

# re-order facets
int_names <- c(
  "0" = "No GLoBI interactions",
  "1" = "GLoBI interactions", 
  "Aves.Aves"="Aves.Aves", 
  "Monocots.Monocots"="Monocots.Monocots", 
  "Insecta.Insecta"="Insecta.Insecta", 
  "Pinopsida.Pinopsida"="Pinopsida.Pinopsida")


fig6 <- new.predict %>% 
  mutate(SERIES.l = unscale(scale.SERIES.l, mean(moddat$SERIES.l), sd(moddat$SERIES.l))) %>% 
  ggplot() +
  geom_line(aes(SERIES.l, `predict(model, newdata = new.data)`, colour=treatment_yn_clean), 
            size=1.2) + 
  # error
  # geom_ribbon(aes(SERIES.l, `predict(model, newdata = new.data)`, ymin = linepred-linepred.error, ymax = linepred+linepred.error,
  #                 group = treatment_yn_clean, fill = treatment_yn_clean), alpha=.2,
  #             show.legend = F) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#5E4141", size=1.1)+
  facet_grid(interaction_present.factor~factor(resolved_taxa_pair, 
                                               levels=c("Aves.Aves","Insecta.Insecta", "Monocots.Monocots", "Pinopsida.Pinopsida")), 
             labeller = as_labeller(int_names)) +
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white")) +   # Hide minor x-axis grid lines
  my.theme +
  labs(x = "Time series length (years)",
       y = "Estimate", 
       color = "Disturbance") +
  scale_color_manual(values=c( "#004A39", "#E69F00")) + 
  scale_fill_manual(values=c( "#004A39", "#E69F00"))   
fig6

ggsave("Revision 1 ecography/output/figures/figure6.pdf", fig6, width=20, height=10, units="in")

# Tables ####
## Table S3 ####
model %>% 
  spread_draws(b_Intercept, b_scale.SERIES.l, b_scale.abs.lat , b_treatment_yn_cleanyes, b_interaction_present.factor1,
               sd_STUDY_ID__Intercept, sd_resolved_taxa_pair__Intercept) %>% 
  #median_qi() %>% 
  summarise_draws() %>% 
  knitr::kable(format="simple", digits=c(6,6,6,6,6,6,6,3,0,0))

## Table S4 ####
# study ID
model %>%
  spread_draws(r_STUDY_ID[STUDY_ID, ]) %>%
  summarise_draws() #%>% 
  #knitr::kable(format="simple", digits=c(6,6,6,6,6,6,6,3,0,0))

# taxonomic category
model %>%
  spread_draws(r_resolved_taxa_pair[resolved_taxa_pair, ]) %>%
  summarise_draws() %>% 
  knitr::kable(format="simple", digits=c(6,6,6,6,6,6,6,3,0,0))

# can also look at taxonomic category per study ID
model %>%
  spread_draws(r_STUDY_ID[STUDY_ID, ], r_resolved_taxa_pair[resolved_taxa_pair, ]) %>%
  summarise_draws() %>% 
  knitr::kable(format="simple", digits=c(6,6,6,6,6,6,6,3,0,0))

#write.csv(tableS4, "Revision 1 ecography/output/figures/tableS4.csv")
  


# Sup mat figures ####

## Barplot genera pair by latitude ######
figS1 <- ggplot(data=moddat, aes(x=abs.lat)) +
  geom_histogram(binwidth = 1.2, bins = 40, fill="grey50", colour="white") +
  labs(x = expression("Absolute latitude"),
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

ggsave("figure s1.pdf", path="Revision 1 ecography/output/figures", width=7.5, height = 8, units="cm")

