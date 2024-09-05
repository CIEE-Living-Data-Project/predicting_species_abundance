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
load(file="Revision 1 ecography/output/meta model/lmer_model_final.Rdata")
model <- lmermod


# Calculate numbers for results ####
### raw correlations ####
# calculate
# use cuttoffs of .167 bc greater than this is a positive correlation as defined by the model
# back transform mu = .167 (mean estimated z-score) from model into pearson correlation
# z = 0.5*log((1+cor)/(1-cor))) #eq 3.11 in met
z=.167
r = (exp(2*z)-1)/(exp(2*z)+1) # pearson correlation is 0.1654646
r_upr <- (exp(2*(z+2*0.0523336911))-1)/(exp(2*(z+2*0.0523336911))+1)
r_lwr <- (exp(2*(z-2*0.0523336911))-1)/(exp(2*(z-2*0.0523336911))+1)

# what % of data falls within estimated 95% CI?
nrow(subset(moddat, cor >= r_lwr & cor <= r_upr)) # 88,543 within 95% CI
88543/nrow(moddat)

# Or how many correlations are weaker than model estimated mean, how many stronger (either negative or positive)
nrow(subset(moddat, cor >= r | cor <= -r)) # 200,861 as strong or stronger
sum(moddat$cor >= r)/nrow(moddat) # strong positive
sum(moddat$cor <= -r)/nrow(moddat) # strong negative

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


##  Figure 3. coef plot  #######
# histogram of pearson cor
hist <- moddat %>% 
  mutate(group=ifelse(cor>=r_lwr & cor<=r_upr, "in \n in", "out \n out")) %>% 
  # mutate(group=ifelse(cor>=r, "Strong \npositive", 
  #                      ifelse(cor<=-r, " Strong \nnegative", "Neutral"))) %>% 
  ggplot(aes(x=cor, fill=group)) + 
  geom_histogram(colour="black", bins=29) +
  theme_bw()+
  xlab("Correlation")+
  ylab("Number of genus pairs")+
  geom_vline(aes(xintercept=r),
             color="black", linetype="dashed", size=1)+
  theme(panel.grid   = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y  = element_text(hjust = 0)) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  scale_fill_manual(values = c("#66C2A5", "#8DA0CB","#FC8D62")) +
  theme(legend.position="top",
        legend.spacing.x = unit(1, 'cm'))+
  guides(fill = guide_legend(label.position = "top", title="Group"))

# time series and interactions error so small can't see, make insets
# better to use confidence intervals?
total <- summary(model)$coefficients %>% 
  data.frame() %>% 
  mutate(term = rownames(.)) %>% 
  ggplot(aes(Estimate, term, xmin=Estimate-2*`Std..Error`, xmax=Estimate+2*`Std..Error`, colour=term)) +
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
      "scale.SERIES.l" = "Time series \nlength",
      "scale.abs.lat" = "Latitude",
      "interaction_present.factor1" = "GLoBI \ninteractions",
      "treatment_yn_cleanyes" = "Disturbance"),
    limits = c("scale.SERIES.l", 
               "treatment_yn_cleanyes",
               "interaction_present.factor1",
               "scale.abs.lat",
               "(Intercept)" 
    )) +
  guides(col = "none") +
  scale_color_manual(values=c("grey30", "#009E73","#0072B2" ,"#CC79A7","#E69F00" ))


time.inset <- summary(model)$coefficients %>% 
  data.frame() %>% 
  mutate(term = rownames(.)) %>% 
  filter(term=="scale.SERIES.l") %>% 
  ggplot(aes(Estimate, term, xmin=Estimate-2*`Std..Error`, xmax=Estimate+2*`Std..Error`, colour=term)) +
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
  scale_x_continuous(breaks=c(.0140, .016)) +
  guides(col = "none") +
  scale_color_manual(values=c("#CC79A7" ))

# GLoBI inset
int.inset <- summary(model)$coefficients %>% 
  data.frame() %>% 
  mutate(term = rownames(.)) %>% 
  filter(term=="interaction_present.factor1") %>% 
  ggplot(aes(Estimate, term, xmin=Estimate-2*`Std..Error`, xmax=Estimate+2*`Std..Error`, colour=term)) +
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
  scale_x_continuous(breaks=c(.006, .011)) 


total2 <- total + 
  annotation_custom(
    ggplotGrob(time.inset), 
    #xmin = .10, xmax = .2, ymin = .6, ymax = 1.6) +
    xmin = .14, xmax = .24, ymin = .6, ymax = 1.6) +
  annotation_custom(
    ggplotGrob(int.inset), 
    #xmin = .05, xmax = .15, ymin = 2.6, ymax = 3.6)
    xmin = .14, xmax = .24, ymin = 2.6, ymax = 3.6) 

fig3 <- plot_grid(hist, total2, rel_widths = c(1,2))

ggsave("Revision 1 ecography/output/figures/figure3_CI.pdf", 
       fig3, width=20, height=8, units="in")

##  Figure 4. ranef plots  #######
str(dd <- as.data.frame(ranef(model)))
ranef_study <- dd %>% 
  filter(grpvar=="STUDY_ID") %>% 
  ggplot(aes(condval+0.1671776, grp,
             xmin=condval+0.1671776 -1.96*condsd, 
             xmax=condval+0.1671776 + 1.96*condsd)) +
  geom_pointrange(linewidth=2, size=1.2) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = 0.1671776, linetype = "dashed", color = "#8DA0CB", size=2)+
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
  ggplot(aes(condval+0.1671776, grp,
             xmin=condval+0.1671776 -1.96*condsd, 
             xmax=condval+0.1671776 +1.96*condsd)) +
  geom_pointrange(linewidth=2, size=1.2) +
  #geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey")+
  geom_vline(xintercept = 0.1671776, linetype = "dashed", color = "#8DA0CB", size=2)+
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

fig4 <- plot_grid(ranef_study, ranef_taxa)
ggsave("Revision 1 ecography/output/figures/figure4.pdf", fig4, width=20, height=20, units="in")

## Figure 5 posterior pred plots ####
library(performance)

# w/o random effects
fig5a_dat <- performance::check_predictions(lmermod, re_formula = NA) %>% 
  pivot_longer(cols=c(sim_1:y)) %>% 
  mutate(Grouping = ifelse(name=="y", "Observed \ndata", "Model-predicted \ndata"))

fig5a <- ggplot(fig5a_dat, 
               aes(value, col=as.factor(Grouping), group=as.factor(name))) + 
  geom_density(size=1.5) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  labs(x = "Value",
       y = "Z-score") +
  scale_colour_manual(name = "", 
                      values=c("#8DA0CB","#FC8D62")) +
  theme(legend.position = "none")


# w random effects
fig5b_dat <- performance::check_predictions(model, re_formula = NULL) %>% 
  pivot_longer(cols=c(sim_1:y)) %>% 
  mutate(Grouping = ifelse(name=="y", "Observed \ndata", "Model-predicted \ndata"))

fig5b<- ggplot(fig5b_dat, 
       aes(value, col=as.factor(Grouping), group=as.factor(name))) + 
  geom_density(size=1.5) +
  theme_bw(base_size = 25) +
  theme(panel.grid.major.x = element_blank(),  # Hide major x-axis grid lines
        panel.grid.minor.x = element_blank()) +   # Hide minor x-axis grid lines
  my.theme +
  labs(x = "Value",
       y = "") +
scale_colour_manual(name = "", 
                    values=c("#8DA0CB","#FC8D62"))

fig5 <- plot_grid(fig5a, fig5b)

ggsave("Revision 1 ecography/output/figures/figure5.pdf", fig5, width=20, height=10, units="in")

# Tables ####
## Table S3 ####
summary(model)$coefficients %>% 
  knitr::kable(format="simple", digits=c(6,6,6,6,6,6,6,3,0,0))

## Table S4 ####
# study ID
# taxonomic category
as.data.frame(ranef(model)) %>% 
  mutate(on.mean = condval+0.1671776) %>% 
  knitr::kable(format="simple", digits=c(6,6,6,6,6,6,6,3,0,0))



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

