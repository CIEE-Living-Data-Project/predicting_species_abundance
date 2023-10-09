# Aims:
# 1. Make a summary plot for the MS fig 1 (Ally requested)
# 2. Make more plots that would be helpful for visualizing the analyses at this stage

# Author: Mia Waters
# Date created: 13 July 2023
# Date updated 

dat <- readRDS('data/preprocessing/log.prop.change.interactions.RDS')
head(dat)


options(scipen = 999) #converts to nice numbers on axes

### Barplot of genera by taxa and realm for MS introduction ####

#need to add line for the path to collated.pairs_standardized_summary

genera_plot <- ggplot(data = collated.pairs_standardized_summary, 
                      aes(x=reorder(TAXA, TAXA,
                      function(x)-length(x)), fill = REALM)) + 
  geom_bar() +
  labs(x = "Taxa", y = "Number of genera") +
  theme_classic() + 
  theme(aspect.ratio = 0.4) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  scale_fill_manual(values = c("#00aeff", "#0057cc", "#66CC66")) #freshwater, marine, terrestrial
#problem for later = what is the "All" group in TAXA? 

ggsave("figures/Genera_overview_plot.png", plot = genera_plot, height = 4, width = 6)

### Positive/Negative interactions within taxa-taxa comparison ####

#### Barplot of genera by taxa in only terrestrial ####
# added by GLL for figure 2 (methods conceptual figure)
# only terrestrial
# only within

dat_terr <- subset(x = dat, subset = REALM1=="Terrestrial" & REALM2=="Terrestrial")
dat_terr <- subset(dat_terr, Metric!="CROSS" &  Type!="Between") %>%
  select(-c("Prop.Change.Gn1", "Prop.Change.Gn2", "YEAR.T","YEAR.T1", "SERIES.start", "SERIES.end")) %>% 
  distinct(.)

dat.terr2 <- dat_terr

dat.terr2$RESOLVED.TAXA.PAIR <- paste0(dat.terr2$RESOLVED.TAXA1, ".",dat.terr2$RESOLVED.TAXA2)

# fix resolved taxa names
#Check the NA values and their corresponding organism values 
dat_na2 <- dat.terr2 %>%
  filter(grepl("\\.NA|NA\\.", RESOLVED.TAXA.PAIR)) 

#double check that the NAs match
unique(dat_na2$ORGANISMS1)
unique(dat_na2$ORGANISMS2)
#Adjust names below 

# taxa 1
dat.terr2$RESOLVED.TAXA1 <- ifelse(
  is.na(dat.terr2$RESOLVED.TAXA1),
  ifelse(
    dat.terr2$ORGANISMS1 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat.terr2$ORGANISMS1 == "birds", "Aves", 
           ifelse(dat.terr2$ORGANISMS1 == "rodents", "Mammalia", NA)
    )
  ),
  dat.terr2$RESOLVED.TAXA1
)


# taxa 2
dat.terr2$RESOLVED.TAXA2 <- ifelse(
  is.na(dat.terr2$RESOLVED.TAXA2)==TRUE,
  ifelse(
    dat.terr2$ORGANISMS2 %in% c("insects", "Grasshoppers", "Acrididae (grasshoppers)"),
    "Insecta",
    ifelse(dat.terr2$ORGANISMS2 == "birds", "Aves", 
           ifelse(dat.terr2$ORGANISMS2 == "rodents", "Mammalia", NA)
    )
  ),
  dat.terr2$RESOLVED.TAXA2
)

sorted_words2 <- apply(dat.terr2[, c('RESOLVED.TAXA1', 'RESOLVED.TAXA2')], 1, function(x) paste(x, collapse = "."))
dat.terr2$RESOLVED.TAXA.PAIR <- sorted_words2
unique(dat.terr2$RESOLVED.TAXA.PAIR)
table(dat.terr2$RESOLVED.TAXA.PAIR)

data.frame(table(dat.terr2$RESOLVED.TAXA.PAIR))


library(forcats)
options(scipen = 999) #converts to nice numbers on axes
ggplot(data.frame(table(dat.terr2$RESOLVED.TAXA1)), 
       aes(x=fct_reorder(Var1, Freq, .desc=TRUE), 
           y=log(Freq))) +
  #geom_bar(stat = "identity", fill="#66CC66") +
  geom_col( fill="#3382BF") +
  labs(x = "Taxanomic category", y = "Number of genera") +
  theme_classic() + 
  theme(aspect.ratio = 0.4) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) 
