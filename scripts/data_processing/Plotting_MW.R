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
  scale_fill_manual(values = c("#00aeff", "#0057cc", "#66CC66"))#freshwater, marine, terrestrial
#problem for later = what is the "All" group in TAXA? 

ggsave("figures/Genera_overview_plot.png", plot = genera_plot, height = 4, width = 6)

### Positive/Negative interactions within taxa-taxa comparison ####


