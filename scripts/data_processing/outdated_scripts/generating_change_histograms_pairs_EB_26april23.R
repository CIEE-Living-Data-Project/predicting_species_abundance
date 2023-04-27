# Date created: 26 Apr 2023
# Date updated: 26 Apr 2023 (NC)

# # LIBRARIES # #
library(tidyverse)
library("dplyr")
library("tibble")
library("readr") 
library("ggplot2")
library("magrittr")

rm(list=ls()) 


log.prop.change.with.meta <- 
  readRDS(" data/log.prop.change.with.meta.RDS")
log_change <- log.prop.change.with.meta

head(log.prop.change.with.meta)

#Investigate the distribution of the change in log % abundance
histogram <- ggplot(log_change, aes(x = Log.prop.change.abun.Gn1)) +
  geom_histogram(aes(y=..density.., fill = "Genus 1"), alpha = 0.5, color = "black", binwidth=0.2) +
  geom_histogram(aes(x = Log.prop.change.abun.Gn2, y=..density.., fill = "Genus 2"), alpha = 0.6, color = "black", 
                 binwidth=0.2) +
  xlab("Log % change in abundance") + ylab("Density") +
  labs(fill="Genus")+
  scale_fill_manual(values = c("Genus 1" = "blue", "Genus 2" = "red")) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
histogram

histogram_bio <- ggplot(log_change, aes(x = Log.prop.change.bio.Gn1)) +
  geom_histogram(aes(y=..density.., fill = "Genus 1"), alpha = 0.5, color = "black", binwidth=0.2) +
  geom_histogram(aes(x = Log.prop.change.bio.Gn2, y=..density.., fill = "Genus 2"), alpha = 0.6, color = "black", 
                 binwidth=0.2) +
  xlab("Log % change in biomass") + ylab("Density") +
  labs(fill="Genus")+
  scale_fill_manual(values = c("Genus 1" = "blue", "Genus 2" = "red")) +
  theme_classic()+
  theme(axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
histogram_bio


