**Project name:** Predicting Species Abundance

**Date created:** 27 April 2023

**Created by:** Sam Straus

**Date updated:**

**Updated by:**

This folder stores intermediate data files created during the data pre-processing steps.

-   `between.studies.overlap.RDS` contains overlap between study pairs set with at least 10 overlapping years, Created by Isaac Eckert

-   `cleaned_collated_standardized_MSF.Rdata` contains version of cleaned_collated_pairs_AD.GLLAF_26April.csv where abundance and biomass measures are standardized at the species level. Species abundance and biomass measures are divided by the sampling effort (EFFORT.YEAR = number of plots sampled per study per year) to produce a standardized measure. Mean, median, min/max, and sd+CoV are then re-calculated from these new standardized measures, Created by Mia Waters, Sophia Fan, and Faraz Khan

-   `dummy.dataset.RDS` a set of random pairs to test modelling framework, Created by Isaac Eckert

-   `genus_interaction_list.csv` list of interactions between genus pairs, Created by Emily Black

-   `log.prop.change.with.meta.RDS` contains the log propotional change for all bio pairs, created by Isaac Eckert
