**Project name:** Predicting Species Abundance

**Date created:** 27 April 2023

**Created by:** Sam Straus

**Date updated:**

**Updated by:**

This folder stores intermediate data files created during the data pre-processing steps.

-   `all.interactions.genus.pairs.RDS` contains pairs of genera included in our main data set and their interactions. Interactions were pulled from GLOBI in the `scripts/Data.prcoessing.R` script, and simplified to four interaction types: dispersal, mutualism, predator_prey, and uncategorized_interaction

-   `between.studies.overlap.RDS` contains overlap between study pairs with at least 10 overlapping years, only accounting for pairs where each species comes from different study ID, Created by Isaac Eckert

-   `class.assignment` contains taxonomic Class information for each genus, pulled using the `taxize` package

-   `collated.pairs_standardized_summary_GLL.Rdata` contains standardized measures of abundance and biomass measures at the species level. Species abundance and biomass measures are divided by the sampling effort (EFFORT.YEAR = number of plots sampled per study per year) to produce a standardized measure. Mean, median, min/max (plus 5% quantiles), and sd+CoV are then re-calculated from these new standardized measures, Created by Mia Waters, Sophia Fan, and Faraz Khan

-   `dummy.dataset.RDS` a set of random pairs to test modelling framework, Created by Isaac Eckert

-   `genus_interaction_list.csv` list of interactions between genus pairs, Created by Emily Black

-   `log.prop.change.full.data.RDS` created in `scripts/data_processing/Portion.Change.R`. This .RDS file contains all the cleaned pair information with resolved taxonomy, does not include interaction data. This is an intermediary file

-   `log.prop.change.interactions.RDS` is created in the `scripts/data_processing/Data.processing.R` script. It integrates interaction information with the `log.prop.change.with.meta.RDS`data file

-   log.prop.change.with.meta.BETWEEN.w.taxa.RDS created in `scripts/data_processing/Portion.Change.R`. This .RDS file contains cleaned pair information, with resolved taxonomy, for genus pairs found only between different studies

-   `log.prop.change.with.meta.RDS` contains the log propotional change for all bio pairs, created by Isaac Eckert

-   log.prop.change.with.meta.WITHIN.w.taxa.RDS created in `scripts/data_processing/Portion.Change.R`. This .RDS file contains cleaned pair information, with resolved taxonomy, for genus pairs found only within the same study

-   `outdated/` folder containing outdated data files that were created during early stages of data processing.

-   `within.studies.overlap.RDS` contains overlap within study ID with at least 10 overlapping years, only accounting for pairs where each species comes from same study ID, Created by Isaac Eckert
