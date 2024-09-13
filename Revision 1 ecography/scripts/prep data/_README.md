**Project name:** Predicting Species Abundance

**Date created:** 23 Feb 2023

**Created by:** Sam Straus

**Date updated:** 18 April 2024

**Updated by:** ENB

This directory stores the scripts used in the data processing and analysis for the working group, generated and run after initial reviewer comments were received. These scripts filter and clean the data and calculate log proportional changes and richness for genera-pairs (results.abundance.data.setup), add the GloBi interactions and taxonomic class with manual corrections being performed (pair_interactions_taxize_all), and make adjustments to studies 221 and 39 which did not use standard binomial species names (throughout, but note pair_interactions_taxize_221_39). Studies 221 and 39 were ultimately removed from analyses, but correction steps are still included. Finally, there is a script which pulls the worldclim summaries for each plot in each study. Note that we did not end up using this data in analyses. 