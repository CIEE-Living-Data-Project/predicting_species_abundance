**Project name:** Predicting Species Abundance

**Date created:** 23 Feb 2023

**Created by:** Sam Straus

**Date updated:** 13 Sept 2024

**Updated by:** GLL


This directory contains the intermediate data the arises from the data cleaning scripts and then gets fed into the modelling script. 

1. `model_data_final.Rdata` - final and cleaned data used to fit the model

2. `results.abundance.csv` - output from `results.abundance.data.setup.R`, called in `setup_model data.R`. Contains log proportional pairwise changes in abundance across genera

3. `within.study.updated.data.050224.csv` - created in `pair_interactions_taxize_all.R` and called in same script

4. `results_abundance_interactions_taxa_032024ENB.RDS` -  created in `pair_interactions_taxize_all.R`, called in `setup_model data.R`. Contains genus pair log proportional changes and whether or not each genus per genera pair has interaction data from Globi.

5. `worldclim.csv` - output from `worldclim summaries.R` script, called in `setup_model data.R`. Contains mean annual temperature and precipitation values for each longitude and latitude coordinate in our cleaned BioTime data.

6. `disturbance cleaning.csv` - manually cleaned BioTime disturbance data, called in `setup_model data.R`
