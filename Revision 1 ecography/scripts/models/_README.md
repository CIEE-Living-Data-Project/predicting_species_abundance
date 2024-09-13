**Project name:** Predicting Species Abundance

**Date created:** 23 Feb 2023

**Created by:** Sam Straus

**Date updated:** 13 September 2024

**Updated by:** GLL

1. `setup_model data` - This file collates all the intermediate data files from the data pre-processing and converts it into the form to fit the model (adds metadata, calculates pearson correlations and z-scores etc)  

2. `lmer_models` - This file script that runs our primary model, calculates the cross-validation, and runs the null model to verify our model results. 