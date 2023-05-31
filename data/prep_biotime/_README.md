**Project name:** Predicting Species Abundance

**Date created:** 23 April 2023

**Created by:** Sam Straus

**Date updated:** 12 May 2023

**Updated by:**

This directory stores the raw data downloaded from BioTIME by Nathalie Chardon on 15 Mar 2023.\
The data are stored as both `.csv` and `.RData` files.

-   The `bio_pairs_1km.*` and `bio_pairs_10km.*` datasets are the records of overlapping species time series at 1km and 10km resolutions, respectively.

-   `collated_pairs.RData` file, which is created using the `scripts/prep_biotime/data_processing.R` script, is the main data file used in the analysis. It contains the abundance and biomass data for each species found at each site/study ID.

-   The `meta_pairs_10km.*` data have the metadata for each study ID found in the `bio_pairs_10km.*` data file.\
    \
    See `_DATA-DICTIONARY.md` for column descriptions.
