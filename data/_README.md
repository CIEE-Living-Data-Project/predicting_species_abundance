**Project name:** Predicting Species Abundance

**Date created:** 23 April 2023

**Created by:** Sam Straus

**Date updated:** 27 April 2023

**Updated by:** Sam Straus

This directory to stores raw data and processed data used in the working group. The file `BioTIMECitations_24_06_2021.csv` is a table with:

`STUDY_ID`: a unique, integer ID number given to each study, and\
`CITATION_LINE`: the full paper citation for each study.

-   There are three subdirectories:\
    `prep_biotime/`: This folder stores raw data queried from BioTIME on 15 Mar 2023 by Nathalie Chardon.

-   `preprocessing/`: This folder stores intermediate data files created during the data pre-processing steps

-   `tidy/`: This folder stores the `collated_pairs.RData` file, which is created using the `scripts/prep_biotime/data_processing.R` script
