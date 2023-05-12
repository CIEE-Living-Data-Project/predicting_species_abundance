**Project name:** Predicting Species Abundance

**Date created:** 27 April 2023

**Created by:** Sam Straus

**Date updated:**

**Updated by:**

This subdirectory of the `scripts/` folder stores the scripts used in the data processing of the BioTIME database. There are two main scripts:

-   `Data.processing.R` is the main data processing file. It calls functions from the `Functions_data.processing.R` script. This script achieves the following data processing steps:

    -   Remove species with fewer than three observation across all studies and years

    -   Standardizing abundances by sampling effort

    -   Find overlapping consecutive years for all pairs

    -   Calculate proportion change for species pairs

-   `Functions_data.processing.R` script contains 5 key functions used in the `Data.processing.R` script:

    -   `set.prog.bar(n_iter)`

    -   `calc.overlap.between(p,data)`

    -   `calc.overlap.within(ID,data)`

    -   `get.log.prop.change(x,data,pairs.keep)`

    -   `make.meta(data,meta)`

-   `/outdated_scripts` is a subdirectory containing scripts that have either been phased out or incorporated into other workflows

-   `Portion.Change.R` Created by Isaac Eckert, this script calculates the proportion change and creates several intermediate data files stored in `data/data_processing/`
