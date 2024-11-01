**Project name:** Predicting Species Abundance

**Date created:** 23 April 2023

**Created by:** Sam Straus

**Date updated:** 12 May 2023

**Updated by:** Sam Straus

**Description of `bio_pairs*`:**\
Both the `bio_pairs_1km.*` and `bio_pairs_10km.*` share the same data dictionary, as follows

| Column Name   | Description                                                              |
|----------------|--------------------------------------------------------|
| dist          | Distance (km) between study sites in ID.1 and ID.2                       |
| ID.1          | Unique study ID number of first study site in pair                       |
| ID.2          | Unique study ID number of second study site in pair                      |
| long.1        | Longitude of first study site in pair                                    |
| long.2        | Longitude of second study site in pair                                   |
| lat.1         | Latitude of first study in pair                                          |
| lat.2         | Latitude of second study in pair                                         |
| habitat.1     | Habitat type of first study in pair                                      |
| habitat.2     | Habitat type of second study in pair                                     |
| taxa.1        | Broadly, the taxa sampled from first study in pair                       |
| taxa.2        | Broadly, the taxa sampled from second study in pair                      |
| organisms.1   | More specifically, the organisms sampled from first study in pair        |
| organisms.2   | More specifically, the organisms sampled from second study in pair       |
| no_species.1  | Number of species sampled from first study in pair                       |
| no_species.2  | Number of species sampled from second study in pair                      |
| start_year.1  | Start year of monitoring of first study in pair                          |
| start_year.2  | Start year of monitoring of second study in pair                         |
| end_year.1    | End year of monitoring of first study in pair                            |
| end_year.2    | End year of monitoring of second study in pair                           |
| overlap.years | Overlap in years of monitoring between studies in pair                   |
| taxa.pairs    | Types of organisms paired in this row, e.g. birds (ID.1) and fish (ID.2) |

**Description of `BioTIMEMetadata_24_06_2021.csv`:** This data file contains detailed study metadata for each unique Study ID. This is the raw data file downloaded from BioTIME by Nathalie Chardon on 15 Mar 2023.

| Column Name       | Description                                                                |
|-----------------|-------------------------------------------------------|
| STUDY_ID          | Unique study ID number of study                                            |
| REALM             | Realm (Terrestrial, Marine, Freshwater) in which the study occurred        |
| CLIMATE           | Climate of study site (Temperate, Tropical, Polar, Global)                 |
| GENERAL_TREAT     | Broadly, treatments of study site, if any                                  |
| TREATMENT         | Spefically, treatments of study site, if any                               |
| TREAT_COMMENTS    | Comments about treatment                                                   |
| TREAT_DATE        | Date of treatment                                                          |
| HABITAT           | Description of study site habitat                                          |
| PROTECTED_AREA    | Study site is in a protected area - TRUE/FALSE                             |
| BIOME_MAP         | Description of study site biome                                            |
| TAXA              | Broadly, taxa monitored in study                                           |
| ORGANISMS         | Specifically, type of organism monitored in study                          |
| TITLE             | Title of study                                                             |
| AB_BIO            | ??                                                                         |
| HAS_PLOT          | ??                                                                         |
| DATA_POINTS       | Number of data points in study                                             |
| START_YEAR        | Start of monitoring period                                                 |
| END_YEAR          | End of monitoring period                                                   |
| CENT_LAT          | Center of latitude of study site                                           |
| CENT_LONG         | Cent of longitude of study site                                            |
| NUMBER_OF_SPECIES | Number of species monitored in study                                       |
| NUMBER_OF_SAMPLES |                                                                            |
| NUMBER_LAT_LONG   |                                                                            |
| TOTAL             |                                                                            |
| GRAIN_SIZE_TEXT   | Text description of spatial extent of study                                |
| GRAIN_SQ_KM       | resolution (aka sampling/plot size)                                        |
| AREA_SQ_KM        | Area covered in study                                                      |
| CONTACT_1         | Name of primary contact for study                                          |
| CONTACT_2         | Name of secondary contact for study                                        |
| CONT_1\_MAIL      | Email of primary contact for study                                         |
| CONT_2\_MAIL      | Email of secondary contact for study                                       |
| LICENSE           | Type of data license                                                       |
| WEB_LINK          | link to data source                                                        |
| DATA_SOURCE       | Name of data source                                                        |
| METHODS           | Description of sampling methods used                                       |
| SUMMARY_METHODS   | Very broadly, type of methods used, e.g. Transects, Quadrats               |
| LINK_ID           |                                                                            |
| COMMENTS          | General notes about the study                                              |
| DATE_STUDY_ADDED  | Date study was added to data sheet                                         |
| ABUNDANCE_TYPE    | Type of abundance measure used: Density, Count, or Presence/Absence        |
| BIOMASS_TYPE      | Type of biomass measure used: Weight or Cover; NOT USING FOR WORKING GROUP |
| SAMPLE_DESC_NAME  |                                                                            |

**Description of `collated.pairs.RData`:**

| Column Name              | Description                                                             |
|--------------------|----------------------------------------------------|
| X                        | Unique record ID                                                        |
| ID                       | Study ID                                                                |
| DAY                      | Day of month                                                            |
| MONTH                    | Month of record, 1-12                                                   |
| YEAR                     | Year of record                                                          |
| SAMPLE_DESC              | unique description code, concatenated latitude, longitude and year      |
| PLOT                     | Plot ID, if applicable                                                  |
| ID_SPECIES               | Unique species ID number                                                |
| LATITUDE                 | Latitude of study site                                                  |
| LONGITUDE                | Longitude of study site                                                 |
| sum.allrawdata.ABUNDANCE | Sum of raw abundance values                                             |
| sum.allrawdata.BIOMASS   | Sun of raw biomass values                                               |
| GENUS                    | Genus of species ID                                                     |
| SPECIES                  | Species epithet of species ID                                           |
| GENUS_SPECIES            | concatenation of Genus and Species columns                              |
| ABUNDANCE_TYPE           | Type of abundance measure (Count, Density, MeanCount, Presence/Absence) |
| BIOMASS_TYPE             | Type of biomass measure (Cover, Weight)                                 |

**Description of `meta_pairs_10km.*`:** This file contains all of the information from `BioTIMEMetadata_24_06_2021.csv` and also includes a few columns added and population by Nathalie Chardon, Courtney Collins, Ryan Langendorf, Haley Branch, and Sam Straus.

| Column Name         | Description                                                                                                                                                                                                                                 |
|-------------|-----------------------------------------------------------|
| STUDY_ID            | Unique study ID number of study                                                                                                                                                                                                             |
| reviewed_by         | \| The name of the person who verified the record                                                                                                                                                                                           |
| REALM               | Realm (Terrestrial, Marine, Freshwater) in which the study occurred                                                                                                                                                                         |
| CLIMATE             | Climate of study site (Temperate, Tropical, Polar, Global)                                                                                                                                                                                  |
| GENERAL_TREAT       | Broadly, treatments of study site, if any                                                                                                                                                                                                   |
| TREATMENT           | Spefically, treatments of study site, if any                                                                                                                                                                                                |
| TREAT_COMMENTS      | Comments about treatment                                                                                                                                                                                                                    |
| TREAT_DATE          | Date of treatment                                                                                                                                                                                                                           |
| HABITAT             | Description of study site habitat                                                                                                                                                                                                           |
| PROTECTED_AREA      | Study site is in a protected area - TRUE/FALSE                                                                                                                                                                                              |
| BIOME_MAP           | Description of study site biome                                                                                                                                                                                                             |
| TAXA                | Broadly, taxa monitored in study                                                                                                                                                                                                            |
| ORGANISMS           | Specifically, type of organism monitored in study                                                                                                                                                                                           |
| TITLE               | Title of study                                                                                                                                                                                                                              |
| AB_BIO              | ??                                                                                                                                                                                                                                          |
| HAS_PLOT            | ??                                                                                                                                                                                                                                          |
| DATA_POINTS         | Number of data points in study                                                                                                                                                                                                              |
| START_YEAR          | Start of monitoring period                                                                                                                                                                                                                  |
| END_YEAR            | End of monitoring period                                                                                                                                                                                                                    |
| CENT_LAT            | Center of latitude of study site                                                                                                                                                                                                            |
| CENT_LONG           | Cent of longitude of study site                                                                                                                                                                                                             |
| NUMBER_OF_SPECIES   | Number of species monitored in study                                                                                                                                                                                                        |
| NUMBER_OF_SAMPLES   |                                                                                                                                                                                                                                             |
| NUMBER_LAT_LONG     |                                                                                                                                                                                                                                             |
| TOTAL               |                                                                                                                                                                                                                                             |
| GRAIN_SIZE_TEXT     | Text description of spatial extent of study                                                                                                                                                                                                 |
| GRAIN_SQ_KM         | resolution (aka sampling/plot size)                                                                                                                                                                                                         |
| AREA_SQ_KM          | Area covered in study                                                                                                                                                                                                                       |
| CONTACT_1           | Name of primary contact for study                                                                                                                                                                                                           |
| CONTACT_2           | Name of secondary contact for study                                                                                                                                                                                                         |
| CONT_1\_MAIL        | Email of primary contact for study                                                                                                                                                                                                          |
| CONT_2\_MAIL        | Email of secondary contact for study                                                                                                                                                                                                        |
| LICENSE             | Type of data license                                                                                                                                                                                                                        |
| WEB_LINK            | link to data source                                                                                                                                                                                                                         |
| DATA_SOURCE         | Name of data source                                                                                                                                                                                                                         |
| METHODS             | Description of sampling methods used                                                                                                                                                                                                        |
| SUMMARY_METHODS     | Very broadly, type of methods used, e.g. Transects, Quadrats                                                                                                                                                                                |
| eco_footprint \| gr | ossly estimated spatial extent of the organism (km)                                                                                                                                                                                         |
| detectability \| hi | gh: captures all individuals with sampling strategy with high precision; medium: captures most individuals with sampling strategy with medium precision; low: does not capture all individuals and likely is at low precision (error prone) |
| LINK_ID             |                                                                                                                                                                                                                                             |
| COMMENTS            | General notes about the study                                                                                                                                                                                                               |
| DATE_STUDY_ADDED    | Date study was added to data sheet                                                                                                                                                                                                          |
| ABUNDANCE_TYPE      | Type of abundance measure used: Density, Count, or Presence/Absence                                                                                                                                                                         |
| BIOMASS_TYPE        | Type of biomass measure used: Weight or Cover; NOT USING FOR WORKING GROUP                                                                                                                                                                  |
| SAMPLE_DESC_NAME    |                                                                                                                                                                                                                                             |
