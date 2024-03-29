**Project name:** Predicting Species Abundance

**Date created:** 23 April 2023

**Created by:** Sam Straus

**Date updated:** 12 May 2023

**Updated by:** Sam Straus

**Description of:** `all.interactions.genus.pairs.RDS`: this RDS file contains information about interactions between genus pairs

| Column Name      | Description                                                                            |
|----------------|--------------------------------------------------------|
| Gn1              | Genus of first species in pair                                                         |
| Gn2              | Genus of second species in pair                                                        |
| interaction_type | Type of interaction: dispersal, mutualism, predator_prey, or uncategorized_interaction |
| num_interactions | The number of times this interaction occurs in our dataset                             |

**Description of: `between.studies.overlap.RDS`:** This RDS file is a spreadsheet containing information about genus pairs and the amount of years (observations?) they overlap between studies. Genus names are aggregated from species-level observations

| Column Name | Description                                                                     |
|-----------------|-------------------------------------------------------|
| Gn1         | Genus 1 of pair                                                                 |
| Gn2         | Genus 2 of pair                                                                 |
| Max.Overlap | Maximum number of overlapping observations per pair                             |
| Type        | Between versus within study pairs, this table contains only between study pairs |
| PAIR.ID     | Unique pair ID, (Gn1 study ID - Gn2 study ID)                                   |

**Description of:** `class.assignment`: RDS object with taxonomic class information pulled using `taxize` package

| Column Name | Description                                       |
|-------------|---------------------------------------------------|
| db          | taxize database that was queried, all from "itis" |
| query       | search query, which is the genus name             |
| class       | Taxonomic class of genus name                     |

**Description of**: `collated.pairs_standardized_summary_GLL.Rdata`: R object name: collated.pairs_standardized_summary. This dataframe is the final collated pairs data used for downstream analysis. It is created using the `scripts/data_processing/Data.processing.R` script. Abundance and biomass measures are standardized at the species level. Species abundance and biomass measures are divided by the sampling effort (EFFORT.YEAR = number of plots sampled per study per year) to produce a standardized measure. Mean, median, min/max, and sd+CoV are then re-calculated from these new standardized measures, Created by Mia Waters, Sophia Fan, and Faraz Khan

| Column Name        | Description                                                                                                                                                                                                                                                               |
|----------------|--------------------------------------------------------|
| ID                 | Study ID                                                                                                                                                                                                                                                                  |
| YEAR               | Year of study                                                                                                                                                                                                                                                             |
| GENUS              | Genus name aggregated from species-level observations                                                                                                                                                                                                                     |
| TAXA               | Organism type: All, Marine plants, Fish, Marine invertebrates, Terrestrial plants, Benthos, Freshwater invertebrates, Freshwater plants, Terrestrial invertebrates, Amphibians, Mammals, Birds, Fungi)                                                                    |
| ORGANISMS          | Further delineation of organism type (e.g. salamanders, tropical algae)                                                                                                                                                                                                   |
| CLIMATE            | Climate of the study site: Temperate/Tropical, Temperate, Polar, Tropical                                                                                                                                                                                                 |
| REALM              | Marine, Terrestrial, Freshwater                                                                                                                                                                                                                                           |
| mean_abun_st       | Standardized mean abundance (if abundance measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                                  |
| median_abun_st     | Standardized median abundance (if abundance measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                                |
| min_abun_st        | Standardized minimum abundance (if abundance measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                               |
| avg_abun_min_5perc | bottom 5 percent quantile of average minimum abundance                                                                                                                                                                                                                    |
| max_abun_st        | Standardized maximum abundance (if abundance measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                               |
| avg_abun_max_5perc | top 5 percent quantile of average maximum abundance                                                                                                                                                                                                                       |
| sd_abun_st         | Standardized standard deviation for abundance (if abundance measure) per genus (for those genera with multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                 |
| CoV_abun_st        | Standardized coefficient of variation (standard deviation/mean) for abundance (if abundance measure) per genus (for those genera with multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year) |
| mean_bio_st        | Standardized mean biomass (if biomass measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                                      |
| median_bio_st      | Standardized median biomass (if biomass measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                                    |
| min_bio_st         | Standardized minimum biomass (if biomass measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                                   |
| avg_bio_min_5perc  | bottom 5 percent quantile of average minimum biomass                                                                                                                                                                                                                      |
| max_bio_st         | Standardized maximum biomass (if biomass measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                                   |
| avg_bio_max_5perc  | top 5 percent quantile of average maximum biomass                                                                                                                                                                                                                         |
| sd_bio_st          | Standardized standard deviation for biomass (if biomass measure) per genus (across multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)                                                    |
| CoV_bio_st         | Standardized coefficient of variation (standard deviation/mean) for biomass (if biomass measure) per genus (for those genera with multiple species), grouped by study ID, year, taxa, organisms, climate, and realm (since these are all the same per study per year)     |

**Description of**: `dummy.dataset.RDS` is a fake dataset using random pairs from `log.prop.change.with.meta.RDS`. This dummy dataset was used to create and test our analysis framework before using our real data.

| Column Name              | Description                                       |
|-------------------------|----------------------------------------------|
| Gn1                      | Genus of first species in pair                    |
| Gn2                      | Genus of second species in pair                   |
| Log.prop.change.abun.Gn1 | Log proportion change in abundance of Gn1         |
| Log.prop.change.abun.Gn2 | Log proportion change in abundance of Gn2         |
| Log.prop.change.bio.Gn1  | Log proportion change in biomass of Gn1           |
| Log.prop.change.bio.Gn2  | Log proportion change in biomass of Gn2           |
| PairID                   | Unique pair ID, (StudyID of Gn1 - StudyID of Gn2) |
| Type                     | Between studies or Within studies                 |
| SERIES.n                 | number of time series in that pair                |
| ID1                      | Study ID of Gn1                                   |
| ID2                      | Study ID of Gn2                                   |
| REALM1                   | Realm of Gn1 (Marine, Terrestrial, Freshwater)    |
| REALM2                   | Realm of Gn2 (Marine, Terrestrial, Freshwater)    |
| TAXA1                    | Taxa of Gn1                                       |
| TAXA2                    | Taxa of Gn2                                       |
| HABITAT1                 | Habitat of Gn1                                    |
| HABITAT2                 | Habitat of Gn2                                    |
| PROTECTED_AREA1          | Is Gn1 observation from a protected area?         |
| PROTECTED_AREA2          | Is Gn2 observation from a protected area?         |
| UNIQUE.PAIR.ID           | A unique pair ID in the format of Gn1-Gn2-ID1-ID2 |
| SERIES.length            | length of series in years                         |

**Description of**: `genus_interaction_list.csv` is a file containing the interaction types between different genera included in our final genus list.

| Column Name | Description                                                                                                                                                                                                                                  |
|----------------|--------------------------------------------------------|
| ...1        | Record number                                                                                                                                                                                                                                |
| Gn1         | Genus 1 of pair                                                                                                                                                                                                                              |
| Gn2         | Genus 2 of pair                                                                                                                                                                                                                              |
| interaction | Type of interaction: adjecentTo, dispersalVectorOf, eatenBy, eats, ecologicallyRelatedTo, flowersVisitedBy, hasDispersalVector, hasVector, interactsWith, mutualistOf, parasiteOf, preyedUponBy, preysOn, visitedBy, visits, visitsFlowersOf |

**Description of:** `log.prop.change.full.data.RDS` is an intermediate RDS file with the log proportional change of abundance, including study metadata and resolved taxonomy

| Column Name     | Description                                       |
|-----------------|---------------------------------------------------|
| Gn1             | Genus of first species in pair                    |
| Gn2             | Genus of second species in pair                   |
| Prop.Change.Gn1 | proportion change of genus 1 (See Metric column)  |
| Prop.Change.Gn2 | proportion change of genus 2 (See Metric column)  |
| PairID          | Unique pair ID, (StudyID of Gn1 - StudyID of Gn2) |
| Type            | Between studies or Within studies                 |
| SERIES.n        | number of time series in that pair                |
| SERIES.start    | Start year of time series for pair                |
| SERIES.end      | End year of time series for pair                  |
| YEAR.T          | year at time T                                    |
| YEAR.T1         | year at time T+1                                  |
| SERIES.l        | Series length in years                            |
| ID1             | Study ID of Gn1                                   |
| ID2             | Study ID of Gn2                                   |
| REALM1          | Realm of Gn1 (Marine, Terrestrial, Freshwater)    |
| REALM2          | Realm of Gn2 (Marine, Terrestrial, Freshwater)    |
| TAXA1           | Taxa of Gn1                                       |
| TAXA2           | Taxa of Gn2                                       |
| ORGANISMS1      | Type of organism for genus 1                      |
| ORGANISMS2      | Type of organism for genus 2                      |
| CLIMATE1        | Climate for ID1                                   |
| CLIMATE2        | Climate for ID2                                   |
| UNIQUE.PAIR.ID  | A unique pair ID in the format of Gn1-Gn2-ID1-ID2 |
| RESOLVED.TAXA1  | Resolved taxonomy for Gn1                         |
| RESOLVED.TAXA2  | Resolved taxonomy for Gn2                         |
| Metric          | Type of metric (abundance, biomass)               |
| dist            | distance between ID1 and ID2                      |

**Description of**: `log.prop.change.interactions.RDS` is an RDS file that integrates interaction information with `log.prop.change.with.meta.RDS`

| Column Name               | Description                                                             |
|-------------------------|----------------------------------------------|
| Gn1                       | Genus of first species in pair                                          |
| Gn2                       | Genus of second species in pair                                         |
| Log.prop.change.abun.Gn1  | log proportional change in abundance of Gn1                             |
| Log.prop.change.abun.Gn2  | log proportional change in abundance of Gn2                             |
| Log.prop.change.bio.Gn1   | log proportional change in biomass of Gn1                               |
| Log.prop.change.bio.Gn2   | log proportional change in biomass of Gn2                               |
| PairID                    | Unique pair ID, (StudyID of Gn1 - StudyID of Gn2)                       |
| Type                      | Between studies or Within studies                                       |
| SERIES.n                  | number of time series in that pair                                      |
| SERIES.start              | Start year of time series for pair                                      |
| SERIES.end                | End year of time series for pair                                        |
| YEAR.T                    | year at time T                                                          |
| YEAR.T1                   | year at time T+1                                                        |
| SERIES.l                  | Series length in years                                                  |
| ID1                       | Study ID of Gn1                                                         |
| ID2                       | Study ID of Gn2                                                         |
| REALM1                    | Realm of Gn1 (Marine, Terrestrial, Freshwater)                          |
| REALM2                    | Realm of Gn2 (Marine, Terrestrial, Freshwater)                          |
| TAXA1                     | Taxa of Gn1                                                             |
| TAXA2                     | Taxa of Gn2                                                             |
| ORGANISMS1                | Type of organism for genus 1                                            |
| ORGANISMS2                | Type of organism for genus 2                                            |
| CLIMATE1                  | Climate for ID1                                                         |
| CLIMATE2                  | Climate for ID2                                                         |
| UNIQUE.PAIR.ID            | A unique pair ID in the format of Gn1-Gn2-ID1-ID2                       |
| interaction_found         | was an interaction found? binary variable. 0 = no, 1 = yes              |
| positive_interaction      | was a positive interaction found? binary variable. 0 = no, 1 = yes      |
| negative_interaction      | was a negative interaction found? binary variable. 0 = no, 1 = yes      |
| neutral_interaction       | was a neutral interaction found? binary variable. 0 = no, 1 = yes       |
| uncategorized_interaction | was a uncategorized interaction found? binary variable. 0 = no, 1 = yes |
| predator_prey_interaction | was a predator-prey interaction found? binary variable. 0 = no, 1 = yes |
| mutualism_interaction     | was a mutualistic interaction found? binary variable. 0 = no, 1 = yes   |
| dispersal_interaction     | was a dispersal interaction found? binary variable. 0 = no, 1 = yes     |

**Description of:** `log.prop.change.with.meta.RDS` is an RDS file with the log proportional change of abundance, including study metadata

| Column Name              | Description                                       |
|-------------------------|----------------------------------------------|
| Gn1                      | Genus of first species in pair                    |
| Gn2                      | Genus of second species in pair                   |
| Log.prop.change.abun.Gn1 | Log proportion change in abundance of Gn1         |
| Log.prop.change.abun.Gn2 | Log proportion change in abundance of Gn2         |
| Log.prop.change.bio.Gn1  | Log proportion change in biomass of Gn1           |
| Log.prop.change.bio.Gn2  | Log proportion change in biomass of Gn2           |
| PairID                   | Unique pair ID, (StudyID of Gn1 - StudyID of Gn2) |
| Type                     | Between studies or Within studies                 |
| SERIES.n                 | number of time series in that pair                |
| SERIES.start             | Start year of time series for pair                |
| SERIES.end               | End year of time series for pair                  |
| YEAR.T                   | year at time T                                    |
| YEAR.T1                  | year at time T+1                                  |
| SERIES.l                 | Series length in years                            |
| ID1                      | Study ID of Gn1                                   |
| ID2                      | Study ID of Gn2                                   |
| REALM1                   | Realm of Gn1 (Marine, Terrestrial, Freshwater)    |
| REALM2                   | Realm of Gn2 (Marine, Terrestrial, Freshwater)    |
| TAXA1                    | Taxa of Gn1                                       |
| TAXA2                    | Taxa of Gn2                                       |
| HABITAT1                 | Habitat of Gn1                                    |
| HABITAT2                 | Habitat of Gn2                                    |
| PROTECTED_AREA1          | Is Gn1 observation from a protected area?         |
| PROTECTED_AREA2          | Is Gn2 observation from a protected area?         |
| UNIQUE.PAIR.ID           | A unique pair ID in the format of Gn1-Gn2-ID1-ID2 |

**Description of:** `log.prop.change.with.meta.WITHIN.w.taxa.RDS` and `log.prop.change.with.meta.WITHIN.w.taxa.RDS` are RDS file with the log proportional change of abundance, including study metadata and resolved taxonomy, either within the same study or between different studies, respectively. Both files have the same column names and descriptions.

| Column Name              | Description                                       |
|-------------------------|----------------------------------------------|
| Gn1                      | Genus of first species in pair                    |
| Gn2                      | Genus of second species in pair                   |
| Log.prop.change.abun.Gn1 | Log proportion change in abundance of Gn1         |
| Log.prop.change.abun.Gn2 | Log proportion change in abundance of Gn2         |
| Log.prop.change.bio.Gn1  | Log proportion change in biomass of Gn1           |
| Log.prop.change.bio.Gn2  | Log proportion change in biomass of Gn2           |
| PairID                   | Unique pair ID, (StudyID of Gn1 - StudyID of Gn2) |
| Type                     | Between studies or Within studies                 |
| SERIES.n                 | number of time series in that pair                |
| SERIES.start             | Start year of time series for pair                |
| SERIES.end               | End year of time series for pair                  |
| YEAR.T                   | year at time T                                    |
| YEAR.T1                  | year at time T+1                                  |
| SERIES.l                 | Series length in years                            |
| ID1                      | Study ID of Gn1                                   |
| ID2                      | Study ID of Gn2                                   |
| REALM1                   | Realm of Gn1 (Marine, Terrestrial, Freshwater)    |
| REALM2                   | Realm of Gn2 (Marine, Terrestrial, Freshwater)    |
| TAXA1                    | Taxa of Gn1                                       |
| TAXA2                    | Taxa of Gn2                                       |
| ORGANISMS1               | Type of organism for genus 1                      |
| ORGANISMS2               | Type of organism for genus 2                      |
| CLIMATE1                 | Climate for ID1                                   |
| CLIMATE2                 | Climate for ID2                                   |
| UNIQUE.PAIR.ID           | A unique pair ID in the format of Gn1-Gn2-ID1-ID2 |
| RESOLVED.TAXA1           | Resolved taxonomy for Gn1                         |
| RESOLVED.TAXA2           | Resolved taxonomy for Gn2                         |

**Description of:** `within.studies.overlap.RDS` this RDS file is contains information about genus pairs and the amount of years (observations?) they overlap within studies. Genus names are aggregated from species-level observations

| Column Name | Description                                                                     |
|------------------|------------------------------------------------------|
| Gn1         | Genus 1 of pair                                                                 |
| Gn2         | Genus 2 of pair                                                                 |
| Max.Overlap | Maximum number of overlapping observations per pair                             |
| Type        | Between versus within study pairs, this table contains only between study pairs |
| PAIR.ID     | Unique pair ID, (Gn1 study ID - Gn2 study ID)                                   |
