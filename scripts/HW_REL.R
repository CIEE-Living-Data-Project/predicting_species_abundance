## Primary Working Directory
pwd <- "/Users/ryan/Windows/Documents/Post UCB/Research/CIEE_Workshop_TrophicInteractionsPredictAbundances/predicting_species_abundance"
# pwd <- "/home/ryan/predicting_species_abundance"

## Packages
library("dplyr")
library("tibble")
library("readr")
library("ggplot2")
library("magrittr")

## Data
data_raw <- readr::read_csv(paste0(pwd, "/data/BioTIMEQuery_24_06_2021.csv"))
bio_pairs_10km <- readr::read_csv(paste0(pwd, "/data/prep_biotime/bio_pairs_10km.csv"))

## Unique ID (1 second resolution) so every run of this script will output without overwriting previous runs
tu_id <- gsub(x = Sys.time(), pattern = "-|:| ", replacement = "")

## Explore a particular pair of time series with lots of overlapping data
pair <- 164

habitat_pair <- bio_pairs_10km$habitat.1[pair]

pair_1_ID <- bio_pairs_10km$ID.1[pair]
pair_2_ID <- bio_pairs_10km$ID.2[pair]

study_1 <- data_raw %>% dplyr::filter(STUDY_ID == pair_1_ID)
study_2 <- data_raw %>% dplyr::filter(STUDY_ID == pair_2_ID)

years_overlap <- unique(study_1$YEAR)[unique(study_1$YEAR) %in% unique(study_2$YEAR)] %>% sort()

## Species in the two studies
study_1_species <- study_1$SPECIES %>% unique()
study_1_species = study_1_species[!is.na(study_1_species)]
study_1_species_length <- sapply(study_1_species, function(x){
    study_1 %>%
        dplyr::filter(SPECIES == x) %>%
        dplyr::select(YEAR) %>%
        unique() %>%
        unlist() %>%
        length()
})

study_2_species <- study_2$SPECIES %>% unique()
study_2_species = study_2_species[!is.na(study_2_species)]
study_2_species_length <- sapply(study_2_species, function(x){
    study_2 %>%
        dplyr::filter(SPECIES == x) %>%
        dplyr::select(YEAR) %>%
        unique() %>%
        unlist() %>%
        length()
})

## Randomly use a species from each study, with at least 10 observations
sp1 = study_1_species[sample(x = which(study_1_species_length > 10), size = 1)]
sp2 = study_2_species[sample(x = which(study_2_species_length > 10), size = 1)]

## Average any spatial replicates within each year
sp1_data <- study_1 %>%
    dplyr::filter(
        SPECIES == sp1,
        YEAR %in% years_overlap
    ) %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarize(
        Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
            mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
            mean(sum.allrawdata.BIOMASS, na.rm = TRUE)
        )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(YEAR)

sp2_data <- study_2 %>%
    dplyr::filter(
        SPECIES == sp2,
        YEAR %in% years_overlap
    ) %>%
    dplyr::group_by(YEAR) %>%
    dplyr::summarize(
        Abundance = sum( ## This assumes there will only ever be either abundance or biomass data
            mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE),
            mean(sum.allrawdata.BIOMASS, na.rm = TRUE)
        )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(YEAR)

## The overlap is for all species in each study, but each pair of species may not overlap across all the years in years_overlap
sp1_data %<>% dplyr::filter(YEAR %in% sp2_data$YEAR)
sp2_data %<>% dplyr::filter(YEAR %in% sp1_data$YEAR)

## Calculate correlation between two species
cor_pair <- cor(sp1_data$Abundance, sp2_data$Abundance)

## Plot
fig_data <- dplyr::bind_rows(
    sp1_data %>% dplyr::mutate(Species = sp1),
    sp2_data %>% dplyr::mutate(Species = sp2)
)

## Rename YEAR to Year
colnames(fig_data)[1] = c("Year")

fig <- fig_data %>%
    ggplot() +
        theme_bw() +
        geom_line(
            aes(
                x = Year,
                y = Abundance,
                color = Species
            ), 
            linewidth = 3
        ) +
        labs(
            title = paste0("Correlation = ", round(cor_pair, 2))
        )

setwd(paste0(pwd, "/figures/explore_biotime_Dec2022"))

grDevices::cairo_pdf(filename = paste0("TimeSeries_Correlation_Pair-", pair, "_", tu_id, ".pdf"), width = 10, height = 5)
    print(fig)
dev.off()
