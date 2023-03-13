pwd <- "/home/ryan/predicting_species_abundance"

library("colorout")
library("tidyverse")
library("magrittr")
library("foreach")
library("doParallel")

## Data
data_raw <- readr::read_csv(paste0(pwd, "/data/BioTIMEQuery_24_06_2021.csv"))
bio_pairs_10km <- readr::read_csv(paste0(pwd, "/data/prep_biotime/bio_pairs_10km.csv"))

## Unique ID (1 second resolution) so every run of this script will output without overwriting previous runs
tu_id <- gsub(x = Sys.time(), pattern = "-|:| ", replacement = "")

## Initialize parallel backend
cluster <- makeCluster(50, outfile = "")
registerDoParallel(cluster)

cor_analysis <- foreach (pair = 1:nrow(bio_pairs_10km), .combine = "rbind", .packages = c("tibble", "dplyr", "magrittr")) %dopar% {
    pair_counter <- pair/nrow(bio_pairs_10km)

    habitat_pair <- bio_pairs_10km$habitat.1[pair]

    pair_1_ID <- bio_pairs_10km$ID.1[pair]
    pair_2_ID <- bio_pairs_10km$ID.2[pair]

    study_1 <- data_raw %>% dplyr::filter(STUDY_ID == pair_1_ID)
    study_2 <- data_raw %>% dplyr::filter(STUDY_ID == pair_2_ID)

    years_overlap <- unique(study_1$YEAR)[unique(study_1$YEAR) %in% unique(study_2$YEAR)] %>% sort()

    output_pair <- tibble(
                Study_1 = numeric(),
                Study_2 = numeric(),
                Species_1 = character(),
                Species_2 = character(),
                Taxa = character(),
                Habitat = character(),
                Correlation = numeric()
            )

    if (length(years_overlap) >= 3) {
        study_1_species <- study_1$SPECIES %>% unique()
        study_1_species = study_1_species[!is.na(study_1_species)]

        study_2_species <- study_2$SPECIES %>% unique()
        study_2_species = study_2_species[!is.na(study_2_species)]

        for (sp1 in study_1_species[1:10]) {
            for (sp2 in study_2_species[1:10]) {
                print(paste0(pair_counter, " | ", sp1, " | ", sp2))

                sp1_data <- study_1 %>%
                    dplyr::filter(
                        SPECIES == sp1,
                        YEAR %in% years_overlap
                    ) %>%
                    dplyr::group_by(YEAR) %>%
                    dplyr::summarize(Abundance = mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE)) %>%
                    dplyr::ungroup() %>%
                    dplyr::arrange(YEAR)
                
                sp2_data <- study_2 %>%
                    dplyr::filter(
                        SPECIES == sp2,
                        YEAR %in% years_overlap
                    ) %>%
                    dplyr::group_by(YEAR) %>%
                    dplyr::summarize(Abundance = mean(sum.allrawdata.ABUNDANCE, na.rm = TRUE)) %>%
                    dplyr::ungroup() %>%
                    dplyr::arrange(YEAR)

                sp1_data %<>% dplyr::filter(YEAR %in% sp2_data$YEAR)
                sp2_data %<>% dplyr::filter(YEAR %in% sp1_data$YEAR)

                if ((nrow(sp1_data) == nrow(sp2_data)) & nrow(sp1_data >= 3)) {
                    # cor_analysis %<>% dplyr::bind_rows(
                    # )
                    cor_pair <- cor(sp1_data$Abundance, sp2_data$Abundance)
                } else {
                    cor_pair <- NA
                }

                output_pair %<>% dplyr::bind_rows(
                    tibble(
                        Study_1 = pair_1_ID,
                        Study_2 = pair_2_ID,
                        Species_1 = sp1,
                        Species_2 = sp2,
                        Taxa = bio_pairs_10km$taxa.pairs[pair],
                        Habitat = habitat_pair,
                        Correlation = cor_pair
                    )
                )

            } ## for (sp2 in study_2_species) {
        } ## for (sp1 in study_1_species) {
    } else { ## if (length(years_overlap >= 2)) {
        output_pair <- tibble(
            Study_1 = pair_1_ID,
            Study_2 = pair_2_ID,
            Species_1 = "Too few overlapping years",
            Species_2 = "Too few overlapping years",
            Taxa = "Too few overlapping years",
            Habitat = habitat_pair,
            Correlation = NA
        )
    }

    return(output_pair)
} ## for (pair in 1:nrow(bio_pairs_10km)) {

parallel::stopCluster(cluster)

setwd(paste0(pwd, "/data/prep_biotime"))
write.table(
    cor_analysis,
    file = paste0("cor_analysis_", tu_id, ".csv"),
    sep = ",",
    col.names = TRUE,
    row.names = FALSE
)


fig <- cor_analysis %>%
    dplyr::filter(
        Correlation != 1,
        Correlation != -1,
        !is.na(Correlation)
    ) %>%
    ggplot() +
        theme_classic() +
        geom_freqpoly(
            aes(
                x = Correlation,
                color = Taxa
            ),
            position = "dodge2"
        ) +
        theme(
            legend.position = "bottom"
        )


setwd(paste0(pwd, "/figures/explore_biotime_Dec2022"))


grDevices::cairo_pdf(filename = paste0("Correlations_Distribution_TaxaPairs_", tu_id, ".pdf"), width = 20, height = 10)
    print(fig)
dev.off()



fig <- cor_analysis %>%
    dplyr::filter(
        Correlation != 1,
        Correlation != -1,
        !is.na(Correlation)
    ) %>%
    ggplot() +
        theme_classic() +
        geom_freqpoly(
            aes(
                x = Correlation,
                color = Habitat
            ),
            position = "dodge2"
        ) +
        theme(
            legend.position = "bottom"
        )


setwd(paste0(pwd, "/figures/explore_biotime_Dec2022"))


grDevices::cairo_pdf(filename = paste0("Correlations_Distribution_Habitat_", tu_id, ".pdf"), width = 20, height = 10)
    print(fig)
dev.off()
