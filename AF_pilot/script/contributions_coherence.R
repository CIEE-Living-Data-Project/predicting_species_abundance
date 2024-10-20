
# PILOT

library(reshape2)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggpubr)

data <- load(here::here("Revision 1 ecography/output/prep_data/model_data_final.Rdata"))

#################################################### FILTER DATA


# Create unique study ID column if not already present
moddat$study_id <- sub("~.*", "", moddat$STUDY_PLOT)

# Get a list of unique study IDs
study_ids <- unique(moddat$study_id)



# Group by STUDY_ID and count the number of unique species (from Gn1 and Gn2)
species_counts <- moddat %>%
  group_by(study_id) %>%
  summarise(
    species_count = length(unique(c(Gn1, Gn2)))  # Ensure we are counting unique species from both columns
  )


# Filter out studies with fewer than 15 species
filtered_studies <- species_counts %>%
  filter(species_count >= 10) %>%
  pull(study_id)

# Filter the main dataset to include only the studies with at least 15 species
filtered_data <- moddat %>%
  filter(study_id %in% filtered_studies & !is.na(cor))



####################### CONTRIBUTIONS TO EXTREME COHERENCE



# Initialize an empty list to store the genus means across studies
genus_means_list <- list()

# Loop through each study
study_ids <- unique(filtered_data$STUDY_ID)

for (study_id in study_ids) {
  
  # Subset the data for the current study
  subset_data <- filtered_data[filtered_data$STUDY_ID == study_id, ]
  
  # Get the list of unique genus (Gn1 + Gn2)
  unique_genus <- unique(c(subset_data$Gn1, subset_data$Gn2))
  
  # Initialize an empty correlation matrix
  cor_matrix <- matrix(NA, nrow = length(unique_genus), ncol = length(unique_genus))
  rownames(cor_matrix) <- unique_genus
  colnames(cor_matrix) <- unique_genus
  
  # Fill the correlation matrix
  for (i in 1:nrow(subset_data)) {
    g1 <- subset_data$Gn1[i]
    g2 <- subset_data$Gn2[i]
    cor_value <- subset_data$cor[i]
    
    cor_matrix[g1, g2] <- cor_value
    cor_matrix[g2, g1] <- cor_value
  }
  
  # Replace the diagonal with NA to exclude self-correlations
  diag(cor_matrix) <- NA
  
  # Replace any other NA values with 0 (if correlations are missing)
 # cor_matrix[is.na(cor_matrix)] <- 0
  
  cor_matrix <- abs(cor_matrix)
  
  # Calculate column means for each genus, excluding NA (i.e., self-correlations)
  genus_means <- colMeans(cor_matrix, na.rm = TRUE)
  
  # Create a data frame with the genus, their mean correlation, and the study ID
  genus_means_df <- data.frame(
    Genus = names(genus_means),
    SumCorrelation = genus_means,
    StudyID = study_id,
    stringsAsFactors = FALSE
  )
  
  # Append to the list
  genus_means_list[[study_id]] <- genus_means_df
}

# Combine all genus means across all studies into a single table
genus_means_all_studies <- do.call(rbind, genus_means_list)



# Create a boxplot to visualize the distribution of SumCorrelation per StudyID
boxplot_meanCorrelation <- ggplot(genus_means_all_studies, aes(x = StudyID, y = SumCorrelation)) +
  geom_boxplot(fill = "lightblue", color = "darkblue") +  # Boxplot aesthetics
  labs(x = "Study ID", y = "Mean Correlations", title = "Distribution of Sum of Correlations per Study") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels for readability
    plot.title = element_text(hjust = 0.5)  # Center the plot title
  )

# Display the plot
print(boxplot_meanCorrelation)


ggsave(here::here("AF_pilot/output/boxplot_contribution_coherence.png"), height = 6, width = 8)
