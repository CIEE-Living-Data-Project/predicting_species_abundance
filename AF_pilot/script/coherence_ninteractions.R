library(rglobi)
library(stringr)
library(reshape2)
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(ggpubr)

data <- load(here::here("Revision 1 ecography/output/prep_data/model_data_final.Rdata"))


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


###############################################################


# Step 1: Get genus names to search for interactions in GLOBI
vec_genus_names_interacting1 <- filtered_data$Gn1
vec_genus_names_interacting2 <- filtered_data$Gn2

# Combine genus names from both columns
vec_genus_names_tosearch <- unique(c(vec_genus_names_interacting1, vec_genus_names_interacting2))

# Create a dataframe to store genus and their number of GLOBI interactions
n_int <- numeric(length(vec_genus_names_tosearch))
df_n_int <- data.frame(Genus = vec_genus_names_tosearch, n_int = n_int)

# Step 2: Obtain interaction data from GLOBI for each genus
for (i in 1:nrow(df_n_int)) {
  # Get interactions for each genus from GLOBI
  df_int <- get_interactions_by_taxa(sourcetaxon = df_n_int$Genus[i])
  
  if (nrow(df_int) == 0) {
    df_n_int$n_int[i] <- NA  # If no interactions are found, set as NA
  } else {
    # Split the genus from species in the source taxon name
    df_int[c('genus_name', 'species_name')] <- str_split_fixed(df_int$source_taxon_name, " ", 2)
    
    # Retain only rows where the genus matches
    df_int <- df_int[df_int$genus_name == df_n_int$Genus[i],]
    
    # Record the number of unique interactions (between genus and target taxon)
    df_n_int$n_int[i] <- nrow(unique(df_int[,c("genus_name", "target_taxon_name")]))
  }
}

# Step 3: Match GLOBI interactions back to the filtered dataset
filtered_data$n_int_g1 <- numeric(nrow(filtered_data))
filtered_data$n_int_g2 <- numeric(nrow(filtered_data))

# For genus Gn1
for (i in 1:nrow(df_n_int)) {
  filtered_data$n_int_g1[which(filtered_data$Gn1 == df_n_int$Genus[i])] <- df_n_int$n_int[i]
  
  # For genus Gn2
  filtered_data$n_int_g2[which(filtered_data$Gn2 == df_n_int$Genus[i])] <- df_n_int$n_int[i]
}

# Step 4: Correlation matrix analysis

# Initialize lists to store the number of interactions, colmeans for each genus, and community sizes
interaction_list <- list()
community_size_list <- list()

# Loop through each study
for (study_id in unique(filtered_data$STUDY_ID)) {
  
  # Subset the data for the current study
  subset_data <- filtered_data[filtered_data$STUDY_ID == study_id, ]
  
  # Get the list of unique genus (Gn1 + Gn2)
  unique_genus <- unique(c(subset_data$Gn1, subset_data$Gn2))
  
  # Record community size (number of unique genus)
  community_size <- length(unique_genus)
  
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
  
  # Calculate the number of interactions for each genus (non-NA values in correlation matrix)
  num_interactions <- apply(cor_matrix, 1, function(x) sum(!is.na(x) & x != 0))
  
  # Calculate the colmeans of the absolute correlation matrix
  mean_correlations <- colMeans(abs(cor_matrix), na.rm = TRUE)
  
  # Create a dataframe for interactions, correlations, and community size for the current study
  interaction_df <- data.frame(
    Genus = names(num_interactions),
    NumInteractions = num_interactions,
    MeanCorrelation = mean_correlations,
    StudyID = study_id,
    CommunitySize = community_size,
    GLOBI_Interactions = df_n_int$n_int[match(names(num_interactions), df_n_int$Genus)],
    stringsAsFactors = FALSE
  )
  
  # Append to the list of results
  interaction_list[[study_id]] <- interaction_df
  
  # Store the community size for plotting later
  community_size_list[[study_id]] <- community_size
}

# Combine all results into a single dataframe
interaction_data <- do.call(rbind, interaction_list)

# Step 5: Plot the relationship between genus interactions and mean correlations, colored by community size
interaction_plot <- ggplot(interaction_data, aes(x = NumInteractions, y = MeanCorrelation)) +
  geom_point(aes(color = CommunitySize), size = 3, alpha = 0.7) +
  scale_color_distiller(palette = "RdBu", direction = -1) +  # Color by community size
  labs(x = "Number of Interactions", y = "Mean Absolute Correlation", 
       title = "Genus Interactions vs Mean Correlation in Study (Colored by Community Size)") +
  theme_classic() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

# Display the plot
print(interaction_plot)

interaction_plot <- ggplot(interaction_data, aes(x = NumInteractions, y = MeanCorrelation)) +
  geom_point(aes(color = StudyID), size = 3, alpha = 0.7) +
  labs(x = "Number of Interactions", y = "Mean Absolute Correlation", 
       title = "Genus Interactions vs Mean Correlation in Study") +
  theme_classic() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

# Display the plot
print(interaction_plot)

interaction_plot <- ggplot(interaction_data, aes(x = log(NumInteractions), y = MeanCorrelation, color = StudyID)) +
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(aes(group = StudyID), method = "gam", se = FALSE) +  # Group by StudyID for separate GAM lines
  labs(x = "Number of Interactions (log scale)", y = "Mean Absolute Correlation", 
       title = "Genus Interactions vs Mean Correlation in Study") +
  theme_classic() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

# Display the plot
print(interaction_plot)


# Optional: Check if the GLOBI interaction data can also be used for analysis
globi_plot <- ggplot(interaction_data, aes(x = GLOBI_Interactions, y = MeanCorrelation)) +
  geom_point(aes(color = CommunitySize), size = 3, alpha = 0.7) +
  scale_color_distiller(palette = "RdBu", direction = -1) +  # Color by community size
  labs(x = "GLOBI Interactions", y = "Mean Absolute Correlation", 
       title = "GLOBI Interactions vs Mean Correlation (Colored by Community Size)") +
  theme_classic() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5))

# Display the GLOBI plot
print(globi_plot)
