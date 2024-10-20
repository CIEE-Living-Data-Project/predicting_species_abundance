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
  filter(species_count >= 8) %>%
  pull(study_id)

# Filter the main dataset to include only the studies with at least 15 species
filtered_data <- moddat %>%
  filter(study_id %in% filtered_studies & !is.na(cor))



###########################################################
# LOOP TO CREATE CORRELATION MATRICES AND DISTRIBUTIONS FOR EACH STUDY

# Initialize lists to store the plots
cor_matrix_plots <- list()
distribution_plots <- list()

# Loop through each study and generate the plots
for (study_id in filtered_studies) {  # Use filtered_studies here
  # Subset the data for the current study ID
  subset_data <- filtered_data[filtered_data$study_id == study_id, ]
  
  # Get the list of unique species (Gn1 + Gn2)
  unique_species <- unique(c(subset_data$Gn1, subset_data$Gn2))
  
  # Initialize an empty matrix to store correlations
  cor_matrix <- matrix(NA, nrow = length(unique_species), ncol = length(unique_species))
  rownames(cor_matrix) <- unique_species
  colnames(cor_matrix) <- unique_species
  
  # Check if there are any species to work with in this study
  if (length(unique_species) < 2) {
    next  # Skip if less than two species are available
  }
  
  # Fill the matrix with correlation values, checking if species exist in the unique_species
  for (i in 1:nrow(subset_data)) {
    sp1 <- subset_data$Gn1[i]
    sp2 <- subset_data$Gn2[i]
    correlation <- subset_data$cor[i]
    
    # Ensure that both species are in the unique_species list and that sp1 and sp2 are not NA or empty
    if (!is.na(sp1) & !is.na(sp2) & sp1 %in% unique_species & sp2 %in% unique_species) {
      # Since it's symmetric, fill both [sp1, sp2] and [sp2, sp1]
      cor_matrix[sp1, sp2] <- correlation
      cor_matrix[sp2, sp1] <- correlation
    } else {
      # Print a message for debugging purposes
      print(paste("Skipping entry: sp1 =", sp1, "sp2 =", sp2))
    }
  }
  
  # Fill diagonal with 1 (since the correlation of a species with itself is 1)
  diag(cor_matrix) <- 1
  
  # Convert the correlation matrix into a format suitable for ggplot
  data_matrix <- melt(cor_matrix, na.rm = TRUE)
  colnames(data_matrix) <- c("x", "y", "value")
  
  # Plot the matrix
  palette <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")
  
  plot_matrix_c <- ggplot(data_matrix, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_gradientn(colors = rev(palette), limits = c(-1, 1), 
                         name = "Correlation") +  # Add the name for the color scale
    labs(x = NULL, y = NULL, title = paste(study_id)) +
    theme_classic() +
    theme(
      legend.position = "right",
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),  
      axis.line = element_blank(),   
      plot.title = element_text(hjust = 0.5)
    ) +
    coord_fixed()
  
  # Store the correlation matrix plot
  cor_matrix_plots[[length(cor_matrix_plots) + 1]] <- plot_matrix_c
  
  # -----------------
  # Now, create the distribution plot for the same study
  # -----------------
  
  # Ungroup the data to avoid issues caused by groupings
  subset_data_ungrouped <- subset_data %>%
    dplyr::ungroup()
  
  # Now try filtering for non-NA 'cor' values
  subset_data_clean <- subset_data_ungrouped %>%
    dplyr::filter(!is.na(cor))
  
  # Bin the correlation values
  data_new <- subset_data_clean %>%
    dplyr::mutate(cor_bin = cut(cor, breaks = seq(-0.9, 0.9, by = 0.1), include.lowest = TRUE, right = FALSE)) %>%
    dplyr::group_by(cor_bin) %>%
    dplyr::summarize(frequency = n())
  
  # Calculate the midpoint for each bin for color mapping
  data_new <- data_new %>%
    dplyr::mutate(cor_bin_mid = (as.numeric(sub("\\[(.+),.+\\)", "\\1", cor_bin)) + 
                                   as.numeric(sub(".+,(.+)\\)", "\\1", cor_bin))) / 2)
  
  # Plot the distribution of correlation values
  plot_distrib_resp <- ggplot(data_new, aes(x = cor_bin, y = frequency, fill = cor_bin_mid)) +
    geom_bar(stat = "identity") +  
    scale_fill_distiller(palette = "RdBu", limits = c(-0.8, 0.8), direction = -1,
                         name = "Correlation") +  
    labs(x = "Co-response", y = "Frequency", title = study_id) +
    theme_classic() +
    scale_y_continuous(breaks = c(0, 1)) +
    scale_x_discrete(labels = c(seq(-0.8, 0.8, by = 0.1), "0.8")) +  
    theme(
      legend.position = "top",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Store the distribution plot
  distribution_plots[[length(distribution_plots) + 1]] <- plot_distrib_resp
}

# Now, arrange the plots using ggarrange (from the ggpubr package)

# Arrange the correlation matrices in a 5x3 grid
ggarrange(plotlist = cor_matrix_plots, ncol = 4, nrow = 6,
                             common.legend = TRUE)
#ggsave(here::here("AF_pilot/output/cormatrix_plots.png"), height = 17, width = 12)


# Arrange the distribution plots in a 5x3 grid
ggarrange(plotlist = distribution_plots, ncol = 4, nrow = 6,
          common.legend = TRUE)
#ggsave(here::here("AF_pilot/output/cordistr_plots.png"), height = 17, width = 12)








####################################################################




######################## PLOT OVERALL COHERENCE


# Combine data for all studies (keep only non-NA correlation values)
combined_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::select(STUDY_ID, cor)

# Calculate the overall mean and SD for all studies combined
overall_mean <- mean(combined_data$cor, na.rm = TRUE)
overall_sd <- sd(combined_data$cor, na.rm = TRUE)

# Plot density curves for each study with transparency
density_plot <- ggplot(combined_data, aes(x = cor, color = STUDY_ID, fill = STUDY_ID)) +
  geom_density(alpha = 0.3, size = 1) +  # Use alpha for transparency and size for line thickness
  scale_color_manual(values = rainbow(length(unique(combined_data$STUDY_ID)))) +  # Different colors for each study
  scale_fill_manual(values = rainbow(length(unique(combined_data$STUDY_ID)))) +  # Matching fill colors
  labs(x = "Interacting species co-response", y = "Density", title = "Ecological Coherence across studies") +
  theme_classic() +
  
  # Add a vertical dashed black line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  # Add overall mean (solid black line) and SD (dashed black lines)
  geom_vline(xintercept = overall_mean, linetype = "solid", color = "black", size = 1.5) +
  geom_vline(xintercept = overall_mean - overall_sd, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = overall_mean + overall_sd, linetype = "dashed", color = "black", size = 1) +
  
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(density_plot)

#ggsave(here::here("AF_pilot/output/overall_coherence.png"), height = 8, width = 10)


######################## PLOT OVERALL COHERENCE - TREATMENT YES/NO


# Calculate mean and SD for each group (treatment yes/no)
summary_stats <- combined_data %>%
  group_by(treatment_yn_clean) %>%
  summarize(
    mean_cor = mean(cor, na.rm = TRUE),
    sd_cor = sd(cor, na.rm = TRUE)
  )

# Plot density curves for each study, colored by treatment status
density_plot <- ggplot(combined_data, aes(x = cor, color = treatment_yn_clean, fill = treatment_yn_clean)) +
  geom_density(alpha = 0.3, size = 1) +  # Use alpha for transparency and size for line thickness
  scale_color_manual(values = c("yes" = "red", "no" = "blue")) +  # Color based on treatment status
  scale_fill_manual(values = c("yes" = "red", "no" = "blue")) +   # Fill colors match line colors
  labs(x = "Interacting species co-response", y = "Density", title = "Density of Correlations: Treatment vs No Treatment") +
  theme_classic() +
  # Add vertical dashed black line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  # Add vertical solid lines for the mean of each group
  geom_vline(data = summary_stats, aes(xintercept = mean_cor, color = treatment_yn_clean), linetype = "solid", size = 1.2) +
  
  # Add dashed lines for mean ± SD of each group
  geom_vline(data = summary_stats, aes(xintercept = mean_cor - sd_cor, color = treatment_yn_clean), linetype = "dashed", size = 1) +
  geom_vline(data = summary_stats, aes(xintercept = mean_cor + sd_cor, color = treatment_yn_clean), linetype = "dashed", size = 1) +
  
  theme(
    legend.position = "right",  # Show legend to distinguish "yes" vs "no"
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(density_plot)


######################## PLOT OVERALL COHERENCE - CLASSES


# Combine RESOLVED.TAXA1 and RESOLVED.TAXA2 to create a single family column
combined_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::mutate(Family = ifelse(RESOLVED.TAXA1 == RESOLVED.TAXA2, RESOLVED.TAXA1, paste(RESOLVED.TAXA1, RESOLVED.TAXA2, sep = "_"))) %>%
  dplyr::select(Family, cor)

# Calculate the overall mean and SD for all families combined
overall_mean <- mean(combined_data$cor, na.rm = TRUE)
overall_sd <- sd(combined_data$cor, na.rm = TRUE)

# Plot density curves for each family with transparency
density_plot <- ggplot(combined_data, aes(x = cor, color = Family, fill = Family)) +
  geom_density(alpha = 0.3, size = 1) +  # Use alpha for transparency and size for line thickness
  scale_color_manual(values = rainbow(length(unique(combined_data$Family)))) +  # Different colors for each family
  scale_fill_manual(values = rainbow(length(unique(combined_data$Family)))) +  # Matching fill colors
  labs(x = "species co-response", y = "Density", title = "Ecological Coherence across Families") +
  theme_classic() +
  
  # Add vertical dashed black line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  # Add overall mean (solid black line) and SD (dashed black lines)
  geom_vline(xintercept = overall_mean, linetype = "solid", color = "black", size = 1.5) +
  geom_vline(xintercept = overall_mean - overall_sd, linetype = "dashed", color = "black", size = 1) +
  geom_vline(xintercept = overall_mean + overall_sd, linetype = "dashed", color = "black", size = 1) +
  
  theme(
    legend.position = "right",  # Show legend to distinguish families
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(density_plot)


######################## COHERENCE - latitude

library(ggridges)

# Filter and prepare the data
combined_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::select(STUDY_ID, cor, LATITUDE)

# Get the range of the correlation values
min_cor <- min(combined_data$cor)
max_cor <- max(combined_data$cor)

# Plot the ridgeline density plot
ridgeline_plot <- ggplot(combined_data, aes(x = cor, y = LATITUDE, group = STUDY_ID, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01) +  # Adjust the scale and height as needed
  scale_fill_distiller(palette = "RdBu", name = "correlation") +  
  labs(x = "species co-response (correlation)", y = "Latitude", 
       title = "Densities of Correlations across Latitudes") +
  theme_classic() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(min_cor, max_cor)  # Set x-axis limits to the actual range of the correlation values

# Display the plot
print(ridgeline_plot)

#ggsave(here::here("AF_pilot/output/coherence_latitude.png"), height = 13, width = 8)

######################## MEAN AND SD COHERENCE - community size

# Step 1: Calculate the mean of the correlation values for each study along with the corrected community size
mean_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::group_by(STUDY_ID) %>%
  dplyr::summarize(
    mean_cor = mean(cor, na.rm = TRUE),
    # Calculate community size as the number of unique genus from both Gn1 and Gn2
    community_size = n_distinct(c(Gn1, Gn2))  
  )

# Step 2: Plot the mean correlation vs community size
mean_community_size_plot <- ggplot(mean_data, aes(x = community_size, y = mean_cor)) +
  geom_point(size = 3, color = "black", alpha = 0.5) +
  labs(x = "Community Size", y = "Mean Correlation") +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

# Display the plot
print(mean_community_size_plot)



# Step 1: Calculate the SD of the correlation values for each study along with the corrected community size
sd_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::group_by(STUDY_ID) %>%
  dplyr::summarize(
    sd_cor = var(cor, na.rm = TRUE),
    # Calculate community size as the number of unique genus from both Gn1 and Gn2
    community_size = n_distinct(c(Gn1, Gn2))  
  )

# Step 2: Plot the SD of the correlation vs community size
var_community_size_plot <- ggplot(sd_data, aes(x = community_size, y = sd_cor)) +
  geom_point(size = 3, color = "black", alpha = 0.5) +
  labs(x = "Community Size", y = "Variance of Correlation") +
  theme_classic()

# Display the plot
print(var_community_size_plot)

ggarrange(mean_community_size_plot,
          var_community_size_plot,
          ncol = 2,
          nrow = 1)


#ggsave(here::here("AF_pilot/output/mean_variance_size.png"), height = 6, width = 9)





######################## SD COHERENCE and treatmenr

# Step 1: Calculate the SD of correlation values for each study along with the treatment status
sd_treatment_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::group_by(STUDY_ID, treatment_yn_clean) %>%
  dplyr::summarize(
    sd_cor = sd(cor, na.rm = TRUE),
    .groups = 'drop'  # Ensure proper ungrouping
  )

# Step 2: Create a boxplot to visualize the relationship between treatment and SD of correlation
sd_treatment_boxplot <- ggplot(sd_treatment_data, aes(x = treatment_yn_clean, y = sd_cor, fill = treatment_yn_clean)) +
  geom_boxplot() +
  scale_fill_manual(values = c("yes" = "red", "no" = "blue")) +  # Custom colors for treatment yes/no
  labs(x = "Treatment Status", y = "SD of Correlation", title = "SD of Correlation vs Treatment Status") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the plot title

# Display the boxplot
print(sd_treatment_boxplot)


######################## MEAN AND SD COHERENCE - latitude

# Step 1: Calculate the mean of the correlation values for each study along with the latitude
mean_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::group_by(STUDY_ID) %>%
  dplyr::summarize(
    mean_cor = mean(cor, na.rm = TRUE),
    LATITUDE = unique(LATITUDE)  # Assuming LATITUDE is constant within each study
  )

# Step 2: Plot the mean correlation vs latitude
mean_latitude_plot <- ggplot(mean_data, aes(x = LATITUDE, y = mean_cor)) +
  geom_point(size = 3, color = "black", alpha = 0.5) +
  labs(x = "Latitude", y = "Mean Correlation", title = "Mean of Correlation Distribution vs Latitude") +
  theme_classic()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")

# Display the plot
print(mean_latitude_plot)

# Step 1: Calculate the SD of the correlation values for each study along with the latitude
sd_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::group_by(STUDY_ID) %>%
  dplyr::summarize(
    sd_cor = sd(cor, na.rm = TRUE),
    LATITUDE = unique(LATITUDE)  # Assuming LATITUDE is constant within each study
  )

# Step 2: Plot the SD of the correlation vs latitude
sd_latitude_plot <- ggplot(sd_data, aes(x = LATITUDE, y = sd_cor)) +
  geom_point(size = 3,color = "black", alpha = 0.5) +
  labs(x = "Latitude", y = "SD of Correlation", title = "SD of Correlation Distribution vs Latitude") +
  theme_classic()

# Display the plot
print(sd_latitude_plot)



########################################### (SIMULATED DATA - GUILDS) AND COHERENCE

combined_data <- filtered_data %>%
  dplyr::filter(!is.na(cor)) %>%
  dplyr::select(study_id, cor)

combined_data <- combined_data %>%
  dplyr::mutate(trophic_guild = sample(c("plant", "herbivore", "predator"), 
                                       size = n(), 
                                       replace = TRUE, 
                                       prob = c(0.3, 0.6, 0.1)))

# Step 1: Simulate trophic guild data (assign randomly to the studies)
set.seed(123)  # For reproducibility
# Ensure the trophic_guild column is a factor to correctly group by categories
unique(combined_data$trophic_guild)

combined_data$trophic_guild <- factor(combined_data$trophic_guild)

table(combined_data$trophic_guild)

summary_stats_guild <- combined_data %>%
  group_by(trophic_guild) %>%
  summarize(
    trophic_guild = first(trophic_guild),  # Ensure the column is retained
    mean_cor = mean(cor, na.rm = TRUE),
    sd_cor = sd(cor, na.rm = TRUE),
    .groups = 'drop'  # Ensure proper ungrouping
  )



# Step 3: Plot density curves for each guild, colored by trophic group
density_plot_guild <- ggplot(combined_data, aes(x = cor, color = trophic_guild, fill = trophic_guild)) +
  geom_density(alpha = 0.3, size = 1) +  # Use alpha for transparency and size for line thickness
  scale_color_manual(values = c("plant" = "green", "herbivore" = "orange", "predator" = "purple")) +  # Custom colors for guilds
  scale_fill_manual(values = c("plant" = "green", "herbivore" = "orange", "predator" = "purple")) +   # Fill colors match line colors
  labs(x = "Interacting species co-response", y = "Density", title = "Density of Correlations by Trophic Guild") +
  theme_classic() +
  
  # Add vertical dashed black line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  
  # Add vertical solid lines for the mean of each guild
  geom_vline(data = summary_stats_guild, aes(xintercept = mean_cor, color = trophic_guild), linetype = "solid", size = 1.2) +
  
  # Add dashed lines for mean ± SD of each guild
  geom_vline(data = summary_stats_guild, aes(xintercept = mean_cor - sd_cor, color = trophic_guild), linetype = "dashed", size = 1) +
  geom_vline(data = summary_stats_guild, aes(xintercept = mean_cor + sd_cor, color = trophic_guild), linetype = "dashed", size = 1) +
  
  theme(
    legend.position = "right",  # Show legend to distinguish guilds
    plot.title = element_text(hjust = 0.5)
  )

# Display the plot
print(density_plot_guild)



########################################### TAXA CONTRIBUTIONS TO COHERENCE


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
  cor_matrix[is.na(cor_matrix)] <- 0
  
  # Calculate column means for each genus, excluding NA (i.e., self-correlations)
  genus_means <- colMeans(cor_matrix, na.rm = TRUE)
  
  # Create a data frame with the genus, their mean correlation, and the study ID
  genus_means_df <- data.frame(
    Genus = names(genus_means),
    MeanCorrelation = genus_means,
    StudyID = study_id,
    stringsAsFactors = FALSE
  )
  
  # Append to the list
  genus_means_list[[study_id]] <- genus_means_df
}

# Combine all genus means across all studies into a single table
genus_means_all_studies <- do.call(rbind, genus_means_list)

# Step 2: Add RESOLVED.TAXA1 for each genus from the original data
genus_means_all_studies <- genus_means_all_studies %>%
  left_join(filtered_data %>% select(Gn1, RESOLVED.TAXA1) %>% distinct(), by = c("Genus" = "Gn1"))

# Remove rows where RESOLVED.TAXA1 is NA
genus_means_all_studies <- genus_means_all_studies %>%
  filter(!is.na(RESOLVED.TAXA1))

# Step 3: Create a boxplot to show mean correlations by RESOLVED.TAXA1 groups

boxplot_taxa <- ggplot(genus_means_all_studies, aes(x = RESOLVED.TAXA1, y = MeanCorrelation)) +
  geom_boxplot() +
  theme_classic() +
  labs(x = "Taxonomic Group (RESOLVED.TAXA1)", y = "Mean Correlation", title = "Mean Genus Correlations by Taxonomic Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the boxplot
print(boxplot_taxa)

# Step 4: Sort the combined table by mean correlation
genus_means_sorted <- genus_means_all_studies %>%
  arrange(desc(MeanCorrelation))

# Display the sorted table
genus_means_sorted




######################################

library(igraph)

# Step 1: Subset the data for a given study (for example, study "18~1")
study_plot_of_interest <- "18~1"
subset_data <- filtered_data[filtered_data$STUDY_PLOT == study_plot_of_interest, ]

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


# Step 1: Add 1 to all correlations to ensure they are positive
cor_matrix <- cor_matrix + 1

# Step 2: Create a network graph from the modified correlation matrix

# Convert the correlation matrix to an edge list for igraph
edge_list <- melt(cor_matrix, na.rm = TRUE)
colnames(edge_list) <- c("from", "to", "weight")

# Create the graph object
graph <- graph_from_data_frame(edge_list, directed = FALSE)

# Perform community detection using the Louvain method
community <- cluster_louvain(graph, weights = E(graph)$weight)

# Calculate modularity
modularity_score <- modularity(community)
cat("Modularity score:", modularity_score, "\n")

# Step 3: Visualize the network and modular groups

# Set node colors based on community membership
V(graph)$community <- membership(community)
plot_colors <- rainbow(length(unique(V(graph)$community)))[V(graph)$community]

# Plot the network with nodes colored by their community
plot(graph, vertex.color = plot_colors, vertex.label = V(graph)$name,
     main = paste("Modularity of Correlation Matrix - Study:", study_plot_of_interest),
     vertex.size = 15, edge.width = E(graph)$weight * 5)

# Step 4: Investigate community structure

# Display the genera and their corresponding modular group
genus_modularity <- data.frame(
  Genus = V(graph)$name,
  Community = V(graph)$community
)
print(genus_modularity)




########################### COHERENT MODULES

pak::pak("evolqg")
library(evolqg)
library(igraph)



# Step 1: Subset the data for a given study (for example, study "18~1")
study_plot_of_interest <- "18~1"
subset_data <- filtered_data[filtered_data$STUDY_PLOT == study_plot_of_interest, ]

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


# Step 1: Add 1 to all correlations to ensure they are positive
cor_matrix <- cor_matrix + 1


# Function to compute and visualize community detection for coherence/incoherence
compute_community_modularity <- function(cor_matrix, coherence = TRUE) {
  
  if (coherence) {
    # Coherence: Keep only positive correlations, set negative to NA
    cor_matrix[cor_matrix < 0] <- NA
  } else {
    # Incoherence: Keep only negative correlations, set positive to NA, and make negative values positive
    cor_matrix[cor_matrix > 0] <- NA
    cor_matrix <- abs(cor_matrix)
    
    # Check if all values have become NA
    if (all(is.na(cor_matrix))) {
      stop("All values in the matrix are NA. Unable to compute incoherent modules.")
    }
  }
  
  # Step 2: Create a network graph from the modified correlation matrix
  edge_list <- melt(cor_matrix, na.rm = TRUE)
  colnames(edge_list) <- c("from", "to", "weight")
  
  # Create the graph object
  graph <- graph_from_data_frame(edge_list, directed = FALSE)
  
  # Check if the graph is connected
  if (!is.connected(graph)) {
    cat("Warning: The graph is disconnected. Adding small weights to ensure connectivity.\n")
    # Add a small value to edges to ensure connectivity (can use different strategies based on need)
    E(graph)$weight <- E(graph)$weight + 0.001
  }
  
  # Perform community detection using the Louvain method
  community <- cluster_louvain(graph, weights = E(graph)$weight)
  
  # Calculate modularity
  modularity_score <- modularity(community)
  cat("Modularity score:", modularity_score, "\n")
  
  # Set node colors based on community membership
  V(graph)$community <- membership(community)
  plot_colors <- rainbow(length(unique(V(graph)$community)))[V(graph)$community]
  
  # Plot the network with nodes colored by their community
  plot(graph, vertex.color = plot_colors, vertex.label = V(graph)$name,
       main = ifelse(coherence, "Coherent Modules", "Incoherent Modules"),
       vertex.size = 15, edge.width = E(graph)$weight * 5)
  
  # Return modularity score and community data
  return(list(modularity_score = modularity_score, community = membership(community)))
}


# Replace diagonal with NA to exclude self-correlations
diag(cor_matrix) <- NA

# Step 3: Compute and visualize coherent modules (positive correlations only)
coherent_results <- compute_community_modularity(cor_matrix, coherence = TRUE)

# Step 4: Compute and visualize incoherent modules (negative correlations only, converted to positive)
incoherent_results <- compute_community_modularity(cor_matrix, coherence = FALSE)



# Function to plot the modules on the correlation matrix
plot_modularity_matrix <- function(cor_matrix, community) {
  
  # Melt the correlation matrix into long format for ggplot
  data_matrix <- melt(cor_matrix, na.rm = TRUE)
  colnames(data_matrix) <- c("x", "y", "value")
  
  # Reorder rows and columns by community membership
  ordered_genus <- names(sort(community))
  data_matrix$x <- factor(data_matrix$x, levels = ordered_genus)
  data_matrix$y <- factor(data_matrix$y, levels = ordered_genus)
  
  # Plot the heatmap with genus names and module color borders
  plot_matrix <- ggplot(data_matrix, aes(x = x, y = y, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradientn(colors = rev(brewer.pal(n = 11, name = "RdBu")), limits = c(-0.8, 0.8)) +
    labs(x = "Genus", y = "Genus", title = "Correlation Matrix with Modules") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.text.y = element_text(size = 10),
      legend.position = "right"
    )
  
  # Display the plot
  print(plot_matrix)
}

# Example usage:
# Assuming `community_structure` is the output from the community detection (Louvain) step
community_structure <- coherent_results$community  # or incoherent_results$community

# Plot the correlation matrix with the modules
plot_modularity_matrix(cor_matrix, community_structure)




################################ TEST COHERENCE - TROPHIC GUILD



  # Coherence: Keep only positive correlations, set negative to NA
  cor_matrix[cor_matrix < 0] <- NA


# Step 2: Create a network graph from the modified correlation matrix
edge_list <- melt(cor_matrix, na.rm = TRUE)
colnames(edge_list) <- c("from", "to", "weight")

# Create the graph object
graph <- graph_from_data_frame(edge_list, directed = FALSE)



# Assuming 'cor_matrix' is your correlation matrix and 'trophic_guilds' is the vector of guilds

# Step 3: Simulate the assignment of trophic guilds
# Randomly assign each genus to one of the three guilds: plant, herbivore, predator
set.seed(123)  # Set seed for reproducibility
trophic_guilds <- sample(c("plant", "herbivore", "predator"), size = length(unique(c(edge_list$from, edge_list$to))), replace = TRUE)
names(trophic_guilds) <- unique(c(edge_list$from, edge_list$to))  # Ensure the names match the genus names

# Step 4: Add trophic guilds as a vertex attribute to the graph
V(graph)$guild <- trophic_guilds[V(graph)$name]

# Step 5: Compute observed modularity based on trophic guilds
observed_modularity <- modularity(graph, as.factor(V(graph)$guild))
cat("Observed Modularity:", observed_modularity, "\n")

# Step 6: Perform a Permutation Test
set.seed(123)  # For reproducibility
num_permutations <- 1000
null_modularities <- numeric(num_permutations)

for (i in 1:num_permutations) {
  # Shuffle trophic guild assignments
  shuffled_guilds <- sample(V(graph)$guild)
  
  # Compute modularity for shuffled guilds
  null_modularities[i] <- modularity(graph, as.factor(shuffled_guilds))
}

# Step 7: Calculate the p-value
p_value <- mean(null_modularities >= observed_modularity)
cat("P-value:", p_value, "\n")

# Step 8: Plot the distribution of null modularities with observed modularity
hist(null_modularities, breaks = 30, main = "Null Distribution of Modularity", xlab = "Modularity")
abline(v = observed_modularity, col = "red", lwd = 2)




####### RDA


library(vegan)

# Step 1: Prepare the correlation matrix
# Keep only positive correlations for the RDA (or use your full matrix with imputed NAs)
cor_matrix[is.na(cor_matrix)] <- 0  # Optionally replace NAs with 0 or another method

# Step 2: Simulate the assignment of trophic guilds
# Randomly assign each genus to one of the three guilds: plant, herbivore, predator
set.seed(123)  # Set seed for reproducibility
genus_names <- rownames(cor_matrix)  # Get the genus names from the correlation matrix
trophic_guilds <- sample(c("plant", "herbivore", "predator"), size = length(genus_names), replace = TRUE)
trophic_guilds <- as.factor(trophic_guilds)

# Step 3: Perform Redundancy Analysis (RDA)
# Convert correlation matrix to a distance matrix (e.g., Euclidean)
#rda_result <- rda(cor_matrix ~ trophic_guilds)

# Step 3: Perform Distance-based Redundancy Analysis (db-RDA)
dist_matrix <- as.dist(1 - cor_matrix)
db_rda_result <- capscale(dist_matrix ~ trophic_guilds)

# Step 4: Summary of the db-RDA results
summary(db_rda_result)

# Step 5: Plot the db-RDA results to visualize genus and guilds
plot(db_rda_result, scaling = 2, main = "db-RDA of Genus Correlation Matrix and Trophic Guilds")

# Step 6: Perform ANOVA to test significance of the db-RDA
anova_db_rda <- anova(db_rda_result, permutations = 999)
cat("ANOVA Results:\n")
print(anova_db_rda)

# Step 7: Extract R-squared value to check explained variance
r_squared <- RsquareAdj(db_rda_result)$r.squared
cat("R-squared:", r_squared, "\n")
