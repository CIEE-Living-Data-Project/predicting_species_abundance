# PILOT

library(reshape2)
library(ggplot2)
library(tidyverse)

data <- load(here::here("Revision 1 ecography/output/prep_data/model_data_final.Rdata"))

####################################################
# LOOP TO CREATE CORRELATION MATRICES AND DISTRIBUTIONS FOR EACH STUDY


library(ggplot2)
library(reshape2)
library(ggpubr)


# Create unique study ID column if not already present
moddat$study_id <- sub("~.*", "", moddat$STUDY_PLOT)

# Get a list of unique study IDs
study_ids <- unique(moddat$study_id)



# Group by STUDY_ID and count the number of unique species (from Gn1 and Gn2)
species_counts <- moddat %>%
  group_by(STUDY_ID) %>%
  summarise(
    species_count = length(unique(c(Gn1, Gn2)))  # Ensure we are counting unique species from both columns
  )


# Filter out studies with fewer than 15 species
filtered_studies <- species_counts %>%
  filter(species_count >= 8) %>%
  pull(STUDY_ID)

# Filter the main dataset to include only the studies with at least 15 species
filtered_data <- moddat %>%
  filter(STUDY_ID %in% filtered_studies & !is.na(cor))



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
    scale_fill_gradientn(colors = palette, limits = c(-1, 1)) +
    labs(x = NULL, y = NULL, title = paste("Co-response matrix - Study", study_id)) +
    theme_classic() +
    theme(
      legend.position = "none",
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
    scale_fill_distiller(palette = "RdBu", limits = c(-0.8, 0.8)) +  
    labs(x = "Co-response", y = "Frequency", title = paste("Distribution - Study", study_id)) +
    theme_classic() +
    scale_y_continuous(breaks = c(0, 1)) +
    scale_x_discrete(labels = c(seq(-0.8, 0.8, by = 0.1), "0.8")) +  
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Store the distribution plot
  distribution_plots[[length(distribution_plots) + 1]] <- plot_distrib_resp
}

# Now, arrange the plots using ggarrange (from the ggpubr package)

# Arrange the correlation matrices in a 5x3 grid
cormatrix_plots <- ggarrange(plotlist = cor_matrix_plots, ncol = 3, nrow = 5)

cormatrix_plots[[1]]
ggsave(here::here("AF_pilot/output/cormatrix_plots1.png"), height = 10, width = 9)
cormatrix_plots[[2]]
ggsave(here::here("AF_pilot/output/cormatrix_plots2.png"), height = 10, width = 9)
cormatrix_plots[[3]]
ggsave(here::here("AF_pilot/output/cormatrix_plots3.png"), height = 10, width = 9)



# Arrange the distribution plots in a 5x3 grid
corditr_plots <- ggarrange(plotlist = distribution_plots, ncol = 3, nrow = 5)

corditr_plots[[1]]
ggsave(here::here("AF_pilot/output/cordistr_plots1.png"), height = 10, width = 9)
corditr_plots[[2]]
ggsave(here::here("AF_pilot/output/cordistr_plots2.png"), height = 10, width = 9)
corditr_plots[[3]]
ggsave(here::here("AF_pilot/output/cordistr_plots3.png"), height = 10, width = 9)





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
  labs(x = "Interacting species co-response", y = "Density", title = "Ecological Coherence across Families") +
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
  scale_fill_distiller(palette = "RdBu") +  
  labs(x = "Interacting species co-response (correlation)", y = "Latitude", 
       title = "Densities of Correlations across Latitudes") +
  theme_classic() +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  xlim(min_cor, max_cor)  # Set x-axis limits to the actual range of the correlation values

# Display the plot
print(ridgeline_plot)


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
  geom_point(size = 3, color = "blue") +
  geom_smooth(method = "lm", color = "darkblue", se = FALSE) +  # Optional: add a trendline
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
  geom_point(size = 3, color = "red") +
  geom_smooth(method = "lm", color = "darkred", se = FALSE) +  # Optional: add a trendline
  labs(x = "Latitude", y = "SD of Correlation", title = "SD of Correlation Distribution vs Latitude") +
  theme_classic()

# Display the plot
print(sd_latitude_plot)


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

# Replace diagonal with NA to exclude self-correlations
diag(cor_matrix) <- NA

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