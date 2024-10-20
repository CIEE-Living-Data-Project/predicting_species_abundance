library(igraph)
library(reshape2)
library(ggplot2)

# Function to compute modules and their sizes
compute_modules <- function(cor_matrix) {
  # Coherent Modules: Keep only positive correlations, set negative to NA
  coherent_matrix <- cor_matrix
  coherent_matrix[coherent_matrix < 0] <- NA
  
  # Incoherent Modules: Keep only negative correlations, set positive to NA and take abs
  incoherent_matrix <- cor_matrix
  incoherent_matrix[incoherent_matrix > 0] <- NA
  incoherent_matrix <- abs(incoherent_matrix)
  
  results <- list(
    coherent_modules = NA,
    coherent_sizes = NA,
    incoherent_modules = NA,
    incoherent_sizes = NA
  )
  
  # Coherent modules computation
  if (sum(!is.na(coherent_matrix)) > 0) {
    edge_list_coherent <- melt(coherent_matrix, na.rm = TRUE)
    colnames(edge_list_coherent) <- c("from", "to", "weight")
    
    graph_coherent <- graph_from_data_frame(edge_list_coherent, directed = FALSE)
    if (vcount(graph_coherent) > 1) {
      community_coherent <- cluster_louvain(graph_coherent, weights = E(graph_coherent)$weight)
      results$coherent_modules <- length(community_coherent)
      results$coherent_sizes <- mean(sizes(community_coherent))
    }
  }
  
  # Incoherent modules computation
  if (sum(!is.na(incoherent_matrix)) > 0) {
    edge_list_incoherent <- melt(incoherent_matrix, na.rm = TRUE)
    colnames(edge_list_incoherent) <- c("from", "to", "weight")
    
    graph_incoherent <- graph_from_data_frame(edge_list_incoherent, directed = FALSE)
    if (vcount(graph_incoherent) > 1) {
      community_incoherent <- cluster_louvain(graph_incoherent, weights = E(graph_incoherent)$weight)
      results$incoherent_modules <- length(community_incoherent)
      results$incoherent_sizes <- mean(sizes(community_incoherent))
    }
  }
  
  return(results)
}


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


################################################### compute modules


# Assuming filtered_data contains your study data and correlation matrices

# Assuming filtered_data contains your study data and correlation matrices

study_ids <- unique(filtered_data$STUDY_ID)

# Initialize a list to store results per study
module_results <- list()

# Loop through each study
for (study_id in study_ids) {
  # Subset the data for the current study
  subset_data <- filtered_data[filtered_data$STUDY_ID == study_id, ]
  
  # Get the list of unique genus (Gn1 + Gn2)
  unique_genus <- unique(c(subset_data$Gn1, subset_data$Gn2))
  
  # Compute the community size as the number of unique genus
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
  diag(cor_matrix) <- 1
  
  # Compute modules for coherent and incoherent cases
  study_result <- compute_modules(cor_matrix)
  
  # Add the study ID and community size to the results
  study_result$study_id <- study_id
  study_result$community_size <- community_size
  
  # Append to the list of results
  module_results[[study_id]] <- study_result
}

# Combine all results into a data frame for easier plotting
module_results_df <- do.call(rbind, lapply(module_results, as.data.frame))




################################################## PLOT BY STUDY ID

# Plot the number of coherent and incoherent modules across studies
plot_modules <- ggplot(module_results_df, aes(x = reorder(study_id, -community_size))) +
  geom_point(aes(y = coherent_modules, color = "Coherent Modules"), size = 3) +
  geom_point(aes(y = incoherent_modules, color = "Incoherent Modules"), size = 3) +
  labs(x = "Study ID (Ordered by Community Size)", y = "Number of Modules") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(
    values = c("Coherent Modules" = "#B2182B",  # Strong Red (RdBu palette)
               "Incoherent Modules" = "#2166AC")  # Strong Blue (RdBu palette)
  )

# Display the plot for the number of modules
print(plot_modules)


# Plot the average module size for coherent and incoherent modules
plot_sizes <- ggplot(module_results_df, aes(x = reorder(study_id, -community_size))) +
  geom_point(aes(y = coherent_sizes, color = "Coherent Module Size"), size = 3) +
  geom_point(aes(y = incoherent_sizes, color = "Incoherent Module Size"), size = 3) +
  labs(x = "Study ID (Ordered by Community Size)", y = "Average Module Size") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(
    values = c("Coherent Module Size" = "#B2182B",  # Strong Red (RdBu palette)
               "Incoherent Module Size" = "#2166AC")  # Strong Blue (RdBu palette)
  )

# Display the plot for module sizes
print(plot_sizes)


ggarrange(plot_modules,
          plot_sizes,
          common.legend = TRUE,
          ncol = 2, 
          nrow = 1)

ggsave(here::here("AF_pilot/output/modules.png"), height = 6, width = 9)


################################################## PLOT BY COMMUNITY SIZE

# Plot the number of coherent and incoherent modules across studies, with community size on the x-axis
plot_modules <- ggplot(module_results_df, aes(x = community_size)) +
  geom_point(aes(y = coherent_modules, color = "Coherent Modules"), size = 3) +
  geom_point(aes(y = incoherent_modules, color = "Incoherent Modules"), size = 3) +
  
  labs(x = "Community Size", y = "Number of Modules", title = "Number of Coherent and Incoherent Modules by Community Size") +
  theme_minimal() +
  scale_color_manual(
    values = c("Coherent Modules" = "#B2182B",  # Strong Red (RdBu palette)
               "Incoherent Modules" = "#2166AC")  # Strong Blue (RdBu palette)
  )

# Display the plot for the number of modules
print(plot_modules)


# Plot the average module size for coherent and incoherent modules, with community size on the x-axis
plot_sizes <- ggplot(module_results_df, aes(x = community_size)) +
  geom_point(aes(y = coherent_sizes, color = "Coherent Module Size"), size = 3) +
  geom_point(aes(y = incoherent_sizes, color = "Incoherent Module Size"), size = 3) +
  
  labs(x = "Community Size", y = "Average Module Size", title = "Average Module Size by Community Size") +
  theme_minimal() +
  scale_color_manual(
    values = c("Coherent Module Size" = "#B2182B",  # Strong Red (RdBu palette)
               "Incoherent Module Size" = "#2166AC")  # Strong Blue (RdBu palette)
  )

# Display the plot for module sizes
print(plot_sizes)




##############################################

#Heterogeneity in modules' coherence

#############################################

# Assuming you have the function to compute coherent and incoherent modules from previous steps

# Initialize a list to store results per study
module_diff_results <- list()

# Loop through each study
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
  
  # Print the correlation matrix for the first study for debugging
  if (study_id == study_ids[1]) {
    print(paste("Correlation matrix for study:", study_id))
    print(cor_matrix)
  }
  
  # Compute modules for the current study
  study_result <- compute_modules(cor_matrix)
  
  # Print module detection result for debugging
  print(paste("Study ID:", study_id, "Num Modules:", length(study_result$modules)))
  
  # Check if any modules were detected
  if (length(study_result$modules) == 0) {
    next  # Skip to the next study if no modules were detected
  }
  
  # Initialize a list to store the means of correlations within each module
  module_means <- c()
  
  # Loop through each module and compute the mean correlation within the module
  for (module in study_result$modules) {
    # Extract the indices corresponding to the module
    module_indices <- match(module, rownames(cor_matrix))
    
    if (length(module_indices) > 1) {  # Ensure module has more than 1 genus
      sub_matrix <- cor_matrix[module_indices, module_indices]
      
      # Compute the mean correlation within the module (excluding NAs)
      module_mean <- mean(sub_matrix, na.rm = TRUE)
      
      # Append the module mean to the list if it's valid
      if (!is.na(module_mean)) {
        module_means <- c(module_means, module_mean)
      }
    }
  }
  
  # Check if we have valid module means before calculating variance
  if (length(module_means) > 1) {
    # Calculate variance of module means to quantify the difference
    module_variance <- var(module_means, na.rm = TRUE)
  } else {
    module_variance <- NA  # Not enough modules to compute variance
  }
  
  # Print intermediate results for debugging
  print(paste("Study ID:", study_id, "Module Means:", paste(module_means, collapse = ","), "Variance:", module_variance))
  
  # Add the results to the list
  module_diff_results[[study_id]] <- data.frame(
    study_id = study_id,
    module_variance = module_variance,
    num_modules = length(study_result$modules)
  )
}

# Combine all results into a data frame for easier plotting
module_diff_df <- do.call(rbind, module_diff_results)

# Plot the module variance across studies
module_diff_plot <- ggplot(module_diff_df, aes(x = study_id, y = module_variance)) +
  geom_point(size = 3, color = "blue") +
  labs(x = "Study ID", y = "Variance of Module Means", title = "Variance of Mean Correlations Across Modules per Study") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(module_diff_plot)


