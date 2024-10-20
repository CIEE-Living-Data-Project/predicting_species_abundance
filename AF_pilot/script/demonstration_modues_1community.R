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



############################ Illustrate modules

# Subset the data for study 333
subset_data <- filtered_data[filtered_data$STUDY_ID == "18", ]

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

# Detect and visualize coherent modules (positive correlations only)
coherent_matrix <- cor_matrix
coherent_matrix[coherent_matrix < 0] <- NA  # Keep only positive correlations

# Plot the coherent correlation matrix
palette <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")
plot_coherent <- ggplot(melt(coherent_matrix, na.rm = TRUE), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = rev(palette), limits = c(-1, 1), na.value = "white", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Coherent Modules for Study 333", x = "Genus", y = "Genus")

# Display the plot
print(plot_coherent)





# Step 1: Subset the data for the current study
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

# --------------------------------------
# Step 2: Coherent Modules (Positive Correlations)
# --------------------------------------

# Keep only positive correlations for coherent modules
coherent_matrix <- cor_matrix
coherent_matrix[coherent_matrix < 0] <- NA  # Set negative values to NA

# Create an igraph object from the correlation matrix
edge_list_coherent <- melt(coherent_matrix, na.rm = TRUE)
colnames(edge_list_coherent) <- c("from", "to", "weight")
graph_coherent <- graph_from_data_frame(edge_list_coherent, directed = FALSE)

# Detect the modules using the Louvain method
coherent_modules <- cluster_louvain(graph_coherent, weights = E(graph_coherent)$weight)

set.seed(40)
# Plot the network with module colors
plot(graph_coherent, 
     vertex.color = membership(coherent_modules), 
     vertex.size = 14, 
     vertex.label.cex = 0.8,
     edge.color = "darkred",  # Set the edge color to dark red
     main = "Coherent Modules in Study 333", 
     edge.width = E(graph_coherent)$weight * 9, 
     layout = layout_with_fr(graph_coherent))



# --------------------------------------
# Step 3: Incoherent Modules (Negative Correlations)
# --------------------------------------

# Keep only negative correlations for incoherent modules, set positive values to NA
incoherent_matrix <- cor_matrix
incoherent_matrix[incoherent_matrix > 0] <- NA  # Set positive values to NA
incoherent_matrix <- abs(incoherent_matrix)  # Convert negative correlations to positive

# Create an igraph object from the incoherent matrix
edge_list_incoherent <- melt(incoherent_matrix, na.rm = TRUE)
colnames(edge_list_incoherent) <- c("from", "to", "weight")
graph_incoherent <- graph_from_data_frame(edge_list_incoherent, directed = FALSE)

# Detect modules for incoherent correlations using the Louvain method
incoherent_modules <- cluster_louvain(graph_incoherent, weights = E(graph_incoherent)$weight)
set.seed(40)
# Plot the network with module colors for incoherent correlations
plot(graph_incoherent, 
     vertex.color = membership(incoherent_modules), 
     vertex.size = 14, 
     vertex.label.cex = 0.8,
     edge.color = "darkblue",  # Set the edge color to dark red
     main = "Incoherent Modules in Study 333", 
     edge.width = E(graph_incoherent)$weight * 9, 
     layout = layout_with_fr(graph_incoherent))



############# incoherent modules


# Detect and visualize incoherent modules (negative correlations only)
incoherent_matrix <- cor_matrix
incoherent_matrix[incoherent_matrix > 0] <- NA  # Keep only negative correlations
incoherent_matrix <- abs(incoherent_matrix)  # Convert negative correlations to positive values for plotting

# Plot the incoherent correlation matrix
plot_incoherent <- ggplot(melt(incoherent_matrix, na.rm = TRUE), aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = palette, limits = c(-1, 1), na.value = "white", name = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Incoherent Modules for Study 333", x = "Genus", y = "Genus")

# Display the plot
print(plot_incoherent)

