library(reshape2)
library(ggplot2)
library(tidyverse)

data <- load(here::here("Revision 1 ecography/output/prep_data/model_data_final.Rdata"))

# Subset the dataset to include only the desired study plot (e.g., "18~1")
study_plot_of_interest <- "18~1"
subset_data <- moddat[moddat$STUDY_PLOT == study_plot_of_interest, ]


# create unique study ID column

moddat$study_id <- sub("~.*", "", moddat$STUDY_PLOT)
# Subset the data for a specific study ID, for example, study ID "18"
study_id_of_interest <- "18"

subset_data <- moddat[moddat$study_id == study_id_of_interest, ]


# CREATE CORRELATION MATRIX

# Get the list of unique species (Gn1 + Gn2)
unique_species <- unique(c(subset_data$Gn1, subset_data$Gn2))

# Initialize an empty matrix to store correlations
cor_matrix <- matrix(NA, nrow = length(unique_species), ncol = length(unique_species))
rownames(cor_matrix) <- unique_species
colnames(cor_matrix) <- unique_species

# Fill the matrix with correlation values
for (i in 1:nrow(subset_data)) {
  sp1 <- subset_data$Gn1[i]
  sp2 <- subset_data$Gn2[i]
  correlation <- subset_data$cor[i]
  
  # Since it's symmetric, fill both [sp1, sp2] and [sp2, sp1]
  cor_matrix[sp1, sp2] <- correlation
  cor_matrix[sp2, sp1] <- correlation
}

# Fill diagonal with 1 (since the correlation of a species with itself is 1)
diag(cor_matrix) <- 1


# PLOT CORRELATION MATRIX

# Convert the correlation matrix into a format suitable for ggplot

data_matrix <- melt(cor_matrix, na.rm = TRUE)
colnames(data_matrix) <- c("x", "y", "value")

# Plot the matrix
# Define a color palette (using the RdBu palette as an example)
palette <- RColorBrewer::brewer.pal(n = 11, name = "RdBu")

plot_matrix_c <- ggplot(data_matrix, aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(colors = palette, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, title = expression("Co-response matrix"~italic("C"))) +
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

# Display the plot
print(plot_matrix_c)


# PLOT DISTRIBUTION



# Ungroup the data to avoid issues caused by groupings
subset_data_ungrouped <- subset_data %>%
  dplyr::ungroup()

# Now try filtering for non-NA 'cor' values
subset_data_clean <- subset_data_ungrouped %>%
  dplyr::filter(!is.na(cor))


# Bin the correlation values, ensuring we include the full range up to 0.8
data_new <- subset_data_clean %>%
  dplyr::mutate(cor_bin = cut(cor, breaks = seq(-0.9, 0.9, by = 0.1), include.lowest = TRUE, right = FALSE)) %>%
  dplyr::group_by(cor_bin) %>%
  dplyr::summarize(frequency = n())

# Calculate the midpoint for each bin for color mapping
data_new <- data_new %>%
  dplyr::mutate(cor_bin_mid = (as.numeric(sub("\\[(.+),.+\\)", "\\1", cor_bin)) + 
                                 as.numeric(sub(".+,(.+)\\)", "\\1", cor_bin))) / 2)

# Plot the distribution of correlation values with a standardized RdBu color scale and x-axis limits
plot_distrib_resp <- ggplot(data_new, aes(x = cor_bin, y = frequency, fill = cor_bin_mid)) +
  geom_bar(stat = "identity") +  
  scale_fill_distiller(palette = "RdBu", limits = c(-0.8, 0.8)) +  # Set fixed limits for color scale
  labs(x = "Interacting species co-response", y = "Frequency", title = "Ecological Network Coherence") +
  theme_classic() +
  scale_y_continuous(breaks = c(0, 1)) +
  scale_x_discrete(labels = c(seq(-0.8, 0.8, by = 0.1), "0.8")) +  # Standardize x-axis labels to show the range -0.8 to 0.8
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels for readability
  )

# Display the plot
print(plot_distrib_resp)