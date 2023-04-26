#Filtering BIOtime (collated pairs)
library(tidyverse) # a suite of tidy packages useful for data manipulation
library("dplyr") # data manipulation
library("tibble") # help to create simple dataframes
library("readr") # flexible package for reading in data files
library("ggplot2") # data visualization
library("magrittr") # aids in syntax and piping


load("data/tidy/collated_pairs.RData") #full dataset with species/abundances

select <- collated.pairs %>%
  select(ID, LATITUDE, LONGITUDE)

unique_data <- select %>%
  distinct(LATITUDE, LONGITUDE, .keep_all=TRUE)

library(RColorBrewer)
n <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
pie(rep(1,n), col=sample(col_vector, n))

pal <- colorFactor(
  palette = col_vector,
  domain = unique_data$ID
)
library(leaflet)

# Create a leaflet map
m <- leaflet() %>%
  
  # Add OpenStreetMap tiles as a base layer
  addTiles() %>%
  
  # Set the initial view to the whole globe
  setView(lng = 0, lat = 0, zoom = 2) %>%
  
  # Add a marker layer for each study, colored by study ID
  addCircleMarkers(data = unique_data, lng = ~LONGITUDE, lat = ~LATITUDE,
                   radius = 5, stroke = FALSE, fillOpacity = 0.7,
                   color = ~pal(ID))

m

library(htmlwidgets
        )

htmlwidgets::saveWidget(m, file = "figures/interactive_sites_map_EB_260423.html", 
                        selfcontained = TRUE)





