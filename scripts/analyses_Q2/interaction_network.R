## Libraries
library("tidyverse")
library("rglobi")
library("igraph")
library("ggrepel")

## Data
## Question: why does using `/data/data_processing/genus_interaction_list.csv` instead result in a different network?
interactions_data <- readRDS("/Users/ryan/Windows/Documents/Post UCB/Research/CIEE_Workshop_TrophicInteractionsPredictAbundances/predicting_species_abundance/data/preprocessing/log.prop.change.interactions.RDS") %>% as_tibble()

## Initialize adjacenty matrix
interactions_data_genus <- c(interactions_data$Gn1, interactions_data$Gn2) %>% unique() %>% sort()
interactions_data_mat <- matrix(0, nrow = length(interactions_data_genus), ncol = length(interactions_data_genus))
rownames(interactions_data_mat) = interactions_data_genus
colnames(interactions_data_mat) = interactions_data_genus

## Loop over rows in the interaction data
for (interaction in 1:nrow(interactions_data)) {

  ## Add an interaction to the adjacency matrix if there is one between these two genera
  ## Note: adding the interaction in both directions = undirected network
  ## Note: multiple instances of the same interaction do not matter = unweighted network
  if (interactions_data$interaction_present[interaction] == 1) {
    interactions_data_mat[interactions_data$Gn2[interaction], interactions_data$Gn1[interaction]] = 1
    interactions_data_mat[interactions_data$Gn1[interaction], interactions_data$Gn2[interaction]] = 1
  }
}

## Only want interspecific interactions
diag(interactions_data_mat)=0

## Convert adjacency matrix into igraph network
net_data <- igraph::graph_from_adjacency_matrix(
  adjmatrix = interactions_data_mat,
  mode = "undirected",
  diag = FALSE
)

## Some plotting options to make the network look nicer
igraph.options(
  plot.layout=layout_with_fr, ## fr = a kind of force-directed layout where edges are modeled as springs
  vertex.size=5,
  vertex.label=NA,
  vertex.color="slateblue",
  edge.curved = 0
)

## Plot the network!
plot(net_data)

## Quantify node importance
centrality_betweenness <- igraph::estimate_betweenness(net_data, cutoff = -1)

## Use the same betweenness measure to find communities
## Hypothesis: genera in the same community will predict each other better than genera in different communities
communities_betweenness <- igraph::cluster_edge_betweenness(net_data)

## Convert the importance of each node into data useful in ggplot
ggdata <- centrality_betweenness %>%
  tibble::enframe() %>%
  dplyr::arrange(desc(value)) %>% ## Sort genera by their importance
  dplyr::mutate(ID = 1:length(centrality_betweenness))

## Plot genera importance (y-axis), ordered by their importance rank (x-axis)
ggdata %>% ggplot() +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  geom_point(aes(x=ID, y=value)) +
  ggrepel::geom_text_repel(
    data = ggdata, ## Not all labels will fit, so might want to filter for only large importances
    aes(x=ID, y=value, label=name)
  ) +
  labs(
    x = "Importance Rank",
    y = "Importance (betweenness centrality)"
  )
