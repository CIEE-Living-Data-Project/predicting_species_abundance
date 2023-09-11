## Libraries
library("dplyr")
library("tibble")
library("magrittr")
library("rglobi")
library("ggplot2")
library("ggrepel")
library("igraph")

### Data
interactions_data <- readRDS("data/data_processing/log.prop.change.interactions.RDS") %>% as_tibble()

## slopes <-
load("outputs/Aug2023/randomslopes_q1model.Rdata")
slopes %<>% as_tibble()

## intmods_terr <-
load("outputs/Aug2023/intercept_only_models.Rdata")


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

## Uncomment if want to plot the network to see what it looks like (eg are there obvious clusters/communities?)
# ## Some plotting options to make the network look nicer
# igraph.options(
#   plot.layout=layout_with_fr, ## fr = a kind of force-directed layout where edges are modeled as springs
#   vertex.size=5,
#   vertex.label=NA,
#   vertex.color="slateblue",
#   edge.curved = 0
# )

# ## Plot the network!
# plot(net_data)

## Quantify node importance
centrality_degree_out <- igraph::degree(net_data, mode = "out") %>% tibble::enframe()
colnames(centrality_degree_out) = c("Name", "Centrality_Degree_Out")

centrality_degree_in <- igraph::degree(net_data, mode = "in") %>% tibble::enframe()
colnames(centrality_degree_in) = c("Name", "Centrality_Degree_In")

centrality_betweenness <- igraph::estimate_betweenness(net_data, cutoff = -1) %>% tibble::enframe() ## cutoff = zero or negative = no limit to length of paths being considered
colnames(centrality_betweenness) = c("Name", "Centrality_Betweenness")

centrality_eigenvector <- igraph::eigen_centrality(net_data, directed = FALSE)$vector %>% tibble::enframe() ## Change directed parameter if make a directed network
colnames(centrality_eigenvector) = c("Name", "Centrality_Eigenvector")

centrality_closeness <- igraph::closeness(net_data) %>% tibble::enframe()
colnames(centrality_closeness) = c("Name", "Centrality_Closeness")

centrality_harmonic <- igraph::harmonic_centrality(net_data, cutoff = -1) %>% tibble::enframe() ## cutoff = zero or negative = no limit to length of paths being considered
colnames(centrality_harmonic) = c("Name", "Centrality_Harmonic")








centrality_nodes <- centrality_degree_out %>%
  dplyr::full_join(centrality_degree_in, by = "Name") %>%
  dplyr::full_join(centrality_betweenness, by = "Name") %>%
  dplyr::full_join(centrality_eigenvector, by = "Name") %>%
  dplyr::full_join(centrality_closeness, by = "Name") %>%
  dplyr::full_join(centrality_harmonic, by = "Name")



## Centrality of edges! Which are easier to compare to the slopes
centrality_betweenness_edges <- cbind(paste0(ends(net_data, es = E(net_data))[,1], "_", ends(net_data, es = E(net_data))[,2]), Betweenness = edge_betweenness(net_data, directed = FALSE, cutoff = -1)) %>% as_tibble()
colnames(centrality_betweenness_edges) = c("Interaction", "Centrality_Betweenness_Edge") ## Since these are undirected networks right now, the from and to nodes can be either Gn1 or Gn2. This is not true if directed networks are used! Use caution in that case to get the directionality correct.


slopes %<>% as_tibble() %>% dplyr::mutate(Gn1 = sapply(slopes$UniquePairID, function(x){strsplit(x, split = "_")[[1]][1]}), Gn2 = sapply(slopes$UniquePairID, function(x){strsplit(x, split = "_")[[1]][2]}))

slopes %<>% as_tibble() %>% dplyr::mutate(Interaction = paste0(Gn1, "_", Gn2))

centrality_edges <- dplyr::full_join(centrality_betweenness_edges, slopes, by = "Interaction") ## Full join so don't accidentally delete any useful indicators of these two datasets not matching up correctly




## Use the same betweenness centrality measure to find communities
## Hypothesis: genera in the same community will predict each other better than genera in different communities
## Note: It seems this way of doing community detection can put the same Genera in multiple communities, which we might want to change
communities_betweenness <- igraph::cluster_edge_betweenness(net_data)

communities <- tibble(
  Name = character(),
  Group = numeric()
)
for (group in seq_along(communities_betweenness)) {
  communities %<>% dplyr::bind_rows(
    tibble(
      Name = communities_betweenness[[1]],
      Group = group
    )
  )
}


## Save all output
setwd("outputs/Aug2023")
saveRDS(object = centrality_nodes, file = "centrality_nodes.RDS")
saveRDS(object = centrality_edges, file = "centrality_edges.RDS")
saveRDS(object = communities, file = "communities.RDS")


## Plot node centrality vs average slope for a taxa
avg_slopes <- slopes %>%
  dplyr::group_by(Gn1) %>%
  dplyr::summarize(Average_Slope = mean(Estimate.Prop.Change.Gn2))
colnames(avg_slopes) = c("Name", "Average_Slope")

dplyr::full_join(centrality_nodes, avg_slopes, by = "Name") %>%
  ggplot(
    aes(
      x = Centrality_Degree_In, ## Can change to other centralities
      y = Average_Slope
    )
  ) +
  theme_minimal() +
  geom_point() +
  geom_smooth(method = "lm")
