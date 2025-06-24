# InfoWalkR
package implementation of the MultiOme Paper

```R
library(InfoWalkR, quietly = TRUE)
```

## Data Import
```R
#read in data
path = '<path_to_your_networks>'
graph_list = read_graphs_from_folder(path)

graph_list
```

All edgelists are structured like this:
```
# Graph Type: undirected
A	B
A1BG	ABCC4
A1BG	ALB
A1BG	APLP2
A1BG	BRPF3
```

The Graph type flag is needed to determine whether the graph is directed or not.
When a third column is present that specifies edge weight, the graph is read as a weighted graph.
As long as the node type is the same in the multiplex, different edge types are allowed across different layers.

```R
#test gene set
path_gene_set = '<path_to_your_test_gene_set>'

# Read the TSV file
data <- read.delim(path_gene_set, sep = "\t", header = TRUE)

# Extract the gene names column as a vector
# in the example, th e
gene_set <- unlist(strsplit(data$all_genes, ";"))

gene_set
```

## Layer Selection

```R
# perform layer selection based on lcc that seed genes form
#alpha sets the threshold for significance
lcc_results = lcc_based_network_selection(graph_list, gene_set,
                                          trial = 100, randomization = 'random',
                                          alpha = 0.05)
#pull results
lcc_results$results

layer_weights = lcc_results$results
sig_layers = lcc_results$significant_networks
```

The lcc_results object contains:
1. the dataframe containing only the significant layers and the associated z_score --> this is for the layer weights (weighted_layer_df)
2. the graph_list containing only the significant layers

## Random Walk
| Parameter       | Description                                                                                                                                                                                                                                                                                                     |
|-----------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| rank_aggregate  | If `TRUE`, aggregates RWR ranking results across all layers and returns summary statistics. If `FALSE`, this step is skipped.                                                                                                                                                                                |
| rank_individual | If `TRUE`, ranks RWR results within each layer and returns results for all layers individually. If `FALSE`, this step is skipped.                                                                                                                                                                            |
| remove_seeds    | If `TRUE`, removes seed genes from ranking. If `FALSE`, seed genes are retained in the ranking process.                                                                                                                                                                                                      |

If both rank_aggregata = FALSE and rank_individual = FALSE, the function returns the raw results for each individual layer


```R
#perform random walk with restart on multiplex network
#the default is equal weights on the seed gene set

infowalk_results = weighted_multiplex_propagation(layer_weights_df = layer_weights,
                                                    graph_list = sig_layers,
                                                    seed_set = gene_set,
                                                    rank_aggregate = TRUE,
                                                    rank_individual = FALSE,
                                                    remove_seeds = FALSE)
head(infowalk_results)
```

#### Description of summary statistics
| Output              | Description                                                                                                                                                                                                                                          |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| mean_arithmetic        | The arithmetic mean of visiting probabilities of a given node (gene-wise average across all layers).                                                                                                                                               |
| mean_geometric         | Geometric mean of visiting probabilities of a given node (gene-wise geometric mean across all layers).                                                                                                                                              |
| rank_all_geometric     | Expresses the visiting probabilities of a node across all layers (aggregated by geometric mean) as a rank. The calculation ranks the visiting probabilities of a node across all layers and then aggregates the ranks by geometric mean.   |
| rank_each_geometric    | Expresses the visiting probabilities of a node across all layers as a rank, but the calculation first ranks the visiting probabilities for the nodes within each layer and then aggregates the result by geometric mean.                     |
| rank_mean_arithmetic   | Expresses the aggregated visiting probabilities for each node across all layers average as rank.                                                                                                                                                   |
| rank_mean_geometric    | Expresses the aggregated visiting probabilities for each node across all layers Geometric Avg as rank.                                                                                                                                              |
#### Results for individual layers

```R
infowalk_results = weighted_multiplex_propagation(layer_weights_df = layer_weights,
                                                    graph_list = sig_layers,
                                                    seed_set = gene_set,
                                                    rank_aggregate = FALSE,
                                                    rank_individual = TRUE,
                                                    remove_seeds = FALSE)
head(infowalk_results)
```

## Visualizations

```R
#plot degree distribution
plot_degree_dist(graph = sig_layers[[1]])
```
![image](https://github.com/user-attachments/assets/d805fe73-3c05-4582-9f71-15436f6e4f5c)


```R
#plot z-score ranking of layers
plot_z_score_ranking(results_df = lcc_results$results, alpha= 0.00001)
```
![image](https://github.com/user-attachments/assets/2f083f0d-ae89-4c83-a215-01e53c12ff15)

```R
#plot histograms of lcc sizes
plot_lcc_overview(results_df = lcc_results$results)
```
![image](https://github.com/user-attachments/assets/cea8235f-8dde-4178-8436-1b4fdc2c8f7b)


#### pre-calculate the supra-adjacency matrix and then perform weighted multiplex propagation
The calculation of the normalized supra-adjacency matrix takes the most time, therefore pre-computing the supra_adjacecncy matrix might be preferable.

```R
supra_adjacency_matrix = construct_supraadjacency_matrix(layer_weights, sig_layers)

infowalk_results = weighted_multiplex_propagation(supra_adjacency_matrix = supra_adjacency_matrix,
                                                    layer_weights_df = layer_weights,
                                                    graph_list = sig_layers,
                                                    seed_set = gene_set,
                                                    rank_aggregate = TRUE,
                                                    rank_individual = FALSE,
                                                    remove_seeds = FALSE)
```



