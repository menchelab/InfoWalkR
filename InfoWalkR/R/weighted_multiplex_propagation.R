# This function performs the informed propagation algorithm. RWR, taking into
# account the weights of the network layers, as calculated proportional to the 
# z-score of the observed to expected LCC size of the given seed gene set

# Parameters:
# supra_adjacency_matrix = Csparse Matrix of the dimensions
# (all_unique_nodes * number of layers) x (all_unique_nodes * number of layers)

# graph_list = list of igraph objects (network layers)

# weighted_layer_df = dataframe containing at minimum the network_name for each layer
# and the corresponding z-score (LCC calculation) to be used as layer weight

# seed_set = seed genes (character vector)

#rank_aggregate = boolean -- determines whether the results are aggregates and summary statistics are returned

#remove_seeds = boolean -- determines if when rank_aggregate == TRUE, seeds are exlcuded from the summary statistics


#' Weighted Multiplex Propagation
#'
#' This function performs weighted multiplex propagation on a set of graphs
#'
#' @param layer_weights_df Data frame containing weights for each layer.
#' @param graph_list List of igraph objects representing graphs for each layer.
#' @param seed_set Vector of seed nodes.
#' @param seed_weights vector of seed, weights
#' @param rank_aggregate Logical indicating whether to rank results (default: FALSE).
#' @param remove_seeds Logical indicating whether to remove seed nodes from results (default: FALSE).
#' @return A dataframe containing the propagation results.
#' @export
#' @examples
#' # Example usage:
#' # weighted_multiplex_propagation(layer_weights_df, graph_list, seed_set)
#' weighted_multiplex_propagation <- function(layer_weights_df, graph_list,
#'                                           seed_set,
#'                                           rank_aggregate = FALSE, remove_seeds = FALSE) {

weighted_multiplex_propagation = function(supra_adjacency_matrix = NULL, graph_list, layer_weights_df,
                                          seed_set, seed_weights = NULL, restart_prop = 0.7,
                                          rank_aggregate = FALSE, rank_individual = FALSE,
                                          remove_seeds = FALSE) {
    
  #needed parameters
  if (is.null(supra_adjacency_matrix)){
    #construct supraadjacency matrix
    supra_adjacency_matrix = construct_supraadjacency_matrix(layer_weights_df,
                                                              graph_list)
  }
  

  el_all_layers = create_complete_edgelist_from_graph_list(graph_list)
  all_unique_nodes = get_all_unique_nodes(el_all_layers)
  number_of_layers = length(el_all_layers)
  names_of_layers = names(el_all_layers)
  
  #----------------------------------------------------------
  # add seed to the computation
  if (is.null(seed_weights)){
    
    #generate equal seed probabilities for all seeds
    seed_prop = generate_seed_probabilities_equal(all_unique_nodes, seed_set)
    #repeat seed vector for each network layers and min-max scale to 1
    seed_prop = rep(seed_prop, number_of_layers)/length(number_of_layers)
  }

  else{
    #check
    stopifnot(length(seed_weights) == length(seed_set),
            "Vectors seed_weights and seed_set must have the same length.")
            
    #produce named seed weight vector
    weighted_seeds = setNames(seed_weights, seed_set)

    #generate weighetd seed propabilities
    seed_prop = generate_seed_probabilities_unequal(all_unique_nodes, weighted_seeds)
     #repeat seed vector for each network layers and min-max scale to 1
    seed_prop = rep(seed_prop, number_of_layers)
    seed_prop = seed_prop/length(seed_prop)

  }

  #----------------------------------------------------------
 # Performing the network-based propagation
  infowalk_result = RWR(M = supra_adjacency_matrix, p_0 = seed_prop, r = restart_prop)
  
  #----------------------------------------------------------
  # annotate the data
  result_df = matrix(infowalk_result, ncol = number_of_layers, byrow = F,
                     dimnames = list(all_unique_nodes, names_of_layers)) %>% as.data.frame()
  
  result_df$gene_name = rownames(result_df)
  seed = all_unique_nodes %in% seed_set
  result_df$seed = seed

  if(rank_aggregate){
    #process results by adding summary statistics
    rank_aggregate_df = combine_rank_values(result_df, number_of_layers, names_of_layers, remove_seeds = remove_seeds)
   
   return(rank_aggregate_df)
  }

  if(rank_individual){
    #rank results within layers
    rank_individual_df = rank_in_layers(result_df, names_of_layers, remove_seeds = remove_seeds)

    return(rank_individual_df)
  }
  
  else{
    return(result_df)
  }
}


#The rank_df this function returns are:
# mean_arithmetic =
#the arithmetic mean of visiting probabilities of a given node (rowwise average across all layers)

# mean_geometric =
# geometric mean of visiting probabilities of a given node (rowise geometric mean across all layers)

# rank_all_geometric =
# expresses the visiting probabilities of a node across all layers (aggregated by geometric mean) as a rank.
# The calculation ranks the visiting probabitiles of a node across all layers and then aggregates the ranks by geometric mean.

# rank_each_geometric =
# expresses the visiting probabilities of a node across all layers as a rank,
# but the calculation first ranks the visiting probabilities for the nodes within each layer and then
# aggregates the result by geometric mean

# rank_mean_arithmetic =
# expresses the aggregated visiting probabilities for each node across all layers avg as rank

# rank_mean_geometric =
# expresses the aggregated visiting probabilities for each node across all layers GeomAvg as rank

combine_rank_values = function(result_df, number_of_layers, names_of_layers, remove_seeds = FALSE){

  if(number_of_layers > 1){
    ## 1. Arithmetic  mean on probability
    ######################################
    result_df$mean_arithmetic = rowMeans(result_df[,1:number_of_layers])
    
    ## 2. Geometric mean on probability
    ######################################
    result_df$mean_geometric = apply(result_df[,1:number_of_layers], 1, function(x) prod(x)^(1/number_of_layers)) 
    
    
    # remove seed from the result list
    if(remove_seeds){
      result_df_noseed = result_df %>% dplyr::filter(!seed)
    }
    
    else{
      result_df_noseed = result_df
    }
   
    
    ## 3. Rank all the nodes across layers together
    ######################################
    rank_alltogether_df = matrix(rank(-result_df_noseed[,1:number_of_layers]), ncol = number_of_layers, byrow = F)
    colnames(rank_alltogether_df) = names_of_layers
    rank_alltogether_df = as.data.frame(rank_alltogether_df)
    rank_alltogether_df$gene_name = result_df_noseed$gene_name
    rank_alltogether_df$rank_all_geometric = apply(rank_alltogether_df[,1:number_of_layers], 1, function(x) prod(x)^(1/number_of_layers)) 
    
    ## 4. Rank each layer separately and aggregate the results by geometric mean
    ######################################
    rank_eachlayer_df = apply(-result_df_noseed[,1:number_of_layers], 2, rank)
    rank_eachlayer_df = as.data.frame(rank_eachlayer_df) %>% mutate(gene_name = result_df_noseed$gene_name)
    rank_eachlayer_df$rank_each_geometric = apply(rank_eachlayer_df[,1:number_of_layers], 1, function(x) prod(x)^(1/number_of_layers)) 
    
    #### combine rank from different sources
    ######################################
    
    rank_df = cbind(result_df_noseed[,-c(1:number_of_layers)], 
                    rank_all_geometric= rank_alltogether_df$rank_all_geometric ,
                    rank_each_geometric =rank_eachlayer_df$rank_each_geometric )
    
    rank_df = rank_df %>% mutate(rank_mean_arithmetic = rank(-mean_arithmetic),
                                 rank_mean_geometric = rank(-mean_geometric),
                                 rank_all_geometric = rank(rank_all_geometric),
                                 rank_each_geometric = rank(rank_each_geometric)
    )
    
  }
  
  else{
   # if there is only one layer, there is no summary statistics
    
    if(remove_seeds){
      result_df = result_df %>% dplyr::filter(!seed)  
      } 
    
    result_df$rank = rank(-result_df[,1])
    rank_df = result_df
  }

  return(rank_df)

}

#############################################################################

rank_in_layers = function(result_df, names_of_layers, remove_seeds = FALSE){

    results_by_layer = list()

    if(remove_seeds){
        df = result_df %>% dplyr::filter(!seed)
    }

    else{
        df = result_df
    }

    for(layer in seq_along(names_of_layers)){
    current_layer = names_of_layers[layer]
    current_layer_result = df %>% dplyr::select(current_layer, gene_name, seed)

    current_layer_result$rank = rank(-current_layer_result[[current_layer]], ties.method = 'min')

    results_by_layer[[current_layer]] = current_layer_result

    
    }

    return(results_by_layer)
}



#############################################################################
generate_seed_probabilities_equal = function(all_unique_nodes, seed_nodes) {
  seed_prop <- ifelse(all_unique_nodes %in% seed_nodes, 1/length(seed_nodes), 0)
  return(seed_prop)
}

#----------------------------------------------------------------------------
#expects one vector with the names of all unique nodes
#expects one named vector with the names and associated weights of the seeds
generate_seed_probabilities_unequal = function(all_unique_nodes, seed_weights) {
    seed_prop <- numeric(length(all_unique_nodes))
    
    for (i in seq_along(all_unique_nodes)) {
        if (all_unique_nodes[i] %in% names(seed_weights)) {
            seed_prop[i] = seed_weights[[all_unique_nodes[i]]]
        } else {
            seed_prop[i] = 0
        }
    }
    
    return(seed_prop)
}