# Function to process the list of graphs
#calculates lcc_size_distribution_under_random expectation
#calcaultes summary statistics and z-score in relation
#to observed lcc size

#' LCC Based Network Selection
#'
#' This function performs selection of significant networks based on the size of 
#' the largest connected component (LCC) compared to random expectations.
#'
#' @param graph_list List of igraph objects representing networks.
#' @param gene_set Vector of genes of interest.
#' @param trial Number of random trials for calculation (default: 1000).
#' @param randomization Method of randomization ('random', 'dpr_exact', or 'dpr_binned') (default:'dpr_binned').
#' @param alpha Significance level for filtering (default: 0.05).
#' @return A list containing results and significant networks.
#' @export
#' @importFrom utils gsub
#' @examples
#' # Example usage:
#' # result <- lcc_based_network_selection(graph_list, gene_set)

lcc_based_network_selection = function(graph_list, gene_set, trial = 1000, randomization = 'dpr_binned', alpha = 0.05) {
  
  results = list()

  for (graph_name in names(graph_list)){
    
    graph = graph_list[[graph_name]]

    observed_lcc_size = calculate_lcc(graph, gene_set)
    lcc_size_distribution = calc_lcc_distribution_random_expectation(graph, gene_set, trial, randomization)
    mean_lcc = mean(lcc_size_distribution)
    sd_lcc = sd(lcc_size_distribution)
    z_score = (observed_lcc_size - mean_lcc)/sd_lcc #is the observed lcc-size larger than mean?
    p_value = 1-pnorm(z_score)

    clean_graph_name = gsub(".tsv$", "", graph_name)
    results[[clean_graph_name]] = list(lcc_size_distribution = lcc_size_distribution,
                                        observed_lcc_size = observed_lcc_size,
                                        mean_lcc = mean_lcc,
                                        sd_lcc = sd_lcc,
                                        z_score = z_score,
                                        p_value = p_value)
  }

  results_df = convert_results_to_df(results)

  #filter results for significant layers only
  results_df = filter_for_sig_layer_results(results_df, alpha)

  #filter graph list for significant layers only
  significant_layers = pull_significant_layers(results_df, alpha, graph_list)

  return(list(results = results_df, significant_networks = significant_layers))
}


#function to convert the named nested list from the lcc_based_network_selection
#function into a data frame
convert_results_to_df = function(results_list){

    results = results_list

    results_df = data.frame()

    for (name in names(results)) {
    sublist <- results[[name]]
    
    row_data <- data.frame(
    network_name = name,
    lcc_size_distribution = I(list(sublist$lcc_size_distribution)),
    observed_lcc_size = sublist$observed_lcc_size,
    mean_lcc_size = sublist$mean_lcc,
    sd_lcc = sublist$sd_lcc,
    z_score = sublist$z_score,
    p_value = sublist$p_value
  )
  results_df <- do.call(rbind, list(results_df, row_data))
}
    return(results_df)
}


#function to retrieve the z-scores of the significant layers only
filter_for_sig_layer_results = function(results_df, alpha){

  #filter_df = results_df %>% select(z_score, p_value, network_name)
  sig_z_scores = results_df %>% filter(p_value < alpha)

  return(sig_z_scores)
}

#############################################################################
#function to filter out significant layers(networks)
# based on a user specified significance threshold
pull_significant_layers = function(results_df, alpha, graph_list){

    filter_df = results_df %>% select(network_name, p_value)
    sig_layers = filter_df %>% filter(p_value < alpha)
    
    sig_layer_names = unique(sig_layers$network_name)
    sig_graphs = intersect(names(graph_list), sig_layer_names)

    sig_graph_list = graph_list[sig_graphs]

    return(sig_graph_list)
}


