
#--------------------------------------------------------------------- #nolint

### MAIN FUNCTION to calculate the LCC size distribution under random expectation
# offers three ways to perform node attribute randomization
# 1. randomly sampled node sets (non-degree preserving)
# 2. exact degree preserving node sets randomly sampled
# 3. degree binned sampling of random node sets
calc_lcc_distribution_random_expectation = function(graph, gene_set, trial, randomization) {
  if (randomization == 'random') {
    # Call the function for 'random' randomization
    lcc_distribution <-calc_lcc_disttribution_rand(graph, gene_set, trial)
  } else if (randomization == 'dpr_exact') {
    # Call the function for 'dpr_exact' randomization
    lcc_distribution <- calc_lcc_disttribution_dpr_exact(graph, gene_set, trial)
  } else if (randomization == 'dpr_binned') {
    # Call the function for 'dpr_binned' randomization
    lcc_distribution <- calc_lcc_disttribution_dpr_binned(graph, gene_set, trial)
  } else {
    stop("Invalid value for 'randomization'. Use 'random', 'dpr_exact', or 'dpr_binned'.")
  }
  
  return(lcc_distribution)
}



#This function  determines the largest connected component (integer)
#that a given node_set (gene_set) forms on a graph

calculate_lcc = function(graph, node_set){ # nolint

  nodes_in_graph = intersect(V(graph)$name, node_set)
  node_set_component = subgraph(graph, nodes_in_graph)
  lcc = components(node_set_component)$csize # nolint
  if(length(lcc) == 0){
    max_lcc = 0
  }
  else{
    max_lcc = max(lcc)
  }

  return(max_lcc)
}

#---------------------------------------------------------------------
### NODE LABEL RANDOMIZATION


#This function samples a random node set of equal size from the graph
sample_random_node_set = function(graph, gene_set){
    
    nodes_full_graph = V(graph)$name
    nodes_in_graph = intersect(nodes_full_graph, gene_set)
    size_of_random_node_set = length(nodes_in_graph)

    random_node_set = sample(nodes_full_graph, size_of_random_node_set)

    return(random_node_set)
}


#This function calculates the expected distribution of lcc sizes
#given a number of trials (min. recommendation 1000) using randomly
#selected node sets of the same size as the gene set of interest
#and calculating the size of the lcc these
#random samples form on the given graph

calc_lcc_disttribution_rand = function(graph, gene_set, trial){

    node_rand_result = sapply(1:trial, function(k)
    calculate_lcc(graph, sample_random_node_set(graph, gene_set)))

    return(node_rand_result)
}

#__________________________________________________________________________
#__________________________________________________________________________

### EXACT DEGREE PRESERVING NODE ATTRIBUTE RANDOMIZATION


#This function calculates the degree sequence of the seed genes on
# the network. This gives the observed degree sequence that should be
# preserved in the sampling of random nodes
calculate_deg_seg_of_seeds_on_original_network = function(graph, gene_set){

    nodes_in_graph = intersect(V(graph)$name, gene_set)
    degrees = degree(graph, nodes_in_graph, mode = "all", loops = 'false')

    return(as.vector(degrees))
}


#
sample_dpr_random_node_set = function(graph, gene_set) {
  degree_sequence <- calculate_deg_seg_of_seeds_on_original_network(graph, gene_set)

  sampled_nodes <- integer(0)

  for (i in 1:length(degree_sequence)) {
    while (TRUE) {
      random_node <- sample(V(graph), 1)
      #if the degree of the sampled node matched the degree at
      # the current index of the degree sequence and is not already in the
      # sampled_nodes set, add node
      if (degree(graph, random_node) == degree_sequence[i] &&
          !(random_node %in% sampled_nodes)) {
        sampled_nodes <- c(sampled_nodes, random_node)
        break  # Exit the while loop once a suitable node is found
      }
      # Continue the while loop until a suitable node is found
    }
  }

  return(names(sampled_nodes))
}

calc_lcc_disttribution_dpr_exact = function(graph, gene_set, trial){

  node_rand_result = sapply(1:trial, function(k)
  calculate_lcc(graph, sample_dpr_random_node_set(graph, gene_set)))

  return(node_rand_result)

}
#-------------------------------------------------------------
#-------------------------------------------------------------
### BINNED DEGREE PRESERVING NODE ATTRIBUTE RANDOMIZATION

# Function to log-bin the degree distribution with nodes in bins
#This function bins the degree distribution with log2 and fills
#the bins with the names of the nodes that have degrees that fall
#into this bin
#returns a list which contains the bin ranges and the node names
#of all the nodes within the bins

log_bin_nodes_by_degree = function(graph) {
  
  degrees = degree(graph)
  max_degree = max(degrees)
  num_bins = ceiling(log2(max_degree))
  
  binned_nodes = vector("list", length = num_bins)
  bin_ranges = character(num_bins)
  
  # fill bin with nodes that fit the degree range of the bins
  for (i in 1:num_bins) {
    lower_bound = 2^(i - 1)
    upper_bound = 2^i
    bin_ranges[i] = paste(lower_bound, upper_bound, sep = "-")
    
    # Find nodes within the current bin's range
    nodes_in_bin = V(graph)[degrees >= lower_bound & degrees < upper_bound]

    # Store nodes in the bin
    binned_nodes[[i]] = names(nodes_in_bin)

    #check if max degree node needs to be included
    max_degree_node = V(graph)$name[which.max(degree(graph))]
    if (!(max_degree_node %in% binned_nodes)){

      #add to last bin
      binned_nodes[[num_bins]] = c(binned_nodes[[num_bins]], max_degree_node)
    }
  }

  # Return the results as list
  binned_data = list(bin_range = bin_ranges, nodes = binned_nodes)

    return(binned_data)
}


#Function to calculate the sampling strategy
#It checks in which degree bins the nodes from the given gene_set
#fall and returns a count from this

get_sampling_count = function(list_of_bins, gene_set) {
  # Initialize a vector to store counts
  counts <- integer(length(list_of_bins))
  
  # Iterate through the list of bins
  for (i in seq_along(list_of_bins)) {
    # Count matches in the current bin
    counts[i] <- sum(gene_set %in% list_of_bins[[i]])
  }
  
  return(counts)
}


# Function to convert character bin ranges to numeric ranges
convert_bin_ranges = function(bin_ranges) {
  numeric_ranges <- lapply(strsplit(bin_ranges, "-"), function(x) as.numeric(x))
  return(numeric_ranges)
}


sample_dpr_bins_random_node_set = function(graph, gene_set){

    binned_data = log_bin_nodes_by_degree(graph)
    degree_bins = convert_bin_ranges(binned_data$bin_range)
    degree_bins_on_original_graph = get_sampling_count(binned_data$nodes, gene_set)

    sampled_nodes <- integer(0)

    for (i in seq_along(degree_bins_on_original_graph)){
        count = 0
        while(count < degree_bins_on_original_graph[i]){
            random_node = sample(V(graph), 1)
            node_check = degree(graph, random_node)

            lower_bound <- degree_bins[[i]][1]
            upper_bound <- degree_bins[[i]][2]
            if (node_check >= lower_bound && node_check < upper_bound) {
              sampled_nodes <- c(sampled_nodes, random_node)
              count = count +1
            }
          }
      }
      return(names(sampled_nodes))
}

calc_lcc_disttribution_dpr_binned = function(graph, gene_set, trial){

    node_rand_result = sapply(1:trial, function(k)
    calculate_lcc(graph, sample_dpr_bins_random_node_set(graph, gene_set)))

    return(node_rand_result)
}