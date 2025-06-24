
### CONSTRUCT SUPRA-ADJACENCY MATRIX
#' Construct supra-adjacency matrix
#'
#' This function constructs a supra-adjacency matrix representing the connectivity 
#' between nodes across multiple network layers.
#'
#' @param weighted_layer_df A data frame containing weighted layer information.
#' @param graph_list A list of igraph objects representing network layers.
#' @return A supra-adjacency matrix representing the connectivity between nodes 
#' across multiple network layers.
#' @export
construct_supraadjacency_matrix = function(weighted_layer_df, graph_list){

  #clean store edge lists for all network layers
  complete_edgelist = create_complete_edgelist_from_graph_list(graph_list)

  #call function to clean graph list
  graph_list = clean_graph_list(graph_list, complete_edgelist)

  all_unique_nodes = get_all_unique_nodes(complete_edgelist)

  number_of_layers = length(complete_edgelist)
  names_of_layers = names(complete_edgelist)
  
  #CALCULATE INTRA TRANSITIONAL PROBABILITIES
  # transitional matrix for each network layer (witin layer transitional probabilities)
  itml = lapply(graph_list, function(x) build_transitional_matrix(x, all_unique_nodes))

  #CALCULATE INTER TRANSITIONAL PROBABILITIES
  # transitional matrix between network layers
  layer_weights = calculate_layer_weights(weighted_layer_df$z_score)
  
  #COMBINED INTRA AND INTRA AND INTER TRANSITIONAL PROBABILITIES
  # supra-transitional matrix
  S = build_supratransitional_matrix(itml, layer_weights)
  
  return(S)
}


#############################################################################
#get unsorted list of all unique node names
#from a network
#or from all networks, if complete_edgelist is given
# a complete edge list is a named list of
# edgelists as tibble
get_all_unique_nodes = function(edgelist) {

 
  if (class(edgelist)== 'list'){
    #extract unique node names from all edgelists in the list
    unique_nodes = unique(unlist(lapply(edgelist, function(edgelist) unique(c(edgelist[,1], edgelist[,2])))))
    
    return(unique_nodes)
  }

  else if(class(edgelist)== 'data.frame'){
    #extract unique nodes from single edgelist
    unique_nodes = unique(c(edgelist[,1], edgelist[,2]))

    return(unique_nodes)
  }
  
}

###########################################################################
#GRAPH LIST CLEANING FUNCTION
#retain the edges which are in the cleaned edgelist
clean_graph_list = function(graph_list, processed_edgelist){

    cleaned_graph_list = list()

    for(network_name in names(graph_list)){
         graph = graph_list[[network_name]]
         clean_edgelist = processed_edgelist[[network_name]]

        #retrieve whether original graph is directed
         ifelse(is_directed(graph), directed_status <- TRUE, directed_status <- FALSE)

         cleaned_graph = graph_from_data_frame(d = clean_edgelist, directed = directed_status)

        cleaned_graph_list[[network_name]] = cleaned_graph
    }

    return(cleaned_graph_list)
}


######################################################
create_complete_edgelist_from_graph_list = function(graph_list){

  complete_edgelist = retrieve_edgelist(graph_list)

  #clean up complete_edgelist by calling process_edgelist()
  complete_edgelist = lapply(complete_edgelist, function(x) process_edgelist(x, distinct = TRUE))

  return(complete_edgelist)
}

#---------------------------------------------------------

### function to retrieve the edgelists from all network layers and
# store them in a named list of edgelists (as tibbles)
retrieve_edgelist = function(graph_list) {
    edgelist <- list()
    
    for (network_name in names(graph_list)) {
      graph = graph_list[[network_name]]
      
        edgelist[[network_name]] = as_data_frame(graph) 
      }
    

    return(edgelist)
}

#-----------------------------------------------------------------

#this function cleans the edgelist of each network by:
# 1. removing unlabeled nodes
# 2. removing loops
# 3. keeping only unique edges to eliminate duplicate edges


process_edgelist = function(edgelist, distinct = TRUE) {
  
  # Remove edges with at least one node without label
  #empty string case
  edgelist = edgelist[!(edgelist[,1] == "" | edgelist[,2] == ""), , drop = FALSE]
  #is.na
  edgelist = edgelist[!(is.na(edgelist[,1]) | is.na(edgelist[,2])), , drop = FALSE]


  # Remove edges where two nodes are the same (self-loops)
  edgelist = edgelist[edgelist[, 1] != edgelist[, 2], , drop = FALSE]

  # Keep only unique rows (unique edges)
  if (distinct) {
    edgelist = edgelist %>% distinct(select(edgelist, 1, 2), .keep_all = TRUE)
  }

  # Standardize column names
  if (ncol(edgelist) == 2) {
    colnames(edgelist) <- c("from", "to")
  } else if (ncol(edgelist) == 3) {
    colnames(edgelist) <- c("from", "to", "weight")
  } else {
    stop("Invalid number of columns in the edge list.")
  }

  return(edgelist)
}


#################################################################################

### INTRA TRANSITIONAL PROBABILITIES

#################################
### FOR UNDIRECTED NETWORKS ####
#################################

#function to build the within-layer transitional matrix
#for a single network
# normalized so that transitional probabilities are
# between 0 and 1

#complete_edge_list = list(edgelist1, edgelist2, edgelist3)
#all_unique_nodes = get_all_unique_nodes(complete_edgelist)

build_transitional_matrix = function(graph, all_unique_nodes){

    #retrieve edgelist
    edgelist = as_data_frame(graph)
    
    
    #positional vectors based on all_unique_nodes vector
    A = as.numeric(factor(edgelist[,1], levels = all_unique_nodes))
    B = as.numeric(factor(edgelist[,2], levels = all_unique_nodes))


    #get dimension of adjacency matrix
    num_nodes = length(all_unique_nodes)

    #intialize
    row_indices = c()
    col_indices = c()

    #record indices for the sparse matrix
    if (is_directed(graph)){
        #do assymmetrical
        for (i in 1:length(A)) {
        row_indices = c(row_indices, A[i])  # Add only A to B for directed network
        col_indices = c(col_indices, B[i])
    }
    }

    else {
        #do symmetrical if network is undirected
        for (i in 1:length(A)) {
        row_indices = c(row_indices, A[i], B[i])
        col_indices = c(col_indices, B[i], A[i])
        }
    }
    
    #create sparse matrix

    #retrieve edge weights
    if(is_weighted(graph)){

        weights = edgelist$weight  

        if(is_directed(graph)){
            weight_indices = weights
        }

        else{            
            weight_indices = rep(weights, each = 2)
            
            
        }
        

      #create transitional matrix for weighted graphs
       A = sparseMatrix(i = row_indices, j = col_indices,
                                    x = weight_indices, dims = c(num_nodes, num_nodes))
    }

    #create transitional matrix for unweighted graphs
    else {
      A = sparseMatrix(i = row_indices, j = col_indices,
                                    x = 1, dims = c(num_nodes, num_nodes))
    }
    
    
    # Set matrix dimension names
    dimnames(A) <- list(all_unique_nodes, all_unique_nodes)


    M = A %*% Matrix::Diagonal(x = 1 / Matrix::colSums(A))
   
    return(M)

}

#################################################################################
### INTER TRANSITIONAL PROBABILITIES

#function to calculate the between_layer transitional matrix
# normalized so that transitional probabilities are
# between 0 and 1

#p[i,j] = the probability of transitioning
#from layer i to layer j
# if z-score of j > z-score of i, p is small
# if z-score of j < z-score of i, p is larger

#input w is a vector containing the z-scores the the layers

calculate_layer_weights = function(w){

  if (any(w < 0)) {
    stop("Input contains values below 0. All input values must be non-negative.")
  }

  L = length(w)
  layer_weights = matrix(0 , nrow = L, ncol = L)
  for(i in 1:L){
    for(j in setdiff(1:L,i)){
      layer_weights[i,j] = (min(1, w[i]/w[j]))/(L)
    }
  }
  diag(layer_weights) = 1 - colSums(layer_weights)
  return(layer_weights)
}

#################################################################################

# COMBINED INTRA AND INTRA AND INTER TRANSITIONAL PROBABILITIES

#function to calculate supra-adjacency matrix
#containing the normalized transition probabilities within layers and
#between layers

#itml - intratrasitional matrix list

build_supratransitional_matrix = function(itml, layer_weights){
  n_gene = nrow(itml[[1]])
  number_layers = length(itml)
  for(i in 1:number_layers){
    
    for(j in 1:number_layers){
      
      # determine the block i,j value if it is the intra- or interlayer transitional
      # interlayer transitional is simply the diagonal
      if(i == j){blockmat = itml[[i]]} else blockmat = Diagonal(n_gene) 
      # adjust the blockmat with the weight
      weighted_blockmat = layer_weights[i,j]*blockmat
      
      # the supra-transitional is adding up each row separately:
      # if j=1 (first column, the beginning of the block), assign it to the block matrix,
      #otherwise append it to the right
      if(j==1){rowblock = weighted_blockmat}
      else rowblock = cbind(rowblock, weighted_blockmat)
    }
    
    # add the rowblock to the existing ones
    if(i==1){S = rowblock}
    else S = rbind(S, rowblock)
   # print(paste0(i, " out of ", nL, " layers completed"))
  }
  
  # correct for the case for nodes don't exist in all layers, leading to ColSums
  #less than 1 for some elements. 
  #Normalising it all solves the problem
  S = S %*% Matrix::Diagonal(x = 1 / Matrix::colSums(S))
  return(S)
}