## export function to write .tsv edfelists from the named graph list
write_graphs_to_tsv <- function(graph_list, output_dir = "./") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  for (network_name in names(graph_list)) {
    graph <- graph_list[[network_name]]
    
    # create edge list from each igraph
    edgelist = as_edgelist(graph, names = TRUE)
    
    # construct the output file path
    output_file <- file.path(output_dir, paste0(network_name, ".tsv"))

    #check for graph type
    if (is_directed(graph)) {
      graph_type <- "directed"
    } else {
      graph_type <- "undirected"
    }

    #insert graph type flag into tsv
    writeLines(paste0("# Graph Type: ", graph_type), output_file)
    
    # write the edge list to tsv
    suppressWarnings(
    write.table(edgelist, file = output_file, append = TRUE, sep = "\t", row.names = FALSE, quote = FALSE)
    )
  }
}