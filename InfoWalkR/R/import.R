#' Read Graphs from Folder
#'
#' This function reads graph data from TSV files in a specified folder and 
#' constructs igraph objects based on the graph type flag.
#' weighting is implicit when weight column is filled
#'
#' @param folder_path The path to the folder containing TSV files.
#' @return A list of igraph objects.
#' @export
#' @importFrom igraph graph_from_data_frame
#' @importFrom utils list.files basename gsub
#' @importFrom stats read.table
#' @importFrom base stop
#' @examples
#' # Example usage:
#' # graph_list <- read_graphs_from_folder("path/to/folder")
read_graphs_from_folder = function(folder_path) {

  tsv_files = list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)
  graph_list = list()

    for (file in tsv_files) {
    edges = read.table(file, header = TRUE, sep = "\t")
    edges = na.omit(edges)
    
    # Read the graph type flag from the TSV file
    graph_type = read_graph_type_flag(file)
    
    # Construct igraph object based on the graph type flag
    if (graph_type == "directed") {
      graph = graph_from_data_frame(edges, directed = TRUE)
    } else if (graph_type == "undirected") {
      graph = graph_from_data_frame(edges, directed = FALSE)
    } else {
      # Handle error if graph type is not recognized
      stop("Invalid graph type specified in the TSV file:", graph_type)
    }
    
    graph_list[[gsub(".tsv$", "", basename(file))]] = graph
  }
  
  return(graph_list)
}


#function to automatically read the graph type (directed or undirected)
read_graph_type_flag <- function(file) {
    
  first_line <- readLines(file, n = 1)
  
  #check graph type flag
  if (grepl("^# Graph Type: directed$", first_line, ignore.case = TRUE)) {
    return("directed")
  } else if (grepl("^# Graph Type: undirected$", first_line, ignore.case = TRUE)) {
    return("undirected")
  } else {
    stop("Graph type flag is not defined in the TSV file.")
  }
}
