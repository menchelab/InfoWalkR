#visualization functions
#---------------------------------------------------------------------

#---------------------------------------------------------------------


#' Plot Degree Distribution
#'
#' This function plots the degree distribution of a graph as a scatterplot on a double logarithmic axis.
#'
#' @param graph The input graph object.
#' @return None
#' @export
plot_degree_dist = function(graph){
# Get the degree distribution of the graph
degree_dist <- igraph::degree_distribution(graph)

# Plot the degree distribution on a double logarithmic scale as a scatterplot
plot(degree_dist, log = "xy", pch = 16, col = "#2e5b99",
     xlab = "Degree", ylab = "Frequency", main = "Degree Distribution")

# Add labels for the axes
axis(side = 1, at = c(1, 10, 100), labels = c(1, 10, 100))
axis(side = 2, at = c(0.1, 1, 10, 100), labels = c(0.1, 1, 10, 100))

# Add a grid
grid(nx = NULL, ny = NULL, lty = 1, col = "lightgray")
}



#function to show z-score based ranking of network modularity
# as lollipop plot
#' Plot Z-Score Ranking
#'
#' This function plots the z-score ranking of observed LCC size to random expectation for each network as lollipop plot.
#'
#' @param results_df The dataframe containing results to be plotted.
#' @param alpha The significance level.
#' @return None
#' @export
plot_z_score_ranking = function(results_df, alpha){

    ordered_data <- results_df[order(results_df$z_score, decreasing = FALSE), ]
    ordered_data$network_name <- factor(ordered_data$network_name,
                                        levels = ordered_data$network_name)

    p = ggplot(ordered_data, aes(x = z_score, y = network_name)) +
        geom_segment(aes(xend = 0, yend = network_name), size = 1, color = "black") +
        geom_point(aes(color = p_value < alpha), size = 2, shape = 21, stroke = 2) +
        scale_x_continuous(expand = c(0.01, 0.01)) +
        scale_color_manual(values = c("TRUE" = "#24ae2b", "FALSE" = "#dc1111")) +  # Define colors
        labs(title = "Z-Scores of Observed LCC Size to Random Expectation",
             x = "z-score",
             y = "Network") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.margin = margin(t = 2, r = 0.5, b = 2, l = 0.2, unit = "cm"))

    print(p)
}


#---------------------------------------------------------------
# function to apply and print single histograms or all in the results 
# data frame, the default is to print all
#' Plot LCC Distribution
#'
#' This function applies and prints single histograms or all in the results dataframe. By default, it prints all histograms, but optionally, you can filter by network name.
#'
#' @param results_df The dataframe containing results to be plotted.
#' @param network_filter An optional parameter to filter by network name. Default is NULL, which prints all histograms.
#' @return None
#' @export
plot_lcc_distribution = function(results_df, network_filter = NULL) {

  # Apply the custom function to each row, optionally filtered by network_name
  if (is.null(network_filter)) {
    histogram_list <- suppressMessages(apply(results_df, 1, create_histograms_detailed))
  } else {
    network_filter <- tolower(network_filter)
    subset_df <- results_df[tolower(results_df$network_name) == network_filter, ]
    histogram_list <- suppressMessages(apply(subset_df, 1, create_histograms_detailed))
  }

  # Print each histogram in the list
  for (i in 1:length(histogram_list)) {
    suppressMessages((print(histogram_list[[i]])))
  }
}


#function to apply and print overview of all histograms
#' Plot LCC Overview
#'
#' This function applies and prints an overview of all histograms for the input dataframe.
#'
#' @param results_df The dataframe containing results to be plotted.
#' @return None
#' @export
plot_lcc_overview = function(results_df){

    histogram_list = suppressMessages(apply(results_df, 1, create_histograms_overview))
    combined_plot = suppressMessages((grid.arrange(grobs = histogram_list, ncol = 3)))
    print(combined_plot)
}


#function to generate histogram for each row of df
# unpacks the vector with lcc size distribution
# returns a list of histograms and adds annotations to each
create_histograms_detailed <- function(row) {
  vector_data <- unlist(row[['lcc_size_distribution']])  # Unlist the list in the df cell
  df_tmp <- data.frame(value = vector_data)
  
  plot_name = row[['network_name']]

  p <- ggplot(df_tmp, aes(x = value)) +
    geom_histogram(fill = "#1c6a228e", color = "white") +
    labs(title = plot_name,
         x = "lcc size",
         y = "frequency") +
    theme_minimal()
  

  # Add a vertical line at the observed LCC size
  observed_lcc_size = row[['observed_lcc_size']]
  p <- p + geom_vline(xintercept = observed_lcc_size,
                                 color = "#a02222",
                                linetype = "solid", size = 0.8)
  
  # Add a text annotation for the z-score
  #z_score = row[['z_score']]
  annotation_text <- round(row[['z_score']], digits = 2)
  pb <- ggplot_build(p) #render plot to pull plot data
  y_pos <- max(range(pb$data[[1]]$y)) #extract fixed annotation position
  p <- p + geom_text(aes(x = observed_lcc_size, y = y_pos -3, label = annotation_text), 
                     size = 4, hjust = -1, vjust = -1)
  
  return(p)
}


create_histograms_overview <- function(row) {
  vector_data <- unlist(row[['lcc_size_distribution']])  # Unlist the list in the df cell
  df_tmp <- data.frame(value = vector_data)
  
  plot_name = row[['network_name']]
  annotation_text <- round(row[['z_score']], digits = 2)

  p <- ggplot(df_tmp, aes(x = value)) +
    geom_histogram(fill = "#1c6a228e", color = "white") +
    labs(title = paste(plot_name, "z_score =", annotation_text),
         x = "lcc size",
         y = "frequency") +
    theme_minimal()
  

  # Add a vertical line at the observed LCC size
  observed_lcc_size = row[['observed_lcc_size']]
  p <- p + geom_vline(xintercept = observed_lcc_size,
                                 color = "#a02222",
                                linetype = "solid", size = 0.5)
  
  
  return(p)
}
