DimMat___GraphMeasures = function(Matrix, threshold = 0.5, path_Export=NULL, visulaize=TRUE, file.name=NULL){
  #=============================================================================
  # Arguments
  #=============================================================================
  # Matrix : only corr mat?
  # threshold : for adjacency matrix





  #=============================================================================
  # Install packages
  #=============================================================================
  if (!require(igraph)) install.packages("igraph")
  library(igraph)






  #=============================================================================
  # Making Adjacency Matrix
  #=============================================================================
  # Convert the correlation matrix to a graph
  graph <- graph_from_adjacency_matrix(Matrix, mode = "undirected", weighted = TRUE, diag = FALSE)






  #=============================================================================
  # Compute graph-measures
  #=============================================================================
  # a list for saving
  Results.list = list()

  # Calculate the diameter of the graph
  Results.list$graph_diameter = graph_diameter <- diameter(graph, directed = FALSE, weights = NA)

  # Calculate the strength (weighted degree) for each node
  Results.list$node_strength = node_strength <- strength(graph, weights = E(graph)$weight)

  # Calculate PageRank for each node
  Results.list$page_rank = page_rank <- page_rank(graph, directed = FALSE)$vector

  # Clustering Coefficient
  Results.list$clustering_coeff = clustering_coeff <- transitivity(graph, type = "average")

  # Global Efficiency
  Results.list$global_efficiency = global_efficiency <- igraph::global_efficiency(graph, weights = E(graph)$weight)

  # Local Efficiency
  Results.list$local_efficiency = local_efficiency <- igraph::local_efficiency(graph)

  # Betweenness Centrality
  Results.list$betweenness = betweenness <- betweenness(graph, normalized = TRUE)

  # Modularity
  communities <- cluster_fast_greedy(graph)
  Results.list$modularity <- modularity(communities)
  print(paste("Modularity:", modularity))

  # Characteristic Path Length
  Results.list$path_length = path_length <- mean_distance(graph, directed = FALSE)


  # Small-Worldness
  Results.list$small_worldness = small_worldness = (clustering_coeff / global_efficiency) / (path_length / global_efficiency)







  #=============================================================================
  # Combine Results
  #=============================================================================
  # Combine the results into a single named vector
  measures_vector <- unlist(Results.list)

  # Now loop over the original list and assign names to the vector elements
  for(measure_name in names(Results.list)) {
    # Skip if the measure is not a vector (has length 1)
    if(length(Results.list[[measure_name]]) == 1) next

    # Generate the names for the vector elements
    measure_names <- paste0(measure_name, "___", seq_along(Results.list[[measure_name]]))

    # Find the indices in the vector that match the measure_name
    measure_indices <- grep(measure_name, names(measures_vector))

    # Assign the new names to those indices
    names(measures_vector)[measure_indices] <- measure_names
  }







  #=============================================================================
  # Visualize
  #=============================================================================
  if(visulaize){
    # Install and load ggraph if you haven't already
    if (!requireNamespace("ggraph", quietly = TRUE)) {
      install.packages("ggraph")
    }
    library(ggraph)


    # Using ggraph to plot
    p =ggraph(graph, layout = "fr") +
      geom_edge_link(aes(edge_alpha = weight), color = "grey") +  # Edges with transparency based on weight
      geom_node_point(color = "blue", size = 3) +                # Nodes as points
      theme_void()                                                # Minimal theme



    # file name
    if(!is.null(file.name)){
      file.name = "Graph"
    }


    # Save image
    if(!is.null(path_Export)){
      # Now use ggsave to save the plot to a file
      ggsave(paste0(path_Export, "/", file.name, ".png"), p)
    }
    # If you want to save as a PDF or other format, just change the extension
    # ggsave("my_network_plot.pdf", plot, width = 10, height = 8, units = "in")
  }






  #=============================================================================
  # Return 2types Results
  #=============================================================================
  list(Vector = measures_vector, List = Results.list) %>% return


}





