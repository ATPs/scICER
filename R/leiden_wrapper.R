#' Perform Leiden clustering on an igraph object
#'
#' @param igraph_obj igraph object to cluster
#' @param resolution Resolution parameter for clustering
#' @param objective_function Objective function ("modularity" or "CPM")
#' @param n_iterations Number of iterations
#' @param beta Beta parameter
#' @param initial_membership Initial cluster membership (optional)
#' @return Vector of cluster assignments (0-based)
#' @keywords internal
leiden_clustering <- function(igraph_obj, resolution, objective_function, 
                             n_iterations, beta, initial_membership = NULL) {
  
  # Check if igraph object has weights
  is_weighted <- is.weighted(igraph_obj)
  
  # Set up parameters for leiden clustering
  weights <- if (is_weighted) E(igraph_obj)$weight else NULL
  
  # Use igraph's leiden clustering
  if (objective_function == "modularity") {
    # For modularity, use the standard leiden algorithm
    result <- cluster_leiden(
      igraph_obj,
      resolution_parameter = resolution,
      weights = weights,
      n_iterations = n_iterations,
      beta = beta
    )
  } else { # CPM
    # For CPM, use the CPM quality function
    result <- cluster_leiden(
      igraph_obj,
      objective_function = "CPM",
      resolution_parameter = resolution,
      weights = weights,
      n_iterations = n_iterations,
      beta = beta
    )
  }
  
  # Return 0-based cluster assignments for consistency with Julia code
  return(as.integer(membership(result)) - 1L)
}

#' Convert Seurat graph to igraph object
#'
#' @param seurat_graph Graph from Seurat object (typically SNN or KNN graph)
#' @return igraph object
#' @keywords internal
graph_to_igraph <- function(seurat_graph) {
  # Convert sparse matrix to igraph
  # Seurat graphs are typically stored as sparse matrices
  if (inherits(seurat_graph, "dgCMatrix") || inherits(seurat_graph, "Matrix")) {
    # Get non-zero entries using Matrix package methods
    indices <- Matrix::which(seurat_graph > 0, arr.ind = TRUE)
    weights <- as.vector(seurat_graph[indices])
    
    # Create edge list (convert to 0-based for igraph)
    edges <- as.vector(t(indices - 1))
    
    # Create igraph object
    n_vertices <- nrow(seurat_graph)
    igraph_obj <- make_empty_graph(n = n_vertices, directed = FALSE)
    
    # Add edges with weights
    if (length(edges) > 0) {
      igraph_obj <- add_edges(igraph_obj, edges)
      E(igraph_obj)$weight <- weights
    }
    
    return(igraph_obj)
  } else {
    stop("Unsupported graph format. Expected sparse matrix.")
  }
}

#' Check if igraph object is weighted
#'
#' @param igraph_obj igraph object
#' @return Logical indicating if graph is weighted
#' @keywords internal
is.weighted.igraph <- function(igraph_obj) {
  return("weight" %in% edge_attr_names(igraph_obj))
}

#' Alternative leiden clustering using leiden package if available
#' 
#' @param adjacency_matrix Adjacency matrix
#' @param resolution Resolution parameter
#' @param objective_function Objective function
#' @param n_iterations Number of iterations
#' @param beta Beta parameter
#' @param initial_membership Initial membership
#' @return Cluster assignments
#' @keywords internal
leiden_clustering_alternative <- function(adjacency_matrix, resolution, objective_function,
                                        n_iterations, beta, initial_membership = NULL) {
  
  # Check if leiden package is available
  if (requireNamespace("leiden", quietly = TRUE)) {
    # Use leiden package directly
    partition_type <- switch(objective_function,
                            "modularity" = "ModularityVertexPartition",
                            "CPM" = "CPMVertexPartition",
                            "RBConfigurationVertexPartition")
    
    result <- leiden::leiden(
      adjacency_matrix,
      partition_type = partition_type,
      resolution_parameter = resolution,
      n_iterations = n_iterations,
      beta = beta,
      initial_membership = initial_membership
    )
    
    return(as.integer(result) - 1L)  # Convert to 0-based
  } else {
    # Fallback to igraph implementation
    igraph_obj <- graph_from_adjacency_matrix(adjacency_matrix, 
                                             mode = "undirected", 
                                             weighted = TRUE)
    return(leiden_clustering(igraph_obj, resolution, objective_function, 
                           n_iterations, beta, initial_membership))
  }
}

#' Run clustering with error handling and retries
#'
#' @param igraph_obj igraph object
#' @param resolution Resolution parameter
#' @param objective_function Objective function
#' @param n_iterations Number of iterations
#' @param beta Beta parameter
#' @param initial_membership Initial membership
#' @param max_retries Maximum number of retries
#' @return Cluster assignments
#' @keywords internal
safe_leiden_clustering <- function(igraph_obj, resolution, objective_function,
                                  n_iterations, beta, initial_membership = NULL,
                                  max_retries = 3) {
  
  for (retry in 1:max_retries) {
    tryCatch({
      result <- leiden_clustering(igraph_obj, resolution, objective_function,
                                n_iterations, beta, initial_membership)
      return(result)
    }, error = function(e) {
      if (retry == max_retries) {
        stop(paste("Leiden clustering failed after", max_retries, "retries:", e$message))
      } else {
        warning(paste("Leiden clustering attempt", retry, "failed, retrying..."))
        Sys.sleep(0.1)  # Brief pause before retry
      }
    })
  }
} 