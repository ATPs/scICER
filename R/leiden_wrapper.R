#' @import igraph
#' @importFrom Matrix which
#' @importFrom methods inherits
#' @importFrom utils head
NULL

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
  weights <- if (is_weighted) igraph::E(igraph_obj)$weight else NULL
  
  # Use igraph's leiden clustering
  if (objective_function == "modularity") {
    # For modularity, use the standard leiden algorithm
    result <- igraph::cluster_leiden(
      igraph_obj,
      resolution_parameter = resolution,
      weights = weights,
      n_iterations = n_iterations,
      beta = beta
    )
  } else { # CPM
    # For CPM, use the CPM quality function
    result <- igraph::cluster_leiden(
      igraph_obj,
      objective_function = "CPM",
      resolution_parameter = resolution,
      weights = weights,
      n_iterations = n_iterations,
      beta = beta
    )
  }
  
  # Return 0-based cluster assignments for consistency with Julia code
  return(as.integer(igraph::membership(result)) - 1L)
}

#' Convert Seurat graph to igraph object
#'
#' @param seurat_graph Graph from Seurat object (typically SNN or KNN graph)
#' @param verbose Whether to print debug messages (default: FALSE)
#' @return igraph object
#' @keywords internal
graph_to_igraph <- function(seurat_graph, verbose = FALSE) {
  # Convert sparse matrix to igraph
  # Seurat graphs are typically stored as sparse matrices
  if (inherits(seurat_graph, "dgCMatrix") || inherits(seurat_graph, "Matrix")) {
    if (verbose) {
      message("GRAPH_TO_IGRAPH: Starting sparse matrix to igraph conversion...")
      message(paste("GRAPH_TO_IGRAPH: Thread context - PID:", Sys.getpid()))
      message(paste("GRAPH_TO_IGRAPH: Matrix class:", paste(class(seurat_graph), collapse = ", ")))
      message(paste("GRAPH_TO_IGRAPH: Matrix dimensions:", nrow(seurat_graph), "x", ncol(seurat_graph)))
      message(paste("GRAPH_TO_IGRAPH: Matrix storage mode:", typeof(seurat_graph)))
    }
    
    # Get non-zero entries using Matrix package methods
    if (verbose) {
      message("GRAPH_TO_IGRAPH: Extracting non-zero entries...")
      matrix_start <- Sys.time()
    }
    
    indices <- Matrix::which(seurat_graph > 0, arr.ind = TRUE)
    weights <- as.vector(seurat_graph[indices])
    
    if (verbose) {
      matrix_time <- as.numeric(difftime(Sys.time(), matrix_start, units = "secs"))
      message(paste("GRAPH_TO_IGRAPH: Extraction completed in", round(matrix_time, 3), "seconds"))
      message(paste("GRAPH_TO_IGRAPH: Found", length(weights), "non-zero entries"))
      total_entries <- nrow(seurat_graph) * ncol(seurat_graph)
      sparsity <- (1 - length(weights) / total_entries) * 100
      message(paste("GRAPH_TO_IGRAPH: Matrix sparsity:", round(sparsity, 2), "%"))
      if (length(weights) > 0) {
        message(paste("GRAPH_TO_IGRAPH: Weight range: [", round(min(weights), 4), ", ", round(max(weights), 4), "]", sep = ""))
        message(paste("GRAPH_TO_IGRAPH: Mean weight:", round(mean(weights), 4)))
        message(paste("GRAPH_TO_IGRAPH: Median weight:", round(median(weights), 4)))
      }
    }
    
    # Create edge list (1-based for igraph)
    # Keep indices as 1-based since igraph in R uses 1-based indexing
    edges_matrix <- indices
    
    if (verbose) {
      message(paste("GRAPH_TO_IGRAPH: Created edge list with", nrow(edges_matrix), "edges"))
      if (nrow(edges_matrix) > 0) {
        message(paste("GRAPH_TO_IGRAPH: Vertex index range: [", min(edges_matrix), ", ", max(edges_matrix), "]", sep = ""))
        message("GRAPH_TO_IGRAPH: Sample edges (from -> to):")
        head_edges <- utils::head(edges_matrix, 5)
        for(i in seq_len(nrow(head_edges))) {
          weight_val <- weights[i]
          message(paste("GRAPH_TO_IGRAPH:   ", head_edges[i,1], " -> ", head_edges[i,2], " (weight: ", round(weight_val, 4), ")", sep = ""))
        }
        if (nrow(edges_matrix) > 5) {
          message(paste("GRAPH_TO_IGRAPH:   ... and", nrow(edges_matrix) - 5, "more edges"))
        }
      }
    }
    
    # Create igraph object
    n_vertices <- nrow(seurat_graph)
    if (verbose) {
      message(paste("GRAPH_TO_IGRAPH: Creating empty graph with", n_vertices, "vertices"))
      graph_creation_start <- Sys.time()
    }
    
    igraph_obj <- igraph::make_empty_graph(n = n_vertices, directed = FALSE)
    
    # Add edges with weights if there are any edges
    if (nrow(edges_matrix) > 0) {
      if (verbose) {
        message("GRAPH_TO_IGRAPH: Adding edges and weights to graph...")
      }
      
      # Add edges from the matrix
      igraph_obj <- igraph::add_edges(igraph_obj, t(edges_matrix))
      igraph::E(igraph_obj)$weight <- weights
      
      if (verbose) {
        graph_creation_time <- as.numeric(difftime(Sys.time(), graph_creation_start, units = "secs"))
        message(paste("GRAPH_TO_IGRAPH: Graph creation completed in", round(graph_creation_time, 3), "seconds"))
        message(paste("GRAPH_TO_IGRAPH: Final graph: vertices =", igraph::vcount(igraph_obj), ", edges =", igraph::ecount(igraph_obj)))
        message(paste("GRAPH_TO_IGRAPH: Graph is directed:", igraph::is_directed(igraph_obj)))
        message(paste("GRAPH_TO_IGRAPH: Graph is weighted:", igraph::is_weighted(igraph_obj)))
        if (igraph::is_weighted(igraph_obj)) {
          edge_weights <- igraph::E(igraph_obj)$weight
          message(paste("GRAPH_TO_IGRAPH: Edge weights summary:"))
          message(paste("GRAPH_TO_IGRAPH:   Min:", round(min(edge_weights), 4)))
          message(paste("GRAPH_TO_IGRAPH:   Max:", round(max(edge_weights), 4)))
          message(paste("GRAPH_TO_IGRAPH:   Mean:", round(mean(edge_weights), 4)))
          message(paste("GRAPH_TO_IGRAPH:   Median:", round(median(edge_weights), 4)))
        }
      }
    } else {
      if (verbose) {
        message("GRAPH_TO_IGRAPH: WARNING - No edges to add! Graph will be empty.")
      }
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
  return("weight" %in% igraph::edge_attr_names(igraph_obj))
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
    igraph_obj <- igraph::graph_from_adjacency_matrix(adjacency_matrix, 
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