#' @importFrom Matrix summary
#' @importFrom utils head
NULL

#' Perform Leiden clustering on an igraph object
#'
#' @description
#' Internal wrapper around \code{igraph::cluster_leiden()} with consistent output.
#'
#' @details
#' This helper standardizes two behaviors used throughout scICER:
#' \itemize{
#'   \item automatically passing edge weights when present;
#'   \item returning integer cluster labels in 0-based indexing.
#' }
#'
#' For \code{objective_function = "CPM"}, CPM mode is explicitly requested.
#' For \code{"modularity"}, the same API is used with modularity semantics.
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
  weights <- if (igraph::is_weighted(igraph_obj)) igraph::E(igraph_obj)$weight else NULL

  if (objective_function == "modularity") {
    result <- igraph::cluster_leiden(
      igraph_obj,
      resolution = resolution,
      weights = weights,
      n_iterations = n_iterations,
      beta = beta,
      initial_membership = initial_membership
    )
  } else {
    result <- igraph::cluster_leiden(
      igraph_obj,
      objective_function = "CPM",
      resolution = resolution,
      weights = weights,
      n_iterations = n_iterations,
      beta = beta,
      initial_membership = initial_membership
    )
  }

  as.integer(igraph::membership(result)) - 1L
}

#' Convert Seurat graph to igraph object
#'
#' @description
#' Converts a sparse Seurat graph matrix into an undirected weighted igraph.
#'
#' @details
#' Conversion uses \code{Matrix::summary()} triplets (\code{i}, \code{j}, \code{x})
#' to avoid dense intermediate objects and to preserve existing edge semantics in
#' the sparse graph representation. This includes self-loops and repeated
#' directional entries when they are present in the source matrix.
#'
#' Verbose mode emits structured timing/statistics diagnostics useful for
#' profiling large graphs.
#'
#' @param seurat_graph Graph from Seurat object (typically SNN or KNN graph)
#' @param verbose Whether to print debug messages (default: FALSE)
#' @return igraph object
#' @keywords internal
graph_to_igraph <- function(seurat_graph, verbose = FALSE) {
  if (!(inherits(seurat_graph, "dgCMatrix") || inherits(seurat_graph, "Matrix"))) {
    stop("Unsupported graph format. Expected sparse matrix.")
  }

  if (verbose) {
    scice_message("GRAPH_TO_IGRAPH: Starting sparse matrix to igraph conversion...")
    scice_message(paste("GRAPH_TO_IGRAPH: Thread context - PID:", Sys.getpid()))
    scice_message(paste("GRAPH_TO_IGRAPH: Matrix class:", paste(class(seurat_graph), collapse = ", ")))
    scice_message(paste("GRAPH_TO_IGRAPH: Matrix dimensions:", nrow(seurat_graph), "x", ncol(seurat_graph)))
    scice_message(paste("GRAPH_TO_IGRAPH: Matrix storage mode:", typeof(seurat_graph)))
    scice_message("GRAPH_TO_IGRAPH: Long step started - extracting sparse triplets...")
    matrix_start <- Sys.time()
  }

  triplets <- Matrix::summary(seurat_graph)
  n_edges <- nrow(triplets)
  weights <- triplets$x
  edges_matrix <- as.matrix(triplets[, c("i", "j"), drop = FALSE])

  if (verbose) {
    matrix_time <- as.numeric(difftime(Sys.time(), matrix_start, units = "secs"))
    scice_message(paste("GRAPH_TO_IGRAPH: Extraction completed in", round(matrix_time, 3), "seconds"))
    scice_message(paste("GRAPH_TO_IGRAPH: Found", length(weights), "non-zero entries"))
    scice_message("GRAPH_TO_IGRAPH: Long step finished - sparse triplet extraction complete.")
    total_entries <- as.numeric(nrow(seurat_graph)) * as.numeric(ncol(seurat_graph))
    sparsity <- if (total_entries > 0) (1 - length(weights) / total_entries) * 100 else NA_real_
    scice_message(paste("GRAPH_TO_IGRAPH: Matrix sparsity:", round(sparsity, 2), "%"))
    if (length(weights) > 0) {
      scice_message(paste("GRAPH_TO_IGRAPH: Weight range: [", round(min(weights), 4), ", ", round(max(weights), 4), "]", sep = ""))
      scice_message(paste("GRAPH_TO_IGRAPH: Mean weight:", round(mean(weights), 4)))
      scice_message(paste("GRAPH_TO_IGRAPH: Median weight:", round(median(weights), 4)))
    }
    scice_message(paste("GRAPH_TO_IGRAPH: Created edge list with", n_edges, "edges"))
    if (n_edges > 0) {
      scice_message(paste("GRAPH_TO_IGRAPH: Vertex index range: [", min(edges_matrix), ", ", max(edges_matrix), "]", sep = ""))
      scice_message("GRAPH_TO_IGRAPH: Sample edges (from -> to):")
      head_edges <- utils::head(triplets, 5)
      for (i in seq_len(nrow(head_edges))) {
        scice_message(
          paste(
            "GRAPH_TO_IGRAPH:   ", head_edges$i[i], " -> ", head_edges$j[i],
            " (weight: ", round(head_edges$x[i], 4), ")",
            sep = ""
          )
        )
      }
      if (n_edges > 5) {
        scice_message(paste("GRAPH_TO_IGRAPH:   ... and", n_edges - 5, "more edges"))
      }
    }
  }

  n_vertices <- nrow(seurat_graph)
  if (verbose) {
    scice_message(paste("GRAPH_TO_IGRAPH: Creating empty graph with", n_vertices, "vertices"))
    graph_creation_start <- Sys.time()
  }

  igraph_obj <- igraph::make_empty_graph(n = n_vertices, directed = FALSE)

  if (n_edges > 0) {
    if (verbose) {
      scice_message(paste("GRAPH_TO_IGRAPH: Long step started - constructing igraph edge structure from", n_edges, "edges..."))
      scice_message("GRAPH_TO_IGRAPH: Adding edges and weights to graph...")
    }

    igraph_obj <- igraph::add_edges(igraph_obj, as.vector(t(edges_matrix)))
    igraph::E(igraph_obj)$weight <- weights

    if (verbose) {
      graph_creation_time <- as.numeric(difftime(Sys.time(), graph_creation_start, units = "secs"))
      scice_message(paste("GRAPH_TO_IGRAPH: Graph creation completed in", round(graph_creation_time, 3), "seconds"))
      scice_message("GRAPH_TO_IGRAPH: Long step finished - igraph edge construction complete.")
      scice_message(paste("GRAPH_TO_IGRAPH: Final graph: vertices =", igraph::vcount(igraph_obj), ", edges =", igraph::ecount(igraph_obj)))
      scice_message(paste("GRAPH_TO_IGRAPH: Graph is directed:", igraph::is_directed(igraph_obj)))
      scice_message(paste("GRAPH_TO_IGRAPH: Graph is weighted:", igraph::is_weighted(igraph_obj)))
      edge_weights <- igraph::E(igraph_obj)$weight
      scice_message("GRAPH_TO_IGRAPH: Edge weights summary:")
      scice_message(paste("GRAPH_TO_IGRAPH:   Min:", round(min(edge_weights), 4)))
      scice_message(paste("GRAPH_TO_IGRAPH:   Max:", round(max(edge_weights), 4)))
      scice_message(paste("GRAPH_TO_IGRAPH:   Mean:", round(mean(edge_weights), 4)))
      scice_message(paste("GRAPH_TO_IGRAPH:   Median:", round(median(edge_weights), 4)))
    }
  } else if (verbose) {
    scice_message("GRAPH_TO_IGRAPH: WARNING - No edges to add! Graph will be empty.")
  }

  igraph_obj
}
