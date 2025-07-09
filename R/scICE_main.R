#' @import Seurat
#' @import SeuratObject
#' @import igraph
#' @import parallel
#' @import foreach
#' @import doParallel
#' @importFrom methods inherits
#' @importFrom stats setNames
#' @importFrom utils getFromNamespace
#' @importFrom foreach getDoParRegistered
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom igraph make_empty_graph add_edges E V vcount ecount cluster_leiden membership graph_from_adjacency_matrix edge_attr_names
#'
#' Single-cell Inconsistency-based Clustering Evaluation (scICE)
#'
#' @description
#' scICE implements a systematic and efficient workflow to evaluate clustering consistency 
#' in single-cell RNA-seq data. It generates multiple cluster labels using the Leiden 
#' algorithm, calculates pairwise similarity with Element-Centric Similarity (ECS), 
#' and provides automated framework for finding consistent clusters.
#'
#' @param object A Seurat object containing single-cell data
#' @param graph_name Name of the graph to use for clustering. If NULL, will use the default SNN graph from the active assay (default: NULL)
#' @param cluster_range Vector of cluster numbers to test (default: 1:20)
#' @param n_workers Number of parallel workers to use (default: 10)
#' @param n_trials Number of clustering trials per resolution (default: 15)
#' @param n_bootstrap Number of bootstrap iterations (default: 100)
#' @param beta Beta parameter for Leiden clustering (default: 0.1)
#' @param n_iterations Number of Leiden iterations (default: 10)
#' @param max_iterations Maximum iterations for optimization (default: 150)
#' @param ic_threshold IC threshold for consistent clustering (default: 1.005)
#' @param objective_function Objective function for Leiden ("modularity" or "CPM", default: "CPM")
#' @param remove_threshold Threshold for removing inconsistent results (default: 1.15)
#' @param resolution_tolerance Tolerance for resolution parameter search (default: 1e-8)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \describe{
#'   \item{gamma}{Resolution parameters for each cluster number}
#'   \item{labels}{Clustering results for each cluster number}
#'   \item{ic}{Inconsistency scores for each cluster number}
#'   \item{ic_vec}{Bootstrap IC distributions}
#'   \item{n_cluster}{Number of clusters tested}
#'   \item{best_labels}{Best clustering labels for each cluster number}
#'   \item{n_iter}{Number of iterations used}
#'   \item{mei}{Mutual Element-wise Information scores}
#'   \item{consistent_clusters}{Cluster numbers meeting consistency threshold}
#' }
#'
#' @examples
#' \dontrun{
#' # Load Seurat object
#' data(pbmc_small)
#' 
#' # Run scICE analysis
#' scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)
#' 
#' # Plot IC results
#' plot_ic(scice_results)
#' 
#' # Extract consistent clustering labels
#' consistent_labels <- get_robust_labels(scice_results, threshold = 1.005)
#' }
#'
#' @export

scICE_clustering <- function(object, 
                            graph_name = NULL,
                            cluster_range = 1:20,
                            n_workers = 10,
                            n_trials = 15,
                            n_bootstrap = 100,
                            beta = 0.1,
                            n_iterations = 10,
                            max_iterations = 150,
                            ic_threshold = Inf,
                            objective_function = "CPM",
                            remove_threshold = Inf,
                            resolution_tolerance = 1e-8,
                            verbose = TRUE) {
  
  # Validate inputs
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Get default graph name if not provided
  if (is.null(graph_name)) {
    default_assay <- DefaultAssay(object)
    graph_name <- paste0(default_assay, "_snn")
  }
  
  # Check if graph exists in object
  if (!(graph_name %in% names(object@graphs))) {
    stop(sprintf("Graph '%s' not found in Seurat object. Available graphs: %s", 
                graph_name, 
                paste(names(object@graphs), collapse = ", ")))
  }
  
  if (verbose) {
    message("Starting scICE clustering analysis...")
    message(paste("Using graph:", graph_name))
    message(paste("Testing cluster range:", paste(range(cluster_range), collapse = "-")))
    message(paste("Using", n_workers, "parallel workers"))
  }
  
  # Set up parallel processing
  if (n_workers > 1) {
    if (!getDoParRegistered()) {
      cl <- parallel::makeCluster(min(n_workers, detectCores() - 1))
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
    }
  }
  
  # Get graph from Seurat object
  graph <- object@graphs[[graph_name]]
  
  if (verbose) {
    message("\nExtracting graph from Seurat object...")
    message(sprintf("Graph class: %s", class(graph)[1]))
    message(sprintf("Graph dimensions: %d x %d", nrow(graph), ncol(graph)))
  }
  
  # Convert to igraph object
  igraph_obj <- graph_to_igraph(graph, verbose = verbose)
  
  if (verbose) {
    message("\nGraph conversion complete")
    message(sprintf("Graph has %d vertices and %d edges", 
                   igraph::vcount(igraph_obj), 
                   igraph::ecount(igraph_obj)))
  }
  
  # Perform clustering analysis
  results <- clustering_main(
    igraph_obj = igraph_obj,
    cluster_range = cluster_range,
    n_workers = n_workers,
    n_trials = n_trials,
    n_bootstrap = n_bootstrap,
    beta = beta,
    n_iterations = n_iterations,
    max_iterations = max_iterations,
    objective_function = objective_function,
    remove_threshold = remove_threshold,
    resolution_tolerance = resolution_tolerance,
    verbose = verbose
  )
  
  # Add cell names to results
  results$cell_names <- Cells(object)
  
  # Store original cluster range for reference
  results$cluster_range_tested <- cluster_range
  
  # Determine consistent clusters using actual cluster numbers, not indices
  if (!is.null(results) && !is.null(results$ic) && length(results$ic) > 0) {
    consistent_indices <- which(results$ic < ic_threshold)
    results$consistent_clusters <- results$n_cluster[consistent_indices]
    
    if (verbose) {
      message(paste("Analysis complete. Found", length(results$consistent_clusters), 
                    "consistent cluster numbers:", paste(results$consistent_clusters, collapse = ", ")))
      
      # Report excluded clusters
      excluded_clusters <- setdiff(cluster_range, results$n_cluster)
      if (length(excluded_clusters) > 0) {
        message(paste("Excluded", length(excluded_clusters), "cluster numbers due to instability:", 
                      paste(excluded_clusters, collapse = ", ")))
      }
      
      # Report tested but inconsistent clusters
      inconsistent_clusters <- setdiff(results$n_cluster, results$consistent_clusters)
      if (length(inconsistent_clusters) > 0) {
        message(paste("Found", length(inconsistent_clusters), "inconsistent cluster numbers:", 
                      paste(inconsistent_clusters, collapse = ", ")))
      }
    }
  } else {
    results$consistent_clusters <- integer(0)
    if (verbose) {
      message("Analysis complete. No consistent cluster numbers found.")
    }
  }
  
  # Add object reference for downstream analysis
  results$seurat_object <- object
  results$graph_name <- graph_name
  
  class(results) <- "scICE"
  return(results)
} 