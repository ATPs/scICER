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
#' @param seed Random seed for reproducibility (default: NULL for random behavior)
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
                            seed = NULL,
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
    start_time <- Sys.time()
    message(paste(rep("=", 80), collapse = ""))
    message("Starting scICE clustering analysis...")
    message(paste("Timestamp:", format(start_time, "%Y-%m-%d %H:%M:%S")))
    message(paste("R Session Info: R", R.version.string))
    message(paste("Platform:", R.version$platform))
    message(paste("OS Type:", .Platform$OS.type))
    message(paste("Process ID:", Sys.getpid()))
    message(paste(rep("-", 80), collapse = ""))
    message("INPUT PARAMETERS:")
    message(paste("  Using graph:", graph_name))
    message(paste("  Testing cluster range:", paste(cluster_range, collapse = ", ")))
    message(paste("  Range: ", min(cluster_range), "-", max(cluster_range), " (", length(cluster_range), " values)"))
    message(paste("  Number of workers:", n_workers))
    message(paste("  Number of trials per resolution:", n_trials))
    message(paste("  Number of bootstrap iterations:", n_bootstrap))
    message(paste("  Random seed:", if(is.null(seed)) "NULL (random)" else seed))
    message(paste("  Beta parameter:", beta))
    message(paste("  Leiden iterations:", n_iterations))
    message(paste("  Maximum optimization iterations:", max_iterations))
    message(paste("  IC threshold:", ic_threshold))
    message(paste("  Objective function:", objective_function))
    message(paste("  Remove threshold:", remove_threshold))
    message(paste("  Resolution tolerance:", resolution_tolerance))
    message(paste(rep("-", 80), collapse = ""))
  }
  
  # Set up parallel processing
  if (verbose) {
    message("PARALLEL PROCESSING SETUP:")
    message(paste("  Detected cores:", detectCores()))
    message(paste("  Requested workers:", n_workers))
    message(paste("  Already registered:", getDoParRegistered()))
  }
  
  if (n_workers > 1) {
    if (!getDoParRegistered()) {
      actual_workers <- min(n_workers, detectCores() - 1)
      if (verbose) {
        message(paste("  Creating cluster with", actual_workers, "workers"))
        message(paste("  Thread setup: Creating parallel cluster"))
        message(paste("  Backend: doParallel with", actual_workers, "cores"))
      }
      cl <- parallel::makeCluster(actual_workers)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      if (verbose) {
        message("  Parallel cluster created successfully")
      }
    } else {
      if (verbose) {
        message("  Using existing parallel backend")
      }
    }
  } else {
    if (verbose) {
      message("  Running in sequential mode (n_workers = 1)")
    }
  }
  
  # Get graph from Seurat object
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("GRAPH EXTRACTION:")
    message(paste("  Accessing graph:", graph_name))
    message(paste("  Available graphs in object:", paste(names(object@graphs), collapse = ", ")))
  }
  
  graph <- object@graphs[[graph_name]]
  
  if (verbose) {
    message(paste("  Graph extraction successful"))
    message(paste("  Graph class:", paste(class(graph), collapse = ", ")))
    message(paste("  Graph dimensions:", nrow(graph), "x", ncol(graph)))
    message(paste("  Graph storage type:", typeof(graph)))
    if (inherits(graph, "dgCMatrix") || inherits(graph, "Matrix")) {
      n_nonzero <- sum(graph > 0)
      message(paste("  Non-zero entries:", n_nonzero))
      message(paste("  Sparsity:", round((1 - n_nonzero / (nrow(graph) * ncol(graph))) * 100, 2), "%"))
      if (n_nonzero > 0) {
        weights <- as.vector(graph[graph > 0])
        message(paste("  Weight range: [", round(min(weights), 4), ", ", round(max(weights), 4), "]", sep = ""))
        message(paste("  Mean weight:", round(mean(weights), 4)))
      }
    }
  }
  
  # Convert to igraph object
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("GRAPH CONVERSION:")
    conversion_start <- Sys.time()
    message(paste("  Starting graph conversion at:", format(conversion_start, "%H:%M:%S")))
  }
  
  igraph_obj <- graph_to_igraph(graph, verbose = verbose)
  
  if (verbose) {
    conversion_end <- Sys.time()
    conversion_time <- as.numeric(difftime(conversion_end, conversion_start, units = "secs"))
    message(paste("  Graph conversion completed in:", round(conversion_time, 3), "seconds"))
    message(paste("  Converted graph vertices:", igraph::vcount(igraph_obj)))
    message(paste("  Converted graph edges:", igraph::ecount(igraph_obj)))
    message(paste("  igraph object class:", paste(class(igraph_obj), collapse = ", ")))
    if (igraph::ecount(igraph_obj) > 0) {
      message(paste("  Graph is weighted:", igraph::is_weighted(igraph_obj)))
      if (igraph::is_weighted(igraph_obj)) {
        edge_weights <- igraph::E(igraph_obj)$weight
        message(paste("  Edge weight range: [", round(min(edge_weights), 4), ", ", round(max(edge_weights), 4), "]", sep = ""))
      }
    }
  }
  
  # Perform clustering analysis
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("CLUSTERING ANALYSIS:")
    clustering_start <- Sys.time()
    message(paste("  Starting clustering analysis at:", format(clustering_start, "%H:%M:%S")))
    message(paste("  Thread context: Main thread (PID:", Sys.getpid(), ")"))
    if (n_workers > 1) {
      message(paste("  Parallel workers will be spawned for sub-tasks"))
    }
  }
  
  results <- clustering_main(
    igraph_obj = igraph_obj,
    cluster_range = cluster_range,
    n_workers = n_workers,
    n_trials = n_trials,
    n_bootstrap = n_bootstrap,
    seed = seed,
    beta = beta,
    n_iterations = n_iterations,
    max_iterations = max_iterations,
    objective_function = objective_function,
    remove_threshold = remove_threshold,
    resolution_tolerance = resolution_tolerance,
    verbose = verbose
  )
  
  if (verbose) {
    clustering_end <- Sys.time()
    clustering_time <- as.numeric(difftime(clustering_end, clustering_start, units = "secs"))
    message(paste("  Clustering analysis completed in:", round(clustering_time, 3), "seconds"))
    message(paste("  Results structure:"))
    message(paste("    - gamma length:", length(results$gamma)))
    message(paste("    - labels length:", length(results$labels)))
    message(paste("    - ic length:", length(results$ic)))
    message(paste("    - n_cluster length:", length(results$n_cluster)))
    if (!is.null(results$excluded)) {
      message(paste("    - excluded info available:", sum(results$excluded), "excluded clusters"))
    }
  }
  
  # Add cell names to results
  results$cell_names <- Cells(object)
  
  # Store original cluster range for reference
  results$cluster_range_tested <- cluster_range
  
  # Determine consistent clusters using actual cluster numbers, not indices
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("RESULTS ANALYSIS:")
    message(paste("  IC threshold for consistency:", ic_threshold))
    message(paste("  Results object is null:", is.null(results)))
    message(paste("  IC vector is null:", is.null(results$ic)))
    message(paste("  IC vector length:", if(!is.null(results$ic)) length(results$ic) else 0))
  }
  
  if (!is.null(results) && !is.null(results$ic) && length(results$ic) > 0) {
    valid_ic <- results$ic[!is.na(results$ic)]
    if (verbose) {
      message(paste("  Valid IC scores:", length(valid_ic)))
      if (length(valid_ic) > 0) {
        message(paste("  IC score range: [", round(min(valid_ic), 4), ", ", round(max(valid_ic), 4), "]", sep = ""))
        message(paste("  IC scores below threshold:", sum(valid_ic < ic_threshold)))
      }
    }
    
    consistent_indices <- which(results$ic < ic_threshold)
    results$consistent_clusters <- results$n_cluster[consistent_indices]
    
    if (verbose) {
      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      message(paste(rep("=", 80), collapse = ""))
      message("ANALYSIS COMPLETE")
      message(paste("  Total execution time:", round(total_time, 3), "seconds"))
      message(paste("  Found", length(results$consistent_clusters), 
                    "consistent cluster numbers:", paste(results$consistent_clusters, collapse = ", ")))
      
      # Report excluded clusters
      excluded_clusters <- setdiff(cluster_range, results$n_cluster)
      if (length(excluded_clusters) > 0) {
        message(paste("  Excluded", length(excluded_clusters), "cluster numbers due to instability:", 
                      paste(excluded_clusters, collapse = ", ")))
      }
      
      # Report tested but inconsistent clusters
      inconsistent_clusters <- setdiff(results$n_cluster, results$consistent_clusters)
      if (length(inconsistent_clusters) > 0) {
        message(paste("  Found", length(inconsistent_clusters), "inconsistent cluster numbers:", 
                      paste(inconsistent_clusters, collapse = ", ")))
      }
      
      # Detailed breakdown
      message("  Detailed breakdown:")
      message(paste("    - Requested clusters:", length(cluster_range)))
      message(paste("    - Excluded clusters:", length(excluded_clusters)))
      message(paste("    - Tested clusters:", length(results$n_cluster)))
      message(paste("    - Consistent clusters:", length(results$consistent_clusters)))
      message(paste("    - Inconsistent clusters:", length(inconsistent_clusters)))
    }
  } else {
    results$consistent_clusters <- integer(0)
    if (verbose) {
      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      message(paste(rep("=", 80), collapse = ""))
      message("ANALYSIS COMPLETE")
      message(paste("  Total execution time:", round(total_time, 3), "seconds"))
      message("  WARNING: No consistent cluster numbers found.")
      message("  Possible reasons:")
      message("    - All clusters were excluded during filtering")
      message("    - IC threshold too strict")
      message("    - Data quality issues")
      message("    - Graph connectivity problems")
      message("  Debugging suggestions:")
      message("    - Try increasing ic_threshold (e.g., to Inf to see all results)")
      message("    - Try increasing remove_threshold")
      message("    - Check graph construction parameters")
      message("    - Verify input data quality")
    }
  }
  
  # Add object reference for downstream analysis
  results$seurat_object <- object
  results$graph_name <- graph_name
  
  class(results) <- "scICE"
  return(results)
} 