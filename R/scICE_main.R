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
#' @param ic_threshold IC threshold for consistent clustering (default: Inf)
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
                            remove_threshold = 1.15,
                            resolution_tolerance = 1e-8,
                            verbose = TRUE) {
  
  # Clear clustering cache at start of new scICE run for optimal performance
  if (exists("clear_clustering_cache") && is.function(clear_clustering_cache)) {
    clear_clustering_cache()
  }
  
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
    scice_message(paste(rep("=", 80), collapse = ""))
    scice_message("Starting scICE clustering analysis...")
    scice_message(paste("Timestamp:", format(start_time, "%Y-%m-%d %H:%M:%S")))
    scice_message(paste("R Session Info: R", R.version.string))
    scice_message(paste("Platform:", R.version$platform))
    scice_message(paste("OS Type:", .Platform$OS.type))
    scice_message(paste("Process ID:", Sys.getpid()))
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("INPUT PARAMETERS:")
    scice_message(paste("  Using graph:", graph_name))
    scice_message(paste("  Testing cluster range:", paste(cluster_range, collapse = ", ")))
    scice_message(paste("  Range: ", min(cluster_range), "-", max(cluster_range), " (", length(cluster_range), " values)"))
    scice_message(paste("  Number of workers:", n_workers))
    scice_message(paste("  Number of trials per resolution:", n_trials))
    scice_message(paste("  Number of bootstrap iterations:", n_bootstrap))
    scice_message(paste("  Random seed:", if(is.null(seed)) "NULL (random)" else seed))
    scice_message(paste("  Beta parameter:", beta))
    scice_message(paste("  Leiden iterations:", n_iterations))
    scice_message(paste("  Maximum optimization iterations:", max_iterations))
    scice_message(paste("  IC threshold:", ic_threshold))
    scice_message(paste("  Objective function:", objective_function))
    scice_message(paste("  Remove threshold:", remove_threshold))
    scice_message(paste("  Resolution tolerance:", resolution_tolerance))
    scice_message(paste(rep("-", 80), collapse = ""))
  }
  
  # Set up parallel processing
  if (verbose) {
    scice_message("PARALLEL PROCESSING SETUP:")
    scice_message(paste("  Detected cores:", detectCores()))
    scice_message(paste("  Requested workers:", n_workers))
    scice_message(paste("  Already registered:", getDoParRegistered()))
  }
  
  if (n_workers > 1) {
    if (!getDoParRegistered()) {
      actual_workers <- min(n_workers, detectCores() - 1)
      if (verbose) {
        scice_message(paste("  Creating cluster with", actual_workers, "workers"))
        scice_message(paste("  Thread setup: Creating parallel cluster"))
        scice_message(paste("  Backend: doParallel with", actual_workers, "cores"))
      }
      cl <- parallel::makeCluster(actual_workers)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      if (verbose) {
        scice_message("  Parallel cluster created successfully")
      }
    } else {
      if (verbose) {
        scice_message("  Using existing parallel backend")
      }
    }
  } else {
    if (verbose) {
      scice_message("  Running in sequential mode (n_workers = 1)")
    }
  }
  
  # Get graph from Seurat object
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("GRAPH EXTRACTION:")
    scice_message(paste("  Accessing graph:", graph_name))
    scice_message(paste("  Available graphs in object:", paste(names(object@graphs), collapse = ", ")))
  }
  
  graph <- object@graphs[[graph_name]]
  
  if (verbose) {
    scice_message(paste("  Graph extraction successful"))
    scice_message(paste("  Graph class:", paste(class(graph), collapse = ", ")))
    scice_message(paste("  Graph dimensions:", nrow(graph), "x", ncol(graph)))
    scice_message(paste("  Graph storage type:", typeof(graph)))
    if (inherits(graph, "dgCMatrix") || inherits(graph, "Matrix")) {
      n_nonzero <- sum(graph > 0)
      scice_message(paste("  Non-zero entries:", n_nonzero))
      scice_message(paste("  Sparsity:", round((1 - n_nonzero / (nrow(graph) * ncol(graph))) * 100, 2), "%"))
      if (n_nonzero > 0) {
        weights <- as.vector(graph[graph > 0])
        scice_message(paste("  Weight range: [", round(min(weights), 4), ", ", round(max(weights), 4), "]", sep = ""))
        scice_message(paste("  Mean weight:", round(mean(weights), 4)))
      }
    }
  }
  
  # Convert to igraph object
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("GRAPH CONVERSION:")
    conversion_start <- Sys.time()
    scice_message(paste("  Starting graph conversion at:", format(conversion_start, "%H:%M:%S")))
  }
  
  igraph_obj <- graph_to_igraph(graph, verbose = verbose)
  
  if (verbose) {
    conversion_end <- Sys.time()
    conversion_time <- as.numeric(difftime(conversion_end, conversion_start, units = "secs"))
    scice_message(paste("  Graph conversion completed in:", round(conversion_time, 3), "seconds"))
    scice_message(paste("  Converted graph vertices:", igraph::vcount(igraph_obj)))
    scice_message(paste("  Converted graph edges:", igraph::ecount(igraph_obj)))
    scice_message(paste("  igraph object class:", paste(class(igraph_obj), collapse = ", ")))
    if (igraph::ecount(igraph_obj) > 0) {
      scice_message(paste("  Graph is weighted:", igraph::is_weighted(igraph_obj)))
      if (igraph::is_weighted(igraph_obj)) {
        edge_weights <- igraph::E(igraph_obj)$weight
        scice_message(paste("  Edge weight range: [", round(min(edge_weights), 4), ", ", round(max(edge_weights), 4), "]", sep = ""))
      }
    }
  }
  
  # Perform clustering analysis
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("CLUSTERING ANALYSIS:")
    clustering_start <- Sys.time()
    scice_message(paste("  Starting clustering analysis at:", format(clustering_start, "%H:%M:%S")))
    scice_message(paste("  Thread context: Main thread (PID:", Sys.getpid(), ")"))
    if (n_workers > 1) {
      scice_message(paste("  Parallel workers will be spawned for sub-tasks"))
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
    verbose = verbose,
    in_parallel_context = FALSE
  )
  
  if (verbose) {
    clustering_end <- Sys.time()
    clustering_time <- as.numeric(difftime(clustering_end, clustering_start, units = "secs"))
    scice_message(paste("  Clustering analysis completed in:", round(clustering_time, 3), "seconds"))
    scice_message(paste("  Results structure:"))
    scice_message(paste("    - gamma length:", length(results$gamma)))
    scice_message(paste("    - labels length:", length(results$labels)))
    scice_message(paste("    - ic length:", length(results$ic)))
    scice_message(paste("    - n_cluster length:", length(results$n_cluster)))
    if (!is.null(results$excluded)) {
      scice_message(paste("    - excluded info available:", sum(results$excluded), "excluded clusters"))
    }
  }
  
  # Add cell names to results
  results$cell_names <- Cells(object)
  
  # Store original cluster range for reference
  results$cluster_range_tested <- cluster_range
  
  # Determine consistent clusters using actual cluster numbers, not indices
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("RESULTS ANALYSIS:")
    scice_message(paste("  IC threshold for consistency:", ic_threshold))
    scice_message(paste("  Results object is null:", is.null(results)))
    scice_message(paste("  IC vector is null:", is.null(results$ic)))
    scice_message(paste("  IC vector length:", if(!is.null(results$ic)) length(results$ic) else 0))
  }
  
  if (!is.null(results) && !is.null(results$ic) && length(results$ic) > 0) {
    valid_ic <- results$ic[!is.na(results$ic)]
    if (verbose) {
      scice_message(paste("  Valid IC scores:", length(valid_ic)))
      if (length(valid_ic) > 0) {
        scice_message(paste("  IC score range: [", round(min(valid_ic), 4), ", ", round(max(valid_ic), 4), "]", sep = ""))
        scice_message(paste("  IC scores below threshold:", sum(valid_ic < ic_threshold)))
      }
    }
    
    consistent_indices <- which(results$ic < ic_threshold)
    results$consistent_clusters <- results$n_cluster[consistent_indices]
    
    if (verbose) {
      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      scice_message(paste(rep("=", 80), collapse = ""))
      scice_message("ANALYSIS COMPLETE")
      scice_message(paste("  Total execution time:", round(total_time, 3), "seconds"))
      scice_message(paste("  Found", length(results$consistent_clusters), 
                    "consistent cluster numbers:", paste(results$consistent_clusters, collapse = ", ")))
      
      # Report excluded clusters
      excluded_clusters <- setdiff(cluster_range, results$n_cluster)
      if (length(excluded_clusters) > 0) {
        scice_message(paste("  Excluded", length(excluded_clusters), "cluster numbers due to instability:", 
                      paste(excluded_clusters, collapse = ", ")))
      }
      
      # Report tested but inconsistent clusters
      inconsistent_clusters <- setdiff(results$n_cluster, results$consistent_clusters)
      if (length(inconsistent_clusters) > 0) {
        scice_message(paste("  Found", length(inconsistent_clusters), "inconsistent cluster numbers:", 
                      paste(inconsistent_clusters, collapse = ", ")))
      }
      
      # Detailed breakdown
      scice_message("  Detailed breakdown:")
      scice_message(paste("    - Requested clusters:", length(cluster_range)))
      scice_message(paste("    - Excluded clusters:", length(excluded_clusters)))
      scice_message(paste("    - Tested clusters:", length(results$n_cluster)))
      scice_message(paste("    - Consistent clusters:", length(results$consistent_clusters)))
      scice_message(paste("    - Inconsistent clusters:", length(inconsistent_clusters)))
    }
  } else {
    results$consistent_clusters <- integer(0)
    if (verbose) {
      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      scice_message(paste(rep("=", 80), collapse = ""))
      scice_message("ANALYSIS COMPLETE")
      scice_message(paste("  Total execution time:", round(total_time, 3), "seconds"))
      scice_message("  WARNING: No consistent cluster numbers found.")
      scice_message("  Possible reasons:")
      scice_message("    - All clusters were excluded during filtering")
      scice_message("    - IC threshold too strict")
      scice_message("    - Data quality issues")
      scice_message("    - Graph connectivity problems")
      scice_message("  Debugging suggestions:")
      scice_message("    - Try increasing ic_threshold (e.g., to Inf to see all results)")
      scice_message("    - Try increasing remove_threshold")
      scice_message("    - Check graph construction parameters")
      scice_message("    - Verify input data quality")
    }
  }
  
  # Add object reference for downstream analysis
  results$seurat_object <- object
  results$graph_name <- graph_name
  
  class(results) <- "scICE"
  return(results)
} 