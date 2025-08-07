#' @import igraph
#' @importFrom stats median
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom parallel detectCores mclapply
#' @importFrom data.table data.table rbindlist
NULL

# Global clustering cache environment for optimization
clustering_cache_env <- new.env(parent = emptyenv())

#' Cached version of leiden clustering to avoid redundant computations
#' @param igraph_obj igraph object to cluster
#' @param resolution Resolution parameter for clustering
#' @param objective_function Objective function ("modularity" or "CPM")
#' @param n_iterations Number of iterations
#' @param beta Beta parameter
#' @param initial_membership Initial cluster membership (optional)
#' @param use_cache Whether to use caching (default: TRUE)
#' @param cache_key_suffix Additional suffix for cache key (default: "")
#' @return Vector of cluster assignments (0-based)
#' @keywords internal
cached_leiden_clustering <- function(igraph_obj, resolution, objective_function, 
                                   n_iterations, beta, initial_membership = NULL,
                                   use_cache = TRUE, cache_key_suffix = "") {
  
  if (!use_cache) {
    return(leiden_clustering(igraph_obj, resolution, objective_function, n_iterations, beta, initial_membership))
  }
  
  # Create cache key based on parameters (excluding initial_membership for broader reuse)
  cache_key <- paste(
    "r", round(resolution, 8),
    "obj", objective_function, 
    "iter", n_iterations,
    "beta", round(beta, 4),
    "suffix", cache_key_suffix,
    sep = "_"
  )
  
  # Check if result exists in cache
  if (exists(cache_key, envir = clustering_cache_env)) {
    return(clustering_cache_env[[cache_key]])
  }
  
  # Compute result and cache it
  result <- leiden_clustering(igraph_obj, resolution, objective_function, n_iterations, beta, initial_membership)
  clustering_cache_env[[cache_key]] <- result
  
  return(result)
}

#' Clear the global clustering cache
#' @keywords internal
clear_clustering_cache <- function() {
  rm(list = ls(envir = clustering_cache_env), envir = clustering_cache_env)
}

#' Get cache statistics
#' @keywords internal  
get_cache_stats <- function() {
  cache_size <- length(ls(envir = clustering_cache_env))
  return(list(cache_entries = cache_size))
}

#' Cross-platform parallel lapply wrapper
#' @param X Vector/list to iterate over
#' @param FUN Function to apply
#' @param mc.cores Number of cores (ignored on Windows)
#' @param ... Additional arguments to FUN
#' @return List of results
#' @keywords internal
cross_platform_mclapply <- function(X, FUN, mc.cores = 1, ...) {
  if (.Platform$OS.type == "windows") {
    # Use regular lapply on Windows
    return(lapply(X, FUN, ...))
  } else {
    # Use mclapply on Unix-like systems
    return(parallel::mclapply(X, FUN, mc.cores = mc.cores, ...))
  }
}

#' Core clustering algorithm implementing binary search and optimization
#'
#' @param igraph_obj igraph object to cluster
#' @param cluster_range Vector of cluster numbers to test
#' @param n_workers Number of parallel workers (default: max(1, parallel::detectCores() - 1))
#' @param n_trials Number of clustering trials per resolution
#' @param n_bootstrap Number of bootstrap iterations
#' @param beta Beta parameter for Leiden clustering
#' @param n_iterations Number of Leiden iterations
#' @param max_iterations Maximum iterations for optimization
#' @param objective_function Objective function ("modularity" or "CPM")
#' @param remove_threshold Threshold for removing inconsistent results
#' @param resolution_tolerance Tolerance for resolution parameter search
#' @param verbose Whether to print progress messages
#' @param in_parallel_context Whether this function is called from within a parallel context (default: FALSE)
#' @return List with clustering results
#' @keywords internal
clustering_main <- function(igraph_obj, cluster_range, n_workers = max(1, parallel::detectCores() - 1), 
                          n_trials, n_bootstrap, seed = NULL, beta, n_iterations, max_iterations, 
                          objective_function, remove_threshold, resolution_tolerance, verbose, 
                          in_parallel_context = FALSE) {
  
  # Clear clustering cache at start for fresh run
  if (verbose) {
    cache_stats_before <- get_cache_stats()
    message(paste("CLUSTERING_MAIN: Cache entries before clearing:", cache_stats_before$cache_entries))
  }
  clear_clustering_cache()
  
  # Set random seed if provided for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    if (verbose) {
      message(paste("CLUSTERING_MAIN: Setting random seed to:", seed))
      message(paste("CLUSTERING_MAIN: Thread context - PID:", Sys.getpid()))
    }
  }
  
  # Initialize results storage using data.table for better performance
  results_dt <- data.table::data.table(
    cluster_number = integer(),
    gamma = numeric(),
    labels = list(),
    ic = numeric(),
    ic_vec = list(),
    best_labels = list(),
    n_iter = integer(),
    mei = list(),
    k = integer(),
    excluded = logical(),
    exclusion_reason = character()
  )
  
  # Determine resolution search bounds
  if (objective_function == "modularity") {
    start_g <- -13
    end_g <- 20  # Increased for higher cluster numbers
  } else { # CPM
    start_g <- log(resolution_tolerance)
    if (start_g < -20) start_g <- -20
    end_g <- 20  # Increased for higher cluster numbers
  }
  
  # Binary search for resolution ranges
  if (verbose) {
    message("CLUSTERING_MAIN: Starting binary search for resolution ranges...")
    message(paste("CLUSTERING_MAIN: Objective function:", objective_function))
    message(paste("CLUSTERING_MAIN: Search bounds: [", start_g, ", ", end_g, "]", sep = ""))
    message(paste("CLUSTERING_MAIN: Target cluster range:", paste(cluster_range, collapse = ", ")))
    resolution_search_start <- Sys.time()
  }
  
  gamma_dict <- find_resolution_ranges(
    igraph_obj, cluster_range, start_g, end_g, objective_function,
    resolution_tolerance, n_workers, verbose, seed, in_parallel_context
  )
  
  if (verbose) {
    resolution_search_time <- as.numeric(difftime(Sys.time(), resolution_search_start, units = "secs"))
    message(paste("CLUSTERING_MAIN: Resolution search completed in", round(resolution_search_time, 3), "seconds"))
    message(paste("CLUSTERING_MAIN: Found resolution ranges for", length(gamma_dict), "cluster numbers"))
    if (length(gamma_dict) > 0) {
      for (i in 1:min(5, length(gamma_dict))) {
        cluster_num <- names(gamma_dict)[i]
        range_vals <- gamma_dict[[cluster_num]]
        message(paste("CLUSTERING_MAIN:   k=", cluster_num, ": gamma ∈ [", 
                     round(range_vals[1], 4), ", ", round(range_vals[2], 4), "]", sep = ""))
      }
      if (length(gamma_dict) > 5) {
        message(paste("CLUSTERING_MAIN:   ... and", length(gamma_dict) - 5, "more"))
      }
    }
  }
  
  # Filter out problematic cluster numbers in parallel
  if (verbose) {
    message("CLUSTERING_MAIN: Starting problematic cluster filtering...")
    message(paste("CLUSTERING_MAIN: Remove threshold:", remove_threshold))
    message(paste("CLUSTERING_MAIN: Platform:", .Platform$OS.type))
    filtering_start <- Sys.time()
  }
  
  # Windows compatibility for parallel processing
  actual_workers <- if (.Platform$OS.type == "windows" && n_workers > 1) 1 else n_workers
  
  if (verbose) {
    message(paste("CLUSTERING_MAIN: Using", actual_workers, "workers for filtering (requested:", n_workers, ")"))
    if (in_parallel_context) {
      message("CLUSTERING_MAIN: Running in parallel context - using 1 worker for nested operations")
    }
  }
  
  cluster_filter_results <- cross_platform_mclapply(cluster_range, function(cluster_num) {
    if (!(as.character(cluster_num) %in% names(gamma_dict))) {
      return(list(excluded = TRUE, cluster_num = cluster_num, reason = "resolution_search_failed"))
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    gamma_test <- seq(gamma_range[1], gamma_range[2], length.out = min(5, abs(diff(gamma_range)) * 100 + 1))
    
    # Test multiple gammas in parallel (distribute workers efficiently when in parallel context)
    nested_workers <- if (in_parallel_context) {
      max(1, as.integer(round(actual_workers / length(cluster_range))))
    } else {
      actual_workers
    }
    
    if (verbose && in_parallel_context) {
      message(paste("CLUSTERING_MAIN: Nested worker optimization - using", nested_workers, 
                   "workers per cluster ( ", actual_workers, "total /", length(cluster_range), "clusters)"))
    }
            ic_scores <- cross_platform_mclapply(gamma_test, function(gamma_val) {
      # Set deterministic seed for filtering if base seed provided
      if (!is.null(seed)) {
        filter_seed <- seed + cluster_num * 100 + as.integer(gamma_val * 1000) %% 10000
        set.seed(filter_seed)
      }
      
      cluster_results <- replicate(10, {
        cached_leiden_clustering(igraph_obj, gamma_val, objective_function, 5, 0.01,
                               cache_key_suffix = paste("filter", cluster_num, sep = "_"))
      }, simplify = TRUE)
      
      extracted_results <- extract_clustering_array(cluster_results)
      ic_result <- calculate_ic_from_extracted(extracted_results)
      return(1 / ic_result)
          }, mc.cores = nested_workers)
      
      ic_scores <- unlist(ic_scores)
      excluded <- min(ic_scores, na.rm = TRUE) >= remove_threshold
    reason <- if (excluded) "high_inconsistency" else "passed_filtering"
    
    return(list(excluded = excluded, cluster_num = cluster_num, reason = reason))
  }, mc.cores = actual_workers)
  
  excluded_numbers <- sapply(cluster_filter_results, function(x) if(x$excluded) x$cluster_num else NULL)
  excluded_numbers <- unlist(excluded_numbers)
  
  if (verbose) {
    filtering_time <- as.numeric(difftime(Sys.time(), filtering_start, units = "secs"))
    message(paste("CLUSTERING_MAIN: Filtering completed in", round(filtering_time, 3), "seconds"))
    message(paste("CLUSTERING_MAIN: Processed", length(cluster_filter_results), "cluster numbers"))
    if (length(excluded_numbers) > 0) {
      message(paste("CLUSTERING_MAIN: Excluded", length(excluded_numbers), "cluster numbers:", paste(excluded_numbers, collapse = ", ")))
      # Report exclusion reasons
      exclusion_reasons <- sapply(cluster_filter_results, function(x) if(x$excluded) paste(x$cluster_num, ":", x$reason) else NULL)
      exclusion_reasons <- unlist(exclusion_reasons)
      if (length(exclusion_reasons) > 0) {
        message("CLUSTERING_MAIN: Exclusion reasons:")
        for (reason in exclusion_reasons) {
          message(paste("CLUSTERING_MAIN:   ", reason))
        }
      }
    } else {
      message("CLUSTERING_MAIN: No clusters excluded during filtering")
    }
  }
  
  valid_clusters <- setdiff(cluster_range, excluded_numbers)
  
  if (verbose) {
    message(paste("CLUSTERING_MAIN: Valid clusters for optimization:", length(valid_clusters)))
    if (length(valid_clusters) > 0) {
      message(paste("CLUSTERING_MAIN: Valid cluster numbers:", paste(valid_clusters, collapse = ", ")))
    }
  }
  
  # Create entries for excluded clusters
  excluded_entries <- lapply(cluster_filter_results, function(x) {
    if (x$excluded) {
      return(data.table::data.table(
        cluster_number = x$cluster_num,
        gamma = NA_real_,
        labels = list(NULL),
        ic = NA_real_,
        ic_vec = list(NULL),
        best_labels = list(NULL),
        n_iter = NA_integer_,
        mei = list(NULL),
        k = NA_integer_,
        excluded = TRUE,
        exclusion_reason = x$reason
      ))
    }
    return(NULL)
  })
  
  # Handle excluded entries properly
  excluded_entries_list <- excluded_entries[!sapply(excluded_entries, is.null)]
  if (length(excluded_entries_list) > 0) {
    excluded_entries <- data.table::rbindlist(excluded_entries_list)
  } else {
    # Create empty data.table with proper structure
    excluded_entries <- data.table::data.table(
      cluster_number = integer(),
      gamma = numeric(),
      labels = list(),
      ic = numeric(),
      ic_vec = list(),
      best_labels = list(),
      n_iter = integer(),
      mei = list(),
      k = integer(),
      excluded = logical(),
      exclusion_reason = character()
    )
  }
  
  if (length(valid_clusters) == 0) {
    if (verbose) {
      message("CLUSTERING_MAIN: WARNING - No valid cluster numbers found!")
      message("CLUSTERING_MAIN: All clusters were excluded during filtering")
      message("CLUSTERING_MAIN: This will result in empty IC results")
      message("CLUSTERING_MAIN: Returning results with only excluded cluster information")
    }
    # Return results with only excluded clusters
    return(list(
      gamma = excluded_entries$gamma,
      labels = excluded_entries$labels,
      ic = excluded_entries$ic,
      ic_vec = excluded_entries$ic_vec,
      n_cluster = excluded_entries$cluster_number,
      best_labels = excluded_entries$best_labels,
      n_iter = excluded_entries$n_iter,
      mei = excluded_entries$mei,
      k = excluded_entries$k,
      excluded = excluded_entries$excluded,
      exclusion_reason = excluded_entries$exclusion_reason
    ))
  }
  
  # Optimize clustering for each valid cluster number in parallel
  if (verbose) {
    message("CLUSTERING_MAIN: Starting clustering optimization...")
    message(paste("CLUSTERING_MAIN: Optimizing", length(valid_clusters), "cluster numbers"))
    optimization_start <- Sys.time()
  }
  
  # Windows compatibility for clustering optimization
  actual_workers_opt <- if (.Platform$OS.type == "windows" && n_workers > 1) 1 else n_workers
  
  if (verbose) {
    message(paste("CLUSTERING_MAIN: Using", actual_workers_opt, "workers for optimization"))
    message(paste("CLUSTERING_MAIN: Progress tracking:"))
    progress_start_time <- Sys.time()
    
    # Pre-calculate estimated time per cluster for progress tracking
    if (length(valid_clusters) > 1) {
      estimated_time_per_cluster <- (n_trials * length(valid_clusters) * n_bootstrap) / (actual_workers_opt * 100)  # rough estimate
      total_estimated_time <- estimated_time_per_cluster * length(valid_clusters)
      message(paste("CLUSTERING_MAIN: Estimated total time:", round(total_estimated_time, 1), "seconds"))
    }
  }
  
  cluster_results <- cross_platform_mclapply(valid_clusters, function(cluster_num) {
    # Thread-local logging for parallel workers
    worker_id <- paste("WORKER", cluster_num)
    
    if (verbose) {
      message(paste(worker_id, ": Starting optimization for k =", cluster_num))
      message(paste(worker_id, ": Thread context - PID:", Sys.getpid()))
    }
    
    if (!(as.character(cluster_num) %in% names(gamma_dict))) {
      if (verbose) {
        message(paste(worker_id, ": ERROR - No gamma range found for k =", cluster_num))
      }
      return(NULL)
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    
    if (verbose) {
      message(paste(worker_id, ": Gamma range [", round(gamma_range[1], 4), ", ", round(gamma_range[2], 4), "]", sep = ""))
      message(paste(worker_id, ": Starting intensive optimization..."))
      opt_start_time <- Sys.time()
    }
    
    # Optimize clustering within this range
    result <- optimize_clustering(
      igraph_obj, cluster_num, gamma_range, objective_function,
      n_trials, n_bootstrap, seed, beta, n_iterations, max_iterations,
      resolution_tolerance, n_workers, verbose, worker_id, in_parallel_context = TRUE
    )
    
    if (verbose) {
      opt_time <- as.numeric(difftime(Sys.time(), opt_start_time, units = "secs"))
      message(paste(worker_id, ": Optimization completed in", round(opt_time, 3), "seconds"))
    }
    
    if (!is.null(result)) {
      # Calculate MEI scores
      mei_scores <- calculate_mei_from_array(result$labels)
      
      # Return as a data.table row
      return(data.table::data.table(
        cluster_number = cluster_num,
        gamma = result$gamma,
        labels = list(result$labels),
        ic = result$ic_median,
        ic_vec = list(result$ic_bootstrap),
        best_labels = list(result$best_labels),
        n_iter = result$n_iterations,
        mei = list(mei_scores),
        k = result$k,
        excluded = FALSE,
        exclusion_reason = "none"
      ))
    }
    return(NULL)
  }, mc.cores = actual_workers_opt)
  
  if (verbose) {
    optimization_time <- as.numeric(difftime(Sys.time(), optimization_start, units = "secs"))
    message(paste("CLUSTERING_MAIN: All optimization workers completed in", round(optimization_time, 3), "seconds"))
    successful_count <- sum(!sapply(cluster_results, is.null))
    message(paste("CLUSTERING_MAIN: Successful optimizations:", successful_count, "/", length(valid_clusters)))
  }
  
  # Combine successful results
  successful_results_list <- cluster_results[!sapply(cluster_results, is.null)]
  if (length(successful_results_list) > 0) {
    successful_results <- data.table::rbindlist(successful_results_list, fill = TRUE)
  } else {
    # Create empty data.table with proper structure
    successful_results <- data.table::data.table(
      cluster_number = integer(),
      gamma = numeric(),
      labels = list(),
      ic = numeric(),
      ic_vec = list(),
      best_labels = list(),
      n_iter = integer(),
      mei = list(),
      k = integer(),
      excluded = logical(),
      exclusion_reason = character()
    )
  }
  
  # Combine all results (excluded + successful) with robust error handling
  tryCatch({
    if (nrow(excluded_entries) > 0 && nrow(successful_results) > 0) {
      # Ensure both data.tables have compatible column structures
      if (all(colnames(excluded_entries) %in% colnames(successful_results)) &&
          all(colnames(successful_results) %in% colnames(excluded_entries))) {
        results_dt <- data.table::rbindlist(list(excluded_entries, successful_results), fill = TRUE)
      } else {
        # If columns don't match, combine manually
        all_cols <- union(colnames(excluded_entries), colnames(successful_results))
        for (col in all_cols) {
          if (!col %in% colnames(excluded_entries)) {
            excluded_entries[[col]] <- switch(col,
              "cluster_number" = integer(nrow(excluded_entries)),
              "gamma" = numeric(nrow(excluded_entries)),
              "labels" = rep(list(NULL), nrow(excluded_entries)),
              "ic" = numeric(nrow(excluded_entries)),
              "ic_vec" = rep(list(NULL), nrow(excluded_entries)),
              "best_labels" = rep(list(NULL), nrow(excluded_entries)),
              "n_iter" = integer(nrow(excluded_entries)),
              "mei" = rep(list(NULL), nrow(excluded_entries)),
              "k" = integer(nrow(excluded_entries)),
              "excluded" = logical(nrow(excluded_entries)),
              "exclusion_reason" = character(nrow(excluded_entries)))
          }
          if (!col %in% colnames(successful_results)) {
            successful_results[[col]] <- switch(col,
              "cluster_number" = integer(nrow(successful_results)),
              "gamma" = numeric(nrow(successful_results)),
              "labels" = rep(list(NULL), nrow(successful_results)),
              "ic" = numeric(nrow(successful_results)),
              "ic_vec" = rep(list(NULL), nrow(successful_results)),
              "best_labels" = rep(list(NULL), nrow(successful_results)),
              "n_iter" = integer(nrow(successful_results)),
              "mei" = rep(list(NULL), nrow(successful_results)),
              "k" = integer(nrow(successful_results)),
              "excluded" = logical(nrow(successful_results)),
              "exclusion_reason" = character(nrow(successful_results)))
          }
        }
        results_dt <- data.table::rbindlist(list(excluded_entries, successful_results), fill = TRUE)
      }
    } else if (nrow(excluded_entries) > 0) {
      results_dt <- excluded_entries
    } else if (nrow(successful_results) > 0) {
      results_dt <- successful_results
    } else {
      # Both are empty - create empty data.table with proper structure
      results_dt <- data.table::data.table(
        cluster_number = integer(),
        gamma = numeric(),
        labels = list(),
        ic = numeric(),
        ic_vec = list(),
        best_labels = list(),
        n_iter = integer(),
        mei = list(),
        k = integer(),
        excluded = logical(),
        exclusion_reason = character()
      )
    }
  }, error = function(e) {
    warning("Error combining results, creating empty results: ", e$message)
    results_dt <<- data.table::data.table(
      cluster_number = integer(),
      gamma = numeric(),
      labels = list(),
      ic = numeric(),
      ic_vec = list(),
      best_labels = list(),
      n_iter = integer(),
      mei = list(),
      k = integer(),
      excluded = logical(),
      exclusion_reason = character()
    )
     })
  
  # Sort by cluster number (only if the column exists and there are rows)
  if (nrow(results_dt) > 0 && "cluster_number" %in% colnames(results_dt)) {
    data.table::setorder(results_dt, cluster_number)
  }
  
  # Report cache statistics for performance monitoring
  if (verbose) {
    cache_stats_final <- get_cache_stats()
    message(paste("CLUSTERING_MAIN: Final cache entries:", cache_stats_final$cache_entries))
    message(paste("CLUSTERING_MAIN: Cache provided", cache_stats_final$cache_entries, "reused clustering results"))
  }
  
  # Convert back to list format for compatibility - with defensive programming
  return(list(
    gamma = if ("gamma" %in% colnames(results_dt)) results_dt$gamma else numeric(0),
    labels = if ("labels" %in% colnames(results_dt)) results_dt$labels else list(),
    ic = if ("ic" %in% colnames(results_dt)) results_dt$ic else numeric(0),
    ic_vec = if ("ic_vec" %in% colnames(results_dt)) results_dt$ic_vec else list(),
    n_cluster = if ("cluster_number" %in% colnames(results_dt)) results_dt$cluster_number else integer(0),
    best_labels = if ("best_labels" %in% colnames(results_dt)) results_dt$best_labels else list(),
    n_iter = if ("n_iter" %in% colnames(results_dt)) results_dt$n_iter else integer(0),
    mei = if ("mei" %in% colnames(results_dt)) results_dt$mei else list(),
    k = if ("k" %in% colnames(results_dt)) results_dt$k else integer(0),
    excluded = if ("excluded" %in% colnames(results_dt)) results_dt$excluded else logical(0),
    exclusion_reason = if ("exclusion_reason" %in% colnames(results_dt)) results_dt$exclusion_reason else character(0)
  ))
}

#' Find resolution parameter ranges for each cluster number using binary search
#' @keywords internal
find_resolution_ranges <- function(igraph_obj, cluster_range, start_g, end_g,
                                  objective_function, resolution_tolerance, n_workers, verbose, seed = NULL, 
                                  in_parallel_context = FALSE) {
  
  # Initialize results storage with data.table
  results_dt <- data.table::data.table(
    cluster_number = integer(),
    left_bound = numeric(),
    right_bound = numeric()
  )
  
  n_preliminary_trials <- 15  # Increased for better precision
  beta_preliminary <- 0.01
  n_iter_preliminary <- 5   # Increased for more stable clustering
  
  # Process cluster numbers in parallel (with Windows compatibility)
  if (.Platform$OS.type == "windows" && n_workers > 1) {
    n_workers <- 1  # Force single worker on Windows
  }
  
  # Use distributed workers if already in parallel context to maximize efficiency
  nested_workers <- if (in_parallel_context) {
    max(1, as.integer(round(n_workers / length(cluster_range))))
  } else {
    n_workers
  }
  
  if (verbose && in_parallel_context) {
    message(paste("RESOLUTION_SEARCH: Nested worker optimization - using", nested_workers, 
                 "workers per cluster ( ", n_workers, "total /", length(cluster_range), "clusters)"))
  }
  
  range_results <- cross_platform_mclapply(cluster_range, function(target_clusters) {
    left <- start_g
    right <- end_g
    max_iterations <- 50  # Limit binary search iterations
    iteration_count <- 0
    
    # Use tighter tolerance for better precision in distinguishing adjacent cluster numbers
    effective_tolerance <- resolution_tolerance / 10
    
    # Binary search for lower bound with fallback
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > effective_tolerance && iteration_count < max_iterations) {
      iteration_count <- iteration_count + 1
      mid <- (left + right) / 2
      gamma_val <- if (objective_function == "modularity") mid else exp(mid)
      
      # Test clustering with current gamma using vectorized operations
      # Set deterministic seed for lower bound search if base seed provided
      if (!is.null(seed)) {
        range_seed <- seed + target_clusters * 10 + iteration_count
        set.seed(range_seed)
      }
      
      cluster_results <- replicate(n_preliminary_trials, {
        cached_leiden_clustering(igraph_obj, gamma_val, objective_function, n_iter_preliminary, beta_preliminary,
                               cache_key_suffix = paste("res_search_lower", target_clusters, sep = "_"))
      }, simplify = TRUE)
      
      # Use apply for efficient calculation - count unique clusters instead of max + 1
      n_clusters_obtained <- stats::median(apply(cluster_results, 2, function(x) length(unique(x))))
      
      if (n_clusters_obtained < target_clusters) {
        left <- mid
      } else {
        right <- mid
      }
    }
    
    left_bound <- right
    
    # Binary search for upper bound with fallback
    left <- left_bound
    right <- end_g
    iteration_count <- 0
    
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > effective_tolerance && iteration_count < max_iterations) {
      iteration_count <- iteration_count + 1
      mid <- (left + right) / 2
      gamma_val <- if (objective_function == "modularity") mid else exp(mid)
      
      # Test clustering with current gamma using vectorized operations
      # Set deterministic seed for upper bound search if base seed provided
      if (!is.null(seed)) {
        range_seed <- seed + target_clusters * 10 + 100 + iteration_count  # +100 to distinguish from lower bound
        set.seed(range_seed)
      }
      
      cluster_results <- replicate(n_preliminary_trials, {
        cached_leiden_clustering(igraph_obj, gamma_val, objective_function, n_iter_preliminary, beta_preliminary,
                               cache_key_suffix = paste("res_search_upper", target_clusters, sep = "_"))
      }, simplify = TRUE)
      
      # Use apply for efficient calculation - count unique clusters instead of max + 1
      n_clusters_obtained <- stats::median(apply(cluster_results, 2, function(x) length(unique(x))))
      
      if (n_clusters_obtained > target_clusters) {
        right <- mid
      } else {
        left <- mid
      }
    }
    
    right_bound <- left
    
    # If bounds are identical or invalid, create a cluster-specific range
    if (identical(left_bound, right_bound) || is.na(left_bound) || is.na(right_bound)) {
      if (objective_function == "CPM") {
        # For CPM, create a range around the found value with cluster-specific adjustment
        center_val <- ifelse(is.na(left_bound), exp((start_g + end_g) / 2), left_bound)
        # Add cluster-specific offset to ensure different ranges for different cluster numbers
        cluster_offset <- (target_clusters - min(cluster_range)) * 0.05  # Small offset based on cluster number
        adjusted_center <- center_val * (1 + cluster_offset)
        left_bound <- max(exp(start_g), adjusted_center * 0.7)
        right_bound <- min(exp(end_g), adjusted_center * 1.3)
      } else {
        # For modularity, create a range around the found value with cluster-specific adjustment
        center_val <- ifelse(is.na(left_bound), (start_g + end_g) / 2, left_bound)
        # Add cluster-specific offset to ensure different ranges for different cluster numbers
        cluster_offset <- (target_clusters - min(cluster_range)) * 0.02  # Small offset based on cluster number
        adjusted_center <- center_val + cluster_offset
        left_bound <- max(start_g, adjusted_center - 0.15)
        right_bound <- min(end_g, adjusted_center + 0.15)
      }
    }
    
    # Return results as a data.table row
    data.table::data.table(
      cluster_number = target_clusters,
      left_bound = if (objective_function == "CPM") exp(left_bound) else left_bound,
      right_bound = if (objective_function == "CPM") exp(right_bound) else right_bound
    )
  }, mc.cores = nested_workers)
  
  # Combine results
  # Filter out NULL results from parallel processing
  valid_results <- range_results[!sapply(range_results, is.null)]
  
  if (length(valid_results) == 0) {
    # If no valid results, return empty gamma_dict
    warning("No valid resolution ranges found for any cluster numbers")
    return(setNames(list(), character(0)))
  }
  
  results_dt <- data.table::rbindlist(valid_results)
  
  # Check if required columns exist
  if (!"left_bound" %in% colnames(results_dt) || !"right_bound" %in% colnames(results_dt)) {
    stop("Missing required columns in resolution range results")
  }
  
  # Sort by cluster number to ensure proper ordering
  data.table::setorder(results_dt, cluster_number)
  
  # Post-process to ensure non-overlapping ranges for adjacent cluster numbers
  if (nrow(results_dt) > 1) {
    for (i in 2:nrow(results_dt)) {
      current_left <- results_dt$left_bound[i]
      current_right <- results_dt$right_bound[i]
      prev_left <- results_dt$left_bound[i-1]
      prev_right <- results_dt$right_bound[i-1]
      
      # Check for overlap or identical ranges
      if (current_left <= prev_right && abs(current_left - prev_right) < resolution_tolerance * 10) {
        if (verbose) {
          message(paste("RESOLUTION_SEARCH: Adjusting overlapping ranges for clusters", 
                       results_dt$cluster_number[i-1], "and", results_dt$cluster_number[i]))
        }
        
        # Calculate a small gap between ranges
        if (objective_function == "CPM") {
          gap <- max(resolution_tolerance * 5, (prev_right - prev_left) * 0.1)
          # Adjust current range to be above previous range
          range_width <- current_right - current_left
          results_dt$left_bound[i] <- prev_right + gap
          results_dt$right_bound[i] <- results_dt$left_bound[i] + range_width
        } else {
          gap <- max(resolution_tolerance * 10, (prev_right - prev_left) * 0.1)
          # Adjust current range to be above previous range
          range_width <- current_right - current_left
          results_dt$left_bound[i] <- prev_right + gap
          results_dt$right_bound[i] <- results_dt$left_bound[i] + range_width
        }
      }
    }
  }
  
  # Convert to named list format - fixing the data.table syntax
  gamma_dict <- lapply(seq_len(nrow(results_dt)), function(i) {
    c(results_dt$left_bound[i], results_dt$right_bound[i])
  })
  names(gamma_dict) <- as.character(results_dt$cluster_number)
  
  return(gamma_dict)
}

#' Optimize clustering within a resolution range
#' @keywords internal
optimize_clustering <- function(igraph_obj, target_clusters, gamma_range, objective_function,
                               n_trials, n_bootstrap, seed = NULL, beta, n_iterations, max_iterations,
                               resolution_tolerance, n_workers, verbose = FALSE, worker_id = "OPTIMIZER", 
                               in_parallel_context = FALSE) {
  
  # Set deterministic seeds for this cluster number if base seed provided
  if (!is.null(seed)) {
    cluster_seed <- seed + target_clusters * 1000  # Different seed per cluster number
    set.seed(cluster_seed)
    if (verbose) {
      message(paste(worker_id, ": Set deterministic seed:", cluster_seed))
    }
  }
  
  if (verbose) {
    message(paste(worker_id, ": Optimization parameters:"))
    message(paste(worker_id, ":   Target clusters:", target_clusters))
    message(paste(worker_id, ":   Trials per gamma:", n_trials))
    message(paste(worker_id, ":   Bootstrap iterations:", n_bootstrap))
    message(paste(worker_id, ":   Max iterations:", max_iterations))
    message(paste(worker_id, ":   Beta:", beta))
    message(paste(worker_id, ":   Leiden iterations:", n_iterations))
  }
  
  n_steps <- 11
  delta_n <- 2
  
  if (verbose) {
    message(paste(worker_id, ": Creating gamma sequence with", n_steps, "steps"))
    message(paste(worker_id, ": Gamma range bounds: [", round(gamma_range[1], 4), ", ", round(gamma_range[2], 4), "]", sep = ""))
  }
  
  # Create gamma sequence using efficient sequence generation
  gamma_sequence <- if (objective_function == "modularity") {
    if (!identical(gamma_range[1], gamma_range[2])) {
      seq(gamma_range[1], gamma_range[2], length.out = n_steps)
    } else {
      delta_g <- resolution_tolerance
      seq(gamma_range[1] - delta_g, gamma_range[1] + delta_g, length.out = n_steps)
    }
  } else { # CPM
    if (!identical(gamma_range[1], gamma_range[2])) {
      exp(seq(log(gamma_range[1]), log(gamma_range[2]), length.out = n_steps))
    } else {
      delta_g <- resolution_tolerance
      exp(seq(log(gamma_range[1]) - delta_g, log(gamma_range[1]) + delta_g, length.out = n_steps))
    }
  }
  
  if (verbose) {
    message(paste(worker_id, ": Generated gamma sequence:"))
    for (i in 1:min(5, length(gamma_sequence))) {
      message(paste(worker_id, ":   γ[", i, "] = ", round(gamma_sequence[i], 4), sep = ""))
    }
    if (length(gamma_sequence) > 5) {
      message(paste(worker_id, ":   ... and", length(gamma_sequence) - 5, "more"))
    }
  }
  
  # Test initial clustering for each gamma in parallel
  if (verbose) {
    message(paste(worker_id, ": Phase 1 - Testing", length(gamma_sequence), "gamma values with", n_trials, "trials each"))
    initial_phase_start <- Sys.time()
    if (in_parallel_context) {
      message(paste(worker_id, ": Running in parallel context - using 1 worker for nested operations"))
    }
  }
  
  # Use distributed workers if already in parallel context to maximize efficiency
  nested_workers <- if (in_parallel_context) {
    max(1, as.integer(round(n_workers / length(gamma_sequence))))
  } else {
    n_workers
  }
  
  if (verbose && in_parallel_context) {
    message(paste(worker_id, ": Nested worker optimization - using", nested_workers, 
                 "workers per gamma ( ", n_workers, "total /", length(gamma_sequence), "gammas)"))
  }
  
  clustering_results <- cross_platform_mclapply(gamma_sequence, function(gamma_val) {
    # Set deterministic seed for this gamma if base seed provided
    if (!is.null(seed)) {
      gamma_seed <- cluster_seed + as.integer(gamma_val * 10000) %% 100000
      set.seed(gamma_seed)
    }
    
    cluster_matrix <- replicate(n_trials, {
      leiden_clustering(igraph_obj, gamma_val, objective_function, n_iterations, beta)
    }, simplify = TRUE)
    
    mean_clusters <- stats::median(apply(cluster_matrix, 2, function(x) length(unique(x))))
    
    # Mini progress report for verbose mode
    if (verbose && length(gamma_sequence) <= 5) {  # Only for small sequences to avoid spam
      message(paste(worker_id, ":   γ =", round(gamma_val, 4), "→ median clusters =", mean_clusters))
    }
    
    list(
      matrix = cluster_matrix,
      mean_clusters = mean_clusters
    )
  }, mc.cores = nested_workers)
  
  if (verbose) {
    initial_phase_time <- as.numeric(difftime(Sys.time(), initial_phase_start, units = "secs"))
    message(paste(worker_id, ": Phase 1 completed in", round(initial_phase_time, 3), "seconds"))
  }
  
  # Extract results
  mean_clusters <- sapply(clustering_results, function(x) x$mean_clusters)
  clustering_matrices <- lapply(clustering_results, function(x) x$matrix)
  
  if (verbose) {
    message(paste(worker_id, ": Phase 2 - Filtering for target cluster count:", target_clusters))
    cluster_counts <- table(mean_clusters)
    for (i in 1:length(cluster_counts)) {
      count_val <- names(cluster_counts)[i]
      freq <- cluster_counts[i]
      message(paste(worker_id, ":   ", freq, "gammas → ", count_val, " clusters", sep = ""))
    }
  }
  
  # Filter for target cluster number
  valid_indices <- which(mean_clusters == target_clusters)
  
  if (length(valid_indices) == 0) {
    if (verbose) {
      message(paste(worker_id, ": ERROR - No gammas produced target cluster count", target_clusters))
    }
    return(NULL)
  }
  
  if (verbose) {
    message(paste(worker_id, ": Found", length(valid_indices), "gammas producing", target_clusters, "clusters"))
  }
  
  gamma_sequence <- gamma_sequence[valid_indices]
  clustering_matrices <- clustering_matrices[valid_indices]
  
  # Calculate IC for each gamma in parallel
  if (verbose) {
    message(paste(worker_id, ": Phase 3 - Calculating IC scores for", length(gamma_sequence), "valid gammas"))
    ic_phase_start <- Sys.time()
  }
  
  ic_scores <- cross_platform_mclapply(clustering_matrices, function(cluster_matrix) {
    extracted <- extract_clustering_array(cluster_matrix)
    ic_result <- calculate_ic_from_extracted(extracted)
    return(1 / ic_result)  # Convert to inconsistency score
  }, mc.cores = nested_workers)
  
  ic_scores <- unlist(ic_scores)
  
  if (verbose) {
    ic_phase_time <- as.numeric(difftime(Sys.time(), ic_phase_start, units = "secs"))
    message(paste(worker_id, ": IC calculation completed in", round(ic_phase_time, 3), "seconds"))
    message(paste(worker_id, ": IC score range: [", round(min(ic_scores), 4), ", ", round(max(ic_scores), 4), "]", sep = ""))
    perfect_scores <- sum(ic_scores == 1)
    if (perfect_scores > 0) {
      message(paste(worker_id, ": Found", perfect_scores, "perfect IC scores (= 1.0)"))
    }
  }
  
  # Find the best gamma
  best_index <- which(ic_scores == 1)[1]
  if (is.na(best_index)) {
    best_index <- which.min(ic_scores)
  }
  
  best_gamma <- gamma_sequence[best_index]
  best_clustering <- clustering_matrices[[best_index]]
  k <- n_iterations
  
  if (verbose) {
    message(paste(worker_id, ": Best gamma:", round(best_gamma, 4), "with IC score:", round(ic_scores[best_index], 4)))
  }
  
  # If no perfect IC found, iterate to improve
  if (ic_scores[best_index] != 1 && length(gamma_sequence) > 1) {
    if (verbose) {
      message(paste(worker_id, ": Phase 4 - Iterative improvement (IC =", round(ic_scores[best_index], 4), "< 1.0)"))
      message(paste(worker_id, ": Starting with", length(gamma_sequence), "gamma values"))
      iterative_start <- Sys.time()
    }
    
    current_results <- clustering_matrices
    current_gammas <- gamma_sequence
    current_ic <- ic_scores
    
    # Track stability using a matrix for efficiency
    ic_history <- matrix(rep(current_ic, 10), nrow = length(current_ic))
    
    iteration_count <- 0
    while (k < max_iterations) {
      k <- k + delta_n
      iteration_count <- iteration_count + 1
      
      if (verbose) {
        message(paste(worker_id, ": Iteration", iteration_count, "(k =", k, ") - Refining", length(current_gammas), "gamma values"))
        iter_start <- Sys.time()
      }
      
      # Update clustering results in parallel
      new_results <- cross_platform_mclapply(seq_along(current_gammas), function(i) {
        gamma_val <- current_gammas[i]
        current_matrix <- current_results[[i]]
        
        # Use previous results as initialization
        # Set deterministic seed for this iteration if base seed provided
        if (!is.null(seed)) {
          iter_seed <- cluster_seed + k * 100 + i
          set.seed(iter_seed)
        }
        
        new_clustering <- replicate(n_trials, {
          init_membership <- current_matrix[, sample.int(ncol(current_matrix), 1)]
          leiden_clustering(igraph_obj, gamma_val, objective_function, delta_n, beta, init_membership)
        }, simplify = TRUE)
        
        # Calculate IC score
        extracted <- extract_clustering_array(new_clustering)
        ic_result <- calculate_ic_from_extracted(extracted)
        
        list(
          matrix = new_clustering,
          ic = 1 / ic_result
        )
      }, mc.cores = nested_workers)
      
      # Extract results
      new_matrices <- lapply(new_results, function(x) x$matrix)
      new_ic <- sapply(new_results, function(x) x$ic)
      
      if (verbose) {
        iter_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
        message(paste(worker_id, ": Iteration", iteration_count, "completed in", round(iter_time, 3), "seconds"))
        message(paste(worker_id, ": New IC range: [", round(min(new_ic), 4), ", ", round(max(new_ic), 4), "]", sep = ""))
      }
      
      # Update IC history
      ic_history <- cbind(ic_history[, -1], new_ic)
      
      # Check for stability and convergence
      stable_indices <- apply(ic_history, 1, function(row) length(unique(row)) == 1)
      perfect_indices <- which(new_ic == 1)
      
      if (verbose) {
        message(paste(worker_id, ": Stable gammas:", sum(stable_indices), "/", length(stable_indices)))
        message(paste(worker_id, ": Perfect IC scores:", length(perfect_indices)))
      }
      
      # Selection criteria
      if (length(perfect_indices) > 0) {
        best_index <- perfect_indices[1]
        best_gamma <- current_gammas[best_index]
        best_clustering <- new_matrices[[best_index]]
        if (verbose) {
          message(paste(worker_id, ": CONVERGED - Found perfect IC score at gamma =", round(best_gamma, 4)))
        }
        break
      } else if (all(stable_indices)) {
        best_index <- which.min(new_ic)
        best_gamma <- current_gammas[best_index]
        best_clustering <- new_matrices[[best_index]]
        if (verbose) {
          message(paste(worker_id, ": CONVERGED - All gammas stable, best IC =", round(new_ic[best_index], 4)))
        }
        break
      } else {
        # Continue with best performing gammas
        keep_indices <- (new_ic <= stats::quantile(new_ic, 0.5)) | stable_indices
        keep_indices[which.min(new_ic)] <- TRUE
        
        if (sum(keep_indices) == 1) {
          best_index <- which(keep_indices)
          best_gamma <- current_gammas[best_index]
          best_clustering <- new_matrices[[best_index]]
          if (verbose) {
            message(paste(worker_id, ": CONVERGED - Single gamma remaining, IC =", round(new_ic[best_index], 4)))
          }
          break
        }
        
        if (verbose) {
          message(paste(worker_id, ": Continuing with", sum(keep_indices), "best gammas"))
        }
        
        current_gammas <- current_gammas[keep_indices]
        current_results <- new_matrices[keep_indices]
        ic_history <- ic_history[keep_indices, , drop = FALSE]
        current_ic <- new_ic[keep_indices]
      }
    }
    
    if (verbose) {
      iterative_time <- as.numeric(difftime(Sys.time(), iterative_start, units = "secs"))
      message(paste(worker_id, ": Phase 4 completed in", round(iterative_time, 3), "seconds after", iteration_count, "iterations"))
    }
  }
  
  # Bootstrap analysis in parallel
  if (verbose) {
    message(paste(worker_id, ": Phase 5 - Bootstrap analysis with", n_bootstrap, "iterations"))
    bootstrap_start <- Sys.time()
  }
  
  ic_bootstrap <- cross_platform_mclapply(seq_len(n_bootstrap), function(i) {
    # Set deterministic seed for this bootstrap iteration if base seed provided
    if (!is.null(seed)) {
      bootstrap_seed <- cluster_seed + 10000 + i
      set.seed(bootstrap_seed)
    }
    
    sample_indices <- sample.int(ncol(best_clustering), ncol(best_clustering), replace = TRUE)
    bootstrap_matrix <- best_clustering[, sample_indices, drop = FALSE]
    extracted <- extract_clustering_array(bootstrap_matrix)
    ic_result <- calculate_ic_from_extracted(extracted)
    return(1 / ic_result)
  }, mc.cores = nested_workers)
  
  ic_bootstrap <- unlist(ic_bootstrap)
  ic_median <- stats::median(ic_bootstrap)
  
  if (verbose) {
    bootstrap_time <- as.numeric(difftime(Sys.time(), bootstrap_start, units = "secs"))
    message(paste(worker_id, ": Bootstrap analysis completed in", round(bootstrap_time, 3), "seconds"))
    message(paste(worker_id, ": Bootstrap IC median:", round(ic_median, 4)))
    message(paste(worker_id, ": Bootstrap IC range: [", round(min(ic_bootstrap), 4), ", ", round(max(ic_bootstrap), 4), "]", sep = ""))
    message(paste(worker_id, ": OPTIMIZATION COMPLETE"))
  }
  
  # Extract best labels
  extracted_best <- extract_clustering_array(best_clustering)
  best_labels <- get_best_clustering(extracted_best)
  
  return(list(
    gamma = best_gamma,
    labels = extracted_best,
    ic_median = ic_median,
    ic_bootstrap = ic_bootstrap,
    best_labels = best_labels,
    n_iterations = k,
    k = k
  ))
} 