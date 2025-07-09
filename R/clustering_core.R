#' @import igraph
#' @importFrom stats median
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom parallel detectCores mclapply
#' @importFrom data.table data.table rbindlist
NULL

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
#' @return List with clustering results
#' @keywords internal
clustering_main <- function(igraph_obj, cluster_range, n_workers = max(1, parallel::detectCores() - 1), 
                          n_trials, n_bootstrap, beta, n_iterations, max_iterations, 
                          objective_function, remove_threshold, resolution_tolerance, verbose) {
  
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
    start_g <- 0
    end_g <- 20  # Increased for higher cluster numbers
  } else { # CPM
    start_g <- log(resolution_tolerance)
    if (start_g < -13) start_g <- -13
    end_g <- 2  # Increased for higher cluster numbers
  }
  
  # Binary search for resolution ranges
  if (verbose) message("Performing binary search for resolution ranges...")
  
  gamma_dict <- find_resolution_ranges(
    igraph_obj, cluster_range, start_g, end_g, objective_function,
    resolution_tolerance, n_workers, verbose
  )
  
  # Filter out problematic cluster numbers in parallel
  if (verbose) message("Filtering problematic clusters...")
  
  # Windows compatibility for parallel processing
  actual_workers <- if (.Platform$OS.type == "windows" && n_workers > 1) 1 else n_workers
  
  cluster_filter_results <- cross_platform_mclapply(cluster_range, function(cluster_num) {
    if (!(as.character(cluster_num) %in% names(gamma_dict))) {
      return(list(excluded = TRUE, cluster_num = cluster_num, reason = "resolution_search_failed"))
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    gamma_test <- seq(gamma_range[1], gamma_range[2], length.out = min(5, diff(gamma_range) * 100 + 1))
    
    # Test multiple gammas in parallel
            ic_scores <- cross_platform_mclapply(gamma_test, function(gamma_val) {
      cluster_results <- replicate(10, {
        leiden_clustering(igraph_obj, gamma_val, objective_function, 5, 0.01)
      }, simplify = TRUE)
      
      extracted_results <- extract_clustering_array(cluster_results)
      ic_result <- calculate_ic_from_extracted(extracted_results)
      return(1 / ic_result)
          }, mc.cores = actual_workers)
      
      ic_scores <- unlist(ic_scores)
      excluded <- min(ic_scores, na.rm = TRUE) >= remove_threshold
    reason <- if (excluded) "high_inconsistency" else "passed_filtering"
    
    return(list(excluded = excluded, cluster_num = cluster_num, reason = reason))
  }, mc.cores = actual_workers)
  
  excluded_numbers <- sapply(cluster_filter_results, function(x) if(x$excluded) x$cluster_num else NULL)
  excluded_numbers <- unlist(excluded_numbers)
  
  if (verbose && length(excluded_numbers) > 0) {
    message(paste("Excluded cluster numbers:", paste(excluded_numbers, collapse = ", ")))
  }
  
  valid_clusters <- setdiff(cluster_range, excluded_numbers)
  
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
    if (verbose) message("No valid cluster numbers found")
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
  if (verbose) message("Optimizing clustering for each valid cluster number...")
  
  # Windows compatibility for clustering optimization
  actual_workers_opt <- if (.Platform$OS.type == "windows" && n_workers > 1) 1 else n_workers
  
  cluster_results <- cross_platform_mclapply(valid_clusters, function(cluster_num) {
    if (!(as.character(cluster_num) %in% names(gamma_dict))) {
      return(NULL)
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    
    # Optimize clustering within this range
    result <- optimize_clustering(
      igraph_obj, cluster_num, gamma_range, objective_function,
      n_trials, n_bootstrap, beta, n_iterations, max_iterations,
      resolution_tolerance, n_workers
    )
    
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
    results_dt <- results_dt[order(results_dt$cluster_number)]
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
                                  objective_function, resolution_tolerance, n_workers, verbose) {
  
  # Initialize results storage with data.table
  results_dt <- data.table::data.table(
    cluster_number = integer(),
    left_bound = numeric(),
    right_bound = numeric()
  )
  
  n_preliminary_trials <- 10
  beta_preliminary <- 0.01
  n_iter_preliminary <- 3
  
  # Process cluster numbers in parallel (with Windows compatibility)
  if (.Platform$OS.type == "windows" && n_workers > 1) {
    n_workers <- 1  # Force single worker on Windows
  }
  
  range_results <- cross_platform_mclapply(cluster_range, function(target_clusters) {
    left <- start_g
    right <- end_g
    max_iterations <- 50  # Limit binary search iterations
    iteration_count <- 0
    
    # Binary search for lower bound with fallback
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > resolution_tolerance && iteration_count < max_iterations) {
      iteration_count <- iteration_count + 1
      mid <- (left + right) / 2
      gamma_val <- if (objective_function == "modularity") mid else exp(mid)
      
      # Test clustering with current gamma using vectorized operations
      cluster_results <- replicate(n_preliminary_trials, {
        leiden_clustering(igraph_obj, gamma_val, objective_function, n_iter_preliminary, beta_preliminary)
      }, simplify = TRUE)
      
      # Use apply for efficient calculation
      n_clusters_obtained <- stats::median(apply(cluster_results, 2, max) + 1)
      
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
    
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > resolution_tolerance && iteration_count < max_iterations) {
      iteration_count <- iteration_count + 1
      mid <- (left + right) / 2
      gamma_val <- if (objective_function == "modularity") mid else exp(mid)
      
      # Test clustering with current gamma using vectorized operations
      cluster_results <- replicate(n_preliminary_trials, {
        leiden_clustering(igraph_obj, gamma_val, objective_function, n_iter_preliminary, beta_preliminary)
      }, simplify = TRUE)
      
      # Use apply for efficient calculation
      n_clusters_obtained <- stats::median(apply(cluster_results, 2, max) + 1)
      
      if (n_clusters_obtained > target_clusters) {
        right <- mid
      } else {
        left <- mid
      }
    }
    
    right_bound <- left
    
    # If bounds are identical or invalid, create a small range around the found value
    if (identical(left_bound, right_bound) || is.na(left_bound) || is.na(right_bound)) {
      if (objective_function == "CPM") {
        # For CPM, create a range around the found value
        center_val <- ifelse(is.na(left_bound), exp((start_g + end_g) / 2), left_bound)
        left_bound <- max(exp(start_g), center_val * 0.8)
        right_bound <- min(exp(end_g), center_val * 1.2)
      } else {
        # For modularity, create a range around the found value
        center_val <- ifelse(is.na(left_bound), (start_g + end_g) / 2, left_bound)
        left_bound <- max(start_g, center_val - 0.1)
        right_bound <- min(end_g, center_val + 0.1)
      }
    }
    
    # Return results as a data.table row
    data.table::data.table(
      cluster_number = target_clusters,
      left_bound = if (objective_function == "CPM") exp(left_bound) else left_bound,
      right_bound = if (objective_function == "CPM") exp(right_bound) else right_bound
    )
  }, mc.cores = n_workers)
  
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
  
  # Convert to named list format - fixing the data.table syntax
  gamma_dict <- lapply(seq_len(nrow(results_dt)), function(i) {
    c(results_dt$left_bound[i], results_dt$right_bound[i])
  })
  names(gamma_dict) <- as.character(results_dt$cluster_number)
  
  return(gamma_dict)
}

#' Filter out cluster numbers that produce inconsistent results
#' @keywords internal
filter_problematic_clusters <- function(gamma_dict, cluster_range, igraph_obj,
                                       objective_function, remove_threshold, n_workers, verbose) {
  
  excluded_numbers <- c()
  
  for (cluster_num in cluster_range) {
    if (!(as.character(cluster_num) %in% names(gamma_dict))) {
      excluded_numbers <- c(excluded_numbers, cluster_num)
      next
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    
    # Test multiple gammas in the range
    gamma_test <- seq(gamma_range[1], gamma_range[2], length.out = min(5, diff(gamma_range) * 100 + 1))
    
    ic_scores <- sapply(gamma_test, function(gamma_val) {
      cluster_results <- replicate(10, {
        leiden_clustering(igraph_obj, gamma_val, objective_function, 5, 0.01)
      }, simplify = TRUE)
      
      extracted_results <- extract_clustering_array(cluster_results)
      ic_result <- calculate_ic_from_extracted(extracted_results)
      return(1 / ic_result)  # Convert to inconsistency score
    })
    
    if (min(ic_scores, na.rm = TRUE) >= remove_threshold) {
      excluded_numbers <- c(excluded_numbers, cluster_num)
    }
  }
  
  return(excluded_numbers)
}

#' Optimize clustering within a resolution range
#' @keywords internal
optimize_clustering <- function(igraph_obj, target_clusters, gamma_range, objective_function,
                               n_trials, n_bootstrap, beta, n_iterations, max_iterations,
                               resolution_tolerance, n_workers) {
  
  n_steps <- 11
  delta_n <- 2
  
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
  
  # Test initial clustering for each gamma in parallel
  clustering_results <- cross_platform_mclapply(gamma_sequence, function(gamma_val) {
    cluster_matrix <- replicate(n_trials, {
      leiden_clustering(igraph_obj, gamma_val, objective_function, n_iterations, beta)
    }, simplify = TRUE)
    
    list(
      matrix = cluster_matrix,
      mean_clusters = stats::median(apply(cluster_matrix, 2, max) + 1)
    )
  }, mc.cores = n_workers)
  
  # Extract results
  mean_clusters <- sapply(clustering_results, function(x) x$mean_clusters)
  clustering_matrices <- lapply(clustering_results, function(x) x$matrix)
  
  # Filter for target cluster number
  valid_indices <- which(mean_clusters == target_clusters)
  
  if (length(valid_indices) == 0) {
    return(NULL)
  }
  
  gamma_sequence <- gamma_sequence[valid_indices]
  clustering_matrices <- clustering_matrices[valid_indices]
  
  # Calculate IC for each gamma in parallel
  ic_scores <- cross_platform_mclapply(clustering_matrices, function(cluster_matrix) {
    extracted <- extract_clustering_array(cluster_matrix)
    ic_result <- calculate_ic_from_extracted(extracted)
    return(1 / ic_result)  # Convert to inconsistency score
  }, mc.cores = n_workers)
  
  ic_scores <- unlist(ic_scores)
  
  # Find the best gamma
  best_index <- which(ic_scores == 1)[1]
  if (is.na(best_index)) {
    best_index <- which.min(ic_scores)
  }
  
  best_gamma <- gamma_sequence[best_index]
  best_clustering <- clustering_matrices[[best_index]]
  k <- n_iterations
  
  # If no perfect IC found, iterate to improve
  if (ic_scores[best_index] != 1 && length(gamma_sequence) > 1) {
    current_results <- clustering_matrices
    current_gammas <- gamma_sequence
    current_ic <- ic_scores
    
    # Track stability using a matrix for efficiency
    ic_history <- matrix(rep(current_ic, 10), nrow = length(current_ic))
    
    while (k < max_iterations) {
      k <- k + delta_n
      
      # Update clustering results in parallel
      new_results <- cross_platform_mclapply(seq_along(current_gammas), function(i) {
        gamma_val <- current_gammas[i]
        current_matrix <- current_results[[i]]
        
        # Use previous results as initialization
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
      }, mc.cores = n_workers)
      
      # Extract results
      new_matrices <- lapply(new_results, function(x) x$matrix)
      new_ic <- sapply(new_results, function(x) x$ic)
      
      # Update IC history
      ic_history <- cbind(ic_history[, -1], new_ic)
      
      # Check for stability and convergence
      stable_indices <- apply(ic_history, 1, function(row) length(unique(row)) == 1)
      perfect_indices <- which(new_ic == 1)
      
      # Selection criteria
      if (length(perfect_indices) > 0) {
        best_index <- perfect_indices[1]
        best_gamma <- current_gammas[best_index]
        best_clustering <- new_matrices[[best_index]]
        break
      } else if (all(stable_indices)) {
        best_index <- which.min(new_ic)
        best_gamma <- current_gammas[best_index]
        best_clustering <- new_matrices[[best_index]]
        break
      } else {
        # Continue with best performing gammas
        keep_indices <- (new_ic <= stats::quantile(new_ic, 0.5)) | stable_indices
        keep_indices[which.min(new_ic)] <- TRUE
        
        if (sum(keep_indices) == 1) {
          best_index <- which(keep_indices)
          best_gamma <- current_gammas[best_index]
          best_clustering <- new_matrices[[best_index]]
          break
        }
        
        current_gammas <- current_gammas[keep_indices]
        current_results <- new_matrices[keep_indices]
        ic_history <- ic_history[keep_indices, , drop = FALSE]
        current_ic <- new_ic[keep_indices]
      }
    }
  }
  
  # Bootstrap analysis in parallel
  ic_bootstrap <- cross_platform_mclapply(seq_len(n_bootstrap), function(i) {
    sample_indices <- sample.int(ncol(best_clustering), ncol(best_clustering), replace = TRUE)
    bootstrap_matrix <- best_clustering[, sample_indices, drop = FALSE]
    extracted <- extract_clustering_array(bootstrap_matrix)
    ic_result <- calculate_ic_from_extracted(extracted)
    return(1 / ic_result)
  }, mc.cores = n_workers)
  
  ic_bootstrap <- unlist(ic_bootstrap)
  ic_median <- stats::median(ic_bootstrap)
  
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