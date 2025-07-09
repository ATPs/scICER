#' @import igraph
#' @importFrom stats median
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom parallel detectCores mclapply
#' @importFrom data.table data.table rbindlist
NULL

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
    end_g <- 10
  } else { # CPM
    start_g <- log(resolution_tolerance)
    if (start_g < -13) start_g <- -13
    end_g <- 0
  }
  
  # Binary search for resolution ranges
  if (verbose) message("Performing binary search for resolution ranges...")
  
  gamma_dict <- find_resolution_ranges(
    igraph_obj, cluster_range, start_g, end_g, objective_function,
    resolution_tolerance, n_workers, verbose
  )
  
  # Filter out problematic cluster numbers in parallel
  if (verbose) message("Filtering problematic clusters...")
  
  cluster_filter_results <- parallel::mclapply(cluster_range, function(cluster_num) {
    if (!(as.character(cluster_num) %in% names(gamma_dict))) {
      return(list(excluded = TRUE, cluster_num = cluster_num, reason = "resolution_search_failed"))
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    gamma_test <- seq(gamma_range[1], gamma_range[2], length.out = min(5, diff(gamma_range) * 100 + 1))
    
    # Test multiple gammas in parallel
    ic_scores <- parallel::mclapply(gamma_test, function(gamma_val) {
      cluster_results <- replicate(10, {
        leiden_clustering(igraph_obj, gamma_val, objective_function, 5, 0.01)
      }, simplify = TRUE)
      
      extracted_results <- extract_clustering_array(cluster_results)
      ic_result <- calculate_ic_from_extracted(extracted_results)
      return(1 / ic_result)
    }, mc.cores = n_workers)
    
    ic_scores <- unlist(ic_scores)
    excluded <- min(ic_scores, na.rm = TRUE) >= remove_threshold
    reason <- if (excluded) "high_inconsistency" else "passed_filtering"
    
    return(list(excluded = excluded, cluster_num = cluster_num, reason = reason))
  }, mc.cores = n_workers)
  
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
  excluded_entries <- data.table::rbindlist(excluded_entries[!sapply(excluded_entries, is.null)])
  
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
  
  cluster_results <- parallel::mclapply(valid_clusters, function(cluster_num) {
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
  }, mc.cores = n_workers)
  
  # Combine successful results
  successful_results <- data.table::rbindlist(cluster_results[!sapply(cluster_results, is.null)])
  
  # Combine all results (excluded + successful)
  if (nrow(excluded_entries) > 0 && nrow(successful_results) > 0) {
    results_dt <- data.table::rbindlist(list(excluded_entries, successful_results))
  } else if (nrow(excluded_entries) > 0) {
    results_dt <- excluded_entries
  } else {
    results_dt <- successful_results
  }
  
  # Sort by cluster number
  results_dt <- results_dt[order(cluster_number)]
  
  # Convert back to list format for compatibility
  return(list(
    gamma = results_dt$gamma,
    labels = results_dt$labels,
    ic = results_dt$ic,
    ic_vec = results_dt$ic_vec,
    n_cluster = results_dt$cluster_number,
    best_labels = results_dt$best_labels,
    n_iter = results_dt$n_iter,
    mei = results_dt$mei,
    k = results_dt$k,
    excluded = results_dt$excluded,
    exclusion_reason = results_dt$exclusion_reason
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
  
  # Process cluster numbers in parallel
  range_results <- parallel::mclapply(cluster_range, function(target_clusters) {
    left <- start_g
    right <- end_g
    
    # Binary search for lower bound
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > resolution_tolerance) {
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
    
    # Binary search for upper bound
    left <- left_bound
    right <- end_g
    
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > resolution_tolerance) {
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
    
    # Return results as a data.table row
    data.table::data.table(
      cluster_number = target_clusters,
      left_bound = if (objective_function == "CPM") exp(left_bound) else left_bound,
      right_bound = if (objective_function == "CPM") exp(right_bound) else right_bound
    )
  }, mc.cores = n_workers)
  
  # Combine results
  results_dt <- data.table::rbindlist(range_results)
  
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
  clustering_results <- parallel::mclapply(gamma_sequence, function(gamma_val) {
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
  ic_scores <- parallel::mclapply(clustering_matrices, function(cluster_matrix) {
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
      new_results <- parallel::mclapply(seq_along(current_gammas), function(i) {
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
  ic_bootstrap <- parallel::mclapply(seq_len(n_bootstrap), function(i) {
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