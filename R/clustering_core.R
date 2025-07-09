#' @import igraph
#' @importFrom stats median
#' @importFrom utils txtProgressBar setTxtProgressBar
NULL

#' Core clustering algorithm implementing binary search and optimization
#'
#' @param igraph_obj igraph object to cluster
#' @param cluster_range Vector of cluster numbers to test
#' @param n_workers Number of parallel workers
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
clustering_main <- function(igraph_obj, cluster_range, n_workers, n_trials, n_bootstrap,
                           beta, n_iterations, max_iterations, objective_function,
                           remove_threshold, resolution_tolerance, verbose) {
  
  # Initialize results storage
  list_gamma <- list()
  list_labels <- list()
  list_incons <- numeric()
  list_all_ic <- list()
  list_best_labels <- list()
  list_n_iter <- numeric()
  list_mn <- numeric()
  list_k <- numeric()
  
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
  
  # Filter out problematic cluster numbers
  excluded_numbers <- filter_problematic_clusters(
    gamma_dict, cluster_range, igraph_obj, objective_function,
    remove_threshold, n_workers, verbose
  )
  
  if (verbose && length(excluded_numbers) > 0) {
    message(paste("Excluded cluster numbers:", paste(excluded_numbers, collapse = ", ")))
  }
  
  valid_clusters <- setdiff(cluster_range, excluded_numbers)
  
  if (length(valid_clusters) == 0) {
    if (verbose) message("No valid cluster numbers found")
    return(NULL)
  }
  
  # Optimize clustering for each valid cluster number
  if (verbose) message("Optimizing clustering for each valid cluster number...")
  
  pb <- NULL
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = length(valid_clusters), style = 3)
  }
  
  for (i in seq_along(valid_clusters)) {
    cluster_num <- valid_clusters[i]
    
    if (verbose) {
      setTxtProgressBar(pb, i)
    }
    
    # Get resolution range for this cluster number
    if (!(cluster_num %in% names(gamma_dict))) {
      next
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    
    # Optimize clustering within this range
    result <- optimize_clustering(
      igraph_obj, cluster_num, gamma_range, objective_function,
      n_trials, n_bootstrap, beta, n_iterations, max_iterations,
      resolution_tolerance, n_workers
    )
    
    if (!is.null(result)) {
      list_gamma <- append(list_gamma, list(result$gamma))
      list_labels <- append(list_labels, list(result$labels))
      list_incons <- c(list_incons, result$ic_median)
      list_all_ic <- append(list_all_ic, list(result$ic_bootstrap))
      list_best_labels <- append(list_best_labels, list(result$best_labels))
      list_n_iter <- c(list_n_iter, result$n_iterations)
      list_mn <- c(list_mn, cluster_num)
      list_k <- c(list_k, result$k)
    }
  }
  
  if (verbose) {
    close(pb)
  }
  
  # Calculate MEI scores
  if (verbose) message("Calculating MEI scores...")
  mei_scores <- lapply(list_labels, calculate_mei_from_array)
  
  # Return results
  return(list(
    gamma = list_gamma,
    labels = list_labels,
    ic = list_incons,
    ic_vec = list_all_ic,
    n_cluster = list_mn,
    best_labels = list_best_labels,
    n_iter = list_n_iter,
    mei = mei_scores,
    k = list_k
  ))
}

#' Find resolution parameter ranges for each cluster number using binary search
#' @keywords internal
find_resolution_ranges <- function(igraph_obj, cluster_range, start_g, end_g,
                                  objective_function, resolution_tolerance, n_workers, verbose) {
  
  gamma_dict <- list()
  n_preliminary_trials <- 10
  beta_preliminary <- 0.01
  n_iter_preliminary <- 3
  
  for (target_clusters in cluster_range) {
    left <- start_g
    right <- end_g
    
    # Binary search for lower bound
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > resolution_tolerance) {
      mid <- (left + right) / 2
      gamma_val <- if (objective_function == "modularity") mid else exp(mid)
      
      # Test clustering with current gamma
      cluster_results <- replicate(n_preliminary_trials, {
        leiden_clustering(igraph_obj, gamma_val, objective_function, n_iter_preliminary, beta_preliminary)
      }, simplify = TRUE)
      
      n_clusters_obtained <- median(apply(cluster_results, 2, max) + 1)
      
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
      
      # Test clustering with current gamma
      cluster_results <- replicate(n_preliminary_trials, {
        leiden_clustering(igraph_obj, gamma_val, objective_function, n_iter_preliminary, beta_preliminary)
      }, simplify = TRUE)
      
      n_clusters_obtained <- median(apply(cluster_results, 2, max) + 1)
      
      if (n_clusters_obtained > target_clusters) {
        right <- mid
      } else {
        left <- mid
      }
    }
    
    right_bound <- left
    
    # Store the range
    if (objective_function == "CPM") {
      gamma_dict[[as.character(target_clusters)]] <- c(exp(left_bound), exp(right_bound))
    } else {
      gamma_dict[[as.character(target_clusters)]] <- c(left_bound, right_bound)
    }
  }
  
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
  
  # Create gamma sequence
  if (objective_function == "modularity") {
    delta_g <- (gamma_range[2] - gamma_range[1]) / n_steps
    if (gamma_range[1] != gamma_range[2]) {
      gamma_sequence <- seq(gamma_range[1], gamma_range[2], by = delta_g)
    } else {
      gamma_sequence <- seq(gamma_range[1] - delta_g, gamma_range[1] + delta_g, by = delta_g / n_steps * 2)
    }
  } else { # CPM
    if (gamma_range[1] != gamma_range[2]) {
      delta_g <- (log(gamma_range[2]) - log(gamma_range[1])) / n_steps
      gamma_sequence <- exp(seq(log(gamma_range[1]), log(gamma_range[2]), by = delta_g))
    } else {
      delta_g <- 0.1
      gamma_sequence <- exp(seq(log(gamma_range[1]) - delta_g, log(gamma_range[1]) + delta_g, by = delta_g / n_steps * 2))
    }
  }
  
  # Test initial clustering for each gamma
  clustering_results <- list()
  mean_clusters <- numeric(length(gamma_sequence))
  
  for (i in seq_along(gamma_sequence)) {
    gamma_val <- gamma_sequence[i]
    cluster_matrix <- replicate(n_trials, {
      leiden_clustering(igraph_obj, gamma_val, objective_function, n_iterations, beta)
    }, simplify = TRUE)
    
    clustering_results[[i]] <- cluster_matrix
    mean_clusters[i] <- median(apply(cluster_matrix, 2, max) + 1)
  }
  
  # Filter for target cluster number
  valid_indices <- which(mean_clusters == target_clusters)
  
  if (length(valid_indices) == 0) {
    return(NULL)
  }
  
  gamma_sequence <- gamma_sequence[valid_indices]
  clustering_results <- clustering_results[valid_indices]
  mean_clusters <- mean_clusters[valid_indices]
  
  # Calculate IC for each gamma
  ic_scores <- sapply(clustering_results, function(cluster_matrix) {
    extracted <- extract_clustering_array(cluster_matrix)
    ic_result <- calculate_ic_from_extracted(extracted)
    return(1 / ic_result)  # Convert to inconsistency score
  })
  
  # Find the best gamma (one with IC = 1 if exists, otherwise lowest IC)
  best_index <- which(ic_scores == 1)[1]
  if (is.na(best_index)) {
    best_index <- which.min(ic_scores)
  }
  
  best_gamma <- gamma_sequence[best_index]
  best_clustering <- clustering_results[[best_index]]
  k <- n_iterations
  
  # If no perfect IC found, iterate to improve
  if (ic_scores[best_index] != 1 && length(gamma_sequence) > 1) {
    # Iterative improvement
    current_results <- clustering_results
    current_gammas <- gamma_sequence
    current_ic <- ic_scores
    
    # Track stability
    ic_history <- matrix(rep(current_ic, 10), nrow = length(current_ic), ncol = 10)
    
    while (k < max_iterations) {
      k <- k + delta_n
      
      # Update clustering results
      for (i in seq_along(current_gammas)) {
        gamma_val <- current_gammas[i]
        
        # Use previous results as initialization
        new_clustering <- replicate(n_trials, {
          init_membership <- current_results[[i]][, sample(ncol(current_results[[i]]), 1)]
          leiden_clustering(igraph_obj, gamma_val, objective_function, delta_n, beta, init_membership)
        }, simplify = TRUE)
        
        current_results[[i]] <- new_clustering
      }
      
      # Recalculate IC scores
      new_ic <- sapply(current_results, function(cluster_matrix) {
        extracted <- extract_clustering_array(cluster_matrix)
        ic_result <- calculate_ic_from_extracted(extracted)
        return(1 / ic_result)
      })
      
      # Update IC history
      ic_history[, 1:(ncol(ic_history)-1)] <- ic_history[, 2:ncol(ic_history)]
      ic_history[, ncol(ic_history)] <- new_ic
      
      # Check for stability and convergence
      stable_indices <- apply(ic_history, 1, function(row) all(diff(row) == 0))
      perfect_indices <- which(new_ic == 1)
      
      # Selection criteria
      if (length(perfect_indices) > 0) {
        best_index <- perfect_indices[1]
        break
      } else if (all(stable_indices)) {
        best_index <- which.min(new_ic)
        break
      } else if (k > 100 && median(new_ic) > 1.1) {
        best_index <- which.min(new_ic)
        break
      } else {
        # Continue with best performing gammas
        keep_indices <- (new_ic <= quantile(new_ic, 0.5)) | stable_indices
        keep_indices[which.min(new_ic)] <- TRUE
        
        if (length(which(keep_indices)) == 1) {
          best_index <- which(keep_indices)
          break
        }
        
        current_gammas <- current_gammas[keep_indices]
        current_results <- current_results[keep_indices]
        ic_history <- ic_history[keep_indices, , drop = FALSE]
        new_ic <- new_ic[keep_indices]
      }
      
      current_ic <- new_ic
    }
    
    best_gamma <- current_gammas[best_index]
    best_clustering <- current_results[[best_index]]
  }
  
  # Bootstrap analysis
  ic_bootstrap <- replicate(n_bootstrap, {
    sample_indices <- sample(ncol(best_clustering), ncol(best_clustering), replace = TRUE)
    bootstrap_matrix <- best_clustering[, sample_indices]
    extracted <- extract_clustering_array(bootstrap_matrix)
    ic_result <- calculate_ic_from_extracted(extracted)
    return(1 / ic_result)
  })
  
  ic_median <- median(ic_bootstrap)
  
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