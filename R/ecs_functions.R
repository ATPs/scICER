#' Calculate Element-Centric Similarity (ECS) between two clustering results
#'
#' @description
#' ECS calculates similarity between two clustering results based on how consistently 
#' each individual cell is grouped with other cells. This is more efficient than 
#' traditional methods that build consensus matrices.
#'
#' @param cluster_a First clustering result (vector of cluster assignments)
#' @param cluster_b Second clustering result (vector of cluster assignments)
#' @param d Damping parameter (default: 0.9)
#' @param return_vector Whether to return similarity scores for each cell (default: FALSE)
#'
#' @return Either mean ECS score (if return_vector=FALSE) or vector of ECS scores for each cell
#' @export
calculate_ecs <- function(cluster_a, cluster_b, d = 0.9, return_vector = FALSE) {
  # Convert to 0-based indexing for consistency with Julia code
  cluster_a <- as.integer(cluster_a) - 1L
  cluster_b <- as.integer(cluster_b) - 1L
  
  n <- length(cluster_a)
  unique_a <- unique(cluster_a)
  unique_b <- unique(cluster_b)
  
  # Create index groups for each cluster
  groups_a <- lapply(unique_a, function(x) which(cluster_a == x))
  groups_b <- lapply(unique_b, function(x) which(cluster_b == x))
  
  # Calculate cluster sizes and damping factors
  cluster_sizes_a <- d / sapply(groups_a, length)
  cluster_sizes_b <- d / sapply(groups_b, length)
  
  # Initialize matrices for memoization
  unique_ecs_vals <- matrix(NaN, nrow = length(unique_a), ncol = length(unique_b))
  ecs_scores <- numeric(n)
  ppr1 <- numeric(n)
  ppr2 <- numeric(n)
  
  for (i in 1:n) {
    # Get cluster indices (1-based for R)
    cluster_idx_a <- cluster_a[i] + 1L
    cluster_idx_b <- cluster_b[i] + 1L
    
    # Find position in unique cluster arrays
    pos_a <- which(unique_a == cluster_a[i])
    pos_b <- which(unique_b == cluster_b[i])
    
    # Check if already calculated
    if (is.nan(unique_ecs_vals[pos_a, pos_b])) {
      # Get neighboring cells
      neighbors_a <- groups_a[[pos_a]]
      neighbors_b <- groups_b[[pos_b]]
      all_neighbors <- unique(c(neighbors_a, neighbors_b))
      
      # Calculate personalized PageRank vectors
      for (idx in neighbors_a) {
        ppr1[idx] <- cluster_sizes_a[pos_a]
      }
      ppr1[i] <- 1.0 - d + cluster_sizes_a[pos_a]
      
      for (idx in neighbors_b) {
        ppr2[idx] <- cluster_sizes_b[pos_b]
      }
      ppr2[i] <- 1.0 - d + cluster_sizes_b[pos_b]
      
      # Calculate L1 distance
      l1_distance <- 0
      for (j in all_neighbors) {
        l1_distance <- l1_distance + abs(ppr2[j] - ppr1[j])
      }
      
      ecs_scores[i] <- l1_distance
      
      # Reset PPR vectors
      for (idx in neighbors_a) {
        ppr1[idx] <- 0.0
      }
      for (idx in neighbors_b) {
        ppr2[idx] <- 0.0
      }
      
      # Store in memoization matrix
      unique_ecs_vals[pos_a, pos_b] <- ecs_scores[i]
    } else {
      ecs_scores[i] <- unique_ecs_vals[pos_a, pos_b]
    }
  }
  
  # Convert to similarity scores
  similarity_scores <- 1 - (1 / (2 * d)) * ecs_scores
  
  if (return_vector) {
    return(similarity_scores)
  } else {
    return(mean(similarity_scores))
  }
}

#' Extract unique clustering arrays and their probabilities
#'
#' @param clustering_matrix Matrix where each column is a clustering result
#' @return List with unique clustering arrays and their probabilities
#' @keywords internal
extract_clustering_array <- function(clustering_matrix) {
  # Convert each column to a vector and find unique patterns
  clustering_vectors <- lapply(1:ncol(clustering_matrix), function(i) {
    clustering_matrix[, i]
  })
  
  # Count occurrences of each unique clustering
  unique_clusterings <- list()
  counts <- numeric()
  
  for (i in seq_along(clustering_vectors)) {
    current_clustering <- clustering_vectors[[i]]
    
    # Check if this clustering already exists
    found <- FALSE
    for (j in seq_along(unique_clusterings)) {
      if (identical(current_clustering, unique_clusterings[[j]])) {
        counts[j] <- counts[j] + 1
        found <- TRUE
        break
      }
    }
    
    if (!found) {
      unique_clusterings <- append(unique_clusterings, list(current_clustering))
      counts <- c(counts, 1)
    }
  }
  
  # Calculate probabilities
  probabilities <- counts / sum(counts)
  
  # Sort by probability (descending)
  sort_indices <- order(probabilities, decreasing = TRUE)
  unique_clusterings <- unique_clusterings[sort_indices]
  probabilities <- probabilities[sort_indices]
  
  # Normalize probabilities
  probabilities <- probabilities / sum(probabilities)
  
  return(list(
    arr = unique_clusterings,
    parr = probabilities
  ))
}

#' Calculate Mutual Element-wise Information (MEI) from clustering array
#'
#' @param clustering_array Result from extract_clustering_array
#' @return Vector of MEI scores
#' @export
calculate_mei_from_array <- function(clustering_array) {
  if (length(clustering_array$arr) == 1) {
    return(rep(1, length(clustering_array$arr[[1]])))
  }
  
  n_clusterings <- length(clustering_array$arr)
  n_cells <- length(clustering_array$arr[[1]])
  
  # Calculate pairwise similarities
  similarity_matrix <- matrix(0, nrow = n_clusterings, ncol = n_clusterings)
  
  for (i in 1:(n_clusterings - 1)) {
    for (j in (i + 1):n_clusterings) {
      sim_scores <- calculate_ecs(
        clustering_array$arr[[i]], 
        clustering_array$arr[[j]], 
        return_vector = TRUE
      )
      weighted_scores <- sim_scores * (clustering_array$parr[i] + clustering_array$parr[j])
      similarity_matrix[i, j] <- mean(weighted_scores)
      similarity_matrix[j, i] <- similarity_matrix[i, j]
    }
  }
  
  # Set diagonal to 1
  diag(similarity_matrix) <- 1
  
  # Calculate MEI as row sums divided by (n_clusterings - 1)
  mei_scores <- rowSums(similarity_matrix) / (n_clusterings - 1)
  
  return(mei_scores)
}

#' Calculate Inconsistency (IC) score from extracted clustering results
#'
#' @param clustering_array Result from extract_clustering_array
#' @return IC score (lower is better, 1 = perfect consistency)
#' @export
calculate_ic <- function(clustering_array) {
  return(calculate_ic_from_extracted(clustering_array))
}

#' Internal function to calculate IC from extracted clustering array
#' @keywords internal
calculate_ic_from_extracted <- function(clustering_array) {
  unique_clusterings <- clustering_array$arr
  probabilities <- clustering_array$parr
  
  n_clusterings <- length(unique_clusterings)
  
  if (n_clusterings == 1) {
    return(1.0)  # Perfect consistency
  }
  
  # Calculate pairwise similarity matrix
  similarity_matrix <- matrix(1, nrow = n_clusterings, ncol = n_clusterings)
  
  for (i in 1:(n_clusterings - 1)) {
    for (j in (i + 1):n_clusterings) {
      similarity <- calculate_ecs(unique_clusterings[[i]], unique_clusterings[[j]])
      similarity_matrix[i, j] <- similarity
      similarity_matrix[j, i] <- similarity
    }
  }
  
  # Calculate weighted IC score
  ic_score <- 0
  for (i in 1:n_clusterings) {
    for (j in 1:n_clusterings) {
      ic_score <- ic_score + probabilities[i] * probabilities[j] * similarity_matrix[i, j]
    }
  }
  
  return(ic_score)
}

#' Get the best clustering from extracted clustering array
#'
#' @param clustering_array Result from extract_clustering_array
#' @return Best clustering labels
#' @keywords internal
get_best_clustering <- function(clustering_array) {
  unique_clusterings <- clustering_array$arr
  n_clusterings <- length(unique_clusterings)
  
  if (n_clusterings == 1) {
    return(unique_clusterings[[1]])
  }
  
  # Calculate similarity matrix
  similarity_matrix <- matrix(1, nrow = n_clusterings, ncol = n_clusterings)
  
  for (i in 1:(n_clusterings - 1)) {
    for (j in (i + 1):n_clusterings) {
      similarity <- calculate_ecs(unique_clusterings[[i]], unique_clusterings[[j]])
      similarity_matrix[i, j] <- similarity
      similarity_matrix[j, i] <- similarity
    }
  }
  
  # Find clustering with highest average similarity to others
  row_sums <- rowSums(similarity_matrix)
  best_index <- which.max(row_sums)
  
  return(unique_clusterings[[best_index]])
} 