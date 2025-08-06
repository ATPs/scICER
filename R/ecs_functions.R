#' @importFrom stats mean
#' @importFrom base diag
NULL

#' Calculate Element-Centric Similarity (ECS) between two clustering results
#'
#' @description
#' ECS calculates similarity between two clustering results based on how consistently 
#' each individual cell is grouped with other cells. This function now uses the 
#' optimized ClustAssess implementation for significantly better performance.
#'
#' @param cluster_a First clustering result (vector of cluster assignments)
#' @param cluster_b Second clustering result (vector of cluster assignments)
#' @param d Damping parameter (default: 0.9)
#' @param return_vector Whether to return similarity scores for each cell (default: FALSE)
#'
#' @return Either mean ECS score (if return_vector=FALSE) or vector of ECS scores for each cell
#' @export
calculate_ecs <- function(cluster_a, cluster_b, d = 0.9, return_vector = FALSE) {
  # Check if ClustAssess is available, if not fall back to our implementation
  if (requireNamespace("ClustAssess", quietly = TRUE)) {
    # Use ClustAssess functions for much better performance
    if (return_vector) {
      return(ClustAssess::element_sim_elscore(cluster_a, cluster_b, alpha = d))
    } else {
      return(ClustAssess::element_sim(cluster_a, cluster_b, alpha = d))
    }
  } else {
    # Fall back to our Julia-like implementation if ClustAssess is not available
    warning("ClustAssess package not found. Using fallback implementation. For better performance, install ClustAssess: install.packages('ClustAssess')")
    
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
    
    # Calculate ECS scores - matching Julia's simmat_v2 behavior
    for (i in 1:n) {
      pos_a <- which(unique_a == cluster_a[i])
      pos_b <- which(unique_b == cluster_b[i])
      
      if (is.nan(unique_ecs_vals[pos_a, pos_b])) {
        # Get all neighbors
        neighbors_a <- groups_a[[pos_a]]
        neighbors_b <- groups_b[[pos_b]]
        all_neighbors <- unique(c(neighbors_a, neighbors_b))
        
        # Initialize PageRank vectors (matching Julia's approach)
        ppr1 <- numeric(n)
        ppr2 <- numeric(n)
        
        # Set cluster members to cluster size (Julia's approach)
        ppr1[neighbors_a] <- cluster_sizes_a[pos_a]
        ppr2[neighbors_b] <- cluster_sizes_b[pos_b]
        
        # Set current node to base score + cluster size (Julia's approach)
        ppr1[i] <- 1.0 - d + cluster_sizes_a[pos_a]
        ppr2[i] <- 1.0 - d + cluster_sizes_b[pos_b]
        
        # Calculate L1 distance
        l1_distance <- sum(abs(ppr2[all_neighbors] - ppr1[all_neighbors]))
        
        ecs_scores[i] <- l1_distance
        unique_ecs_vals[pos_a, pos_b] <- l1_distance
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
}

#' Extract unique clustering arrays and their probabilities
#'
#' @param clustering_matrix Matrix where each column is a clustering result
#' @return List with unique clustering arrays and their probabilities
#' @keywords internal
extract_clustering_array <- function(clustering_matrix) {
  # Convert each column to a single string for fast unique identification
  clustering_strings <- apply(clustering_matrix, 2, paste, collapse = ",")
  
  # Count frequencies of each unique clustering string
  clustering_counts <- table(clustering_strings)
  
  # Get the unique clusterings as vectors
  unique_strings <- names(clustering_counts)
  unique_clusterings <- lapply(strsplit(unique_strings, ","), as.integer)
  
  # Get probabilities
  probabilities <- as.vector(clustering_counts) / ncol(clustering_matrix)
  
  # Sort by probability (descending)
  sort_indices <- order(probabilities, decreasing = TRUE)
  
  # Return sorted results
  return(list(
    arr = unique_clusterings[sort_indices],
    parr = probabilities[sort_indices]
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
  
  # Calculate pairwise similarity matrix using vectorized operations
  similarity_matrix <- matrix(1, nrow = n_clusterings, ncol = n_clusterings)
  
  # Convert list of clusterings to a matrix for faster operations
  clustering_matrix <- do.call(cbind, unique_clusterings)
  
  # Calculate similarities for upper triangle only
  for (i in 1:(n_clusterings - 1)) {
    # Calculate similarities for all j > i at once
    similarities <- sapply((i + 1):n_clusterings, function(j) {
      calculate_ecs(clustering_matrix[, i], clustering_matrix[, j])
    })
    
    # Fill both upper and lower triangle
    similarity_matrix[i, (i + 1):n_clusterings] <- similarities
    similarity_matrix[(i + 1):n_clusterings, i] <- similarities
  }
  
  # Calculate weighted IC score using matrix operations
  ic_score <- sum(outer(probabilities, probabilities) * similarity_matrix)
  
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