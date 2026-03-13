NULL

#' Calculate Element-Centric Similarity (ECS) between two clustering results
#'
#' @description
#' ECS calculates similarity between two clustering results based on how consistently 
#' each individual cell is grouped with other cells. This function uses the 
#' optimized ClustAssess implementation for significantly better performance.
#'
#' @details
#' This function is the core pairwise similarity primitive used across scICER.
#' It performs strict input validation before calling C++ routines from
#' \pkg{ClustAssess}:
#' \itemize{
#'   \item cluster vectors must be non-empty and have identical length;
#'   \item NA labels are not allowed;
#'   \item if names are present, both vectors must have exactly the same names.
#' }
#'
#' Labels are internally remapped to compact integer IDs to avoid numeric range
#' issues in downstream routines while preserving partition structure.
#'
#' When \code{return_vector = FALSE}, the returned scalar ECS is computed as
#' the mean of element-level scores (\code{element_sim_elscore}) for numerical
#' stability on large cell counts.
#'
#' @param cluster_a First clustering result (vector of cluster assignments)
#' @param cluster_b Second clustering result (vector of cluster assignments)
#' @param d Damping parameter (default: 0.9)
#' @param return_vector Whether to return similarity scores for each cell (default: FALSE)
#'
#' @return Either mean ECS score (if return_vector=FALSE) or vector of ECS scores for each cell
#' @export
calculate_ecs <- function(cluster_a, cluster_b, d = 0.9, return_vector = FALSE) {
  valid_label_type <- c("numeric", "integer", "factor", "character")
  if (!inherits(cluster_a, valid_label_type) || !inherits(cluster_b, valid_label_type)) {
    stop("cluster_a and cluster_b must be vectors of numeric/integer/factor/character labels.")
  }
  if (length(cluster_a) == 0 || length(cluster_b) == 0) {
    stop("cluster_a and cluster_b must be non-empty.")
  }
  if (length(cluster_a) != length(cluster_b)) {
    stop("cluster_a and cluster_b must have the same length.")
  }
  if (!is.numeric(d) || length(d) != 1 || is.na(d)) {
    stop("d must be a single numeric value.")
  }
  if (anyNA(cluster_a) || anyNA(cluster_b)) {
    stop("cluster_a and cluster_b cannot contain NA values.")
  }

  names_a <- names(cluster_a)
  names_b <- names(cluster_b)
  if (!is.null(names_a) || !is.null(names_b)) {
    if (is.null(names_a) || is.null(names_b) || any(names_a != names_b)) {
      stop("cluster_a and cluster_b must have identical names when names are provided.")
    }
  }

  # Check if ClustAssess is available
  if (!requireNamespace("ClustAssess", quietly = TRUE)) {
    stop(
      "\nClustAssess package is required for ECS calculations but is not installed.\n",
      "Please install it using:\n\n",
      "  install.packages('ClustAssess')\n\n",
      "ClustAssess provides ~150x faster performance for clustering similarity calculations.\n",
      "More information: https://github.com/Core-Bioinformatics/ClustAssess"
    )
  }

  # Compact relabeling avoids integer range issues in ClustAssess C++ routines.
  compact_labels <- function(x) {
    out <- as.integer(factor(x, levels = unique(x)))
    names(out) <- names(x)
    out
  }

  cluster_a <- compact_labels(cluster_a)
  cluster_b <- compact_labels(cluster_b)

  # Scalar ECS must be computed from element-level scores for numerical stability
  # on large cell counts.
  if (return_vector) {
    return(ClustAssess::element_sim_elscore(cluster_a, cluster_b, alpha = d))
  } else {
    return(mean(ClustAssess::element_sim_elscore(cluster_a, cluster_b, alpha = d)))
  }
}

#' Extract unique clustering arrays and their probabilities
#'
#' @description
#' Deduplicates repeated clustering outcomes and computes their empirical weights.
#'
#' @details
#' This helper converts each clustering column to a compact string signature,
#' counts frequency of each unique signature, and converts back to integer label
#' vectors. The result is sorted by decreasing probability, which improves
#' downstream IC/MEI calculations that iterate over unique outcomes.
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
#' @description
#' Computes per-cell stability scores from a set of repeated clustering outcomes.
#'
#' @details
#' The input should be the output of \code{extract_clustering_array()}, where
#' \code{arr} contains unique clustering label vectors and \code{parr} contains
#' their empirical probabilities.
#'
#' For each pair of unique clustering outcomes, MEI uses element-level ECS scores,
#' weighted by the pair's probabilities, then aggregates across all pairs.
#' Higher values indicate that a cell is assigned consistently across repeated runs.
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
#' @description
#' Calculates a global consistency score for repeated clustering outcomes.
#'
#' @details
#' IC is computed from a probability-weighted pairwise ECS similarity matrix
#' over unique clustering outcomes.
#'
#' Interpretation:
#' \itemize{
#'   \item \code{1.0}: perfect consistency (all effective outcomes agree);
#'   \item values near \code{1.0}: highly stable clustering;
#'   \item larger values: lower stability / higher inconsistency.
#' }
#'
#' @param clustering_array Result from extract_clustering_array
#' @return IC score (lower is better, 1 = perfect consistency)
#' @export
calculate_ic <- function(clustering_array) {
  return(calculate_ic_from_extracted(clustering_array))
}

#' Internal function to calculate IC from extracted clustering array
#'
#' @description
#' Low-level implementation of IC scoring used by public wrappers.
#'
#' @details
#' Given unique clustering outcomes and their probabilities, this function builds
#' a symmetric pairwise ECS matrix and computes the probability-weighted sum.
#' It is separated from \code{calculate_ic()} so internal workflows can reuse
#' the same scoring logic without extra checks or conversions.
#'
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
#' @description
#' Selects a representative clustering from a set of unique outcomes.
#'
#' @details
#' The selected solution is the one with the largest total similarity to all
#' other unique outcomes (row-sum criterion on the pairwise ECS matrix).
#' This heuristic favors the most central/stable partition among candidates.
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
