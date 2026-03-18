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
#' When \code{return_vector = FALSE}, the returned scalar ECS is the mean of the
#' element-level scores from \code{ClustAssess::element_sim_elscore()}.
#' This avoids the known negative-value bug in \code{ClustAssess::element_sim()}
#' on very large flat partitions.
#'
#' @param cluster_a First clustering result (vector of cluster assignments)
#' @param cluster_b Second clustering result (vector of cluster assignments)
#' @param d Damping parameter (default: 0.9)
#' @param return_vector Whether to return similarity scores for each cell (default: FALSE)
#'
#' @return Either scalar ECS score (if \code{return_vector = FALSE}) or vector
#'   of ECS scores for each cell
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

  if (return_vector) {
    return(ClustAssess::element_sim_elscore(cluster_a, cluster_b, alpha = d))
  } else {
    return(mean(ClustAssess::element_sim_elscore(cluster_a, cluster_b, alpha = d)))
  }
}

normalize_clustering_array_input <- function(clustering_array) {
  if (!is.list(clustering_array) ||
      is.null(clustering_array$arr) ||
      is.null(clustering_array$parr)) {
    stop("clustering_array must be a list with `arr` and `parr` entries.")
  }

  clusterings <- clustering_array$arr
  probabilities <- as.numeric(clustering_array$parr)

  if (length(clusterings) == 0L) {
    stop("clustering_array$arr must contain at least one clustering.")
  }
  if (length(probabilities) != length(clusterings)) {
    stop("clustering_array$parr must have the same length as clustering_array$arr.")
  }
  if (any(!is.finite(probabilities)) || any(probabilities < 0)) {
    stop("clustering_array$parr must contain finite non-negative weights.")
  }

  n_cells_each <- vapply(clusterings, length, integer(1))
  if (any(n_cells_each == 0L)) {
    stop("Each clustering in clustering_array$arr must be non-empty.")
  }
  if (length(unique(n_cells_each)) != 1L) {
    stop("All clusterings in clustering_array$arr must have the same length.")
  }

  total_weight <- sum(probabilities)
  if (!is.finite(total_weight) || total_weight <= 0) {
    stop("clustering_array$parr must sum to a positive value.")
  }

  list(
    arr = clusterings,
    parr = probabilities / total_weight,
    n_cells = n_cells_each[[1]]
  )
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
#' MEI is computed as the expected element-level ECS between two independent
#' draws from the empirical clustering distribution. Diagonal self-match
#' probability contributes a value of 1 for every cell, while off-diagonal pairs
#' contribute the pairwise element-level ECS weighted by the product of their
#' normalized probabilities.
#'
#' Higher values indicate that a cell is assigned consistently across repeated
#' runs.
#'
#' @param clustering_array Result from extract_clustering_array
#' @return Vector of MEI scores
#' @export
calculate_mei_from_array <- function(clustering_array) {
  normalized <- normalize_clustering_array_input(clustering_array)
  clusterings <- normalized$arr
  probabilities <- normalized$parr
  n_clusterings <- length(clusterings)
  n_cells <- normalized$n_cells

  if (n_clusterings == 1L) {
    return(rep(1, n_cells))
  }

  mei_scores <- rep(sum(probabilities^2), n_cells)

  for (i in 1:(n_clusterings - 1)) {
    for (j in (i + 1):n_clusterings) {
      sim_scores <- calculate_ecs(clusterings[[i]], clusterings[[j]], return_vector = TRUE)
      mei_scores <- mei_scores + (2 * probabilities[[i]] * probabilities[[j]] * sim_scores)
    }
  }

  mei_scores
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
  normalized <- normalize_clustering_array_input(clustering_array)
  unique_clusterings <- normalized$arr
  probabilities <- normalized$parr
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
  normalized <- normalize_clustering_array_input(clustering_array)
  unique_clusterings <- normalized$arr
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
