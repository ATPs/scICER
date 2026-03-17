#' @import ggplot2
#' @importFrom stats quantile median
NULL

#' Plot Inconsistency (IC) scores across different cluster numbers
#'
#' @description
#' Creates a boxplot showing the distribution of IC scores for each cluster number 
#' tested. Lower IC scores indicate more consistent clustering results.
#'
#' @details
#' \code{plot_ic()} is typically the first diagnostic plot after
#' \code{scICE_clustering()}:
#' \itemize{
#'   \item each x-position is a tested cluster number;
#'   \item box/jitter values are bootstrap IC scores (\code{ic_vec});
#'   \item the threshold line marks your chosen consistency cutoff.
#' }
#'
#' If no valid IC values are available (for example all targets excluded),
#' the function returns an informative empty plot instead of failing.
#'
#' @param scice_results Results object from scICE_clustering function
#' @param threshold IC threshold line to display (default: 1.005)
#' @param figure_size Figure size as c(width, height) (default: c(10, 6))
#' @param title Plot title (default: "Clustering Consistency Analysis")
#' @param show_threshold Whether to show threshold line (default: TRUE)
#' @param show_gamma Whether to add selected gamma values below each x-axis
#'   cluster label with a 45 degree rotation (default: TRUE)
#'
#' @return ggplot object
#' @export
#' @examples
#' \dontrun{
#' # Run scICE analysis
#' scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)
#' 
#' # Plot IC results
#' plot_ic(scice_results)
#' 
#' # Customize plot
#' plot_ic(scice_results, threshold = 1.01, title = "My Clustering Analysis")
#' }
plot_ic <- function(scice_results, threshold = 1.005, figure_size = c(10, 6),
                   title = "Clustering Consistency Analysis", show_threshold = TRUE,
                   show_gamma = TRUE) {

  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }

  # Check if there are any IC results to plot
  if (length(scice_results$ic_vec) == 0) {
    # Create empty plot with message when no IC results are found
    empty_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = title) +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = "NO consistent cluster numbers found", 
                       size = 8, color = "red", fontface = "bold") +
      ggplot2::xlim(0, 1) +
      ggplot2::ylim(0, 1)
    return(empty_plot)
  }

  # Handle both old results (without excluded info) and new results (with excluded info)
  has_exclusion_info <- !is.null(scice_results$excluded)
  
  # Prepare data for plotting - only include non-excluded clusters with IC data
  if (has_exclusion_info) {
    valid_indices <- which(!scice_results$excluded & !is.na(scice_results$ic))
  } else {
    valid_indices <- which(!is.na(scice_results$ic))
  }
  
  if (length(valid_indices) == 0) {
    # Create empty plot with message when no valid IC results are found
    empty_plot <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = title) +
      ggplot2::annotate("text", x = 0.5, y = 0.5, 
                       label = "NO consistent cluster numbers found", 
                       size = 8, color = "red", fontface = "bold") +
      ggplot2::xlim(0, 1) +
      ggplot2::ylim(0, 1)
    return(empty_plot)
  }
  
  # Create plot data for valid clusters
  plot_data <- data.frame(
    cluster_number = rep(scice_results$n_cluster[valid_indices],
                        sapply(scice_results$ic_vec[valid_indices], length)),
    ic_score = unlist(scice_results$ic_vec[valid_indices])
  )
  
  # Determine consistent clusters
  consistent_cluster_numbers <- character()
  if (has_exclusion_info) {
    consistent_indices <- which(!scice_results$excluded & scice_results$ic < threshold)
    consistent_cluster_numbers <- scice_results$n_cluster[consistent_indices]
  } else {
    consistent_indices <- which(scice_results$ic < threshold)
    consistent_cluster_numbers <- scice_results$n_cluster[consistent_indices]
  }
  
  # Add consistency status to plot data
  plot_data$is_consistent <- plot_data$cluster_number %in% consistent_cluster_numbers
  
  # Build subtitle text
  subtitle_parts <- c(
    paste("Lower IC scores indicate more consistent clustering. Threshold:", threshold)
  )
  
  if (has_exclusion_info) {
    excluded_clusters <- scice_results$n_cluster[scice_results$excluded]
    if (length(excluded_clusters) > 0) {
      subtitle_parts <- c(
        subtitle_parts,
        paste("Excluded clusters:", paste(excluded_clusters, collapse = ", "))
      )
    }
  }
  
  cluster_levels <- as.character(sort(unique(plot_data$cluster_number)))
  x_axis_labels <- stats::setNames(cluster_levels, cluster_levels)
  x_axis_text <- ggplot2::element_text(size = 10)
  x_axis_title_margin <- ggplot2::margin(t = 0)
  bottom_plot_margin <- 5.5

  if (isTRUE(show_gamma)) {
    format_gamma <- function(x) {
      if (!is.finite(x)) {
        return("NA")
      }
      sprintf("%.2e", x)
    }

    gamma_values <- stats::setNames(
      rep("NA", length(cluster_levels)),
      cluster_levels
    )
    gamma_indices <- valid_indices[!is.na(scice_results$gamma[valid_indices])]

    if (length(gamma_indices) > 0) {
      gamma_values[as.character(scice_results$n_cluster[gamma_indices])] <-
        vapply(scice_results$gamma[gamma_indices], format_gamma, character(1))
    }

    x_axis_labels <- stats::setNames(
      paste0(cluster_levels, "\n", gamma_values[cluster_levels]),
      cluster_levels
    )
    x_axis_text <- ggplot2::element_text(
      size = 10,
      angle = 45,
      hjust = 1,
      vjust = 1
    )
    x_axis_title_margin <- ggplot2::margin(t = 12)
    bottom_plot_margin <- 18
  }
  subtitle_text <- paste(subtitle_parts, collapse = "\n")

  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = factor(cluster_number), y = ic_score)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = cluster_number, fill = is_consistent),
                 outlier.alpha = 0.6,
                 alpha = 0.7) +
    ggplot2::scale_fill_manual(
      name = "Consistency",
      values = c("TRUE" = "lightgreen", "FALSE" = "lightgray"),
      labels = c("TRUE" = paste("Consistent (IC <", threshold, ")"), 
                "FALSE" = paste("Inconsistent (IC >=", threshold, ")"))
    ) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.4, size = 0.8, height = 0) +
    ggplot2::scale_x_discrete(
      name = "Number of Clusters",
      labels = x_axis_labels
    ) +
    ggplot2::scale_y_continuous(name = "Inconsistency (IC) Score",
                      limits = c(0.99, max(plot_data$ic_score) * 1.05)) +
    ggplot2::labs(title = title,
         subtitle = subtitle_text) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.title.x = ggplot2::element_text(size = 12, margin = x_axis_title_margin),
      axis.text = ggplot2::element_text(size = 10),
      axis.text.x = x_axis_text,
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(5.5, 5.5, bottom_plot_margin, 5.5)
    )

  # Add threshold line if requested
  if (show_threshold) {
    p <- p + ggplot2::geom_hline(yintercept = threshold,
                       color = "red",
                       linetype = "dashed",
                       linewidth = 1) +
             ggplot2::annotate("text",
                     x = length(unique(plot_data$cluster_number)) * 0.8,
                     y = threshold + 0.001,
                     label = paste("Threshold =", threshold),
                     color = "red",
                     size = 3.5)
  }

  return(p)
}

#' Extract robust clustering labels from scICE results
#'
#' @description
#' Extracts clustering labels for cluster numbers that meet the consistency threshold.
#' Returns a data frame with cell names and cluster assignments for each consistent
#' cluster number.
#'
#' @details
#' Returned labels are converted to 1-based indexing for standard R usage.
#' Output columns are:
#' \itemize{
#'   \item \code{cell_id}: cell names from the original Seurat object;
#'   \item \code{clusters_<k>}: robust labels for consistent cluster number \code{k}.
#' }
#'
#' If \code{return_seurat = TRUE}, these \code{clusters_<k>} columns are added to
#' \code{meta.data} of \code{object}. The function validates that
#' \code{scice_results$cell_names} and \code{Cells(object)} refer to the same cells.
#' Cell order differences are handled by cell-name matching.
#'
#' @param scice_results Results object from scICE_clustering function
#' @param threshold IC threshold for determining consistent clusters (default: 1.005)
#' @param return_seurat Whether to return results as Seurat metadata (default: FALSE)
#' @param object Optional Seurat object used when \code{return_seurat = TRUE}.
#'   Pass the original input Seurat object from \code{scICE_clustering()}.
#'
#' @return Data frame with cell names and cluster assignments, or Seurat object with updated metadata
#' @export
#' @examples
#' \dontrun{
#' # Extract consistent clustering labels
#' consistent_labels <- get_robust_labels(scice_results, threshold = 1.005)
#' 
#' # View the results
#' head(consistent_labels)
#' 
#' # Add to Seurat object directly
#' seurat_obj <- get_robust_labels(scice_results, return_seurat = TRUE, object = seurat_obj)
#' }
get_robust_labels <- function(scice_results, threshold = 1.005, return_seurat = FALSE, object = NULL) {
  
  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }
  
  # Handle both old results (without excluded info) and new results (with excluded info)
  has_exclusion_info <- !is.null(scice_results$excluded)
  
  # Find consistent cluster numbers, excluding already-excluded clusters
  if (has_exclusion_info) {
    # For new format: only consider non-excluded clusters with valid IC scores
    valid_indices <- which(!scice_results$excluded & !is.na(scice_results$ic))
    consistent_indices <- valid_indices[scice_results$ic[valid_indices] < threshold]
  } else {
    # For old format: use all available clusters
    consistent_indices <- which(!is.na(scice_results$ic) & scice_results$ic < threshold)
  }
  
  if (length(consistent_indices) == 0) {
    if (has_exclusion_info) {
      excluded_clusters <- scice_results$n_cluster[scice_results$excluded]
      if (length(excluded_clusters) > 0) {
        warning(paste("No cluster numbers meet the consistency threshold.",
                     "Excluded clusters:", paste(excluded_clusters, collapse = ", ")))
      } else {
        warning("No cluster numbers meet the consistency threshold")
      }
    } else {
      warning("No cluster numbers meet the consistency threshold")
    }
    return(NULL)
  }
  
  consistent_clusters <- scice_results$n_cluster[consistent_indices]
  
  # Create output data frame
  output_df <- data.frame(
    cell_id = scice_results$cell_names,
    stringsAsFactors = FALSE
  )
  
  # Add clustering results for each consistent cluster number
  for (i in seq_along(consistent_indices)) {
    cluster_num <- consistent_clusters[i]
    cluster_labels <- scice_results$best_labels[[consistent_indices[i]]]
    
    # Skip if labels are NULL (shouldn't happen for consistent clusters, but be safe)
    if (is.null(cluster_labels)) {
      warning(paste("No labels found for cluster number", cluster_num))
      next
    }
    
    # Convert to 1-based indexing for R
    cluster_labels <- cluster_labels + 1
    
    # Add to data frame
    col_name <- paste0("clusters_", cluster_num)
    output_df[[col_name]] <- cluster_labels
  }
  
  # Check if we actually added any clustering columns
  if (ncol(output_df) == 1) {
    warning("No valid clustering results could be extracted")
    return(NULL)
  }
  
  if (return_seurat) {
    seurat_obj <- object
    
    if (is.null(seurat_obj)) {
      if (!is.null(scice_results$seurat_object)) {
        warning(
          "`return_seurat = TRUE` without `object` is deprecated. ",
          "Using `scice_results$seurat_object` for backward compatibility."
        )
        seurat_obj <- scice_results$seurat_object
      } else {
        stop("When `return_seurat = TRUE`, please provide `object` (the original Seurat object).")
      }
    }
    
    if (!inherits(seurat_obj, "Seurat")) {
      stop("`object` must be a Seurat object when `return_seurat = TRUE`.")
    }
    
    if (is.null(scice_results$cell_names) || length(scice_results$cell_names) == 0) {
      stop("`scice_results$cell_names` is missing. Cannot map labels back to Seurat metadata.")
    }
    
    result_cells <- as.character(scice_results$cell_names)
    object_cells <- SeuratObject::Cells(seurat_obj)
    
    if (length(object_cells) != length(result_cells)) {
      stop(
        "Cell count mismatch between `object` (", length(object_cells),
        ") and `scice_results` (", length(result_cells), ")."
      )
    }
    
    if (!setequal(object_cells, result_cells)) {
      stop("Cell identifiers in `object` and `scice_results$cell_names` do not match.")
    }
    
    for (i in 2:ncol(output_df)) {
      col_name <- names(output_df)[i]
      label_map <- stats::setNames(output_df[[i]], output_df$cell_id)
      aligned_values <- unname(label_map[object_cells])
      if (any(is.na(aligned_values))) {
        stop("Failed to align cluster labels to Seurat metadata for column: ", col_name)
      }
      seurat_obj@meta.data[[col_name]] <- aligned_values
    }
    
    return(seurat_obj)
  } else {
    return(output_df)
  }
}

#' Extract consistent clusters summary
#'
#' @description
#' Provides a summary of consistent clusters found by scICE analysis, including
#' cluster numbers, IC scores, and other relevant statistics.
#'
#' @details
#' The summary includes per-cluster median IC and optimization metadata, and
#' adds bootstrap confidence intervals when \code{ic_vec} is available.
#'
#' For newer result objects containing exclusion metadata, aggregate counts
#' (tested/excluded/inconsistent/consistent) are attached as attributes to the
#' returned data frame.
#'
#' @param scice_results Results object from scICE_clustering function
#' @param threshold IC threshold for determining consistent clusters (default: 1.005)
#'
#' @return Data frame summarizing consistent clusters
#' @export
extract_consistent_clusters <- function(scice_results, threshold = 1.005) {
  
  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }
  
  # Handle both old results (without excluded info) and new results (with excluded info)
  has_exclusion_info <- !is.null(scice_results$excluded)
  
  # Find consistent cluster numbers, excluding already-excluded clusters
  if (has_exclusion_info) {
    # For new format: only consider non-excluded clusters with valid IC scores
    valid_indices <- which(!scice_results$excluded & !is.na(scice_results$ic))
    consistent_indices <- valid_indices[scice_results$ic[valid_indices] < threshold]
  } else {
    # For old format: use all available clusters
    consistent_indices <- which(!is.na(scice_results$ic) & scice_results$ic < threshold)
  }
  
  if (length(consistent_indices) == 0) {
    if (has_exclusion_info) {
      excluded_clusters <- scice_results$n_cluster[scice_results$excluded]
      if (length(excluded_clusters) > 0) {
        warning(paste("No cluster numbers meet the consistency threshold.",
                     "Excluded clusters:", paste(excluded_clusters, collapse = ", ")))
      } else {
        warning("No cluster numbers meet the consistency threshold")
      }
    } else {
      warning("No cluster numbers meet the consistency threshold")
    }
    return(data.frame())
  }
  
  # Create summary data frame
  summary_df <- data.frame(
    cluster_number = scice_results$n_cluster[consistent_indices],
    ic_median = scice_results$ic[consistent_indices],
    resolution_parameter = sapply(scice_results$gamma[consistent_indices], function(x) x),
    n_iterations = scice_results$n_iter[consistent_indices],
    stringsAsFactors = FALSE
  )
  
  # Add IC confidence intervals if bootstrap data available
  if (!is.null(scice_results$ic_vec) && length(scice_results$ic_vec) > 0) {
    ic_ci_lower <- sapply(consistent_indices, function(i) {
      quantile(scice_results$ic_vec[[i]], 0.025)
    })
    ic_ci_upper <- sapply(consistent_indices, function(i) {
      quantile(scice_results$ic_vec[[i]], 0.975)
    })
    
    summary_df$ic_ci_lower <- ic_ci_lower
    summary_df$ic_ci_upper <- ic_ci_upper
  }
  
  # Sort by cluster number
  summary_df <- summary_df[order(summary_df$cluster_number), ]
  rownames(summary_df) <- NULL
  
  # Add summary information as attributes
  if (has_exclusion_info) {
    excluded_count <- sum(scice_results$excluded)
    tested_count <- length(scice_results$n_cluster)
    consistent_count <- nrow(summary_df)
    
    attr(summary_df, "summary_info") <- list(
      total_tested = tested_count,
      excluded = excluded_count,
      inconsistent = tested_count - excluded_count - consistent_count,
      consistent = consistent_count,
      threshold_used = threshold
    )
  }
  
  return(summary_df)
}

#' Plot clustering stability across resolution parameters
#'
#' @description
#' Visualizes stability diagnostics from an \code{scICE} result object.
#'
#' @details
#' Two modes are supported:
#' \itemize{
#'   \item \code{cluster_number = NULL}: scatter/line view of IC vs cluster number
#'   across all tested solutions;
#'   \item \code{cluster_number = k}: histogram of bootstrap IC scores for one
#'   specific cluster number.
#' }
#'
#' @param scice_results Results object from scICE_clustering function
#' @param cluster_number Specific cluster number to plot (optional)
#' @return ggplot object
#' @export
plot_stability <- function(scice_results, cluster_number = NULL) {
  
  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }
  
  if (is.null(cluster_number)) {
    # Plot all cluster numbers
    plot_data <- data.frame(
      cluster_number = scice_results$n_cluster,
      ic_score = scice_results$ic,
      resolution = sapply(scice_results$gamma, function(x) x),
      stringsAsFactors = FALSE
    )
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cluster_number, y = ic_score)) +
      ggplot2::geom_point(ggplot2::aes(color = resolution), size = 3) +
      ggplot2::geom_line(alpha = 0.6) +
      ggplot2::scale_color_gradient(name = "Resolution", low = "blue", high = "red") +
      ggplot2::scale_x_continuous(name = "Number of Clusters", 
                                 breaks = unique(plot_data$cluster_number)) +
      ggplot2::scale_y_continuous(name = "IC Score") +
      ggplot2::labs(title = "Clustering Stability Across Cluster Numbers") +
      ggplot2::theme_minimal()
      
  } else {
    # Plot specific cluster number
    if (!(cluster_number %in% scice_results$n_cluster)) {
      stop(paste("Cluster number", cluster_number, "not found in results"))
    }
    
    idx <- which(scice_results$n_cluster == cluster_number)
    ic_bootstrap <- scice_results$ic_vec[[idx]]
    
    plot_data <- data.frame(
      ic_score = ic_bootstrap,
      bootstrap_sample = 1:length(ic_bootstrap)
    )
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = ic_score)) +
      ggplot2::geom_histogram(bins = 30, fill = "lightblue", alpha = 0.7) +
      ggplot2::geom_vline(ggplot2::aes(xintercept = median(ic_bootstrap)), 
                        color = "red", linetype = "dashed") +
      ggplot2::scale_x_continuous(name = "IC Score") +
      ggplot2::scale_y_continuous(name = "Frequency") +
      ggplot2::labs(title = paste("IC Score Distribution for", cluster_number, "Clusters"),
                   subtitle = paste("Median IC:", round(median(ic_bootstrap), 4))) +
      ggplot2::theme_minimal()
  }
  
  return(p)
} 
