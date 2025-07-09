#' @import ggplot2
#' @importFrom stats quantile median
#' @importFrom methods inherits
NULL

#' Plot Inconsistency (IC) scores across different cluster numbers
#'
#' @description
#' Creates a boxplot showing the distribution of IC scores for each cluster number 
#' tested. Lower IC scores indicate more consistent clustering results.
#'
#' @param scice_results Results object from scICE_clustering function
#' @param threshold IC threshold line to display (default: 1.005)
#' @param figure_size Figure size as c(width, height) (default: c(10, 6))
#' @param title Plot title (default: "Clustering Consistency Analysis")
#' @param show_threshold Whether to show threshold line (default: TRUE)
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
                   title = "Clustering Consistency Analysis", show_threshold = TRUE) {

  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }

  if (length(scice_results$ic_vec) == 0) {
    stop("No IC results found in scICE object")
  }

  # Prepare data for plotting
  plot_data <- data.frame(
    cluster_number = rep(scice_results$n_cluster,
                        sapply(scice_results$ic_vec, length)),
    ic_score = unlist(scice_results$ic_vec)
  )

  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = factor(cluster_number), y = ic_score)) +
    ggplot2::geom_boxplot(ggplot2::aes(group = cluster_number),
                 outlier.alpha = 0.6,
                 fill = "lightblue",
                 alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.4, size = 0.8, height = 0) +
    ggplot2::scale_x_discrete(name = "Number of Clusters") +
    ggplot2::scale_y_continuous(name = "Inconsistency (IC) Score",
                      limits = c(0.99, max(plot_data$ic_score) * 1.05)) +
    ggplot2::labs(title = title,
         subtitle = paste("Lower IC scores indicate more consistent clustering.",
                         "Threshold:", threshold)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 14, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
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

  # Highlight consistent clusters
  consistent_clusters <- scice_results$consistent_clusters
  if (length(consistent_clusters) > 0) {
    consistent_data <- plot_data[plot_data$cluster_number %in% consistent_clusters, ]
    p <- p + ggplot2::geom_boxplot(data = consistent_data,
                         ggplot2::aes(x = factor(cluster_number), y = ic_score),
                         fill = "lightgreen",
                         alpha = 0.8)
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
#' @param scice_results Results object from scICE_clustering function
#' @param threshold IC threshold for determining consistent clusters (default: 1.005)
#' @param return_seurat Whether to return results as Seurat metadata (default: FALSE)
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
#' seurat_obj <- get_robust_labels(scice_results, return_seurat = TRUE)
#' }
get_robust_labels <- function(scice_results, threshold = 1.005, return_seurat = FALSE) {
  
  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }
  
  # Find consistent cluster numbers
  consistent_indices <- which(scice_results$ic < threshold)
  
  if (length(consistent_indices) == 0) {
    warning("No cluster numbers meet the consistency threshold")
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
    
    # Convert to 1-based indexing for R
    cluster_labels <- cluster_labels + 1
    
    # Add to data frame
    col_name <- paste0("clusters_", cluster_num)
    output_df[[col_name]] <- cluster_labels
  }
  
  if (return_seurat) {
    if (is.null(scice_results$seurat_object)) {
      stop("No Seurat object found in scICE results")
    }
    
    # Add metadata to Seurat object
    seurat_obj <- scice_results$seurat_object
    
    for (i in 2:ncol(output_df)) {
      col_name <- names(output_df)[i]
      seurat_obj@meta.data[[col_name]] <- output_df[[i]]
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
#' @param scice_results Results object from scICE_clustering function
#' @param threshold IC threshold for determining consistent clusters (default: 1.005)
#'
#' @return Data frame summarizing consistent clusters
#' @export
extract_consistent_clusters <- function(scice_results, threshold = 1.005) {
  
  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }
  
  # Find consistent cluster numbers
  consistent_indices <- which(scice_results$ic < threshold)
  
  if (length(consistent_indices) == 0) {
    warning("No cluster numbers meet the consistency threshold")
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
  
  return(summary_df)
}

#' Plot clustering stability across resolution parameters
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