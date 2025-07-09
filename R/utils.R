#' Check if Seurat object is ready for scICER analysis
#'
#' @description
#' Validates that a Seurat object has the necessary preprocessing steps completed
#' for scICER analysis. Checks for normalization, variable features, scaling, PCA,
#' and neighbor graphs.
#'
#' @param seurat_obj A Seurat object to validate
#' @param graph_name Name of the graph to check for (default: "snn")
#' @return List with validation results and suggestions
#' @export
#' @examples
#' \dontrun{
#' library(Seurat)
#' data("pbmc_small")
#' 
#' # Check readiness
#' validation <- check_seurat_ready(pbmc_small)
#' print(validation$summary)
#' 
#' if (!validation$ready) {
#'   cat("Suggestions:\n")
#'   cat(paste(validation$suggestions, collapse = "\n"))
#' }
#' }
check_seurat_ready <- function(seurat_obj, graph_name = "snn") {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  issues <- character()
  suggestions <- character()
  
  # Check for normalization
  if (!"data" %in% names(seurat_obj@assays$RNA@layers) && 
      is.null(seurat_obj@assays$RNA@data) || 
      all(seurat_obj@assays$RNA@data == seurat_obj@assays$RNA@counts)) {
    issues <- c(issues, "Data not normalized")
    suggestions <- c(suggestions, "Run: seurat_obj <- NormalizeData(seurat_obj)")
  }
  
  # Check for variable features
  if (length(VariableFeatures(seurat_obj)) == 0) {
    issues <- c(issues, "No variable features found")
    suggestions <- c(suggestions, "Run: seurat_obj <- FindVariableFeatures(seurat_obj)")
  }
  
  # Check for scaling
  if (is.null(seurat_obj@assays$RNA@scale.data) || nrow(seurat_obj@assays$RNA@scale.data) == 0) {
    issues <- c(issues, "Data not scaled")
    suggestions <- c(suggestions, "Run: seurat_obj <- ScaleData(seurat_obj)")
  }
  
  # Check for PCA
  if (!"pca" %in% names(seurat_obj@reductions)) {
    issues <- c(issues, "PCA not computed")
    suggestions <- c(suggestions, "Run: seurat_obj <- RunPCA(seurat_obj)")
  }
  
  # Check for neighbor graph
  if (!(graph_name %in% names(seurat_obj@graphs))) {
    issues <- c(issues, paste("Graph", graph_name, "not found"))
    suggestions <- c(suggestions, "Run: seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)")
  }
  
  ready <- length(issues) == 0
  
  summary_msg <- if (ready) {
    "✓ Seurat object is ready for scICER analysis"
  } else {
    paste("✗ Seurat object needs preprocessing:", paste(issues, collapse = ", "))
  }
  
  return(list(
    ready = ready,
    issues = issues,
    suggestions = suggestions,
    summary = summary_msg
  ))
}

#' Get recommended parameters for scICER analysis based on dataset size
#'
#' @description
#' Provides parameter recommendations for scICER analysis based on the number of cells
#' in the dataset. Helps optimize performance and accuracy for different dataset sizes.
#'
#' @param n_cells Number of cells in the dataset
#' @param analysis_type Type of analysis: "quick", "standard", or "thorough" (default: "standard")
#' @return List of recommended parameters
#' @export
#' @examples
#' \dontrun{
#' # Get recommendations for different dataset sizes
#' small_params <- get_recommended_parameters(1000)
#' medium_params <- get_recommended_parameters(10000)
#' large_params <- get_recommended_parameters(100000)
#' 
#' # Use recommendations
#' results <- scICE_clustering(
#'   seurat_obj, 
#'   cluster_range = small_params$cluster_range,
#'   n_workers = small_params$n_workers,
#'   n_trials = small_params$n_trials,
#'   n_bootstrap = small_params$n_bootstrap
#' )
#' }
get_recommended_parameters <- function(n_cells, analysis_type = "standard") {
  
  if (!analysis_type %in% c("quick", "standard", "thorough")) {
    stop("analysis_type must be 'quick', 'standard', or 'thorough'")
  }
  
  # Base parameters
  if (n_cells < 1000) {
    category <- "small"
    base_cluster_range <- 2:8
    base_n_workers <- 2
    base_n_trials <- 10
    base_n_bootstrap <- 50
    objective_function <- "modularity"
  } else if (n_cells < 10000) {
    category <- "medium"
    base_cluster_range <- 2:12
    base_n_workers <- 4
    base_n_trials <- 12
    base_n_bootstrap <- 75
    objective_function <- "CPM"
  } else if (n_cells < 50000) {
    category <- "large"
    base_cluster_range <- 3:20
    base_n_workers <- 6
    base_n_trials <- 10
    base_n_bootstrap <- 50
    objective_function <- "CPM"
  } else {
    category <- "very_large"
    base_cluster_range <- 5:25
    base_n_workers <- 8
    base_n_trials <- 8
    base_n_bootstrap <- 30
    objective_function <- "CPM"
  }
  
  # Adjust based on analysis type
  multiplier <- switch(analysis_type,
                      "quick" = 0.5,
                      "standard" = 1.0,
                      "thorough" = 1.5)
  
  params <- list(
    dataset_category = category,
    n_cells = n_cells,
    cluster_range = base_cluster_range,
    n_workers = min(base_n_workers, parallel::detectCores() - 1),
    n_trials = max(5, round(base_n_trials * multiplier)),
    n_bootstrap = max(20, round(base_n_bootstrap * multiplier)),
    objective_function = objective_function,
    ic_threshold = 1.005,
    analysis_type = analysis_type
  )
  
  # Add performance notes
  estimated_time <- switch(category,
                          "small" = "2-5 minutes",
                          "medium" = "5-15 minutes", 
                          "large" = "15-45 minutes",
                          "very_large" = "30-120 minutes")
  
  params$estimated_runtime <- paste(estimated_time, "for", analysis_type, "analysis")
  
  # Add recommendations
  recommendations <- c(
    paste("Dataset category:", category, "(", n_cells, "cells)"),
    paste("Recommended analysis type:", analysis_type),
    paste("Estimated runtime:", params$estimated_runtime)
  )
  
  if (category %in% c("large", "very_large")) {
    recommendations <- c(recommendations,
                        "Consider installing 'leiden' package for better performance",
                        "Use SSD storage for faster I/O if possible")
  }
  
  params$recommendations <- recommendations
  
  return(params)
}

#' Create a summary report of scICER results
#'
#' @description
#' Generates a comprehensive text summary of scICER analysis results,
#' including consistent clusters, parameter settings, and interpretation guidance.
#'
#' @param scice_results Results object from scICE_clustering function
#' @param threshold IC threshold used for consistency (default: 1.005)
#' @return Character vector with formatted summary
#' @export
#' @examples
#' \dontrun{
#' # Run analysis and generate summary
#' results <- scICE_clustering(seurat_obj, cluster_range = 2:10)
#' summary_report <- create_results_summary(results)
#' cat(paste(summary_report, collapse = "\n"))
#' 
#' # Save to file
#' writeLines(summary_report, "scICER_analysis_summary.txt")
#' }
create_results_summary <- function(scice_results, threshold = 1.005) {
  
  if (!inherits(scice_results, "scICE")) {
    stop("Input must be a scICE results object")
  }
  
  # Header
  summary <- c(
    "=====================================",
    "       scICER Analysis Summary       ",
    "=====================================",
    "",
    paste("Analysis Date:", Sys.time()),
    paste("Package Version: scICER 1.0.0"),
    ""
  )
  
  # Analysis parameters
  summary <- c(summary,
              "ANALYSIS PARAMETERS:",
              paste("- Cluster range tested:", paste(range(scice_results$n_cluster), collapse = "-")),
              paste("- Total cluster numbers:", length(scice_results$n_cluster)),
              paste("- IC threshold:", threshold),
              "")
  
  # Results overview
  consistent_clusters <- which(scice_results$ic < threshold)
  n_consistent <- length(consistent_clusters)
  
  summary <- c(summary,
              "RESULTS OVERVIEW:",
              paste("- Consistent cluster numbers found:", n_consistent),
              if (n_consistent > 0) {
                paste("- Consistent cluster numbers:", paste(scice_results$n_cluster[consistent_clusters], collapse = ", "))
              } else {
                "- No cluster numbers meet the consistency threshold"
              },
              "")
  
  # IC scores summary
  ic_stats <- c(
    paste("- Best IC score:", round(min(scice_results$ic), 4)),
    paste("- Worst IC score:", round(max(scice_results$ic), 4)),
    paste("- Mean IC score:", round(mean(scice_results$ic), 4)),
    paste("- Median IC score:", round(median(scice_results$ic), 4))
  )
  
  summary <- c(summary,
              "IC SCORE STATISTICS:",
              ic_stats,
              "")
  
  # Detailed results for consistent clusters
  if (n_consistent > 0) {
    summary <- c(summary, "CONSISTENT CLUSTERS DETAIL:")
    
    for (i in consistent_clusters) {
      cluster_num <- scice_results$n_cluster[i]
      ic_score <- scice_results$ic[i]
      resolution <- scice_results$gamma[[i]]
      iterations <- scice_results$n_iter[i]
      
      summary <- c(summary,
                  paste("  ", cluster_num, "clusters:"),
                  paste("    - IC score:", round(ic_score, 4)),
                  paste("    - Resolution:", round(resolution, 6)),
                  paste("    - Iterations:", iterations))
    }
    summary <- c(summary, "")
  }
  
  # Interpretation guidance
  summary <- c(summary,
              "INTERPRETATION GUIDANCE:",
              "- IC < 1.005: Highly consistent (recommended)",
              "- IC 1.005-1.01: Moderately consistent",
              "- IC 1.01-1.02: Low consistency (use with caution)",
              "- IC > 1.02: Poor consistency (likely overfitting)",
              "")
  
  # Recommendations
  recommendations <- character()
  
  if (n_consistent == 0) {
    recommendations <- c(
      "No clusters meet the strict threshold. Consider:",
      "- Using a more lenient threshold (e.g., 1.01)",
      "- Expanding the cluster range",
      "- Checking data preprocessing steps"
    )
  } else if (n_consistent == 1) {
    best_cluster <- scice_results$n_cluster[consistent_clusters[1]]
    recommendations <- c(
      paste("Single consistent clustering found with", best_cluster, "clusters."),
      "This suggests a clear optimal cluster number for your data."
    )
  } else {
    best_cluster <- scice_results$n_cluster[consistent_clusters[which.min(scice_results$ic[consistent_clusters])]]
    recommendations <- c(
      paste("Multiple consistent clusterings found."),
      paste("Most consistent:", best_cluster, "clusters"),
      "Consider biological relevance when choosing final clustering."
    )
  }
  
  summary <- c(summary,
              "RECOMMENDATIONS:",
              paste("-", recommendations),
              "")
  
  # Footer
  summary <- c(summary,
              "=====================================",
              "For detailed visualization:",
              "  plot_ic(results)",
              "  plot_stability(results)",
              "",
              "To extract clustering labels:",
              "  seurat_obj <- get_robust_labels(results, return_seurat = TRUE)",
              "=====================================")
  
  return(summary)
} 