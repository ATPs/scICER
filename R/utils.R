#' @import Seurat
#' @import SeuratObject
#' @importFrom parallel detectCores
NULL

#' Check if Seurat object is ready for scICER analysis
#'
#' @description
#' Validates that a Seurat object has the necessary preprocessing steps completed
#' for scICER analysis. Checks for normalization, variable features, scaling, PCA,
#' and neighbor graphs.
#'
#' @details
#' This helper is intended as a preflight check before expensive clustering runs.
#' It returns a structured report rather than throwing for missing preprocessing,
#' so users can decide whether to proceed or follow suggested commands.
#'
#' @param seurat_obj A Seurat object to validate
#' @param graph_name Name of the graph to check for. If \code{NULL}, uses the
#'   default assay SNN graph name \code{<DefaultAssay>_snn}. For backward
#'   compatibility, \code{"snn"} also falls back to that assay-specific graph
#'   when present.
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
check_seurat_ready <- function(seurat_obj, graph_name = NULL) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  if (!nzchar(assay_name) || !(assay_name %in% names(seurat_obj@assays))) {
    assay_name <- names(seurat_obj@assays)[[1]]
  }

  get_assay_data_safe <- function(layer_name) {
    layer_data <- tryCatch(
      SeuratObject::GetAssayData(seurat_obj, assay = assay_name, layer = layer_name),
      error = function(e) NULL
    )
    if (!is.null(layer_data)) {
      return(layer_data)
    }

    slot_data <- tryCatch(
      SeuratObject::GetAssayData(seurat_obj, assay = assay_name, slot = layer_name),
      error = function(e) NULL
    )
    if (!is.null(slot_data)) {
      return(slot_data)
    }

    assay <- seurat_obj@assays[[assay_name]]
    if (!is.null(assay) && layer_name %in% methods::slotNames(assay)) {
      return(methods::slot(assay, layer_name))
    }

    if (!is.null(assay) && "layers" %in% methods::slotNames(assay)) {
      assay_layers <- methods::slot(assay, "layers")
      if (layer_name %in% names(assay_layers)) {
        return(assay_layers[[layer_name]])
      }
    }

    NULL
  }

  default_graph_name <- paste0(assay_name, "_snn")
  resolved_graph_name <- graph_name
  if (is.null(resolved_graph_name)) {
    resolved_graph_name <- default_graph_name
  } else if (identical(resolved_graph_name, "snn") &&
             !(resolved_graph_name %in% names(seurat_obj@graphs)) &&
             default_graph_name %in% names(seurat_obj@graphs)) {
    resolved_graph_name <- default_graph_name
  }
  
  issues <- character()
  suggestions <- character()
  data_layer <- get_assay_data_safe("data")
  counts_layer <- get_assay_data_safe("counts")
  scale_layer <- get_assay_data_safe("scale.data")
  
  # Check for normalization
  has_normalized_data <- !is.null(data_layer) && length(data_layer) > 0
  data_matches_counts <- FALSE
  if (has_normalized_data && !is.null(counts_layer) && identical(dim(data_layer), dim(counts_layer))) {
    data_matches_counts <- isTRUE(all.equal(data_layer, counts_layer))
  }

  if (!has_normalized_data || data_matches_counts) {
    issues <- c(issues, "Data not normalized")
    suggestions <- c(suggestions, "Run: seurat_obj <- NormalizeData(seurat_obj)")
  }
  
  # Check for variable features
  if (length(SeuratObject::VariableFeatures(seurat_obj)) == 0) {
    issues <- c(issues, "No variable features found")
    suggestions <- c(suggestions, "Run: seurat_obj <- FindVariableFeatures(seurat_obj)")
  }
  
  # Check for scaling
  if (is.null(scale_layer) || nrow(scale_layer) == 0) {
    issues <- c(issues, "Data not scaled")
    suggestions <- c(suggestions, "Run: seurat_obj <- ScaleData(seurat_obj)")
  }
  
  # Check for PCA
  if (!"pca" %in% names(seurat_obj@reductions)) {
    issues <- c(issues, "PCA not computed")
    suggestions <- c(suggestions, "Run: seurat_obj <- RunPCA(seurat_obj)")
  }
  
  # Check for neighbor graph
  if (!(resolved_graph_name %in% names(seurat_obj@graphs))) {
    issues <- c(issues, paste("Graph", resolved_graph_name, "not found"))
    suggestions <- c(suggestions, "Run: seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)")
  }
  
  ready <- length(issues) == 0
  
  summary_msg <- if (ready) {
    "OK: Seurat object is ready for scICER analysis"
  } else {
    paste("ERROR: Seurat object needs preprocessing:", paste(issues, collapse = ", "))
  }
  
  return(list(
    ready = ready,
    issues = issues,
    suggestions = suggestions,
    summary = summary_msg,
    assay = assay_name,
    graph_name = resolved_graph_name
  ))
}

#' Get recommended parameters for scICER analysis based on dataset size
#'
#' @description
#' Provides parameter recommendations for scICER analysis based on the number of cells
#' in the dataset. Helps optimize performance and accuracy for different dataset sizes.
#'
#' @details
#' Recommendations are generated from dataset-size buckets and analysis depth:
#' \itemize{
#'   \item \code{quick}: faster screening with fewer trials/bootstraps;
#'   \item \code{standard}: balanced defaults;
#'   \item \code{thorough}: higher robustness with longer runtime.
#' }
#'
#' Returned values are intended as starting points and can be adjusted based on
#' graph density, hardware limits, and desired precision.
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
#' @details
#' The report is formatted as a character vector so it can be printed directly
#' (\code{cat(paste(summary, collapse = "\\n"))}) or written to file with
#' \code{writeLines()}.
#'
#' It includes:
#' \itemize{
#'   \item overall run metadata and threshold used;
#'   \item IC distribution summary;
#'   \item per-cluster detail for consistent solutions;
#'   \item interpretation notes for downstream reporting.
#' }
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

  package_version <- tryCatch(
    as.character(utils::packageVersion("scICER")),
    error = function(e) NA_character_
  )
  if (!is.character(package_version) || length(package_version) != 1L || is.na(package_version)) {
    package_version <- "unknown"
  }

  format_cluster_values <- function(values) {
    values <- sort(unique(as.integer(values[is.finite(values)])))
    if (length(values) == 0L) {
      return("none")
    }
    if (length(values) == 1L) {
      return(as.character(values))
    }
    if (identical(values, seq.int(min(values), max(values)))) {
      return(paste0(min(values), "-", max(values), " (", length(values), " values)"))
    }
    paste0(paste(values, collapse = ", "), " (", length(values), " values)")
  }

  analysis_mode <- if (!is.null(scice_results$analysis_mode) &&
                       length(scice_results$analysis_mode) == 1L &&
                       nzchar(scice_results$analysis_mode)) {
    scice_results$analysis_mode
  } else {
    "cluster_range"
  }

  target_diagnostics <- scice_results$target_diagnostics
  has_target_diagnostics <- is.data.frame(target_diagnostics)
  search_diagnostics <- scice_results$resolution_search_diagnostics
  has_search_diagnostics <- is.data.frame(search_diagnostics)
  excluded_targets <- if (has_target_diagnostics &&
                          all(c("searched_target_cluster", "excluded") %in% colnames(target_diagnostics))) {
    as.integer(target_diagnostics$searched_target_cluster[target_diagnostics$excluded])
  } else {
    integer(0)
  }
  
  # Header
  summary <- c(
    "=====================================",
    "       scICER Analysis Summary       ",
    "=====================================",
    "",
    paste("Analysis Date:", Sys.time()),
    paste("Package Version: scICER", package_version),
    ""
  )
  
  # Analysis parameters
  parameter_lines <- c(
    "ANALYSIS PARAMETERS:",
    paste("- Analysis mode:", analysis_mode)
  )
  if (identical(analysis_mode, "resolution")) {
    manual_resolutions <- if (!is.null(scice_results$resolution_input) &&
                              length(scice_results$resolution_input) > 0) {
      paste(signif(scice_results$resolution_input, 6), collapse = ", ")
    } else {
      "none"
    }
    parameter_lines <- c(
      parameter_lines,
      paste("- Manual resolutions evaluated:", manual_resolutions),
      paste("- Returned cluster numbers:", paste(scice_results$n_cluster, collapse = ", "))
    )
  } else {
    requested_range_label <- if (!is.null(scice_results$requested_cluster_range)) {
      format_cluster_values(scice_results$requested_cluster_range)
    } else {
      "none"
    }
    searched_range_label <- if (!is.null(scice_results$searched_target_cluster_range)) {
      format_cluster_values(scice_results$searched_target_cluster_range)
    } else if (has_target_diagnostics && "searched_target_cluster" %in% colnames(target_diagnostics)) {
      format_cluster_values(target_diagnostics$searched_target_cluster)
    } else {
      "none"
    }
    returned_final_label <- format_cluster_values(scice_results$n_cluster)
    parameter_lines <- c(
      parameter_lines,
      paste("- Requested final cluster range:", requested_range_label),
      paste("- Searched target cluster range:", searched_range_label),
      paste("- Returned final cluster numbers:", returned_final_label),
      paste("- Total searched targets:", if (has_target_diagnostics) nrow(target_diagnostics) else length(scice_results$n_cluster)),
      paste("- Total returned final clusters:", length(scice_results$n_cluster)),
      paste("- Shared gamma probes:", if (has_search_diagnostics) nrow(search_diagnostics) else 0),
      paste("- Shared gamma sweep coverage complete:", if (!is.null(scice_results$search_coverage_complete)) isTRUE(scice_results$search_coverage_complete) else "unknown"),
      paste("- Coverage complete:", if (!is.null(scice_results$coverage_complete)) isTRUE(scice_results$coverage_complete) else "unknown"),
      paste("- Excluded searched targets:", length(excluded_targets))
    )
    if (!is.null(scice_results$uncovered_targets) && length(scice_results$uncovered_targets) > 0L) {
      parameter_lines <- c(
        parameter_lines,
        paste("- Uncovered requested final targets:", paste(scice_results$uncovered_targets, collapse = ", "))
      )
    }
    if (!is.null(scice_results$plateau_stop)) {
      parameter_lines <- c(
        parameter_lines,
        paste("- Plateau stop triggered:", isTRUE(scice_results$plateau_stop))
      )
    }
  }
  if (!is.null(scice_results$best_cluster) && length(scice_results$best_cluster) == 1L &&
      !is.na(scice_results$best_cluster)) {
    parameter_lines <- c(
      parameter_lines,
      paste(
        "- Lowest-IC solution:",
        scice_results$best_cluster,
        "clusters at resolution",
        signif(scice_results$best_resolution, 6)
      )
    )
  }
  summary <- c(summary, parameter_lines, paste("- IC threshold:", threshold), "")
  
  # Results overview
  consistent_clusters <- which(is.finite(scice_results$ic) & scice_results$ic < threshold)
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
  finite_ic <- scice_results$ic[is.finite(scice_results$ic)]
  if (length(finite_ic) > 0) {
    ic_stats <- c(
      paste("- Best IC score:", round(min(finite_ic), 4)),
      paste("- Worst IC score:", round(max(finite_ic), 4)),
      paste("- Mean IC score:", round(mean(finite_ic), 4)),
      paste("- Median IC score:", round(median(finite_ic), 4))
    )
  } else {
    ic_stats <- c(
      "- Best IC score: NA",
      "- Worst IC score: NA",
      "- Mean IC score: NA",
      "- Median IC score: NA"
    )
  }
  
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
      if (identical(analysis_mode, "resolution")) {
        "- Trying additional manual resolution values"
      } else {
        "- Expanding the requested final cluster range"
      },
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
              "  seurat_obj <- get_robust_labels(results, return_seurat = TRUE, object = seurat_obj)",
              "=====================================")
  
  return(summary)
} 
