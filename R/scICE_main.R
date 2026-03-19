#' @import Seurat
#' @import SeuratObject
#' @import parallel
#' @importFrom igraph make_empty_graph add_edges E vcount ecount cluster_leiden membership
NULL

emit_scice_notice <- function(message_text, verbose = FALSE) {
  if (isTRUE(verbose)) {
    scice_message(message_text)
  } else {
    message(message_text)
  }
}

append_target_results_dt <- function(existing_dt, new_dt) {
  if (is.null(new_dt) || nrow(new_dt) == 0L) {
    return(if (is.null(existing_dt)) empty_target_results_dt() else existing_dt)
  }
  if (is.null(existing_dt) || nrow(existing_dt) == 0L) {
    return(data.table::copy(new_dt))
  }
  data.table::rbindlist(list(existing_dt, new_dt), fill = TRUE)
}

build_target_diagnostics_df <- function(target_results_dt, requested_cluster_range, gamma_dict = NULL) {
  requested_cluster_range <- as.integer(requested_cluster_range)
  searched_targets <- requested_cluster_range

  lookup_bounds <- function(target_cluster) {
    if (is.null(gamma_dict) || !(as.character(target_cluster) %in% names(gamma_dict))) {
      return(c(NA_real_, NA_real_))
    }
    bounds <- as.numeric(gamma_dict[[as.character(target_cluster)]])
    if (length(bounds) != 2L) {
      return(c(NA_real_, NA_real_))
    }
    bounds
  }
  gamma_bounds <- t(vapply(
    searched_targets,
    lookup_bounds,
    numeric(2)
  ))

  diagnostics_df <- data.frame(
    requested_target_cluster = searched_targets,
    searched_target_cluster = searched_targets,
    requested_by_user = TRUE,
    gamma_left = as.numeric(gamma_bounds[, 1]),
    gamma_right = as.numeric(gamma_bounds[, 2]),
    gamma = rep(NA_real_, length(searched_targets)),
    ic = rep(NA_real_, length(searched_targets)),
    effective_cluster_median = rep(NA_real_, length(searched_targets)),
    raw_cluster_median = rep(NA_real_, length(searched_targets)),
    final_cluster_median = rep(NA_real_, length(searched_targets)),
    best_labels_raw_cluster_count = rep(NA_integer_, length(searched_targets)),
    best_labels_final_cluster_count = rep(NA_integer_, length(searched_targets)),
    excluded = rep(FALSE, length(searched_targets)),
    exclusion_reason = rep("missing_target_result", length(searched_targets)),
    selected_main_result = rep(FALSE, length(searched_targets)),
    result_status = rep("missing_target_result", length(searched_targets)),
    stringsAsFactors = FALSE
  )

  if (is.null(target_results_dt) || nrow(target_results_dt) == 0L) {
    rownames(diagnostics_df) <- NULL
    return(diagnostics_df)
  }

  target_results_dt <- data.table::as.data.table(data.table::copy(target_results_dt))
  target_results_dt[, source_target_cluster := as.integer(source_target_cluster)]
  if (!"result_status" %in% colnames(target_results_dt)) {
    target_results_dt[, result_status := data.table::fifelse(
      selected_main_result,
      "selected_main_result",
      data.table::fifelse(excluded, as.character(exclusion_reason), "deduplicated")
    )]
  }
  target_results_dt[, diagnostics_priority := data.table::fcase(
    result_status == "selected_main_result", 1L,
    !excluded, 2L,
    excluded, 3L,
    default = 4L
  )]
  data.table::setorder(target_results_dt, source_target_cluster, diagnostics_priority, ic, gamma)
  target_results_dt <- target_results_dt[, .SD[1], by = source_target_cluster]

  matched_idx <- match(diagnostics_df$searched_target_cluster, target_results_dt$source_target_cluster)
  has_match <- !is.na(matched_idx)

  diagnostics_df$gamma[has_match] <- as.numeric(target_results_dt$gamma[matched_idx[has_match]])
  diagnostics_df$ic[has_match] <- as.numeric(target_results_dt$ic[matched_idx[has_match]])
  diagnostics_df$effective_cluster_median[has_match] <- as.numeric(target_results_dt$effective_cluster_median[matched_idx[has_match]])
  diagnostics_df$raw_cluster_median[has_match] <- as.numeric(target_results_dt$raw_cluster_median[matched_idx[has_match]])
  if ("final_cluster_median" %in% colnames(target_results_dt)) {
    diagnostics_df$final_cluster_median[has_match] <- as.numeric(target_results_dt$final_cluster_median[matched_idx[has_match]])
  }
  diagnostics_df$best_labels_raw_cluster_count[has_match] <- as.integer(target_results_dt$best_labels_raw_cluster_count[matched_idx[has_match]])
  diagnostics_df$best_labels_final_cluster_count[has_match] <- as.integer(target_results_dt$best_labels_final_cluster_count[matched_idx[has_match]])
  diagnostics_df$excluded[has_match] <- as.logical(target_results_dt$excluded[matched_idx[has_match]])
  diagnostics_df$exclusion_reason[has_match] <- as.character(target_results_dt$exclusion_reason[matched_idx[has_match]])
  diagnostics_df$selected_main_result[has_match] <- as.logical(target_results_dt$selected_main_result[matched_idx[has_match]])
  diagnostics_df$result_status[has_match] <- as.character(target_results_dt$result_status[matched_idx[has_match]])

  missing_match <- !has_match
  diagnostics_df$excluded[missing_match] <- TRUE
  diagnostics_df$exclusion_reason[missing_match] <- ifelse(
    is.na(diagnostics_df$gamma_left[missing_match]),
    "resolution_search_failed",
    "optimization_admission_failed"
  )
  diagnostics_df$result_status[missing_match] <- diagnostics_df$exclusion_reason[missing_match]

  diagnostics_df <- data.frame(
    requested_target_cluster = as.integer(diagnostics_df$requested_target_cluster),
    searched_target_cluster = as.integer(diagnostics_df$searched_target_cluster),
    requested_by_user = as.logical(diagnostics_df$requested_by_user),
    gamma_left = as.numeric(diagnostics_df$gamma_left),
    gamma_right = as.numeric(diagnostics_df$gamma_right),
    gamma = as.numeric(diagnostics_df$gamma),
    ic = as.numeric(diagnostics_df$ic),
    effective_cluster_median = as.numeric(diagnostics_df$effective_cluster_median),
    raw_cluster_median = as.numeric(diagnostics_df$raw_cluster_median),
    final_cluster_median = as.numeric(diagnostics_df$final_cluster_median),
    best_labels_raw_cluster_count = as.integer(diagnostics_df$best_labels_raw_cluster_count),
    best_labels_final_cluster_count = as.integer(diagnostics_df$best_labels_final_cluster_count),
    excluded = as.logical(diagnostics_df$excluded),
    exclusion_reason = as.character(diagnostics_df$exclusion_reason),
    selected_main_result = as.logical(diagnostics_df$selected_main_result),
    result_status = as.character(diagnostics_df$result_status),
    stringsAsFactors = FALSE
  )
  rownames(diagnostics_df) <- NULL
  diagnostics_df
}

finalize_cluster_range_results <- function(target_results_dt, requested_cluster_range,
                                           searched_target_cluster_range,
                                           search_coverage_complete,
                                           gamma_dict = NULL,
                                           resolution_search_diagnostics = NULL,
                                           plateau_stop = FALSE,
                                           search_uncovered_targets = integer(0),
                                           discovered_upper_gamma = NULL,
                                           upper_cap_stop_reason = NULL,
                                           coarse_probe_count = NULL) {
  rekeyed_results <- rekey_target_results_by_final_cluster(target_results_dt)
  result <- cluster_results_dt_to_list(rekeyed_results$main_results_dt)
  result$requested_cluster_range <- as.integer(requested_cluster_range)
  result$searched_target_cluster_range <- as.integer(searched_target_cluster_range)
  result$search_coverage_complete <- isTRUE(search_coverage_complete)
  result$resolution_search_diagnostics <- resolution_search_diagnostics
  result$plateau_stop <- isTRUE(plateau_stop)
  result$search_uncovered_targets <- as.integer(search_uncovered_targets)
  result$discovered_upper_gamma <- if (length(discovered_upper_gamma) == 1L) as.numeric(discovered_upper_gamma) else NULL
  result$upper_cap_stop_reason <- if (length(upper_cap_stop_reason) == 1L) as.character(upper_cap_stop_reason) else NULL
  result$coarse_probe_count <- if (length(coarse_probe_count) == 1L) as.integer(coarse_probe_count) else NULL
  result$uncovered_targets <- as.integer(setdiff(requested_cluster_range, as.integer(result$n_cluster)))
  result$coverage_complete <- length(result$uncovered_targets) == 0L
  result$target_diagnostics <- build_target_diagnostics_df(
    rekeyed_results$target_results_dt,
    requested_cluster_range = requested_cluster_range,
    gamma_dict = gamma_dict
  )
  result
}

run_cluster_range_mode <- function(igraph_obj, requested_cluster_range, n_workers,
                                   n_trials, n_bootstrap, seed, beta,
                                   n_iterations, max_iterations,
                                   objective_function, remove_threshold,
                                   snn_graph, min_cluster_size,
                                   resolution_tolerance, verbose,
                                   runtime_context = NULL) {
  requested_cluster_range <- sort(unique(as.integer(requested_cluster_range)))
  searched_target_cluster_range <- requested_cluster_range
  if (exists("clear_clustering_cache") && is.function(clear_clustering_cache)) {
    clear_clustering_cache()
  }

  search_start_g <- if (objective_function == "modularity") {
    -13
  } else {
    lower <- log(resolution_tolerance)
    if (lower < -20) lower <- -20
    lower
  }
  search_end_g <- 20

  if (verbose) {
    scice_message(
      paste(
        "  Cluster-range mode uses a shared gamma sweep for requested final targets:",
        paste(requested_cluster_range, collapse = ", ")
      )
    )
  }

  gamma_dict <- find_resolution_ranges(
    igraph_obj = igraph_obj,
    cluster_range = requested_cluster_range,
    start_g = search_start_g,
    end_g = search_end_g,
    objective_function = objective_function,
    resolution_tolerance = resolution_tolerance,
    n_workers = n_workers,
    verbose = verbose,
    seed = seed,
    snn_graph = snn_graph,
    min_cluster_size = min_cluster_size,
    in_parallel_context = FALSE,
    runtime_context = runtime_context
  )
  resolution_search_diagnostics <- attr(gamma_dict, "resolution_search_diagnostics")
  coverage_complete <- isTRUE(attr(gamma_dict, "coverage_complete"))
  plateau_stop <- isTRUE(attr(gamma_dict, "plateau_stop"))
  uncovered_targets <- attr(gamma_dict, "uncovered_targets")
  discovered_upper_gamma <- attr(gamma_dict, "discovered_upper_gamma")
  upper_cap_stop_reason <- attr(gamma_dict, "upper_cap_stop_reason")
  coarse_probe_count <- attr(gamma_dict, "coarse_probe_count")
  if (is.null(uncovered_targets)) {
    uncovered_targets <- setdiff(requested_cluster_range, as.integer(names(gamma_dict)))
  }

  batch_results <- clustering_main(
    igraph_obj = igraph_obj,
    cluster_range = requested_cluster_range,
    n_workers = n_workers,
    n_trials = n_trials,
    n_bootstrap = n_bootstrap,
    seed = seed,
    beta = beta,
    n_iterations = n_iterations,
    max_iterations = max_iterations,
    objective_function = objective_function,
    remove_threshold = remove_threshold,
    snn_graph = snn_graph,
    min_cluster_size = min_cluster_size,
    resolution_tolerance = resolution_tolerance,
    verbose = verbose,
    in_parallel_context = FALSE,
    runtime_context = runtime_context,
    precomputed_gamma_dict = gamma_dict,
    precomputed_resolution_search_diagnostics = resolution_search_diagnostics,
    precomputed_coverage_complete = coverage_complete
  )

  final_results <- finalize_cluster_range_results(
    target_results_dt = batch_results$target_results,
    requested_cluster_range = requested_cluster_range,
    searched_target_cluster_range = searched_target_cluster_range,
    search_coverage_complete = coverage_complete,
    gamma_dict = gamma_dict,
    resolution_search_diagnostics = resolution_search_diagnostics,
    plateau_stop = plateau_stop,
    search_uncovered_targets = uncovered_targets,
    discovered_upper_gamma = discovered_upper_gamma,
    upper_cap_stop_reason = upper_cap_stop_reason,
    coarse_probe_count = coarse_probe_count
  )

  if (!isTRUE(final_results$coverage_complete) || length(final_results$uncovered_targets) > 0L) {
    warning(
      paste(
        "scICE_clustering: requested final targets were not all recovered as returned final merged clusters;",
        "highest requested target =", max(requested_cluster_range),
        "| uncovered targets =",
        if (length(final_results$uncovered_targets) > 0L) paste(final_results$uncovered_targets, collapse = ",") else "none"
      ),
      call. = FALSE
    )
  }

  final_results
}

build_manual_resolution_results <- function(igraph_obj, resolution_values, n_workers,
                                            n_trials, n_bootstrap, seed, beta,
                                            n_iterations, objective_function,
                                            snn_graph = NULL, min_cluster_size = 1L,
                                            verbose = FALSE, runtime_context = NULL) {
  resolution_values <- as.numeric(resolution_values)
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  n_vertices <- igraph::vcount(igraph_obj)

  actual_workers <- if (.Platform$OS.type == "windows" && n_workers > 1L) 1L else as.integer(n_workers)
  active_resolution_workers <- max(1L, min(as.integer(length(resolution_values)), actual_workers))
  per_resolution_worker_budget <- max(
    1L,
    as.integer(floor(as.double(actual_workers) / as.double(active_resolution_workers)))
  )
  per_resolution_worker_budget <- cap_workers_by_memory(
    per_resolution_worker_budget,
    estimate_trial_matrix_bytes(n_vertices, n_trials, 1L),
    runtime_context
  )

  if (verbose) {
    scice_message("CLUSTERING_MAIN: Using manual resolution mode.")
    scice_message(
      paste(
        "CLUSTERING_MAIN: Resolution worker layout -",
        active_resolution_workers, "resolution workers x",
        per_resolution_worker_budget, "trial/bootstrap workers"
      )
    )
  }

  resolution_results <- cross_platform_mclapply(
    seq_along(resolution_values),
    function(resolution_idx) {
      resolution_value <- resolution_values[[resolution_idx]]
      worker_id <- paste("RESOLUTION", resolution_idx)
      if (verbose) {
        scice_message(
          paste(worker_id, ": Starting manual resolution evaluation for gamma =", signif(resolution_value, 6))
        )
      }

      evaluate_fixed_resolution(
        igraph_obj = igraph_obj,
        resolution = resolution_value,
        objective_function = objective_function,
        n_trials = n_trials,
        n_bootstrap = n_bootstrap,
        seed = seed,
        beta = beta,
        n_iterations = n_iterations,
        n_workers = per_resolution_worker_budget,
        snn_graph = snn_graph,
        min_cluster_size = min_cluster_size,
        verbose = verbose,
        worker_id = worker_id,
        in_parallel_context = TRUE,
        runtime_context = runtime_context
      )
    },
    mc.cores = active_resolution_workers,
    mc.preschedule = FALSE
  )

  cluster_numbers <- vapply(
    resolution_results,
    function(x) as.integer(x$best_labels_final_cluster_count),
    integer(1)
  )
  ic_scores <- vapply(
    resolution_results,
    function(x) as.numeric(x$ic_median),
    numeric(1)
  )

  pick_best_resolution <- function(indices) {
    finite_local <- which(is.finite(ic_scores[indices]))
    if (length(finite_local) == 0L) {
      return(indices[[1]])
    }
    indices[[finite_local[[which.min(ic_scores[indices][finite_local])]]]]
  }

  grouped_indices <- split(seq_along(cluster_numbers), cluster_numbers)
  selected_indices <- vapply(grouped_indices, pick_best_resolution, integer(1))
  selected_indices <- selected_indices[order(cluster_numbers[selected_indices], resolution_values[selected_indices])]
  selected_mask <- rep(FALSE, length(resolution_results))
  selected_mask[selected_indices] <- TRUE

  extract_final_cluster_median <- function(x) {
    if (!is.null(x$final_cluster_median) && length(x$final_cluster_median) == 1L) {
      return(as.numeric(x$final_cluster_median))
    }
    if (!is.null(x$best_labels_final_cluster_count) && length(x$best_labels_final_cluster_count) == 1L) {
      return(as.numeric(x$best_labels_final_cluster_count))
    }
    NA_real_
  }

  resolution_diagnostics <- data.frame(
    resolution = resolution_values,
    cluster_number = cluster_numbers,
    ic = ic_scores,
    effective_cluster_median = vapply(
      resolution_results,
      function(x) as.numeric(x$effective_cluster_median),
      numeric(1)
    ),
    raw_cluster_median = vapply(
      resolution_results,
      function(x) as.numeric(x$raw_cluster_median),
      numeric(1)
    ),
    final_cluster_median = vapply(
      resolution_results,
      extract_final_cluster_median,
      numeric(1)
    ),
    best_labels_raw_cluster_count = vapply(
      resolution_results,
      function(x) as.integer(x$best_labels_raw_cluster_count),
      integer(1)
    ),
    best_labels_final_cluster_count = vapply(
      resolution_results,
      function(x) as.integer(x$best_labels_final_cluster_count),
      integer(1)
    ),
    n_iter = vapply(
      resolution_results,
      function(x) as.integer(x$n_iterations),
      integer(1)
    ),
    selected = selected_mask,
    stringsAsFactors = FALSE
  )
  rownames(resolution_diagnostics) <- NULL

  selected_results <- resolution_results[selected_indices]

  list(
    gamma = vapply(selected_results, function(x) as.numeric(x$gamma), numeric(1)),
    labels = lapply(selected_results, function(x) x$labels),
    ic = vapply(selected_results, function(x) as.numeric(x$ic_median), numeric(1)),
    ic_vec = lapply(selected_results, function(x) x$ic_bootstrap),
    n_cluster = vapply(
      selected_results,
      function(x) as.integer(x$best_labels_final_cluster_count),
      integer(1)
    ),
    best_labels = lapply(selected_results, function(x) x$best_labels),
    effective_cluster_median = vapply(
      selected_results,
      function(x) as.numeric(x$effective_cluster_median),
      numeric(1)
    ),
    raw_cluster_median = vapply(
      selected_results,
      function(x) as.numeric(x$raw_cluster_median),
      numeric(1)
    ),
    final_cluster_median = vapply(
      selected_results,
      extract_final_cluster_median,
      numeric(1)
    ),
    admission_mode = vapply(selected_results, function(x) as.character(x$admission_mode), character(1)),
    best_labels_raw_cluster_count = vapply(
      selected_results,
      function(x) as.integer(x$best_labels_raw_cluster_count),
      integer(1)
    ),
    best_labels_final_cluster_count = vapply(
      selected_results,
      function(x) as.integer(x$best_labels_final_cluster_count),
      integer(1)
    ),
    n_iter = vapply(selected_results, function(x) as.integer(x$n_iterations), integer(1)),
    mei = lapply(selected_results, function(x) calculate_mei_from_array(x$labels)),
    k = vapply(selected_results, function(x) as.integer(x$k), integer(1)),
    source_target_cluster = rep(NA_integer_, length(selected_results)),
    excluded = rep(FALSE, length(selected_results)),
    exclusion_reason = rep("none", length(selected_results)),
    resolution_diagnostics = resolution_diagnostics
  )
}

#' Single-cell Inconsistency-based Clustering Evaluation (scICE)
#'
#' @description
#' scICE implements a systematic and efficient workflow to evaluate clustering consistency 
#' in single-cell RNA-seq data. It generates multiple cluster labels using the Leiden 
#' algorithm, calculates pairwise similarity with Element-Centric Similarity (ECS), 
#' and provides automated framework for finding consistent clusters.
#'
#' @details
#' Workflow overview:
#' \enumerate{
#'   \item extract a Seurat graph (for example \code{RNA_snn}) and convert to \pkg{igraph};
#'   \item in \code{cluster_range} mode, treat \code{cluster_range} as requested
#'   final merged cluster counts, run a shared global gamma coarse-to-refine
#'   sweep, and derive per-target gamma intervals from the observed
#'   effective/raw/final count curves without expanding to surrogate target
#'   cluster numbers;
#'   \item run repeated Leiden trials and bootstrap consistency scoring;
#'   \item compute IC/MEI metrics, merge undersized final clusters if requested,
#'   and key the main result object by the true final merged cluster counts.
#' }
#'
#' Practical guidance:
#' \itemize{
#'   \item use \code{objective_function = "CPM"} for most medium/large datasets;
#'   \item set \code{seed} for reproducible runs;
#'   \item set \code{remove_threshold = Inf} to keep all searched cluster numbers
#'   and skip inconsistency-based pre-filtering;
#'   \item use \code{plot_ic()} and \code{extract_consistent_clusters()} to inspect
#'   and report stable solutions.
#' }
#'
#' @param object A Seurat object containing single-cell data
#' @param graph_name Name of the graph to use for clustering. If NULL, will use the default SNN graph from the active assay (default: NULL)
#' @param cluster_range Vector of requested final merged cluster numbers to test
#'   (default: 1:20). In \code{cluster_range} mode, scICER keeps this set fixed
#'   and expands coverage only along the gamma axis via a shared global sweep;
#'   it no longer searches higher surrogate target cluster numbers.
#' @param n_workers Number of parallel workers to use (default: 10)
#' @param n_trials Number of clustering trials per resolution (default: 15)
#' @param n_bootstrap Number of bootstrap iterations (default: 100)
#' @param seed Random seed for reproducibility (default: NULL for random behavior)
#' @param beta Beta parameter for Leiden clustering (default: 0.1)
#' @param n_iterations Number of Leiden iterations (default: 10)
#' @param max_iterations Maximum iterations for optimization (default: 150)
#' @param ic_threshold IC threshold for consistent clustering (default: Inf)
#' @param objective_function Objective function for Leiden ("modularity" or "CPM", default: "CPM")
#' @param remove_threshold Threshold for removing inconsistent results (default: 1.15)
#' @param min_cluster_size Minimum number of cells required per cluster. Clusters
#'   smaller than this threshold are excluded from target-cluster counting during
#'   resolution search and optimization. IC/MEI are always computed from raw
#'   (unmerged) trial labels, while the final \code{best_labels} are merged so
#'   undersized clusters are not returned as separate output groups. Set to
#'   \code{1} to disable effective-cluster matching and final best-label merge.
#'   Default: \code{2}, matching Seurat singleton behavior.
#' @param resolution_tolerance Tolerance for resolution parameter search (default: 1e-8)
#' @param verbose Whether to print progress messages (default: TRUE)
#' @param resolution Optional manual Leiden resolution (gamma) value or vector.
#'   When provided, scICER skips \code{cluster_range}-based resolution search,
#'   evaluates the supplied gamma values directly, and keeps the lowest-IC
#'   solution for each final cluster number. If both \code{resolution} and
#'   \code{cluster_range} are provided, \code{cluster_range} is ignored with a
#'   user-facing message.
#'
#' @return A list containing:
#' \describe{
#'   \item{gamma}{Resolution parameters for each returned final cluster number}
#'   \item{labels}{Raw clustering arrays used for IC/MEI calculations (no small-cluster merge)}
#'   \item{ic}{Inconsistency scores for each returned final cluster number}
#'   \item{ic_vec}{Bootstrap IC distributions}
#'   \item{n_cluster}{Returned final merged cluster counts in the main result object}
#'   \item{best_labels}{Final best labels after applying the \code{min_cluster_size} merge rule}
#'   \item{effective_cluster_median}{Median effective cluster count at the selected gamma}
#'   \item{raw_cluster_median}{Median raw cluster count at the selected gamma}
#'   \item{admission_mode}{Phase-2 admission family used for the selected gamma}
#'   \item{best_labels_raw_cluster_count}{Raw cluster count of the selected best trial before final merge}
#'   \item{best_labels_final_cluster_count}{True final merged cluster count of the selected best trial}
#'   \item{source_target_cluster}{Requested target cluster number that produced each returned final result}
#'   \item{n_iter}{Number of iterations used}
#'   \item{mei}{Mutual Element-wise Information scores}
#'   \item{consistent_clusters}{Cluster numbers meeting consistency threshold}
#'   \item{min_cluster_size}{Minimum cluster size used for effective-count matching}
#'   \item{analysis_mode}{Whether the run used \code{"cluster_range"} or \code{"resolution"} mode}
#'   \item{resolution_input}{Manual gamma values requested by the user (if any)}
#'   \item{resolution_diagnostics}{Per-gamma evaluation summary for manual resolution mode}
#'   \item{requested_cluster_range}{User-requested final merged cluster numbers in \code{cluster_range} mode}
#'   \item{searched_target_cluster_range}{Target cluster numbers searched in \code{cluster_range} mode; after the shared gamma-sweep refactor this matches the requested final targets rather than an expanded surrogate target set}
#'   \item{coverage_complete}{Whether every requested final merged cluster target is present in the returned main result object}
#'   \item{search_coverage_complete}{Whether the shared gamma sweep resolved all requested targets tightly enough to hand off usable gamma intervals to optimization}
#'   \item{resolution_search_diagnostics}{Per-gamma shared sweep diagnostic table with probe stage, discovery-round labels, effective/raw/final cluster counts, high-gamma-degeneracy markers, coarse-sweep metadata, refinement-interval widths and allocations, refinement flags, and plateau-round markers}
#'   \item{discovered_upper_gamma}{Shared gamma sweep upper bound discovered before the first coarse sweep in \code{cluster_range} mode}
#'   \item{upper_cap_stop_reason}{Reason the shared gamma upper-cap discovery stopped; typically \code{target_covered}, \code{high_gamma_degenerate}, or \code{hard_cap}}
#'   \item{coarse_probe_count}{Number of probe gamma values used in the first shared coarse sweep}
#'   \item{target_diagnostics}{One-row-per-requested-target diagnostic table with requested/searched targets, gamma interval bounds, selected gamma, IC, exclusion status, effective/raw/final medians, and result-status labels such as \code{selected_main_result}, \code{deduplicated}, \code{resolution_search_failed}, or \code{optimization_admission_failed}}
#'   \item{plateau_stop}{Whether the shared gamma sweep stopped after consecutive plateau rounds without new coverage}
#'   \item{uncovered_targets}{Requested final cluster targets that are missing from the returned main result object}
#'   \item{search_uncovered_targets}{Requested final cluster targets that the shared gamma sweep did not resolve into optimization-ready gamma intervals}
#'   \item{best_cluster}{Cluster number with the lowest IC score among returned results}
#'   \item{best_resolution}{Gamma corresponding to \code{best_cluster}}
#' }
#'
#' @examples
#' \dontrun{
#' # Load Seurat object
#' data(pbmc_small)
#' 
#' # Run scICE analysis
#' scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)
#' 
#' # Plot IC results
#' plot_ic(scice_results)
#' 
#' # Extract consistent clustering labels
#' consistent_labels <- get_robust_labels(scice_results, threshold = 1.005)
#'
#' # Evaluate manually chosen gamma values
#' scice_manual <- scICE_clustering(
#'   pbmc_small,
#'   resolution = c(0.2, 0.4, 0.8)
#' )
#' }
#'
#' @export

scICE_clustering <- function(object, 
                            graph_name = NULL,
                            cluster_range = 1:20,
                            n_workers = 10,
                            n_trials = 15,
                            n_bootstrap = 100,
                            seed = NULL,
                            beta = 0.1,
                            n_iterations = 10,
                            max_iterations = 150,
                            ic_threshold = Inf,
                            objective_function = "CPM",
                            remove_threshold = 1.15,
                            min_cluster_size = 2,
                            resolution_tolerance = 1e-8,
                            verbose = TRUE,
                            resolution = NULL) {
  cluster_range_was_supplied <- !missing(cluster_range)
  resolution_mode <- !missing(resolution) && !is.null(resolution)
  resolution_input <- NULL
  if (resolution_mode) {
    if (!is.numeric(resolution) || length(resolution) == 0L) {
      stop("resolution must be a non-empty numeric value or vector.")
    }
    if (anyNA(resolution) || any(!is.finite(resolution))) {
      stop("resolution must contain only finite numeric values.")
    }
    resolution_input <- as.numeric(resolution)
  }
  
  # Clear clustering cache at start of new scICE run for optimal performance
  if (exists("clear_clustering_cache") && is.function(clear_clustering_cache)) {
    clear_clustering_cache()
  }
  
  # Validate inputs
  if (!inherits(object, "Seurat")) {
    stop("object must be a Seurat object")
  }
  
  # Get default graph name if not provided
  if (is.null(graph_name)) {
    default_assay <- DefaultAssay(object)
    graph_name <- paste0(default_assay, "_snn")
  }
  
  # Check if graph exists in object
  if (!(graph_name %in% names(object@graphs))) {
    stop(sprintf("Graph '%s' not found in Seurat object. Available graphs: %s", 
                graph_name, 
                paste(names(object@graphs), collapse = ", ")))
  }

  if (!is.numeric(n_workers) || length(n_workers) != 1 || is.na(n_workers)) {
    stop("n_workers must be a single numeric value.")
  }
  if (n_workers < 1) {
    stop("n_workers must be >= 1.")
  }

  if (!is.numeric(min_cluster_size) || length(min_cluster_size) != 1 || is.na(min_cluster_size)) {
    stop("min_cluster_size must be a single numeric value.")
  }
  if (abs(min_cluster_size - round(min_cluster_size)) > .Machine$double.eps^0.5) {
    stop("min_cluster_size must be an integer >= 1.")
  }
  min_cluster_size <- as.integer(round(min_cluster_size))
  if (min_cluster_size < 1L) {
    stop("min_cluster_size must be >= 1.")
  }

  if (!resolution_mode) {
    if (!is.numeric(cluster_range) || length(cluster_range) == 0L) {
      stop("cluster_range must be a non-empty numeric vector.")
    }
    if (anyNA(cluster_range) || any(!is.finite(cluster_range))) {
      stop("cluster_range must contain only finite numeric values.")
    }
    if (any(abs(cluster_range - round(cluster_range)) > .Machine$double.eps^0.5)) {
      stop("cluster_range must contain only integers >= 1.")
    }
    cluster_range <- sort(unique(as.integer(round(cluster_range))))
    if (any(cluster_range < 1L)) {
      stop("cluster_range must contain only integers >= 1.")
    }
  }

  if (resolution_mode && cluster_range_was_supplied) {
    emit_scice_notice(
      "scICE_clustering: `resolution` provided; ignoring `cluster_range` and evaluating manual gamma values directly.",
      verbose = verbose
    )
  }
  if (resolution_mode && anyDuplicated(resolution_input)) {
    duplicate_count <- sum(duplicated(resolution_input))
    resolution_input <- resolution_input[!duplicated(resolution_input)]
    emit_scice_notice(
      paste(
        "scICE_clustering: removed", duplicate_count,
        "duplicated manual resolution value(s) before evaluation."
      ),
      verbose = verbose
    )
  }

  requested_workers <- as.integer(n_workers)
  detected_cores <- parallel::detectCores()
  if (is.na(detected_cores) || detected_cores < 1) {
    detected_cores <- 1L
  }
  max_workers <- if (detected_cores > 1L) detected_cores - 1L else 1L
  effective_workers <- min(requested_workers, max_workers)
  if (.Platform$OS.type == "windows") {
    effective_workers <- 1L
  }
  n_workers <- max(1L, effective_workers)

  runtime_context <- create_runtime_context()
  on.exit(cleanup_runtime_spill(runtime_context), add = TRUE)
  
  if (verbose) {
    start_time <- Sys.time()
    scice_message(paste(rep("=", 80), collapse = ""))
    scice_message("Starting scICE clustering analysis...")
    scice_message(paste("Timestamp:", format(start_time, "%Y-%m-%d %H:%M:%S")))
    scice_message(paste("R Session Info: R", R.version.string))
    scice_message(paste("Platform:", R.version$platform))
    scice_message(paste("OS Type:", .Platform$OS.type))
    scice_message(paste("Process ID:", Sys.getpid()))
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("INPUT PARAMETERS:")
    scice_message(paste("  Using graph:", graph_name))
    if (resolution_mode) {
      scice_message(paste("  Analysis mode: manual resolution"))
      scice_message(paste("  Manual resolutions:", paste(signif(resolution_input, 6), collapse = ", ")))
      scice_message(paste("  Resolution count:", length(resolution_input)))
      if (cluster_range_was_supplied) {
        scice_message("  cluster_range: ignored because `resolution` was provided")
      }
    } else {
      scice_message(paste("  Analysis mode: cluster range search"))
      scice_message(paste("  Testing cluster range:", paste(cluster_range, collapse = ", ")))
      scice_message(paste("  Range: ", min(cluster_range), "-", max(cluster_range), " (", length(cluster_range), " values)"))
    }
    scice_message(paste("  Requested workers:", requested_workers))
    scice_message(paste("  Effective workers:", n_workers))
    scice_message(paste("  Internal memory budget (bytes):", format(runtime_context$memory_budget_bytes, scientific = FALSE)))
    scice_message(paste("  Number of trials per resolution:", n_trials))
    scice_message(paste("  Number of bootstrap iterations:", n_bootstrap))
    scice_message(paste("  Random seed:", if(is.null(seed)) "NULL (random)" else seed))
    scice_message(paste("  Beta parameter:", beta))
    scice_message(paste("  Leiden iterations:", n_iterations))
    scice_message(paste("  Maximum optimization iterations:", max_iterations))
    scice_message(paste("  IC threshold:", ic_threshold))
    scice_message(paste("  Objective function:", objective_function))
    if (resolution_mode) {
      scice_message("  Remove threshold: ignored in manual resolution mode")
    } else {
      scice_message(paste("  Remove threshold:", remove_threshold))
    }
    scice_message(paste("  Minimum cluster size:", min_cluster_size))
    if (min_cluster_size > 1L) {
      scice_message("  min_cluster_size semantics: counting uses effective clusters; final best_labels are merged")
    }
    if (resolution_mode) {
      scice_message("  Resolution tolerance: not used in manual resolution mode")
    } else {
      scice_message(paste("  Resolution tolerance:", resolution_tolerance))
    }
    scice_message(paste(rep("-", 80), collapse = ""))
  }
  
  # Parallel setup (mclapply backend in clustering_runtime.R)
  if (verbose) {
    scice_message("PARALLEL PROCESSING SETUP:")
    scice_message(paste("  Detected cores:", detected_cores))
    scice_message(paste("  Requested workers:", requested_workers))
    scice_message(paste("  Effective workers:", n_workers))
    if (n_workers == 1L) {
      scice_message("  Running in sequential mode (effective n_workers = 1)")
    } else {
      scice_message("  Using parallel::mclapply backend")
    }
  }
  
  # Get graph from Seurat object
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("GRAPH EXTRACTION:")
    scice_message(paste("  Accessing graph:", graph_name))
    scice_message(paste("  Available graphs in object:", paste(names(object@graphs), collapse = ", ")))
  }
  
  graph <- object@graphs[[graph_name]]
  
  if (verbose) {
    scice_message(paste("  Graph extraction successful"))
    scice_message(paste("  Graph class:", paste(class(graph), collapse = ", ")))
    scice_message(paste("  Graph dimensions:", nrow(graph), "x", ncol(graph)))
    scice_message(paste("  Graph storage type:", typeof(graph)))
    if (inherits(graph, "sparseMatrix")) {
      n_nonzero <- length(graph@x)
      scice_message(paste("  Non-zero entries:", n_nonzero))
      total_entries <- as.numeric(nrow(graph)) * as.numeric(ncol(graph))
      sparsity <- if (total_entries > 0) {
        (1 - n_nonzero / total_entries) * 100
      } else {
        NA_real_
      }
      scice_message(paste("  Sparsity:", round(sparsity, 2), "%"))
      if (n_nonzero > 0) {
        weights <- graph@x
        scice_message(paste("  Weight range: [", round(min(weights), 4), ", ", round(max(weights), 4), "]", sep = ""))
        scice_message(paste("  Mean weight:", round(mean(weights), 4)))
      }
    }
  }
  
  # Convert to igraph object
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("GRAPH CONVERSION:")
    conversion_start <- Sys.time()
    scice_message(paste("  Starting graph conversion at:", format(conversion_start, "%H:%M:%S")))
  }
  
  igraph_obj <- graph_to_igraph(graph, verbose = verbose)
  
  if (verbose) {
    conversion_end <- Sys.time()
    conversion_time <- as.numeric(difftime(conversion_end, conversion_start, units = "secs"))
    scice_message(paste("  Graph conversion completed in:", round(conversion_time, 3), "seconds"))
    scice_message(paste("  Converted graph vertices:", igraph::vcount(igraph_obj)))
    scice_message(paste("  Converted graph edges:", igraph::ecount(igraph_obj)))
    scice_message(paste("  igraph object class:", paste(class(igraph_obj), collapse = ", ")))
    if (igraph::ecount(igraph_obj) > 0) {
      scice_message(paste("  Graph is weighted:", igraph::is_weighted(igraph_obj)))
      if (igraph::is_weighted(igraph_obj)) {
        edge_weights <- igraph::E(igraph_obj)$weight
        scice_message(paste("  Edge weight range: [", round(min(edge_weights), 4), ", ", round(max(edge_weights), 4), "]", sep = ""))
      }
    }
  }
  
  # Perform clustering analysis
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("CLUSTERING ANALYSIS:")
    clustering_start <- Sys.time()
    scice_message(paste("  Starting clustering analysis at:", format(clustering_start, "%H:%M:%S")))
    scice_message(paste("  Thread context: Main thread (PID:", Sys.getpid(), ")"))
    if (n_workers > 1) {
      scice_message(paste("  Parallel workers will be spawned for sub-tasks"))
    }
  }

  if (resolution_mode) {
    if (verbose) {
      scice_message("  Manual resolution mode selected - skipping cluster_range search and evaluating supplied gamma values directly.")
    }
    results <- build_manual_resolution_results(
      igraph_obj = igraph_obj,
      resolution_values = resolution_input,
      n_workers = n_workers,
      n_trials = n_trials,
      n_bootstrap = n_bootstrap,
      seed = seed,
      beta = beta,
      n_iterations = n_iterations,
      objective_function = objective_function,
      snn_graph = graph,
      min_cluster_size = min_cluster_size,
      verbose = verbose,
      runtime_context = runtime_context
    )
  } else {
    results <- run_cluster_range_mode(
      igraph_obj = igraph_obj,
      requested_cluster_range = cluster_range,
      n_workers = n_workers,
      n_trials = n_trials,
      n_bootstrap = n_bootstrap,
      seed = seed,
      beta = beta,
      n_iterations = n_iterations,
      max_iterations = max_iterations,
      objective_function = objective_function,
      remove_threshold = remove_threshold,
      snn_graph = graph,
      min_cluster_size = min_cluster_size,
      resolution_tolerance = resolution_tolerance,
      verbose = verbose,
      runtime_context = runtime_context
    )
  }
  
  if (verbose) {
    clustering_end <- Sys.time()
    clustering_time <- as.numeric(difftime(clustering_end, clustering_start, units = "secs"))
    scice_message(paste("  Clustering analysis completed in:", round(clustering_time, 3), "seconds"))
    scice_message(paste("  Results structure:"))
    scice_message(paste("    - gamma length:", length(results$gamma)))
    scice_message(paste("    - labels length:", length(results$labels)))
    scice_message(paste("    - ic length:", length(results$ic)))
    scice_message(paste("    - n_cluster length:", length(results$n_cluster)))
    scice_message(paste("    - effective_cluster_median length:", length(results$effective_cluster_median)))
    scice_message(paste("    - raw_cluster_median length:", length(results$raw_cluster_median)))
    scice_message(paste("    - final_cluster_median length:", length(results$final_cluster_median)))
    scice_message(paste("    - admission_mode length:", length(results$admission_mode)))
    scice_message(paste("    - best_labels_raw_cluster_count length:", length(results$best_labels_raw_cluster_count)))
    scice_message(paste("    - best_labels_final_cluster_count length:", length(results$best_labels_final_cluster_count)))
    scice_message(paste("    - source_target_cluster length:", length(results$source_target_cluster)))
    if (!is.null(results$target_diagnostics)) {
      scice_message(paste("    - target diagnostics rows:", nrow(results$target_diagnostics)))
    }
    if (!is.null(results$resolution_search_diagnostics)) {
      scice_message(paste("    - resolution search diagnostics rows:", nrow(results$resolution_search_diagnostics)))
    }
    if (!is.null(results$search_coverage_complete) && !is.na(results$search_coverage_complete)) {
      scice_message(paste("    - search coverage complete:", results$search_coverage_complete))
    }
    if (!is.null(results$coverage_complete) && !is.na(results$coverage_complete)) {
      scice_message(paste("    - coverage complete:", results$coverage_complete))
    }
  }
  
  # Add cell names to results
  results$cell_names <- Cells(object)
  
  results$analysis_mode <- if (resolution_mode) "resolution" else "cluster_range"
  results$resolution_input <- if (resolution_mode) resolution_input else NULL
  results$resolution_diagnostics <- if (resolution_mode) results$resolution_diagnostics else NULL
  results$resolution_search_diagnostics <- if (resolution_mode) NULL else results$resolution_search_diagnostics
  results$requested_cluster_range <- if (resolution_mode) NULL else results$requested_cluster_range
  results$searched_target_cluster_range <- if (resolution_mode) NULL else results$searched_target_cluster_range
  results$search_coverage_complete <- if (resolution_mode) TRUE else isTRUE(results$search_coverage_complete)
  results$coverage_complete <- if (resolution_mode) TRUE else isTRUE(results$coverage_complete)
  results$target_diagnostics <- if (resolution_mode) NULL else results$target_diagnostics
  results$cluster_range_tested <- results$n_cluster
  
  # Determine consistent clusters using actual cluster numbers, not indices
  if (verbose) {
    scice_message(paste(rep("-", 80), collapse = ""))
    scice_message("RESULTS ANALYSIS:")
    scice_message(paste("  IC threshold for consistency:", ic_threshold))
    scice_message(paste("  Results object is null:", is.null(results)))
    scice_message(paste("  IC vector is null:", is.null(results$ic)))
    scice_message(paste("  IC vector length:", if(!is.null(results$ic)) length(results$ic) else 0))
  }
  
  if (!is.null(results) && !is.null(results$ic) && length(results$ic) > 0) {
    valid_ic <- results$ic[!is.na(results$ic)]
    if (verbose) {
      scice_message(paste("  Valid IC scores:", length(valid_ic)))
      if (length(valid_ic) > 0) {
        scice_message(paste("  IC score range: [", round(min(valid_ic), 4), ", ", round(max(valid_ic), 4), "]", sep = ""))
        scice_message(paste("  IC scores below threshold:", sum(valid_ic < ic_threshold)))
      }
    }
    
    consistent_indices <- which(results$ic < ic_threshold)
    results$consistent_clusters <- results$n_cluster[consistent_indices]
    
    if (verbose) {
      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      scice_message(paste(rep("=", 80), collapse = ""))
      scice_message("ANALYSIS COMPLETE")
      scice_message(paste("  Total execution time:", round(total_time, 3), "seconds"))
      scice_message(paste("  Found", length(results$consistent_clusters), 
                    "consistent cluster numbers:", paste(results$consistent_clusters, collapse = ", ")))
      
      target_diagnostics <- if (!resolution_mode) results$target_diagnostics else NULL
      excluded_targets <- if (!is.null(target_diagnostics)) {
        target_diagnostics$searched_target_cluster[target_diagnostics$excluded]
      } else {
        integer(0)
      }
      if (!resolution_mode && length(excluded_targets) > 0) {
        scice_message(paste(
          "  Excluded", length(excluded_targets),
          "searched target cluster numbers:",
          paste(excluded_targets, collapse = ", ")
        ))
      }
      if (!resolution_mode && !isTRUE(results$coverage_complete)) {
        requested_max <- max(results$requested_cluster_range)
        returned_max <- if (length(results$n_cluster) > 0L) max(results$n_cluster) else NA_integer_
        scice_message(paste(
          "  Coverage incomplete:",
          "requested final max =", requested_max,
          "| searched target max =", max(results$searched_target_cluster_range),
          "| returned final max =", ifelse(is.na(returned_max), "NA", as.character(returned_max))
        ))
        if (!is.null(results$search_coverage_complete)) {
          scice_message(paste(
            "  Shared gamma sweep coverage complete:",
            isTRUE(results$search_coverage_complete)
          ))
        }
        if (!is.null(results$uncovered_targets) && length(results$uncovered_targets) > 0L) {
          scice_message(paste(
            "  Uncovered requested final targets:",
            paste(results$uncovered_targets, collapse = ", ")
          ))
        }
      }
      if (resolution_mode && !is.null(results$resolution_diagnostics)) {
        deduplicated_count <- sum(!results$resolution_diagnostics$selected)
        scice_message(paste("  Manual resolutions evaluated:", nrow(results$resolution_diagnostics)))
        scice_message(paste("  Manual resolutions retained after per-cluster IC selection:", length(results$n_cluster)))
        if (deduplicated_count > 0) {
          scice_message(paste("  Manual resolutions superseded by lower-IC matches:", deduplicated_count))
        }
      }
      
      # Report tested but inconsistent clusters
      inconsistent_clusters <- setdiff(results$n_cluster, results$consistent_clusters)
      if (length(inconsistent_clusters) > 0) {
        scice_message(paste("  Found", length(inconsistent_clusters), "inconsistent cluster numbers:", 
                      paste(inconsistent_clusters, collapse = ", ")))
      }
      
      # Detailed breakdown
      scice_message("  Detailed breakdown:")
      if (resolution_mode) {
        scice_message(paste("    - Requested clusters:", length(resolution_input)))
      } else {
        scice_message(paste("    - Requested final cluster targets:", length(results$requested_cluster_range)))
        scice_message(paste("    - Searched target clusters:", length(results$searched_target_cluster_range)))
        if (!is.null(results$resolution_search_diagnostics)) {
          scice_message(paste("    - Shared gamma probes:", nrow(results$resolution_search_diagnostics)))
        }
        if (!is.null(results$search_coverage_complete)) {
          scice_message(paste("    - Shared gamma sweep coverage complete:", isTRUE(results$search_coverage_complete)))
        }
      }
      scice_message(paste("    - Excluded searched targets:", if (resolution_mode) 0 else length(excluded_targets)))
      scice_message(paste("    - Returned final clusters:", length(results$n_cluster)))
      scice_message(paste("    - Consistent clusters:", length(results$consistent_clusters)))
      scice_message(paste("    - Inconsistent clusters:", length(inconsistent_clusters)))
    }
  } else {
    results$consistent_clusters <- integer(0)
    if (verbose) {
      total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      scice_message(paste(rep("=", 80), collapse = ""))
      scice_message("ANALYSIS COMPLETE")
      scice_message(paste("  Total execution time:", round(total_time, 3), "seconds"))
      scice_message("  WARNING: No consistent cluster numbers found.")
      scice_message("  Possible reasons:")
      scice_message("    - All clusters were excluded during filtering")
      scice_message("    - IC threshold too strict")
      scice_message("    - Data quality issues")
      scice_message("    - Graph connectivity problems")
      scice_message("  Debugging suggestions:")
      scice_message("    - Try increasing ic_threshold (e.g., to Inf to see all results)")
      scice_message("    - Try increasing remove_threshold")
      scice_message("    - Check graph construction parameters")
      scice_message("    - Verify input data quality")
    }
  }
  
  # Keep lightweight result object (do not store full Seurat object)
  results$graph_name <- graph_name
  results$min_cluster_size <- min_cluster_size
  valid_best_indices <- which(!is.na(results$ic))
  if (length(valid_best_indices) > 0) {
    best_index <- valid_best_indices[[which.min(results$ic[valid_best_indices])]]
    results$best_cluster <- results$n_cluster[[best_index]]
    results$best_resolution <- results$gamma[[best_index]]
  } else {
    results$best_cluster <- NA_integer_
    results$best_resolution <- NA_real_
  }
  
  class(results) <- "scICE"
  return(results)
} 
