#' @importFrom stats median setNames
#' @importFrom parallel detectCores mclapply
#' @importFrom data.table data.table rbindlist
NULL

empty_cluster_results_dt <- function() {
  data.table::data.table(
    cluster_number = integer(),
    gamma = numeric(),
    labels = list(),
    ic = numeric(),
    ic_vec = list(),
    best_labels = list(),
    effective_cluster_median = numeric(),
    raw_cluster_median = numeric(),
    final_cluster_median = numeric(),
    admission_mode = character(),
    best_labels_raw_cluster_count = integer(),
    best_labels_final_cluster_count = integer(),
    n_iter = integer(),
    mei = list(),
    k = integer(),
    source_target_cluster = integer(),
    excluded = logical(),
    exclusion_reason = character(),
    result_status = character(),
    phase1_primary_gamma_count = integer(),
    phase1_secondary_gamma_count = integer(),
    phase1_total_gamma_count = integer(),
    phase1_elapsed_sec = numeric(),
    phase1_leiden_runs = integer(),
    secondary_phase1_used = logical(),
    exact_hit_gamma_count = integer(),
    phase4_iterations = integer(),
    phase4_elapsed_sec = numeric(),
    phase5_elapsed_sec = numeric(),
    optimization_elapsed_sec = numeric()
  )
}

empty_target_results_dt <- function() {
  dt <- empty_cluster_results_dt()
  dt$selected_main_result <- logical()
  dt
}

select_lowest_ic_index <- function(ic_values) {
  ic_values <- as.numeric(ic_values)
  finite_indices <- which(is.finite(ic_values))
  if (length(finite_indices) == 0L) {
    return(1L)
  }
  finite_indices[[which.min(ic_values[finite_indices])]]
}

rekey_target_results_by_final_cluster <- function(target_results_dt) {
  if (is.null(target_results_dt) || nrow(target_results_dt) == 0L) {
    return(list(
      main_results_dt = empty_cluster_results_dt(),
      target_results_dt = empty_target_results_dt()
    ))
  }

  target_results_dt <- data.table::copy(target_results_dt)
  if (!"source_target_cluster" %in% colnames(target_results_dt)) {
    target_results_dt[, source_target_cluster := as.integer(cluster_number)]
  }
  target_results_dt[, selected_main_result := FALSE]
  target_results_dt[, target_row_id := seq_len(.N)]

  candidate_rows <- target_results_dt[
    !excluded & !is.na(best_labels_final_cluster_count)
  ]

  if (nrow(candidate_rows) == 0L) {
    target_results_dt[, result_status := as.character(exclusion_reason)]
    target_results_dt[, target_row_id := NULL]
    return(list(
      main_results_dt = empty_cluster_results_dt(),
      target_results_dt = target_results_dt
    ))
  }

  data.table::setorder(
    candidate_rows,
    best_labels_final_cluster_count,
    source_target_cluster,
    gamma
  )

  selected_row_ids <- candidate_rows[, {
    chosen_idx <- select_lowest_ic_index(ic)
    .(target_row_id = target_row_id[[chosen_idx]])
  }, by = best_labels_final_cluster_count]$target_row_id

  target_results_dt[target_row_id %in% selected_row_ids, selected_main_result := TRUE]
  target_results_dt[
    ,
    result_status := data.table::fifelse(
      selected_main_result,
      "selected_main_result",
      data.table::fifelse(excluded, as.character(exclusion_reason), "deduplicated")
    )
  ]

  main_results_dt <- data.table::copy(
    target_results_dt[selected_main_result == TRUE]
  )
  main_results_dt[, cluster_number := as.integer(best_labels_final_cluster_count)]
  main_results_dt[, excluded := FALSE]
  main_results_dt[, exclusion_reason := "none"]
  main_results_dt[, result_status := "selected_main_result"]
  main_results_dt[, selected_main_result := NULL]
  main_results_dt[, target_row_id := NULL]
  data.table::setorder(main_results_dt, cluster_number, source_target_cluster, gamma)

  target_results_dt[, target_row_id := NULL]
  data.table::setorder(target_results_dt, source_target_cluster)

  list(
    main_results_dt = main_results_dt,
    target_results_dt = target_results_dt
  )
}

cluster_results_dt_to_list <- function(results_dt, target_results_dt = NULL) {
  if (is.null(results_dt) || nrow(results_dt) == 0L) {
    results_dt <- empty_cluster_results_dt()
  }

  result <- list(
    gamma = if ("gamma" %in% colnames(results_dt)) results_dt$gamma else numeric(0),
    labels = if ("labels" %in% colnames(results_dt)) results_dt$labels else list(),
    ic = if ("ic" %in% colnames(results_dt)) results_dt$ic else numeric(0),
    ic_vec = if ("ic_vec" %in% colnames(results_dt)) results_dt$ic_vec else list(),
    n_cluster = if ("cluster_number" %in% colnames(results_dt)) results_dt$cluster_number else integer(0),
    best_labels = if ("best_labels" %in% colnames(results_dt)) results_dt$best_labels else list(),
    effective_cluster_median = if ("effective_cluster_median" %in% colnames(results_dt)) results_dt$effective_cluster_median else numeric(0),
    raw_cluster_median = if ("raw_cluster_median" %in% colnames(results_dt)) results_dt$raw_cluster_median else numeric(0),
    final_cluster_median = if ("final_cluster_median" %in% colnames(results_dt)) results_dt$final_cluster_median else numeric(0),
    admission_mode = if ("admission_mode" %in% colnames(results_dt)) results_dt$admission_mode else character(0),
    best_labels_raw_cluster_count = if ("best_labels_raw_cluster_count" %in% colnames(results_dt)) results_dt$best_labels_raw_cluster_count else integer(0),
    best_labels_final_cluster_count = if ("best_labels_final_cluster_count" %in% colnames(results_dt)) results_dt$best_labels_final_cluster_count else integer(0),
    n_iter = if ("n_iter" %in% colnames(results_dt)) results_dt$n_iter else integer(0),
    mei = if ("mei" %in% colnames(results_dt)) results_dt$mei else list(),
    k = if ("k" %in% colnames(results_dt)) results_dt$k else integer(0),
    source_target_cluster = if ("source_target_cluster" %in% colnames(results_dt)) results_dt$source_target_cluster else integer(0),
    excluded = if ("excluded" %in% colnames(results_dt)) results_dt$excluded else logical(0),
    exclusion_reason = if ("exclusion_reason" %in% colnames(results_dt)) results_dt$exclusion_reason else character(0),
    result_status = if ("result_status" %in% colnames(results_dt)) results_dt$result_status else character(0),
    phase1_primary_gamma_count = if ("phase1_primary_gamma_count" %in% colnames(results_dt)) results_dt$phase1_primary_gamma_count else integer(0),
    phase1_secondary_gamma_count = if ("phase1_secondary_gamma_count" %in% colnames(results_dt)) results_dt$phase1_secondary_gamma_count else integer(0),
    phase1_total_gamma_count = if ("phase1_total_gamma_count" %in% colnames(results_dt)) results_dt$phase1_total_gamma_count else integer(0),
    phase1_elapsed_sec = if ("phase1_elapsed_sec" %in% colnames(results_dt)) results_dt$phase1_elapsed_sec else numeric(0),
    phase1_leiden_runs = if ("phase1_leiden_runs" %in% colnames(results_dt)) results_dt$phase1_leiden_runs else integer(0),
    secondary_phase1_used = if ("secondary_phase1_used" %in% colnames(results_dt)) results_dt$secondary_phase1_used else logical(0),
    exact_hit_gamma_count = if ("exact_hit_gamma_count" %in% colnames(results_dt)) results_dt$exact_hit_gamma_count else integer(0),
    phase4_iterations = if ("phase4_iterations" %in% colnames(results_dt)) results_dt$phase4_iterations else integer(0),
    phase4_elapsed_sec = if ("phase4_elapsed_sec" %in% colnames(results_dt)) results_dt$phase4_elapsed_sec else numeric(0),
    phase5_elapsed_sec = if ("phase5_elapsed_sec" %in% colnames(results_dt)) results_dt$phase5_elapsed_sec else numeric(0),
    optimization_elapsed_sec = if ("optimization_elapsed_sec" %in% colnames(results_dt)) results_dt$optimization_elapsed_sec else numeric(0)
  )

  if (!is.null(target_results_dt)) {
    result$target_results <- target_results_dt
  }

  result
}

build_target_gamma_seed_table <- function(target_cluster, gamma_dict,
                                          target_gamma_seeds = NULL,
                                          target_interval_details = NULL,
                                          resolution_search_diagnostics = NULL) {
  target_key <- as.character(target_cluster)
  interval_detail <- NULL
  if (!is.null(target_interval_details) && target_key %in% names(target_interval_details)) {
    interval_detail <- target_interval_details[[target_key]]
  }

  gamma_left <- if (!is.null(interval_detail) && isTRUE(is.finite(interval_detail$gamma_left))) {
    as.numeric(interval_detail$gamma_left)
  } else if (!is.null(gamma_dict) && target_key %in% names(gamma_dict)) {
    as.numeric(gamma_dict[[target_key]][1])
  } else {
    NA_real_
  }
  gamma_right <- if (!is.null(interval_detail) && isTRUE(is.finite(interval_detail$gamma_right))) {
    as.numeric(interval_detail$gamma_right)
  } else if (!is.null(gamma_dict) && target_key %in% names(gamma_dict)) {
    as.numeric(gamma_dict[[target_key]][2])
  } else {
    NA_real_
  }

  seed_values <- if (!is.null(target_gamma_seeds) && target_key %in% names(target_gamma_seeds)) {
    as.numeric(target_gamma_seeds[[target_key]])
  } else {
    numeric(0)
  }
  seed_values <- sort(unique(seed_values[is.finite(seed_values)]))

  exact_probe_values <- if (!is.null(interval_detail) && !is.null(interval_detail$exact_probe_values)) {
    sort(unique(as.numeric(interval_detail$exact_probe_values)))
  } else {
    numeric(0)
  }
  near_probe_values <- if (!is.null(interval_detail) && !is.null(interval_detail$near_probe_values)) {
    sort(unique(as.numeric(interval_detail$near_probe_values)))
  } else {
    numeric(0)
  }
  near_probe_values <- setdiff(near_probe_values, exact_probe_values)

  selected_gamma <- NA_real_
  if (length(seed_values) > 0L && is.finite(gamma_left) && is.finite(gamma_right)) {
    interval_midpoint <- mean(c(gamma_left, gamma_right))
    selected_gamma <- seed_values[[which.min(abs(seed_values - interval_midpoint))]]
  } else if (length(seed_values) > 0L) {
    selected_gamma <- seed_values[[1]]
  }

  role_rows <- list(
    data.frame(gamma = gamma_left, seed_role = "left", stringsAsFactors = FALSE),
    data.frame(gamma = gamma_right, seed_role = "right", stringsAsFactors = FALSE),
    data.frame(gamma = selected_gamma, seed_role = "selected", stringsAsFactors = FALSE),
    data.frame(gamma = exact_probe_values, seed_role = "exact", stringsAsFactors = FALSE),
    data.frame(gamma = near_probe_values, seed_role = "near", stringsAsFactors = FALSE)
  )
  covered_values <- unique(c(gamma_left, gamma_right, selected_gamma, exact_probe_values, near_probe_values))
  generic_seed_values <- setdiff(seed_values, covered_values)
  role_rows[[length(role_rows) + 1L]] <- data.frame(
    gamma = generic_seed_values,
    seed_role = rep("seed", length(generic_seed_values)),
    stringsAsFactors = FALSE
  )

  seed_table <- data.table::rbindlist(role_rows, fill = TRUE)
  if (nrow(seed_table) == 0L) {
    return(data.frame(
      gamma = numeric(0),
      seed_role = character(0),
      final_cluster_count = numeric(0),
      raw_cluster_count = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  seed_table <- seed_table[is.finite(gamma)]
  seed_table[, final_cluster_count := NA_real_]
  seed_table[, raw_cluster_count := NA_real_]

  if (!is.null(resolution_search_diagnostics) &&
      is.data.frame(resolution_search_diagnostics) &&
      nrow(resolution_search_diagnostics) > 0L &&
      all(c("gamma", "final_cluster_count", "raw_cluster_count") %in% colnames(resolution_search_diagnostics))) {
    diagnostics_dt <- data.table::as.data.table(resolution_search_diagnostics)[
      ,
      .(gamma = as.numeric(gamma),
        final_cluster_count = as.numeric(final_cluster_count),
        raw_cluster_count = as.numeric(raw_cluster_count))
    ]
    tol <- max(.Machine$double.eps^0.5, abs(gamma_right - gamma_left) * 1e-8, 1e-12)
    matched_counts <- lapply(seed_table$gamma, function(gamma_value) {
      delta <- abs(diagnostics_dt$gamma - gamma_value)
      if (length(delta) == 0L || all(!is.finite(delta))) {
        return(c(NA_real_, NA_real_))
      }
      idx <- which.min(delta)
      if (!is.finite(delta[[idx]]) || delta[[idx]] > tol) {
        return(c(NA_real_, NA_real_))
      }
      c(
        diagnostics_dt$final_cluster_count[[idx]],
        diagnostics_dt$raw_cluster_count[[idx]]
      )
    })
    matched_counts <- do.call(rbind, matched_counts)
    if (!is.null(matched_counts)) {
      seed_table[, final_cluster_count := as.numeric(matched_counts[, 1])]
      seed_table[, raw_cluster_count := as.numeric(matched_counts[, 2])]
    }
  }

  seed_table <- unique(seed_table, by = c("gamma", "seed_role"))
  seed_table <- seed_table[order(gamma, seed_role)]
  as.data.frame(seed_table)
}

#' Core clustering algorithm implementing binary search and optimization
#'
#' @description
#' Internal engine that executes the end-to-end scICER clustering workflow.
#'
#' @details
#' Major stages:
#' \enumerate{
#'   \item clear cache and initialize run state;
#'   \item search feasible resolution ranges per target cluster number;
#'   \item optionally filter unstable targets (skipped when threshold is \code{Inf});
#'   \item optimize each valid target via repeated Leiden trials and bootstrap IC;
#'   \item merge successful/excluded entries into a compatibility-preserving result list.
#' }
#'
#' The function includes extensive defensive handling for empty/failed branches,
#' making it safe for large parameter sweeps where some targets may fail to
#' converge.
#'
#' @param igraph_obj igraph object to cluster
#' @param cluster_range Vector of cluster numbers to test
#' @param n_workers Number of parallel workers (default: max(1, parallel::detectCores() - 1))
#' @param n_trials Number of clustering trials per resolution
#' @param n_bootstrap Number of bootstrap iterations
#' @param seed Optional deterministic seed
#' @param beta Beta parameter for Leiden clustering
#' @param n_iterations Number of Leiden iterations
#' @param max_iterations Maximum iterations for optimization
#' @param objective_function Objective function ("modularity" or "CPM")
#' @param remove_threshold Threshold for removing inconsistent results
#' @param snn_graph Seurat SNN graph matrix for final small-cluster merge
#' @param min_cluster_size Minimum required cells per effective cluster count
#' @param resolution_tolerance Tolerance for resolution parameter search
#' @param verbose Whether to print progress messages
#' @param in_parallel_context Whether this function is called from within a parallel context (default: FALSE)
#' @param runtime_context Internal mutable context for spill/memory controls
#' @return List with clustering results
#' @keywords internal
clustering_main <- function(igraph_obj, cluster_range, n_workers = max(1, parallel::detectCores() - 1), 
                          n_trials, n_bootstrap, seed = NULL, beta, n_iterations, max_iterations, 
                          objective_function, remove_threshold, snn_graph = NULL,
                          min_cluster_size = 1L, resolution_tolerance, verbose, 
                          in_parallel_context = FALSE, runtime_context = NULL,
                          precomputed_gamma_dict = NULL,
                          precomputed_resolution_search_diagnostics = NULL,
                          precomputed_coverage_complete = NULL) {
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  if (min_cluster_size > 1L && is.null(snn_graph)) {
    stop("snn_graph must be provided when min_cluster_size > 1.")
  }
  if (!is.null(snn_graph) && igraph::vcount(igraph_obj) != nrow(snn_graph)) {
    stop("snn_graph row count must match the number of graph vertices.")
  }
  if (!is.null(snn_graph) && nrow(snn_graph) != ncol(snn_graph)) {
    stop("snn_graph must be a square matrix.")
  }
  
  # Clear clustering cache at start for fresh run
  if (verbose) {
    cache_stats_before <- get_cache_stats()
    scice_message(paste("CLUSTERING_MAIN: Cache entries before clearing:", cache_stats_before$cache_entries))
  }
  clear_clustering_cache()
  
  # Set random seed if provided for reproducibility
  if (!is.null(seed)) {
    set.seed(seed)
    if (verbose) {
      scice_message(paste("CLUSTERING_MAIN: Setting random seed to:", seed))
      scice_message(paste("CLUSTERING_MAIN: Thread context - PID:", Sys.getpid()))
    }
  }
  
  # Initialize results storage using data.table for better performance
  results_dt <- empty_target_results_dt()
  
  # Determine resolution search bounds
  if (objective_function == "modularity") {
    start_g <- -13
    end_g <- 20  # Increased for higher cluster numbers
  } else { # CPM
    start_g <- log(resolution_tolerance)
    if (start_g < -20) start_g <- -20
    end_g <- 20  # Increased for higher cluster numbers
  }
  
  if (verbose) {
    scice_message("CLUSTERING_MAIN: Starting resolution search...")
    scice_message(paste("CLUSTERING_MAIN: Objective function:", objective_function))
    scice_message(paste("CLUSTERING_MAIN: Search bounds: [", start_g, ", ", end_g, "]", sep = ""))
    scice_message(paste("CLUSTERING_MAIN: Target cluster range:", paste(cluster_range, collapse = ", ")))
    resolution_search_start <- Sys.time()
  }

  if (is.null(precomputed_gamma_dict)) {
    gamma_dict <- find_resolution_ranges(
      igraph_obj, cluster_range, start_g, end_g, objective_function,
      resolution_tolerance, n_workers, verbose, seed, snn_graph, min_cluster_size,
      in_parallel_context, runtime_context
    )
    resolution_search_diagnostics <- attr(gamma_dict, "resolution_search_diagnostics")
    search_coverage_complete <- attr(gamma_dict, "coverage_complete")
  } else {
    gamma_dict <- precomputed_gamma_dict
    resolution_search_diagnostics <- precomputed_resolution_search_diagnostics
    search_coverage_complete <- precomputed_coverage_complete
  }
  target_gamma_seeds <- attr(gamma_dict, "target_gamma_seeds")
  if (is.null(target_gamma_seeds)) {
    target_gamma_seeds <- setNames(vector("list", length(cluster_range)), as.character(cluster_range))
  }
  target_interval_details <- attr(gamma_dict, "target_interval_details")
  if (is.null(target_interval_details)) {
    target_interval_details <- setNames(vector("list", length(cluster_range)), as.character(cluster_range))
  }

  if (verbose) {
    resolution_search_time <- as.numeric(difftime(Sys.time(), resolution_search_start, units = "secs"))
    scice_message(paste("CLUSTERING_MAIN: Resolution search completed in", round(resolution_search_time, 3), "seconds"))
    scice_message(paste("CLUSTERING_MAIN: Found resolution ranges for", length(gamma_dict), "cluster numbers"))
    if (!is.null(resolution_search_diagnostics)) {
      scice_message(paste("CLUSTERING_MAIN: Shared gamma probes evaluated:", nrow(resolution_search_diagnostics)))
      if (!is.null(search_coverage_complete)) {
        scice_message(paste("CLUSTERING_MAIN: Search coverage complete:", isTRUE(search_coverage_complete)))
      }
    }
    if (length(gamma_dict) > 0) {
      for (i in 1:min(5, length(gamma_dict))) {
        cluster_num <- names(gamma_dict)[i]
        range_vals <- gamma_dict[[cluster_num]]
        scice_message(paste("CLUSTERING_MAIN:   k=", cluster_num, ": gamma in [",
                     signif(range_vals[1], 6), ", ", signif(range_vals[2], 6), "]", sep = ""))
      }
      if (length(gamma_dict) > 5) {
        scice_message(paste("CLUSTERING_MAIN:   ... and", length(gamma_dict) - 5, "more"))
      }
    }
  }
  
  # Filter out problematic cluster numbers in parallel
  if (verbose) {
    scice_message("CLUSTERING_MAIN: Starting problematic cluster filtering...")
    scice_message(paste("CLUSTERING_MAIN: Remove threshold:", remove_threshold))
    scice_message(paste("CLUSTERING_MAIN: Platform:", .Platform$OS.type))
    filtering_start <- Sys.time()
  }
  
  if (is.infinite(remove_threshold)) {
    cluster_filter_results <- lapply(cluster_range, function(cluster_num) {
      if (!(as.character(cluster_num) %in% names(gamma_dict))) {
        list(excluded = TRUE, cluster_num = cluster_num, reason = "resolution_search_failed")
      } else {
        list(excluded = FALSE, cluster_num = cluster_num, reason = "filtering_skipped_inf_threshold")
      }
    })
    excluded_numbers <- sapply(cluster_filter_results, function(x) if (x$excluded) x$cluster_num else NULL)
    excluded_numbers <- unlist(excluded_numbers)
    if (verbose) {
      scice_message("CLUSTERING_MAIN: remove_threshold is Inf - skipping filtering step.")
      filtering_time <- as.numeric(difftime(Sys.time(), filtering_start, units = "secs"))
      scice_message(paste("CLUSTERING_MAIN: Filtering skipped in", round(filtering_time, 3), "seconds"))
    }
  } else {
    # Windows compatibility for parallel processing
    actual_workers <- if (.Platform$OS.type == "windows" && n_workers > 1) 1 else n_workers
    active_filter_workers <- max(1L, min(as.integer(length(cluster_range)), as.integer(actual_workers)))
    filter_gamma_workers <- max(1L, as.integer(floor(as.double(actual_workers) / active_filter_workers)))
    filter_gamma_workers <- cap_workers_by_memory(
      filter_gamma_workers,
      estimate_trial_matrix_bytes(igraph::vcount(igraph_obj), 10L, 1L),
      runtime_context
    )
    
    if (verbose) {
      scice_message(paste("CLUSTERING_MAIN: Using", actual_workers, "workers for filtering (requested:", n_workers, ")"))
      scice_message(
        paste(
          "CLUSTERING_MAIN: Filtering worker layout -",
          active_filter_workers, "cluster workers x",
          filter_gamma_workers, "gamma workers"
        )
      )
      if (in_parallel_context) {
        scice_message("CLUSTERING_MAIN: Running in parallel context - using 1 worker for nested operations")
      }
    }
    
    cluster_filter_results <- cross_platform_mclapply(cluster_range, function(cluster_num) {
      if (!(as.character(cluster_num) %in% names(gamma_dict))) {
        return(list(excluded = TRUE, cluster_num = cluster_num, reason = "resolution_search_failed"))
      }
      
      gamma_range <- gamma_dict[[as.character(cluster_num)]]
      gamma_test_len <- min(5L, max(1L, as.integer(ceiling(abs(diff(gamma_range)) * 100)) + 1L))
      gamma_test <- seq(gamma_range[1], gamma_range[2], length.out = gamma_test_len)
      
      # Test multiple gammas in parallel (distribute workers efficiently when in parallel context)
      nested_workers <- if (in_parallel_context) {
        max(1L, as.integer(round(actual_workers / max(1L, length(cluster_range)))))
      } else {
        filter_gamma_workers
      }
      nested_workers <- cap_workers_by_memory(
        nested_workers,
        estimate_trial_matrix_bytes(igraph::vcount(igraph_obj), 10L, 1L),
        runtime_context
      )
      
      if (verbose && in_parallel_context) {
        scice_message(paste("CLUSTERING_MAIN: Nested worker optimization - using", nested_workers, 
                     "workers per cluster ( ", actual_workers, "total /", length(cluster_range), "clusters)"))
      }
      ic_scores <- cross_platform_mclapply(gamma_test, function(gamma_val) {
        # Set deterministic seed for filtering if base seed provided
        if (!is.null(seed)) {
          gamma_component <- if (is.finite(gamma_val)) {
            as.integer(floor(abs(gamma_val * 1000) %% 10000))
          } else {
            0L
          }
          filter_seed <- as.integer((as.double(seed) + cluster_num * 100 + gamma_component) %% .Machine$integer.max)
          if (is.na(filter_seed) || filter_seed <= 0L) {
            filter_seed <- 1L
          }
          set.seed(filter_seed)
        }
        
        cluster_results <- replicate(10, {
          cached_leiden_clustering(igraph_obj, gamma_val, objective_function, 5, 0.01,
                                 cache_key_suffix = paste("filter", cluster_num, sep = "_"),
                                 snn_graph = snn_graph, min_cluster_size = min_cluster_size)
        }, simplify = TRUE)
        
        extracted_results <- extract_clustering_array(cluster_results)
        ic_result <- calculate_ic_from_extracted(extracted_results)
        return(1 / ic_result)
      }, mc.cores = nested_workers)
      
      ic_scores <- unlist(ic_scores)
      excluded <- min(ic_scores, na.rm = TRUE) >= remove_threshold
      reason <- if (excluded) "high_inconsistency" else "passed_filtering"
      
      list(excluded = excluded, cluster_num = cluster_num, reason = reason)
    }, mc.cores = actual_workers)
    
    excluded_numbers <- sapply(cluster_filter_results, function(x) if (x$excluded) x$cluster_num else NULL)
    excluded_numbers <- unlist(excluded_numbers)
    
    if (verbose) {
      filtering_time <- as.numeric(difftime(Sys.time(), filtering_start, units = "secs"))
      scice_message(paste("CLUSTERING_MAIN: Filtering completed in", round(filtering_time, 3), "seconds"))
      scice_message(paste("CLUSTERING_MAIN: Processed", length(cluster_filter_results), "cluster numbers"))
      if (length(excluded_numbers) > 0) {
        scice_message(paste("CLUSTERING_MAIN: Excluded", length(excluded_numbers), "cluster numbers:", paste(excluded_numbers, collapse = ", ")))
        # Report exclusion reasons
        exclusion_reasons <- sapply(cluster_filter_results, function(x) if (x$excluded) paste(x$cluster_num, ":", x$reason) else NULL)
        exclusion_reasons <- unlist(exclusion_reasons)
        if (length(exclusion_reasons) > 0) {
          scice_message("CLUSTERING_MAIN: Exclusion reasons:")
          for (reason in exclusion_reasons) {
            scice_message(paste("CLUSTERING_MAIN:   ", reason))
          }
        }
      } else {
        scice_message("CLUSTERING_MAIN: No clusters excluded during filtering")
      }
    }
  }
  
  valid_clusters <- setdiff(cluster_range, excluded_numbers)
  
  if (verbose) {
    scice_message(paste("CLUSTERING_MAIN: Valid clusters for optimization:", length(valid_clusters)))
    if (length(valid_clusters) > 0) {
      scice_message(paste("CLUSTERING_MAIN: Valid cluster numbers:", paste(valid_clusters, collapse = ", ")))
    }
  }
  
  # Create entries for excluded clusters
  excluded_entries <- lapply(cluster_filter_results, function(x) {
    if (x$excluded) {
      return(data.table::data.table(
        cluster_number = x$cluster_num,
        gamma = NA_real_,
        labels = list(NULL),
        ic = NA_real_,
        ic_vec = list(NULL),
        best_labels = list(NULL),
        effective_cluster_median = NA_real_,
        raw_cluster_median = NA_real_,
        final_cluster_median = NA_real_,
        admission_mode = NA_character_,
        best_labels_raw_cluster_count = NA_integer_,
        best_labels_final_cluster_count = NA_integer_,
        n_iter = NA_integer_,
        mei = list(NULL),
        k = NA_integer_,
        source_target_cluster = x$cluster_num,
        excluded = TRUE,
        exclusion_reason = x$reason,
        selected_main_result = FALSE,
        result_status = x$reason,
        phase1_primary_gamma_count = NA_integer_,
        phase1_secondary_gamma_count = NA_integer_,
        phase1_total_gamma_count = NA_integer_,
        phase1_elapsed_sec = NA_real_,
        phase1_leiden_runs = NA_integer_,
        secondary_phase1_used = NA,
        exact_hit_gamma_count = NA_integer_,
        phase4_iterations = NA_integer_,
        phase4_elapsed_sec = NA_real_,
        phase5_elapsed_sec = NA_real_,
        optimization_elapsed_sec = NA_real_
      ))
    }
    return(NULL)
  })
  
  # Handle excluded entries properly
  excluded_entries_list <- excluded_entries[!sapply(excluded_entries, is.null)]
  if (length(excluded_entries_list) > 0) {
    excluded_entries <- data.table::rbindlist(excluded_entries_list)
  } else {
    excluded_entries <- empty_target_results_dt()
  }
  
  if (length(valid_clusters) == 0) {
    if (verbose) {
      scice_message("CLUSTERING_MAIN: WARNING - No valid cluster numbers found!")
      scice_message("CLUSTERING_MAIN: All clusters were excluded during filtering")
      scice_message("CLUSTERING_MAIN: This will result in empty IC results")
      scice_message("CLUSTERING_MAIN: Returning results with only excluded cluster information")
    }
    rekeyed_results <- rekey_target_results_by_final_cluster(excluded_entries)
    result <- cluster_results_dt_to_list(
      rekeyed_results$main_results_dt,
      target_results_dt = rekeyed_results$target_results_dt
    )
    result$resolution_search_diagnostics <- resolution_search_diagnostics
    result$search_coverage_complete <- search_coverage_complete
    return(result)
  }
  
  # Optimize clustering for each valid cluster number in parallel
  if (verbose) {
    scice_message("CLUSTERING_MAIN: Starting clustering optimization...")
    scice_message(paste("CLUSTERING_MAIN: Optimizing", length(valid_clusters), "cluster numbers"))
    optimization_start <- Sys.time()
  }
  
  # Windows compatibility for clustering optimization
  actual_workers_opt <- if (.Platform$OS.type == "windows" && n_workers > 1) 1 else n_workers
  active_cluster_workers_opt <- 1L
  cluster_worker_budget <- 1L

  if (length(valid_clusters) > 0) {
    n_vertices <- igraph::vcount(igraph_obj)
    gamma_steps_estimates <- vapply(valid_clusters, function(cluster_num) {
      gamma_range <- gamma_dict[[as.character(cluster_num)]]
      range_width <- abs(gamma_range[2] - gamma_range[1])
      if (n_vertices >= 200000 &&
          range_width <= max(resolution_tolerance * 10, .Machine$double.eps)) {
        5L
      } else {
        11L
      }
    }, integer(1))

    max_cluster_bytes <- max(estimate_trial_matrix_bytes(
      n_cells = n_vertices,
      n_trials = n_trials,
      n_gamma = gamma_steps_estimates
    ))
    actual_workers_opt <- cap_workers_by_memory(actual_workers_opt, max_cluster_bytes, runtime_context)
    active_cluster_workers_opt <- max(1L, min(as.integer(length(valid_clusters)), as.integer(actual_workers_opt)))
    cluster_worker_budget <- max(1L, as.integer(actual_workers_opt))

    if (should_enable_spill(runtime_context, max_cluster_bytes)) {
      activate_runtime_spill(runtime_context, estimated_bytes = max_cluster_bytes)
    }
  }
  
  if (verbose) {
    scice_message(paste("CLUSTERING_MAIN: Using", actual_workers_opt, "workers for optimization"))
    scice_message(paste("CLUSTERING_MAIN: Active cluster workers:", active_cluster_workers_opt))
    scice_message(paste("CLUSTERING_MAIN: Per-cluster worker budget:", cluster_worker_budget))
    scice_message("CLUSTERING_MAIN: Optimization scheduling - dynamic worker queue (mc.preschedule = FALSE)")
    scice_message(paste("CLUSTERING_MAIN: Progress tracking:"))
    progress_start_time <- Sys.time()
    
    # Pre-calculate estimated time per cluster for progress tracking
    if (length(valid_clusters) > 1) {
      estimated_time_per_cluster <- (n_trials * length(valid_clusters) * n_bootstrap) / (actual_workers_opt * 100)  # rough estimate
      total_estimated_time <- estimated_time_per_cluster * length(valid_clusters)
      scice_message(paste("CLUSTERING_MAIN: Estimated total time:", round(total_estimated_time, 1), "seconds"))
    }
  }
  
  cluster_results <- cross_platform_mclapply(valid_clusters, function(cluster_num) {
    # Thread-local logging for parallel workers
    worker_id <- paste("WORKER", cluster_num)
    
    if (verbose) {
      scice_message(paste(worker_id, ": Starting optimization for k =", cluster_num))
      scice_message(paste(worker_id, ": Thread context - PID:", Sys.getpid()))
    }
    
    if (!(as.character(cluster_num) %in% names(gamma_dict))) {
      if (verbose) {
        scice_message(paste(worker_id, ": ERROR - No gamma range found for k =", cluster_num))
      }
      return(NULL)
    }
    
    gamma_range <- gamma_dict[[as.character(cluster_num)]]
    
    if (verbose) {
      scice_message(paste(worker_id, ": Gamma range [", signif(gamma_range[1], 6), ", ", signif(gamma_range[2], 6), "]", sep = ""))
      scice_message(paste(worker_id, ": Starting intensive optimization..."))
      opt_start_time <- Sys.time()
    }
    
    # Optimize clustering within this range
    gamma_seed_table <- build_target_gamma_seed_table(
      target_cluster = cluster_num,
      gamma_dict = gamma_dict,
      target_gamma_seeds = target_gamma_seeds,
      target_interval_details = target_interval_details,
      resolution_search_diagnostics = resolution_search_diagnostics
    )
    result <- optimize_clustering(
      igraph_obj, cluster_num, gamma_range, objective_function,
      n_trials, n_bootstrap, seed, beta, n_iterations, max_iterations,
      resolution_tolerance, cluster_worker_budget, snn_graph,
      gamma_seed_table, min_cluster_size, verbose,
      worker_id, in_parallel_context = TRUE,
      runtime_context = runtime_context
    )
    
    if (verbose) {
      opt_time <- as.numeric(difftime(Sys.time(), opt_start_time, units = "secs"))
      scice_message(paste(worker_id, ": Optimization completed in", round(opt_time, 3), "seconds"))
    }
    
    if (!is.null(result)) {
      if (!isTRUE(result$success)) {
        return(data.table::data.table(
          cluster_number = cluster_num,
          gamma = result$gamma,
          labels = list(NULL),
          ic = NA_real_,
          ic_vec = list(NULL),
          best_labels = list(NULL),
          effective_cluster_median = result$effective_cluster_median,
          raw_cluster_median = result$raw_cluster_median,
          final_cluster_median = result$final_cluster_median,
          admission_mode = result$admission_mode,
          best_labels_raw_cluster_count = result$best_labels_raw_cluster_count,
          best_labels_final_cluster_count = result$best_labels_final_cluster_count,
          n_iter = result$n_iterations,
          mei = list(NULL),
          k = result$k,
          source_target_cluster = cluster_num,
          excluded = TRUE,
          exclusion_reason = result$failure_reason,
          selected_main_result = FALSE,
          result_status = result$failure_reason,
          phase1_primary_gamma_count = result$phase1_primary_gamma_count,
          phase1_secondary_gamma_count = result$phase1_secondary_gamma_count,
          phase1_total_gamma_count = result$phase1_total_gamma_count,
          phase1_elapsed_sec = result$phase1_elapsed_sec,
          phase1_leiden_runs = result$phase1_leiden_runs,
          secondary_phase1_used = result$secondary_phase1_used,
          exact_hit_gamma_count = result$exact_hit_gamma_count,
          phase4_iterations = result$phase4_iterations,
          phase4_elapsed_sec = result$phase4_elapsed_sec,
          phase5_elapsed_sec = result$phase5_elapsed_sec,
          optimization_elapsed_sec = result$optimization_elapsed_sec
        ))
      }
      # Calculate MEI scores
      mei_scores <- calculate_mei_from_array(result$labels)
      
      # Return as a data.table row
      return(data.table::data.table(
        cluster_number = cluster_num,
        gamma = result$gamma,
        labels = list(result$labels),
        ic = result$ic_median,
        ic_vec = list(result$ic_bootstrap),
        best_labels = list(result$best_labels),
        effective_cluster_median = result$effective_cluster_median,
        raw_cluster_median = result$raw_cluster_median,
        final_cluster_median = result$final_cluster_median,
        admission_mode = result$admission_mode,
        best_labels_raw_cluster_count = result$best_labels_raw_cluster_count,
        best_labels_final_cluster_count = result$best_labels_final_cluster_count,
        n_iter = result$n_iterations,
        mei = list(mei_scores),
        k = result$k,
        source_target_cluster = cluster_num,
        excluded = FALSE,
        exclusion_reason = "none",
        selected_main_result = FALSE,
        result_status = "candidate",
        phase1_primary_gamma_count = result$phase1_primary_gamma_count,
        phase1_secondary_gamma_count = result$phase1_secondary_gamma_count,
        phase1_total_gamma_count = result$phase1_total_gamma_count,
        phase1_elapsed_sec = result$phase1_elapsed_sec,
        phase1_leiden_runs = result$phase1_leiden_runs,
        secondary_phase1_used = result$secondary_phase1_used,
        exact_hit_gamma_count = result$exact_hit_gamma_count,
        phase4_iterations = result$phase4_iterations,
        phase4_elapsed_sec = result$phase4_elapsed_sec,
        phase5_elapsed_sec = result$phase5_elapsed_sec,
        optimization_elapsed_sec = result$optimization_elapsed_sec
      ))
    }
    return(NULL)
  }, mc.cores = active_cluster_workers_opt, mc.preschedule = FALSE)
  
  if (verbose) {
    optimization_time <- as.numeric(difftime(Sys.time(), optimization_start, units = "secs"))
    scice_message(paste("CLUSTERING_MAIN: All optimization workers completed in", round(optimization_time, 3), "seconds"))
    successful_count <- sum(!sapply(cluster_results, is.null))
    scice_message(paste("CLUSTERING_MAIN: Successful optimizations:", successful_count, "/", length(valid_clusters)))
  }
  
  # Combine successful results
  successful_results_list <- cluster_results[!sapply(cluster_results, is.null)]
  if (length(successful_results_list) > 0) {
    successful_results <- data.table::rbindlist(successful_results_list, fill = TRUE)
  } else {
    successful_results <- empty_target_results_dt()
  }
  
  # Combine all target-level results (excluded + successful)
  tryCatch({
    result_parts <- list(excluded_entries, successful_results)
    result_parts <- result_parts[vapply(result_parts, nrow, integer(1)) > 0L]
    if (length(result_parts) > 0L) {
      results_dt <- data.table::rbindlist(result_parts, fill = TRUE)
    } else {
      results_dt <- empty_target_results_dt()
    }
  }, error = function(e) {
    warning("Error combining results, creating empty results: ", e$message)
    results_dt <<- empty_target_results_dt()
  })
  
  # Sort by cluster number (only if the column exists and there are rows)
  if (nrow(results_dt) > 0 && "cluster_number" %in% colnames(results_dt)) {
    data.table::setorder(results_dt, cluster_number)
  }

  if (verbose &&
      nrow(results_dt) > 0 &&
      all(c(
        "cluster_number", "gamma", "effective_cluster_median",
        "raw_cluster_median", "admission_mode",
        "best_labels_raw_cluster_count", "best_labels_final_cluster_count",
        "excluded", "selected_main_result",
        "phase1_total_gamma_count", "phase4_iterations",
        "optimization_elapsed_sec"
      ) %in% colnames(results_dt))) {
    format_diag_value <- function(x, digits = 6L) {
      if (length(x) == 0L || is.na(x)) {
        return("NA")
      }
      as.character(signif(x, digits))
    }
    scice_message("CLUSTERING_MAIN: Per-k selection diagnostics:")
    for (row_idx in seq_len(nrow(results_dt))) {
      cluster_num <- results_dt$cluster_number[[row_idx]]
      if (isTRUE(results_dt$excluded[[row_idx]])) {
        scice_message(
          paste(
            "CLUSTERING_MAIN:   k =", cluster_num,
            "- excluded = TRUE",
            "- reason =", results_dt$exclusion_reason[[row_idx]]
          )
        )
        next
      }
      scice_message(
        paste(
          "CLUSTERING_MAIN:   k =", cluster_num,
          "- gamma =", format_diag_value(results_dt$gamma[[row_idx]]),
          "- effective_median =", format_diag_value(results_dt$effective_cluster_median[[row_idx]]),
          "- raw_median =", format_diag_value(results_dt$raw_cluster_median[[row_idx]]),
          "- admission_mode =", ifelse(
            is.na(results_dt$admission_mode[[row_idx]]),
            "NA",
            results_dt$admission_mode[[row_idx]]
          ),
          "- best_labels_raw_clusters =", ifelse(
            is.na(results_dt$best_labels_raw_cluster_count[[row_idx]]),
            "NA",
            as.character(results_dt$best_labels_raw_cluster_count[[row_idx]])
          ),
          "- best_labels_final_clusters =", ifelse(
            is.na(results_dt$best_labels_final_cluster_count[[row_idx]]),
            "NA",
            as.character(results_dt$best_labels_final_cluster_count[[row_idx]])
          ),
          "- phase1_gammas =", ifelse(
            is.na(results_dt$phase1_total_gamma_count[[row_idx]]),
            "NA",
            as.character(results_dt$phase1_total_gamma_count[[row_idx]])
          ),
          "- phase4_iterations =", ifelse(
            is.na(results_dt$phase4_iterations[[row_idx]]),
            "NA",
            as.character(results_dt$phase4_iterations[[row_idx]])
          ),
          "- optimization_elapsed_sec =", ifelse(
            is.na(results_dt$optimization_elapsed_sec[[row_idx]]),
            "NA",
            as.character(round(results_dt$optimization_elapsed_sec[[row_idx]], 3))
          ),
          "- selected_main_result =", results_dt$selected_main_result[[row_idx]]
        )
      )
    }
  }

  rekeyed_results <- rekey_target_results_by_final_cluster(results_dt)

  if (verbose &&
      nrow(rekeyed_results$main_results_dt) > 0 &&
      all(c("cluster_number", "source_target_cluster", "ic") %in% colnames(rekeyed_results$main_results_dt))) {
    scice_message("CLUSTERING_MAIN: Final-count keyed main results:")
    for (row_idx in seq_len(nrow(rekeyed_results$main_results_dt))) {
      scice_message(
        paste(
          "CLUSTERING_MAIN:   final k =", rekeyed_results$main_results_dt$cluster_number[[row_idx]],
          "- source target =", rekeyed_results$main_results_dt$source_target_cluster[[row_idx]],
          "- IC =", round(rekeyed_results$main_results_dt$ic[[row_idx]], 4),
          "- raw clusters =", rekeyed_results$main_results_dt$best_labels_raw_cluster_count[[row_idx]],
          "- final clusters =", rekeyed_results$main_results_dt$best_labels_final_cluster_count[[row_idx]]
        )
      )
    }
  }
  
  # Report cache statistics for performance monitoring
  if (verbose) {
    cache_stats_final <- get_cache_stats()
    scice_message(paste("CLUSTERING_MAIN: Final cache entries:", cache_stats_final$cache_entries))
    scice_message(paste("CLUSTERING_MAIN: Cache provided", cache_stats_final$cache_entries, "reused clustering results"))
  }
  
  result <- cluster_results_dt_to_list(
    rekeyed_results$main_results_dt,
    target_results_dt = rekeyed_results$target_results_dt
  )
  result$resolution_search_diagnostics <- resolution_search_diagnostics
  result$search_coverage_complete <- search_coverage_complete
  return(result)
}
