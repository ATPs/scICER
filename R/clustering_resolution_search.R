#' Reindex cluster labels to contiguous 0-based integers
#'
#' @description
#' Normalizes arbitrary integer cluster identifiers to \code{0:(k-1)}.
#'
#' @param labels Integer cluster labels
#' @return Integer labels reindexed to contiguous 0-based values
#' @keywords internal
reindex_cluster_labels <- function(labels) {
  labels <- as.integer(labels)
  if (length(labels) == 0L) {
    return(labels)
  }

  unique_ids <- sort(unique(labels))
  if (length(unique_ids) == 0L) {
    return(integer(0))
  }

  expected_ids <- seq.int(0L, length(unique_ids) - 1L)
  if (identical(unique_ids, expected_ids)) {
    return(labels)
  }

  id_map <- stats::setNames(expected_ids, as.character(unique_ids))
  as.integer(unname(id_map[as.character(labels)]))
}

#' Count effective clusters under a minimum-size rule
#'
#' @description
#' Counts clusters whose size is at least \code{min_cluster_size}.
#'
#' @details
#' This helper is used for target-cluster matching during resolution search and
#' optimization. Clusters smaller than the threshold are ignored. When all
#' clusters are undersized, the function returns \code{0}.
#'
#' @param labels Integer cluster labels
#' @param min_cluster_size Minimum required cells per counted cluster
#' @return Integer count of effective clusters
#' @keywords internal
count_effective_clusters <- function(labels, min_cluster_size = 1L) {
  labels <- as.integer(labels)
  min_cluster_size <- max(1L, as.integer(min_cluster_size))

  if (length(labels) == 0L) {
    return(0L)
  }

  if (min_cluster_size <= 1L) {
    return(as.integer(length(unique(labels))))
  }

  cluster_sizes <- tabulate(reindex_cluster_labels(labels) + 1L)
  as.integer(sum(cluster_sizes >= min_cluster_size))
}

# Raw-cluster guard thresholds for effective-cluster matching.
raw_cluster_guard_limits <- function(target_clusters) {
  target_clusters <- max(1L, as.integer(target_clusters))
  list(
    soft = as.integer(max(target_clusters + 3L, ceiling(target_clusters * 1.1))),
    hard = as.integer(max(target_clusters + 5L, ceiling(target_clusters * 1.5)))
  )
}

raw_cluster_search_upper <- function(target_clusters) {
  target_clusters <- max(1L, as.integer(target_clusters))
  if (target_clusters <= 10L) {
    return(as.integer(target_clusters + 1L))
  }
  as.integer(max(target_clusters + 2L, ceiling(target_clusters * 1.05)))
}

build_gamma_sequence_for_range <- function(gamma_range, objective_function,
                                           resolution_tolerance,
                                           n_vertices = NA_integer_) {
  gamma_range <- as.numeric(gamma_range)
  if (length(gamma_range) != 2L || anyNA(gamma_range) || any(!is.finite(gamma_range))) {
    stop("gamma_range must contain two finite numeric bounds.")
  }

  range_width <- abs(gamma_range[2] - gamma_range[1])
  n_steps <- if (!is.na(n_vertices) &&
                 n_vertices >= 200000 &&
                 range_width <= max(resolution_tolerance * 10, .Machine$double.eps)) {
    5L
  } else {
    11L
  }

  if (objective_function == "modularity") {
    if (abs(gamma_range[2] - gamma_range[1]) > resolution_tolerance) {
      return(seq(gamma_range[1], gamma_range[2], length.out = n_steps))
    }
    delta_g <- resolution_tolerance
    return(seq(gamma_range[1] - delta_g, gamma_range[1] + delta_g, length.out = n_steps))
  }

  lower <- max(min(gamma_range), .Machine$double.xmin)
  upper <- max(max(gamma_range), .Machine$double.xmin)
  if (abs(upper - lower) > max(resolution_tolerance, lower * 1e-6)) {
    return(exp(seq(log(lower), log(upper), length.out = n_steps)))
  }

  delta_log <- max(resolution_tolerance, 1e-4)
  exp(seq(log(lower) - delta_log, log(lower) + delta_log, length.out = n_steps))
}

stabilize_probe_raw_medians <- function(raw_cluster_medians) {
  raw_cluster_medians <- as.numeric(raw_cluster_medians)
  finite_indices <- which(is.finite(raw_cluster_medians))
  if (length(finite_indices) <= 1L) {
    return(raw_cluster_medians)
  }

  raw_cluster_medians[finite_indices] <- cummax(raw_cluster_medians[finite_indices])
  raw_cluster_medians
}

clamp_gamma_range_to_raw_plateau <- function(gamma_sequence, raw_cluster_medians,
                                             target_clusters) {
  gamma_sequence <- as.numeric(gamma_sequence)
  raw_cluster_medians <- as.numeric(raw_cluster_medians)
  target_clusters <- as.integer(target_clusters)

  if (length(gamma_sequence) != length(raw_cluster_medians) || length(gamma_sequence) == 0L) {
    stop("gamma_sequence and raw_cluster_medians must have the same non-zero length.")
  }

  stabilized_raw_medians <- stabilize_probe_raw_medians(raw_cluster_medians)
  coarse_bounds <- range(gamma_sequence)
  exact_indices <- which(
    !is.na(stabilized_raw_medians) & stabilized_raw_medians == target_clusters
  )
  if (length(exact_indices) > 0L) {
    return(list(
      bounds = range(gamma_sequence[exact_indices]),
      mode = "raw_exact",
      indices = exact_indices,
      stabilized_raw_medians = stabilized_raw_medians
    ))
  }

  left_raw <- stabilized_raw_medians[-length(stabilized_raw_medians)]
  right_raw <- stabilized_raw_medians[-1L]
  crossing_indices <- which(
    !is.na(left_raw) &
      !is.na(right_raw) &
      ((left_raw < target_clusters & right_raw > target_clusters) |
         (left_raw > target_clusters & right_raw < target_clusters))
  )
  if (length(crossing_indices) > 0L) {
    widths <- abs(gamma_sequence[crossing_indices + 1L] - gamma_sequence[crossing_indices])
    best_pair_idx <- crossing_indices[[which.min(widths)]]
    return(list(
      bounds = sort(gamma_sequence[c(best_pair_idx, best_pair_idx + 1L)]),
      mode = "raw_bracket",
      indices = c(best_pair_idx, best_pair_idx + 1L),
      stabilized_raw_medians = stabilized_raw_medians
    ))
  }

  near_target_indices <- which(
    !is.na(stabilized_raw_medians) & abs(stabilized_raw_medians - target_clusters) <= 1
  )
  if (length(near_target_indices) > 0L) {
    return(list(
      bounds = range(gamma_sequence[near_target_indices]),
      mode = "raw_near_target",
      indices = near_target_indices,
      stabilized_raw_medians = stabilized_raw_medians
    ))
  }

  list(
    bounds = coarse_bounds,
    mode = "coarse",
    indices = seq_along(gamma_sequence),
    stabilized_raw_medians = stabilized_raw_medians
  )
}

# Reject pathological effective matches with absurd raw cluster counts.
passes_raw_cluster_guard <- function(raw_cluster_median, target_clusters,
                                     min_cluster_size = 1L,
                                     level = c("soft", "hard")) {
  level <- match.arg(level)
  target_clusters <- as.integer(target_clusters)
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  if (min_cluster_size <= 1L || is.na(target_clusters) || target_clusters <= 1L) {
    return(rep(TRUE, length(raw_cluster_median)))
  }

  raw_cluster_median <- as.numeric(raw_cluster_median)
  limits <- raw_cluster_guard_limits(target_clusters)
  upper_bound <- if (identical(level, "soft")) limits$soft else limits$hard

  !is.na(raw_cluster_median) & raw_cluster_median <= upper_bound
}

# Classify one resolution-search gamma using raw/effective median counts.
classify_resolution_search_state <- function(raw_cluster_median, effective_cluster_median,
                                             target_clusters, min_cluster_size = 1L) {
  target_clusters <- as.integer(target_clusters)
  raw_cluster_median <- as.numeric(raw_cluster_median)
  effective_cluster_median <- as.numeric(effective_cluster_median)
  raw_guard_soft <- passes_raw_cluster_guard(
    raw_cluster_median,
    target_clusters,
    min_cluster_size = min_cluster_size,
    level = "soft"
  )
  raw_guard_search <- if (min_cluster_size <= 1L || is.na(target_clusters) || target_clusters <= 1L) {
    TRUE
  } else {
    !is.na(raw_cluster_median) && raw_cluster_median <= raw_cluster_search_upper(target_clusters)
  }
  raw_below <- !is.na(raw_cluster_median) && raw_cluster_median < target_clusters
  raw_above_soft <- !raw_below && !isTRUE(raw_guard_search)
  raw_in_band <- !raw_below && !raw_above_soft
  effective_meets_target <- !is.na(effective_cluster_median) &&
    effective_cluster_median >= target_clusters
  over_fragmented <- !raw_below && !effective_meets_target
  raw_class <- if (raw_below) {
    "raw_below"
  } else if (raw_above_soft) {
    "raw_above_soft"
  } else {
    "raw_in_band"
  }
  lower_action <- if (raw_below) "increase_gamma" else "decrease_gamma"
  upper_action <- if (raw_below || (raw_in_band && effective_meets_target)) {
    "increase_gamma"
  } else {
    "decrease_gamma"
  }
  list(
    raw_class = raw_class,
    raw_below = raw_below,
    raw_in_band = raw_in_band,
    raw_above_soft = raw_above_soft,
    raw_guard_soft = raw_guard_soft,
    raw_guard_search = raw_guard_search,
    effective_meets_target = effective_meets_target,
    over_fragmented = over_fragmented,
    lower_action = lower_action,
    upper_action = upper_action
  )
}

#' Find resolution parameter ranges for each cluster number using binary search
#'
#' @description
#' Estimates lower/upper resolution bounds that reproduce each target effective cluster count.
#'
#' @details
#' For each requested \code{k}, this helper runs two binary searches:
#' one for the smallest gamma yielding at least \code{k} effective clusters and
#' one for the largest gamma yielding at most \code{k} effective clusters
#' (CPM/modularity aware).
#' Effective clusters are those with size \code{>= min_cluster_size}; all-small
#' trials are counted as \code{0}.
#'
#' Large graphs use lighter preliminary settings to reduce runtime pressure.
#' After per-k search, adjacent ranges are lightly adjusted when nearly touching
#' to reduce ambiguity between neighboring target counts.
#'
#' @param igraph_obj igraph object to cluster
#' @param cluster_range Target cluster numbers to search
#' @param start_g Lower bound of search interval (log scale for CPM path)
#' @param end_g Upper bound of search interval (log scale for CPM path)
#' @param objective_function Objective function ("modularity" or "CPM")
#' @param resolution_tolerance Resolution convergence tolerance
#' @param n_workers Requested worker count
#' @param verbose Whether to print progress diagnostics
#' @param seed Optional deterministic seed
#' @param snn_graph Seurat SNN graph matrix for final best-label merging
#' @param min_cluster_size Minimum required cells per effective cluster count
#' @param in_parallel_context Whether called from an outer parallel worker
#' @param runtime_context Internal mutable context for spill/memory controls
#' @return Named list mapping cluster number to \code{c(left_bound, right_bound)}
#' @keywords internal
find_resolution_ranges <- function(igraph_obj, cluster_range, start_g, end_g,
                                  objective_function, resolution_tolerance, n_workers, verbose, seed = NULL,
                                  snn_graph = NULL, min_cluster_size = 1L,
                                  in_parallel_context = FALSE, runtime_context = NULL) {
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  if (min_cluster_size > 1L && is.null(snn_graph)) {
    stop("snn_graph must be provided when min_cluster_size > 1.")
  }
  
  # Initialize results storage with data.table
  results_dt <- data.table::data.table(
    cluster_number = integer(),
    left_bound = numeric(),
    right_bound = numeric()
  )
  
  n_vertices <- igraph::vcount(igraph_obj)
  n_preliminary_trials <- if (n_vertices >= 200000) {
    3L
  } else if (n_vertices >= 100000) {
    5L
  } else {
    15L
  }
  beta_preliminary <- 0.01
  n_iter_preliminary <- if (n_vertices >= 200000) 3L else 5L
  
  if (verbose) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Using", n_preliminary_trials,
        "preliminary trials per step (graph vertices:", n_vertices, ")"
      )
    )
  }
  
  max_search_iterations <- if (n_vertices >= 200000) 30L else 50L
  effective_tolerance <- if (n_vertices >= 200000) {
    max(resolution_tolerance, .Machine$double.eps * 100)
  } else {
    resolution_tolerance / 10
  }
  
  # Process cluster numbers in parallel (with Windows compatibility)
  requested_search_workers <- max(1L, as.integer(n_workers))
  if (.Platform$OS.type == "windows" && requested_search_workers > 1L) {
    requested_search_workers <- 1L  # Force single worker on Windows
  }

  n_cluster_targets <- max(1L, as.integer(length(cluster_range)))
  # In nested context, split the provided budget across cluster targets.
  search_worker_budget <- if (in_parallel_context) {
    max(1L, as.integer(round(as.double(requested_search_workers) / n_cluster_targets)))
  } else {
    requested_search_workers
  }
  search_worker_budget <- cap_workers_by_memory(
    search_worker_budget,
    estimate_trial_matrix_bytes(n_vertices, n_preliminary_trials, 1L),
    runtime_context
  )

  # Outer RESOLUTION_SEARCH workers should never exceed the number of targets.
  active_cluster_workers <- max(1L, min(n_cluster_targets, as.integer(search_worker_budget)))
  preliminary_trial_worker_capacity <- if (in_parallel_context) {
    1L
  } else {
    max(1L, as.integer(floor(as.double(search_worker_budget) / active_cluster_workers)))
  }
  preliminary_trial_worker_capacity <- min(preliminary_trial_worker_capacity, n_preliminary_trials)
  preliminary_trial_worker_capacity <- cap_workers_by_memory(
    preliminary_trial_worker_capacity,
    estimate_trial_matrix_bytes(n_vertices, 1L, 1L),
    runtime_context
  )

  min_parallel_preliminary_workers <- getOption("scICER.internal_preliminary_parallel_min_workers", 3L)
  if (!is.numeric(min_parallel_preliminary_workers) ||
      length(min_parallel_preliminary_workers) != 1L ||
      !is.finite(min_parallel_preliminary_workers) ||
      min_parallel_preliminary_workers < 2) {
    min_parallel_preliminary_workers <- 3L
  } else {
    min_parallel_preliminary_workers <- as.integer(round(min_parallel_preliminary_workers))
  }
  max_parallel_preliminary_workers <- getOption("scICER.internal_preliminary_parallel_max_workers", 4L)
  if (!is.numeric(max_parallel_preliminary_workers) ||
      length(max_parallel_preliminary_workers) != 1L ||
      !is.finite(max_parallel_preliminary_workers) ||
      max_parallel_preliminary_workers < min_parallel_preliminary_workers) {
    max_parallel_preliminary_workers <- max(4L, min_parallel_preliminary_workers)
  } else {
    max_parallel_preliminary_workers <- as.integer(round(max_parallel_preliminary_workers))
  }

  use_parallel_preliminary_trials <- !in_parallel_context &&
    preliminary_trial_worker_capacity >= min_parallel_preliminary_workers &&
    n_preliminary_trials >= min_parallel_preliminary_workers

  preliminary_trial_workers <- if (use_parallel_preliminary_trials) {
    min(preliminary_trial_worker_capacity, n_preliminary_trials, max_parallel_preliminary_workers)
  } else {
    1L
  }
  
  if (verbose) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Worker allocation - requested:", as.integer(n_workers),
        "| budget after context/memory:", search_worker_budget,
        "| active cluster workers:", active_cluster_workers,
        "| preliminary worker capacity per gamma:", preliminary_trial_worker_capacity,
        "| preliminary workers per gamma:", preliminary_trial_workers
      )
    )
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Preliminary trial workers per gamma:",
        preliminary_trial_workers
      )
    )
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Preliminary trial mode:",
        if (preliminary_trial_workers > 1L) "parallel" else "serial",
        "(parallel requires worker capacity >=",
        min_parallel_preliminary_workers, ")"
      )
    )
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Preliminary trial strategy - reuse one representative preliminary clustering per gamma step and summarize nominal median effective/raw counts."
      )
    )
  }
  
  range_results <- cross_platform_mclapply(cluster_range, function(target_clusters) {
    trial_steps <- 0L
    total_preliminary_elapsed <- 0
    heartbeat <- create_heartbeat_logger(
      verbose = verbose,
      context = paste("RESOLUTION_SEARCH: k =", target_clusters)
    )
    
    derive_trial_seed <- function(seed_base, trial_idx) {
      if (is.null(seed_base)) {
        return(NULL)
      }
      trial_seed <- as.integer((as.double(seed_base) + trial_idx) %% .Machine$integer.max)
      if (is.na(trial_seed) || trial_seed <= 0L) {
        trial_seed <- 1L
      }
      trial_seed
    }
    
    run_single_trial_count <- function(trial_idx, gamma_val, cache_suffix_base, seed_base,
                                       use_cache = FALSE) {
      trial_seed <- derive_trial_seed(seed_base, trial_idx)
      if (!is.null(trial_seed)) {
        set.seed(trial_seed)
      }
      cache_key_suffix <- if (isTRUE(use_cache)) {
        cache_suffix_base
      } else {
        paste(cache_suffix_base, "trial", trial_idx, sep = "_")
      }
      labels <- cached_leiden_clustering(
        igraph_obj,
        gamma_val,
        objective_function,
        n_iter_preliminary,
        beta_preliminary,
        use_cache = use_cache,
        cache_key_suffix = cache_key_suffix,
        snn_graph = snn_graph,
        min_cluster_size = min_cluster_size
      )
      effective_count <- count_effective_clusters(labels, min_cluster_size = min_cluster_size)
      raw_count <- as.integer(length(unique(labels)))
      c(
        effective_count = as.integer(effective_count),
        raw_count = raw_count
      )
    }
    
    run_preliminary_trials <- function(gamma_val, cache_suffix_base, seed_base = NULL,
                                       phase_label, iteration_idx) {
      step_start_time <- Sys.time()
      trial_counts <- integer(n_preliminary_trials)
      raw_trial_counts <- integer(n_preliminary_trials)
      unpack_trial_counts <- function(trial_value, trial_idx_for_error = NA_integer_) {
        parsed <- suppressWarnings(as.numeric(trial_value))
        if (length(parsed) == 1L && !is.na(parsed[1])) {
          return(c(effective = parsed[1], raw = parsed[1]))
        }
        if (length(parsed) >= 2L && !anyNA(parsed[1:2])) {
          return(c(effective = parsed[1], raw = parsed[2]))
        }
        stop(
          paste(
            "Invalid preliminary trial output for k =", target_clusters,
            "phase =", phase_label,
            "iteration =", iteration_idx,
            "trial =", trial_idx_for_error
          )
        )
      }
      heartbeat(function() {
        paste(
          phase_label, "step", iteration_idx,
          "- representative preliminary evaluation started for gamma", signif(gamma_val, 6)
        )
      })
      trial_result <- run_single_trial_count(
        trial_idx = 0L,
        gamma_val = gamma_val,
        cache_suffix_base = cache_suffix_base,
        seed_base = seed_base,
        use_cache = TRUE
      )
      trial_pair <- unpack_trial_counts(trial_result, trial_idx_for_error = 0L)
      trial_counts[] <- as.integer(trial_pair[["effective"]])
      raw_trial_counts[] <- as.integer(trial_pair[["raw"]])
      
      if (length(trial_counts) == 0L ||
          anyNA(trial_counts) ||
          anyNA(raw_trial_counts)) {
        stop(
          paste(
            "No preliminary trial results for k =", target_clusters,
            "phase =", phase_label,
            "iteration =", iteration_idx
          )
        )
      }
      
      n_clusters_obtained <- as.numeric(stats::median(trial_counts))
      raw_clusters_obtained <- as.numeric(stats::median(raw_trial_counts))
      elapsed_seconds <- as.numeric(difftime(Sys.time(), step_start_time, units = "secs"))
      
      if (verbose) {
        search_state <- classify_resolution_search_state(
          raw_cluster_median = raw_clusters_obtained,
          effective_cluster_median = n_clusters_obtained,
          target_clusters = target_clusters,
          min_cluster_size = min_cluster_size
        )
        scice_message(
          paste(
            "RESOLUTION_SEARCH: k =", target_clusters,
            phase_label, "step", iteration_idx,
            "- representative preliminary evaluation completed",
            "(", n_preliminary_trials, "/", n_preliminary_trials, ")",
            "- median effective clusters =", signif(n_clusters_obtained, 6),
            "- median raw clusters =", signif(raw_clusters_obtained, 6),
            "- raw class =", search_state$raw_class,
            "- over_fragmented =", search_state$over_fragmented,
            "- elapsed", round(elapsed_seconds, 3), "s"
          )
        )
      }
      
      list(
        n_clusters_obtained = n_clusters_obtained,
        raw_clusters_obtained = raw_clusters_obtained,
        completed_trials = n_preliminary_trials,
        elapsed_seconds = elapsed_seconds
      )
    }

    left <- start_g
    right <- end_g
    max_iterations <- max_search_iterations
    iteration_count <- 0
    
    if (verbose) {
      scice_message(paste("RESOLUTION_SEARCH: k =", target_clusters, "- starting lower bound search"))
    }
    
    # Binary search for lower bound with fallback
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > effective_tolerance && iteration_count < max_iterations) {
      iteration_count <- iteration_count + 1
      mid <- (left + right) / 2
      gamma_val <- if (objective_function == "modularity") mid else exp(mid)
      heartbeat(function() {
        paste(
          "lower-bound loop iteration", iteration_count, "/", max_iterations,
          "- gamma", signif(gamma_val, 6),
          "- current interval [", signif(left, 6), ",", signif(right, 6), "]"
        )
      })
      
      # Test clustering with current gamma using vectorized operations
      # Set deterministic seed for lower bound search if base seed provided
      range_seed <- NULL
      if (!is.null(seed)) {
        range_seed <- as.integer((as.double(seed) + target_clusters * 10 + iteration_count) %% .Machine$integer.max)
        if (is.na(range_seed) || range_seed <= 0L) {
          range_seed <- 1L
        }
      }

      cluster_results <- run_preliminary_trials(
        gamma_val,
        cache_suffix_base = paste("res_search_lower", target_clusters, sep = "_"),
        seed_base = range_seed,
        phase_label = "lower",
        iteration_idx = iteration_count
      )
      
      n_clusters_obtained <- cluster_results$n_clusters_obtained
      raw_clusters_obtained <- cluster_results$raw_clusters_obtained
      search_state <- classify_resolution_search_state(
        raw_cluster_median = raw_clusters_obtained,
        effective_cluster_median = n_clusters_obtained,
        target_clusters = target_clusters,
        min_cluster_size = min_cluster_size
      )
      trial_steps <- trial_steps + 1L
      total_preliminary_elapsed <- total_preliminary_elapsed + cluster_results$elapsed_seconds
      if (verbose && (iteration_count == 1 || iteration_count %% 5 == 0)) {
        interval_span <- if (objective_function == "modularity") {
          abs(right - left)
        } else {
          abs(exp(right) - exp(left))
        }
        scice_message(
          paste(
            "RESOLUTION_SEARCH: k =", target_clusters,
            "lower-bound progress", iteration_count, "/", max_iterations,
            "- gamma =", signif(gamma_val, 6),
            "- median effective clusters =", n_clusters_obtained,
            "- median raw clusters =", raw_clusters_obtained,
            "- raw class =", search_state$raw_class,
            "- over_fragmented =", search_state$over_fragmented,
            "- interval span =", signif(interval_span, 6)
          )
        )
      }
      
      if (identical(search_state$lower_action, "increase_gamma")) {
        left <- mid
      } else {
        right <- mid
      }
    }
    
    left_bound <- right
    
    # Binary search for upper bound with fallback
    left <- left_bound
    right <- end_g
    iteration_count <- 0
    
    if (verbose) {
      scice_message(paste("RESOLUTION_SEARCH: k =", target_clusters, "- starting upper bound search"))
    }
    
    while (abs(if (objective_function == "modularity") left - right else exp(left) - exp(right)) > effective_tolerance && iteration_count < max_iterations) {
      iteration_count <- iteration_count + 1
      mid <- (left + right) / 2
      gamma_val <- if (objective_function == "modularity") mid else exp(mid)
      heartbeat(function() {
        paste(
          "upper-bound loop iteration", iteration_count, "/", max_iterations,
          "- gamma", signif(gamma_val, 6),
          "- current interval [", signif(left, 6), ",", signif(right, 6), "]"
        )
      })
      
      # Test clustering with current gamma using vectorized operations
      # Set deterministic seed for upper bound search if base seed provided
      range_seed <- NULL
      if (!is.null(seed)) {
        range_seed <- as.integer((as.double(seed) + target_clusters * 10 + 100 + iteration_count) %% .Machine$integer.max)
        if (is.na(range_seed) || range_seed <= 0L) {
          range_seed <- 1L
        }
      }

      cluster_results <- run_preliminary_trials(
        gamma_val,
        cache_suffix_base = paste("res_search_upper", target_clusters, sep = "_"),
        seed_base = range_seed,
        phase_label = "upper",
        iteration_idx = iteration_count
      )
      
      n_clusters_obtained <- cluster_results$n_clusters_obtained
      raw_clusters_obtained <- cluster_results$raw_clusters_obtained
      search_state <- classify_resolution_search_state(
        raw_cluster_median = raw_clusters_obtained,
        effective_cluster_median = n_clusters_obtained,
        target_clusters = target_clusters,
        min_cluster_size = min_cluster_size
      )
      trial_steps <- trial_steps + 1L
      total_preliminary_elapsed <- total_preliminary_elapsed + cluster_results$elapsed_seconds
      if (verbose && (iteration_count == 1 || iteration_count %% 5 == 0)) {
        interval_span <- if (objective_function == "modularity") {
          abs(right - left)
        } else {
          abs(exp(right) - exp(left))
        }
        scice_message(
          paste(
            "RESOLUTION_SEARCH: k =", target_clusters,
            "upper-bound progress", iteration_count, "/", max_iterations,
            "- gamma =", signif(gamma_val, 6),
            "- median effective clusters =", n_clusters_obtained,
            "- median raw clusters =", raw_clusters_obtained,
            "- raw class =", search_state$raw_class,
            "- over_fragmented =", search_state$over_fragmented,
            "- interval span =", signif(interval_span, 6)
          )
        )
      }
      
      if (identical(search_state$upper_action, "increase_gamma")) {
        left <- mid
      } else {
        right <- mid
      }
    }
    
    right_bound <- left
    
    # If bounds are identical or invalid, create a cluster-specific range
    if (identical(left_bound, right_bound) || is.na(left_bound) || is.na(right_bound)) {
      if (objective_function == "CPM") {
        # CPM searches in log-gamma space, so convert through gamma scale before
        # widening an exact/degenerate interval.
        center_gamma <- if (is.na(left_bound)) {
          exp((start_g + end_g) / 2)
        } else {
          exp(left_bound)
        }
        # When an exact CPM bound is already identified, keep the widened fallback
        # centered on that bound instead of shifting by cluster number.
        left_gamma <- max(exp(start_g), center_gamma * 0.7)
        right_gamma <- min(exp(end_g), center_gamma * 1.3)
        left_bound <- log(left_gamma)
        right_bound <- log(right_gamma)
      } else {
        # For modularity, create a range around the found value with cluster-specific adjustment
        center_val <- ifelse(is.na(left_bound), (start_g + end_g) / 2, left_bound)
        # Add cluster-specific offset to ensure different ranges for different cluster numbers
        cluster_offset <- (target_clusters - min(cluster_range)) * 0.02  # Small offset based on cluster number
        adjusted_center <- center_val + cluster_offset
        left_bound <- max(start_g, adjusted_center - 0.15)
        right_bound <- min(end_g, adjusted_center + 0.15)
      }
    }

    coarse_gamma_range <- sort(
      if (objective_function == "CPM") {
        exp(c(left_bound, right_bound))
      } else {
        c(left_bound, right_bound)
      }
    )
    final_gamma_range <- coarse_gamma_range
    plateau_clamp_mode <- "coarse"
    plateau_probe_indices <- integer(0)

    if (min_cluster_size > 1L && all(is.finite(coarse_gamma_range))) {
      probe_gamma_sequence <- build_gamma_sequence_for_range(
        coarse_gamma_range,
        objective_function = objective_function,
        resolution_tolerance = resolution_tolerance,
        n_vertices = igraph::vcount(igraph_obj)
      )
      probe_raw_medians <- numeric(length(probe_gamma_sequence))

      for (probe_idx in seq_along(probe_gamma_sequence)) {
        plateau_seed <- NULL
        if (!is.null(seed)) {
          plateau_seed <- as.integer(
            (as.double(seed) + target_clusters * 10 + 1000 + probe_idx) %% .Machine$integer.max
          )
          if (is.na(plateau_seed) || plateau_seed <= 0L) {
            plateau_seed <- 1L
          }
        }
        probe_result <- run_preliminary_trials(
          probe_gamma_sequence[[probe_idx]],
          cache_suffix_base = paste("res_search_clamp", target_clusters, sep = "_"),
          seed_base = plateau_seed,
          phase_label = "clamp",
          iteration_idx = probe_idx
        )
        trial_steps <- trial_steps + 1L
        total_preliminary_elapsed <- total_preliminary_elapsed + probe_result$elapsed_seconds
        probe_raw_medians[[probe_idx]] <- probe_result$raw_clusters_obtained
      }

      clamp_result <- clamp_gamma_range_to_raw_plateau(
        probe_gamma_sequence,
        probe_raw_medians,
        target_clusters = target_clusters
      )
      final_gamma_range <- sort(as.numeric(clamp_result$bounds))
      plateau_clamp_mode <- clamp_result$mode
      plateau_probe_indices <- as.integer(clamp_result$indices)
    }

    if (verbose) {
      scice_message(
        paste0(
          "RESOLUTION_SEARCH: k = ", target_clusters,
          " - coarse bounds [", signif(coarse_gamma_range[1], 6), ", ", signif(coarse_gamma_range[2], 6), "]"
        )
      )
      if (min_cluster_size > 1L) {
        scice_message(
          paste0(
            "RESOLUTION_SEARCH: k = ", target_clusters,
            " - plateau clamp mode = ", plateau_clamp_mode,
            " -> final bounds [", signif(final_gamma_range[1], 6), ", ", signif(final_gamma_range[2], 6), "]"
          )
        )
        if (!is.null(clamp_result$stabilized_raw_medians) &&
            !isTRUE(all.equal(
              as.numeric(probe_raw_medians),
              as.numeric(clamp_result$stabilized_raw_medians),
              tolerance = 0
            ))) {
          scice_message(
            paste(
              "RESOLUTION_SEARCH: k =", target_clusters,
              "- plateau clamp stabilized raw medians:",
              paste(signif(as.numeric(clamp_result$stabilized_raw_medians), 6), collapse = ", ")
            )
          )
        }
        if (length(plateau_probe_indices) > 0L && plateau_clamp_mode != "coarse") {
          scice_message(
            paste(
              "RESOLUTION_SEARCH: k =", target_clusters,
              "- plateau clamp probe indices:",
              paste(plateau_probe_indices, collapse = ", ")
            )
          )
        }
      } else {
        scice_message(
          paste0(
            "RESOLUTION_SEARCH: k = ", target_clusters,
            " - bounds identified [", signif(final_gamma_range[1], 6), ", ", signif(final_gamma_range[2], 6), "]"
          )
        )
      }
      scice_message(
        paste(
          "RESOLUTION_SEARCH: k =", target_clusters,
          "trial summary - full preliminary trial sets:", trial_steps,
          "- trials per step:", n_preliminary_trials,
          "- cumulative preliminary elapsed:", round(total_preliminary_elapsed, 3), "s"
        )
      )
    }
    
    # Return results as a data.table row
    data.table::data.table(
      cluster_number = target_clusters,
      left_bound = final_gamma_range[1],
      right_bound = final_gamma_range[2]
    )
  }, mc.cores = active_cluster_workers)
  
  # Combine results
  # Filter out NULL results from parallel processing
  valid_results <- range_results[!sapply(range_results, is.null)]
  
  if (length(valid_results) == 0) {
    # If no valid results, return empty gamma_dict
    warning("No valid resolution ranges found for any cluster numbers")
    return(setNames(list(), character(0)))
  }
  
  results_dt <- data.table::rbindlist(valid_results)
  
  # Check if required columns exist
  if (!"left_bound" %in% colnames(results_dt) || !"right_bound" %in% colnames(results_dt)) {
    stop("Missing required columns in resolution range results")
  }
  
  # Sort by cluster number to ensure proper ordering
  data.table::setorder(results_dt, cluster_number)
  
  # Post-process to ensure non-overlapping ranges for adjacent cluster numbers
  if (nrow(results_dt) > 1) {
    for (i in 2:nrow(results_dt)) {
      current_left <- results_dt$left_bound[i]
      current_right <- results_dt$right_bound[i]
      prev_left <- results_dt$left_bound[i-1]
      prev_right <- results_dt$right_bound[i-1]
      
      # Check for overlap or identical ranges
      if (current_left <= prev_right && abs(current_left - prev_right) < resolution_tolerance * 10) {
        if (verbose) {
          scice_message(paste("RESOLUTION_SEARCH: Adjusting overlapping ranges for clusters", 
                       results_dt$cluster_number[i-1], "and", results_dt$cluster_number[i]))
        }
        
        # Calculate a small gap between ranges
        if (objective_function == "CPM") {
          gap <- max(resolution_tolerance * 5, (prev_right - prev_left) * 0.1)
          # Adjust current range to be above previous range
          range_width <- current_right - current_left
          results_dt$left_bound[i] <- prev_right + gap
          results_dt$right_bound[i] <- results_dt$left_bound[i] + range_width
        } else {
          gap <- max(resolution_tolerance * 10, (prev_right - prev_left) * 0.1)
          # Adjust current range to be above previous range
          range_width <- current_right - current_left
          results_dt$left_bound[i] <- prev_right + gap
          results_dt$right_bound[i] <- results_dt$left_bound[i] + range_width
        }
      }
    }
  }
  
  # Convert to named list format - fixing the data.table syntax
  gamma_dict <- lapply(seq_len(nrow(results_dt)), function(i) {
    c(results_dt$left_bound[i], results_dt$right_bound[i])
  })
  names(gamma_dict) <- as.character(results_dt$cluster_number)
  
  return(gamma_dict)
}
