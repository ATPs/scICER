# Shared global gamma sweep utilities for cluster-range mode.

empty_resolution_search_diagnostics_df <- function() {
  data.frame(
    sweep_round = integer(),
    gamma = numeric(),
    effective_cluster_count = numeric(),
    raw_cluster_count = numeric(),
    final_cluster_count = numeric(),
    raw_class = character(),
    over_fragmented = logical(),
    selected_for_refinement = logical(),
    selected_for_target_interval = logical(),
    plateau_round = integer(),
    stringsAsFactors = FALSE
  )
}

stabilize_monotone_probe_counts <- function(values) {
  values <- as.numeric(values)
  finite_indices <- which(is.finite(values))
  if (length(finite_indices) <= 1L) {
    return(values)
  }
  values[finite_indices] <- cummax(values[finite_indices])
  values
}

global_resolution_search_midpoint <- function(left, right, objective_function) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  if (!is.finite(left) || !is.finite(right)) {
    return(NA_real_)
  }
  if (identical(objective_function, "CPM")) {
    exp((log(left) + log(right)) / 2)
  } else {
    (left + right) / 2
  }
}

global_resolution_search_interval_small <- function(left, right, objective_function,
                                                    resolution_tolerance) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  if (!is.finite(left) || !is.finite(right) || left >= right) {
    return(TRUE)
  }
  tolerance <- max(resolution_tolerance, .Machine$double.eps * 100)
  if (identical(objective_function, "CPM")) {
    abs(log(right) - log(left)) <= tolerance
  } else {
    abs(right - left) <= tolerance
  }
}

global_resolution_search_probe_batch <- function(igraph_obj, gamma_values, sweep_round,
                                                 objective_function, n_iter_preliminary,
                                                 beta_preliminary, requested_max,
                                                 min_cluster_size, snn_graph,
                                                 active_probe_workers, verbose, seed) {
  gamma_values <- sort(unique(as.numeric(gamma_values[is.finite(gamma_values)])))
  if (length(gamma_values) == 0L) {
    return(data.table::data.table())
  }

  if (verbose) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Sweep round", sweep_round,
        "- evaluating", length(gamma_values), "shared gamma probes with",
        min(length(gamma_values), active_probe_workers), "worker(s)"
      )
    )
  }

  probe_results <- cross_platform_mclapply(
    seq_along(gamma_values),
    function(idx) {
      gamma_val <- gamma_values[[idx]]
      probe_seed <- NULL
      if (!is.null(seed)) {
        gamma_component <- if (is.finite(gamma_val)) {
          as.integer(floor(abs(gamma_val * 1000000) %% 100000000))
        } else {
          0L
        }
        probe_seed <- as.integer(
          (as.double(seed) + as.double(gamma_component) + as.double(sweep_round) + as.double(idx)) %%
            .Machine$integer.max
        )
        if (is.na(probe_seed) || probe_seed <= 0L) {
          probe_seed <- 1L
        }
        set.seed(probe_seed)
      }

      labels <- cached_leiden_clustering(
        igraph_obj,
        gamma_val,
        objective_function,
        n_iter_preliminary,
        beta_preliminary,
        use_cache = TRUE,
        cache_key_suffix = "global_resolution_search",
        snn_graph = snn_graph,
        min_cluster_size = min_cluster_size
      )
      effective_count <- count_effective_clusters(labels, min_cluster_size = min_cluster_size)
      raw_count <- as.integer(length(unique(labels)))
      final_count <- if (min_cluster_size > 1L) {
        merged_labels <- merge_small_clusters_to_neighbors(
          labels = labels,
          snn_graph = snn_graph,
          min_cluster_size = min_cluster_size
        )
        as.integer(length(unique(merged_labels)))
      } else {
        raw_count
      }
      state <- classify_resolution_search_state(
        raw_cluster_median = raw_count,
        effective_cluster_median = effective_count,
        target_clusters = requested_max,
        min_cluster_size = min_cluster_size
      )
      data.table::data.table(
        sweep_round = as.integer(sweep_round),
        gamma = as.numeric(gamma_val),
        effective_cluster_count = as.numeric(effective_count),
        raw_cluster_count = as.numeric(raw_count),
        final_cluster_count = as.numeric(final_count),
        raw_class = as.character(state$raw_class),
        over_fragmented = isTRUE(state$over_fragmented),
        selected_for_refinement = sweep_round > 1L,
        selected_for_target_interval = FALSE,
        plateau_round = 0L
      )
    },
    mc.cores = min(length(gamma_values), active_probe_workers),
    mc.preschedule = FALSE
  )

  data.table::rbindlist(probe_results, fill = TRUE)
}

derive_shared_gamma_intervals <- function(probes_dt, cluster_range, gamma_bounds) {
  if (is.null(probes_dt) || nrow(probes_dt) == 0L) {
    target_gamma_seeds <- setNames(vector("list", length(cluster_range)), as.character(cluster_range))
    target_interval_details <- setNames(vector("list", length(cluster_range)), as.character(cluster_range))
    for (target_clusters in as.integer(cluster_range)) {
      target_interval_details[[as.character(target_clusters)]] <- list(
        target_cluster = as.integer(target_clusters),
        bracketed = FALSE,
        optimization_ready = FALSE,
        has_exact_probe = FALSE,
        gamma_left = NA_real_,
        gamma_right = NA_real_,
        seed_gamma_values = numeric(0),
        exact_probe_values = numeric(0),
        near_probe_values = numeric(0)
      )
    }
    return(list(
      gamma_dict = setNames(list(), character(0)),
      bracketed_targets = integer(0),
      optimization_ready_targets = integer(0),
      unresolved_targets = as.integer(cluster_range),
      selected_gamma_values = numeric(0),
      unresolved_intervals = list(),
      target_gamma_seeds = target_gamma_seeds,
      target_interval_details = target_interval_details
    ))
  }

  probes_dt <- data.table::copy(probes_dt)
  data.table::setorder(probes_dt, gamma)
  probes_dt[, stabilized_final_cluster_count := stabilize_monotone_probe_counts(final_cluster_count)]

  gamma_dict <- setNames(vector("list", length(cluster_range)), as.character(cluster_range))
  bracketed_targets <- integer(0)
  optimization_ready_targets <- integer(0)
  unresolved_targets <- integer(0)
  selected_gamma_values <- numeric(0)
  unresolved_intervals <- list()
  target_gamma_seeds <- setNames(vector("list", length(cluster_range)), as.character(cluster_range))
  target_interval_details <- setNames(vector("list", length(cluster_range)), as.character(cluster_range))

  gamma_values <- probes_dt$gamma
  stabilized_final <- probes_dt$stabilized_final_cluster_count

  for (target_clusters in as.integer(cluster_range)) {
    exact_indices <- which(stabilized_final == target_clusters)
    near_indices <- which(abs(stabilized_final - target_clusters) <= 1)
    bounds <- NULL
    interval_indices <- integer(0)
    bracketed <- FALSE
    optimization_ready <- FALSE

    if (length(exact_indices) > 0L) {
      interval_indices <- exact_indices
      bounds <- range(gamma_values[exact_indices])
      if (identical(bounds[1], bounds[2])) {
        left_neighbor <- suppressWarnings(max(which(gamma_values < bounds[1])))
        right_neighbor <- suppressWarnings(min(which(gamma_values > bounds[2])))
        neighbor_indices <- integer(0)
        if (length(left_neighbor) == 1L && is.finite(left_neighbor)) {
          neighbor_indices <- c(neighbor_indices, left_neighbor)
        }
        if (length(right_neighbor) == 1L && is.finite(right_neighbor)) {
          neighbor_indices <- c(neighbor_indices, right_neighbor)
        }
        if (length(neighbor_indices) > 0L) {
          interval_indices <- sort(unique(c(interval_indices, neighbor_indices)))
          bounds <- range(gamma_values[interval_indices])
        }
      }
      bracketed <- length(bounds) == 2L && all(is.finite(bounds))
      optimization_ready <- isTRUE(bracketed)
    } else {
      below_indices <- which(stabilized_final < target_clusters)
      above_indices <- which(stabilized_final > target_clusters)
      if (length(below_indices) > 0L && length(above_indices) > 0L) {
        left_idx <- max(below_indices)
        right_idx <- min(above_indices)
        if (left_idx < right_idx) {
          interval_indices <- c(left_idx, right_idx)
          bounds <- sort(gamma_values[interval_indices])
          bracketed <- TRUE
        }
      }
    }

    seed_indices <- sort(unique(c(interval_indices, exact_indices, near_indices)))
    target_gamma_seeds[[as.character(target_clusters)]] <- sort(unique(gamma_values[seed_indices]))

    if (isTRUE(bracketed)) {
      gamma_dict[[as.character(target_clusters)]] <- as.numeric(bounds)
      bracketed_targets <- c(bracketed_targets, target_clusters)
      if (isTRUE(optimization_ready)) {
        optimization_ready_targets <- c(optimization_ready_targets, target_clusters)
      }
      selected_gamma_values <- c(selected_gamma_values, gamma_values[seed_indices])
    } else {
      unresolved_targets <- c(unresolved_targets, target_clusters)
      below_indices <- which(stabilized_final < target_clusters)
      above_indices <- which(stabilized_final > target_clusters)
      interval <- NULL
      if (length(below_indices) > 0L && length(above_indices) > 0L) {
        left_idx <- max(below_indices)
        right_idx <- min(above_indices)
        if (left_idx < right_idx) {
          interval <- sort(gamma_values[c(left_idx, right_idx)])
        }
      } else if (length(below_indices) > 0L) {
        interval <- c(max(gamma_values[below_indices]), gamma_bounds[2])
      } else if (length(above_indices) > 0L) {
        interval <- c(gamma_bounds[1], min(gamma_values[above_indices]))
      } else {
        interval <- gamma_bounds
      }
      unresolved_intervals[[as.character(target_clusters)]] <- as.numeric(sort(interval))
    }

    if (isTRUE(bracketed) && !isTRUE(optimization_ready)) {
      unresolved_targets <- c(unresolved_targets, target_clusters)
      unresolved_intervals[[as.character(target_clusters)]] <- as.numeric(sort(bounds))
    }

    target_interval_details[[as.character(target_clusters)]] <- list(
      target_cluster = as.integer(target_clusters),
      bracketed = isTRUE(bracketed),
      optimization_ready = isTRUE(optimization_ready),
      has_exact_probe = length(exact_indices) > 0L,
      gamma_left = if (isTRUE(bracketed)) as.numeric(bounds[[1]]) else NA_real_,
      gamma_right = if (isTRUE(bracketed)) as.numeric(bounds[[2]]) else NA_real_,
      seed_gamma_values = sort(unique(gamma_values[seed_indices])),
      exact_probe_values = sort(unique(gamma_values[exact_indices])),
      near_probe_values = sort(unique(gamma_values[near_indices]))
    )
  }

  gamma_dict <- gamma_dict[vapply(gamma_dict, function(x) length(x) == 2L, logical(1))]
  list(
    gamma_dict = gamma_dict,
    bracketed_targets = sort(unique(as.integer(bracketed_targets))),
    optimization_ready_targets = sort(unique(as.integer(optimization_ready_targets))),
    unresolved_targets = sort(unique(as.integer(unresolved_targets))),
    selected_gamma_values = sort(unique(as.numeric(selected_gamma_values))),
    unresolved_intervals = unresolved_intervals,
    target_gamma_seeds = target_gamma_seeds,
    target_interval_details = target_interval_details
  )
}

# Override the legacy per-target binary search with a shared global gamma sweep.
find_resolution_ranges <- function(igraph_obj, cluster_range, start_g, end_g,
                                  objective_function, resolution_tolerance, n_workers, verbose, seed = NULL,
                                  snn_graph = NULL, min_cluster_size = 1L,
                                  in_parallel_context = FALSE, runtime_context = NULL) {
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  if (min_cluster_size > 1L && is.null(snn_graph)) {
    stop("snn_graph must be provided when min_cluster_size > 1.")
  }

  cluster_range <- sort(unique(as.integer(cluster_range)))
  if (length(cluster_range) == 0L) {
    gamma_dict <- setNames(list(), character(0))
    attr(gamma_dict, "resolution_search_diagnostics") <- empty_resolution_search_diagnostics_df()
    attr(gamma_dict, "coverage_complete") <- TRUE
    attr(gamma_dict, "plateau_stop") <- FALSE
    attr(gamma_dict, "uncovered_targets") <- integer(0)
    attr(gamma_dict, "target_gamma_seeds") <- setNames(list(), character(0))
    attr(gamma_dict, "target_interval_details") <- setNames(list(), character(0))
    return(gamma_dict)
  }

  n_vertices <- igraph::vcount(igraph_obj)
  n_preliminary_trials <- if (n_vertices >= 200000) {
    3L
  } else if (n_vertices >= 100000) {
    5L
  } else {
    15L
  }
  n_iter_preliminary <- if (n_vertices >= 200000) 3L else 5L
  beta_preliminary <- 0.01
  max_search_iterations <- if (n_vertices >= 200000) 30L else 50L
  requested_search_workers <- max(1L, as.integer(n_workers))
  if (.Platform$OS.type == "windows" && requested_search_workers > 1L) {
    requested_search_workers <- 1L
  }
  search_worker_budget <- if (in_parallel_context) 1L else requested_search_workers
  search_worker_budget <- cap_workers_by_memory(
    search_worker_budget,
    estimate_trial_matrix_bytes(n_vertices, n_preliminary_trials, 1L),
    runtime_context
  )
  active_probe_workers <- max(1L, as.integer(search_worker_budget))
  gamma_bounds <- if (identical(objective_function, "CPM")) {
    c(exp(start_g), exp(end_g))
  } else {
    c(start_g, end_g)
  }
  requested_max <- max(cluster_range)

  if (verbose) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Using", n_preliminary_trials,
        "preliminary trials per step (graph vertices:", n_vertices, ")"
      )
    )
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Worker allocation - requested:", as.integer(n_workers),
        "| budget after context/memory:", search_worker_budget,
        "| active probe workers:", active_probe_workers
      )
    )
    scice_message("RESOLUTION_SEARCH: Preliminary trial workers per gamma: 1")
    scice_message("RESOLUTION_SEARCH: Preliminary trial mode: serial (shared one-clustering-per-gamma sweep)")
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Preliminary trial strategy - reuse one representative preliminary clustering per shared gamma probe and derive target intervals from final merged cluster counts."
      )
    )
  }

  all_probe_results <- global_resolution_search_probe_batch(
    igraph_obj = igraph_obj,
    gamma_values = build_gamma_sequence_for_range(
      gamma_range = gamma_bounds,
      objective_function = objective_function,
      resolution_tolerance = resolution_tolerance,
      n_vertices = n_vertices
    ),
    sweep_round = 1L,
    objective_function = objective_function,
    n_iter_preliminary = n_iter_preliminary,
    beta_preliminary = beta_preliminary,
    requested_max = requested_max,
    min_cluster_size = min_cluster_size,
    snn_graph = snn_graph,
    active_probe_workers = active_probe_workers,
    verbose = verbose,
    seed = seed
  )
  data.table::setorder(all_probe_results, gamma)

  interval_state <- derive_shared_gamma_intervals(all_probe_results, cluster_range, gamma_bounds)
  previous_max_final <- if (nrow(all_probe_results) > 0L) {
    max(stabilize_monotone_probe_counts(all_probe_results$final_cluster_count), na.rm = TRUE)
  } else {
    NA_real_
  }
  previous_ready_count <- length(interval_state$optimization_ready_targets)
  plateau_count <- 0L
  plateau_stop <- FALSE

  if (verbose) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Initial shared sweep covered final merged cluster counts up to",
        ifelse(is.finite(previous_max_final), format(signif(previous_max_final, 6)), "NA"),
        "| optimization-ready targets =",
        previous_ready_count, "/", length(cluster_range)
      )
    )
  }

  for (sweep_round in seq_len(max_search_iterations) + 1L) {
    if (length(interval_state$unresolved_targets) == 0L) {
      break
    }

    next_probe_values <- numeric(0)
    for (target_clusters in interval_state$unresolved_targets) {
      interval <- interval_state$unresolved_intervals[[as.character(target_clusters)]]
      if (is.null(interval) || length(interval) != 2L || any(!is.finite(interval))) {
        next
      }
      interval <- sort(as.numeric(interval))
      if (interval[1] >= interval[2] ||
          global_resolution_search_interval_small(
            interval[1], interval[2], objective_function, resolution_tolerance
          )) {
        next
      }
      next_probe_values <- c(
        next_probe_values,
        global_resolution_search_midpoint(interval[1], interval[2], objective_function)
      )
    }
    next_probe_values <- sort(unique(next_probe_values))
    next_probe_values <- next_probe_values[!next_probe_values %in% all_probe_results$gamma]

    if (length(next_probe_values) == 0L) {
      break
    }

    new_probe_results <- global_resolution_search_probe_batch(
      igraph_obj = igraph_obj,
      gamma_values = next_probe_values,
      sweep_round = sweep_round,
      objective_function = objective_function,
      n_iter_preliminary = n_iter_preliminary,
      beta_preliminary = beta_preliminary,
      requested_max = requested_max,
      min_cluster_size = min_cluster_size,
      snn_graph = snn_graph,
      active_probe_workers = active_probe_workers,
      verbose = verbose,
      seed = seed
    )
    all_probe_results <- data.table::rbindlist(list(all_probe_results, new_probe_results), fill = TRUE)
    all_probe_results <- unique(all_probe_results, by = "gamma")
    data.table::setorder(all_probe_results, gamma)

    interval_state <- derive_shared_gamma_intervals(all_probe_results, cluster_range, gamma_bounds)
    current_max_final <- if (nrow(all_probe_results) > 0L) {
      max(stabilize_monotone_probe_counts(all_probe_results$final_cluster_count), na.rm = TRUE)
    } else {
      NA_real_
    }
    current_ready_count <- length(interval_state$optimization_ready_targets)
    no_growth <- (!is.finite(previous_max_final) && !is.finite(current_max_final)) ||
      (is.finite(previous_max_final) && is.finite(current_max_final) && current_max_final <= previous_max_final)
    no_new_ready <- current_ready_count <= previous_ready_count

    if (no_growth && no_new_ready) {
      plateau_count <- plateau_count + 1L
      current_round <- as.integer(sweep_round)
      all_probe_results[sweep_round == current_round, plateau_round := plateau_count]
    } else {
      plateau_count <- 0L
    }

    if (verbose) {
      scice_message(
        paste(
          "RESOLUTION_SEARCH: Sweep round", sweep_round,
          "- probed", length(next_probe_values), "new gamma values",
          "- current final max =", ifelse(is.finite(current_max_final), format(signif(current_max_final, 6)), "NA"),
          "- optimization-ready targets =", current_ready_count, "/", length(cluster_range),
          "- plateau count =", plateau_count
        )
      )
    }

    previous_max_final <- current_max_final
    previous_ready_count <- current_ready_count
    if (plateau_count >= 2L) {
      plateau_stop <- TRUE
      break
    }
  }

  all_probe_results[, selected_for_target_interval := gamma %in% interval_state$selected_gamma_values]
  stabilized_final <- stabilize_monotone_probe_counts(all_probe_results$final_cluster_count)
  coverage_complete <- length(interval_state$unresolved_targets) == 0L &&
    any(is.finite(stabilized_final)) &&
    max(stabilized_final, na.rm = TRUE) >= requested_max
  uncovered_targets <- setdiff(cluster_range, interval_state$optimization_ready_targets)

  if (verbose) {
    if (length(interval_state$gamma_dict) > 0L) {
      for (cluster_num in names(interval_state$gamma_dict)) {
        bounds <- interval_state$gamma_dict[[cluster_num]]
        scice_message(
          paste0(
            "RESOLUTION_SEARCH: k = ", cluster_num,
            " - shared-sweep bounds [", signif(bounds[1], 6), ", ", signif(bounds[2], 6), "]"
          )
        )
      }
    }
    if (length(uncovered_targets) > 0L) {
      scice_message(
        paste("RESOLUTION_SEARCH: Uncovered requested final targets:", paste(uncovered_targets, collapse = ", "))
      )
    }
    if (plateau_stop) {
      scice_message("RESOLUTION_SEARCH: Stopped after two consecutive plateau rounds without new coverage.")
    }
  }

  diagnostics_df <- as.data.frame(all_probe_results)
  rownames(diagnostics_df) <- NULL
  attr(interval_state$gamma_dict, "resolution_search_diagnostics") <- diagnostics_df
  attr(interval_state$gamma_dict, "coverage_complete") <- isTRUE(coverage_complete)
  attr(interval_state$gamma_dict, "plateau_stop") <- isTRUE(plateau_stop)
  attr(interval_state$gamma_dict, "uncovered_targets") <- as.integer(uncovered_targets)
  attr(interval_state$gamma_dict, "target_gamma_seeds") <- interval_state$target_gamma_seeds
  attr(interval_state$gamma_dict, "target_interval_details") <- interval_state$target_interval_details
  interval_state$gamma_dict
}
