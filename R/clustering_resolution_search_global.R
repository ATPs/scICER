# Shared global gamma sweep utilities for cluster-range mode.

empty_resolution_search_diagnostics_df <- function() {
  data.frame(
    sweep_round = integer(),
    discovery_round = integer(),
    probe_stage = character(),
    probe_index = integer(),
    probe_pid = integer(),
    probe_elapsed_sec = numeric(),
    gamma = numeric(),
    upper_cap_discovery_gamma = numeric(),
    degenerate_high_gamma = logical(),
    scheduled_probe_workers = integer(),
    coarse_probe_count = integer(),
    discovered_upper_gamma = numeric(),
    upper_cap_stop_reason = character(),
    refinement_interval_width = numeric(),
    refinement_interval_id = integer(),
    refinement_points_per_interval = integer(),
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

is_high_gamma_degenerate_probe <- function(effective_cluster_count, raw_cluster_count,
                                           final_cluster_count, n_vertices) {
  isTRUE(
    is.finite(effective_cluster_count) &&
      is.finite(raw_cluster_count) &&
      is.finite(final_cluster_count) &&
      effective_cluster_count == 0 &&
      final_cluster_count == 1 &&
      raw_cluster_count >= (0.98 * as.numeric(n_vertices))
  )
}

build_cpm_discovery_batch_gamma_values <- function(current_gamma, hard_cap_gamma,
                                                   batch_size, step_ratio = 4) {
  current_gamma <- max(as.numeric(current_gamma), .Machine$double.xmin)
  hard_cap_gamma <- max(as.numeric(hard_cap_gamma), current_gamma)
  batch_size <- max(0L, as.integer(batch_size))
  step_ratio <- max(1 + .Machine$double.eps, as.numeric(step_ratio))
  if (batch_size <= 0L) {
    return(numeric(0))
  }
  sort(unique(pmin(
    current_gamma * (step_ratio ^ seq.int(0L, batch_size - 1L)),
    hard_cap_gamma
  )))
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

global_resolution_search_interval_width <- function(left, right, objective_function) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  if (!is.finite(left) || !is.finite(right) || left >= right) {
    return(NA_real_)
  }
  if (identical(objective_function, "CPM")) {
    return(log(right) - log(left))
  }
  right - left
}

global_resolution_search_internal_points <- function(left, right, objective_function, n_points) {
  left <- as.numeric(left)
  right <- as.numeric(right)
  n_points <- max(0L, as.integer(n_points))
  if (n_points <= 0L || !is.finite(left) || !is.finite(right) || left >= right) {
    return(numeric(0))
  }
  fractions <- seq_len(n_points) / (n_points + 1L)
  if (identical(objective_function, "CPM")) {
    return(exp(log(left) + fractions * (log(right) - log(left))))
  }
  left + fractions * (right - left)
}

build_refinement_probe_plan <- function(unresolved_intervals, objective_function,
                                        resolution_tolerance, active_probe_workers,
                                        existing_gamma_values = numeric(0)) {
  interval_rows <- lapply(names(unresolved_intervals), function(name) {
    interval <- unresolved_intervals[[name]]
    if (is.null(interval) || length(interval) != 2L || any(!is.finite(interval))) {
      return(NULL)
    }
    interval <- sort(as.numeric(interval))
    if (interval[1] >= interval[2] ||
        global_resolution_search_interval_small(
          interval[1], interval[2], objective_function, resolution_tolerance
        )) {
      return(NULL)
    }
    data.table::data.table(
      refinement_interval_id = as.integer(name),
      gamma_left = interval[1],
      gamma_right = interval[2],
      refinement_interval_width = global_resolution_search_interval_width(
        interval[1], interval[2], objective_function
      )
    )
  })
  intervals_dt <- data.table::rbindlist(interval_rows, fill = TRUE)
  if (nrow(intervals_dt) == 0L) {
    return(list(
      probe_metadata = data.table::data.table(),
      interval_summary = data.table::data.table()
    ))
  }

  data.table::setorder(intervals_dt, -refinement_interval_width, refinement_interval_id)
  n_intervals <- nrow(intervals_dt)
  intervals_dt[, refinement_points_per_interval := 1L]

  if (n_intervals < active_probe_workers) {
    max_points_per_interval <- min(8L, max(1L, floor(active_probe_workers / n_intervals)))
    max_extra_per_interval <- max(0L, max_points_per_interval - 1L)
    remaining_points <- min(
      max(0L, active_probe_workers - n_intervals),
      n_intervals * max_extra_per_interval
    )
    while (remaining_points > 0L && any(intervals_dt$refinement_points_per_interval < max_points_per_interval)) {
      for (idx in seq_len(nrow(intervals_dt))) {
        if (remaining_points <= 0L) {
          break
        }
        if (intervals_dt$refinement_points_per_interval[[idx]] < max_points_per_interval) {
          data.table::set(intervals_dt, i = idx, j = "refinement_points_per_interval",
                          value = intervals_dt$refinement_points_per_interval[[idx]] + 1L)
          remaining_points <- remaining_points - 1L
        }
      }
    }
  }

  probe_rows <- lapply(seq_len(nrow(intervals_dt)), function(idx) {
    row <- intervals_dt[idx]
    gammas <- global_resolution_search_internal_points(
      left = row$gamma_left,
      right = row$gamma_right,
      objective_function = objective_function,
      n_points = row$refinement_points_per_interval
    )
    if (length(gammas) == 0L) {
      return(NULL)
    }
    data.table::data.table(
      gamma = as.numeric(gammas),
      refinement_interval_id = as.integer(row$refinement_interval_id),
      refinement_interval_width = as.numeric(row$refinement_interval_width),
      refinement_points_per_interval = as.integer(row$refinement_points_per_interval)
    )
  })
  probe_metadata <- data.table::rbindlist(probe_rows, fill = TRUE)
  if (nrow(probe_metadata) == 0L) {
    return(list(
      probe_metadata = data.table::data.table(),
      interval_summary = intervals_dt
    ))
  }

  data.table::setorder(probe_metadata, -refinement_interval_width, refinement_interval_id, gamma)
  probe_metadata <- unique(probe_metadata, by = "gamma")
  probe_metadata <- probe_metadata[!gamma %in% as.numeric(existing_gamma_values)]
  data.table::setorder(probe_metadata, gamma)

  list(
    probe_metadata = probe_metadata,
    interval_summary = intervals_dt
  )
}

global_resolution_search_probe_batch <- function(igraph_obj, gamma_values, sweep_round,
                                                 objective_function, n_iter_preliminary,
                                                 beta_preliminary, requested_max,
                                                 min_cluster_size, snn_graph,
                                                 active_probe_workers, verbose, seed,
                                                 probe_stage = c("coarse", "refinement", "upper_cap_discovery"),
                                                 coarse_probe_count = NA_integer_,
                                                 probe_metadata = NULL) {
  probe_stage <- match.arg(probe_stage)
  gamma_values <- sort(unique(as.numeric(gamma_values[is.finite(gamma_values)])))
  if (length(gamma_values) == 0L) {
    return(data.table::data.table())
  }
  probe_metadata_dt <- if (is.null(probe_metadata)) {
    data.table::data.table(gamma = gamma_values)
  } else {
    probe_metadata_dt <- data.table::as.data.table(probe_metadata)
    if (!("gamma" %in% names(probe_metadata_dt))) {
      stop("probe_metadata must contain a gamma column.")
    }
    probe_metadata_dt <- probe_metadata_dt[is.finite(gamma)]
    probe_metadata_dt <- unique(probe_metadata_dt, by = "gamma")
    missing_gamma_values <- setdiff(gamma_values, probe_metadata_dt$gamma)
    if (length(missing_gamma_values) > 0L) {
      probe_metadata_dt <- data.table::rbindlist(
        list(
          probe_metadata_dt,
          data.table::data.table(gamma = as.numeric(missing_gamma_values))
        ),
        fill = TRUE
      )
    }
    probe_metadata_dt
  }
  data.table::setorder(probe_metadata_dt, gamma)
  probe_metadata_dt <- probe_metadata_dt[gamma %in% gamma_values]
  gamma_values <- probe_metadata_dt$gamma
  scheduled_probe_workers <- min(length(gamma_values), active_probe_workers)

  if (verbose) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH:",
        if (identical(probe_stage, "upper_cap_discovery")) "Upper-cap discovery" else paste("Sweep round", sweep_round),
        "- evaluating", length(gamma_values), paste0(probe_stage, " probe(s) with"),
        scheduled_probe_workers, "worker(s)"
      )
    )
  }

  probe_results <- cross_platform_mclapply(
    seq_along(gamma_values),
    function(idx) {
      gamma_val <- gamma_values[[idx]]
      metadata_row <- probe_metadata_dt[idx]
      probe_pid <- Sys.getpid()
      probe_start <- Sys.time()
      if (verbose) {
        scice_message(
          paste(
            "RESOLUTION_SEARCH: Probe start - stage", probe_stage,
            "| round", sweep_round,
            "| probe", idx, "/", length(gamma_values),
            "| gamma =", signif(gamma_val, 6),
            "| pid =", probe_pid
          )
        )
      }
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
      degenerate_high_gamma <- is_high_gamma_degenerate_probe(
        effective_cluster_count = effective_count,
        raw_cluster_count = raw_count,
        final_cluster_count = final_count,
        n_vertices = igraph::vcount(igraph_obj)
      )
      probe_elapsed_sec <- as.numeric(difftime(Sys.time(), probe_start, units = "secs"))
      if (verbose) {
        scice_message(
          paste(
            "RESOLUTION_SEARCH: Probe finish - stage", probe_stage,
            "| round", sweep_round,
            "| probe", idx, "/", length(gamma_values),
            "| gamma =", signif(gamma_val, 6),
            "| pid =", probe_pid,
            "| elapsed =", round(probe_elapsed_sec, 3), "sec",
            "| effective =", effective_count,
            "| raw =", raw_count,
            "| final =", final_count,
            "| degenerate_high_gamma =", degenerate_high_gamma
          )
        )
      }
      data.table::data.table(
        sweep_round = as.integer(sweep_round),
        discovery_round = if ("discovery_round" %in% names(metadata_row)) as.integer(metadata_row$discovery_round) else NA_integer_,
        probe_stage = as.character(probe_stage),
        probe_index = as.integer(idx),
        probe_pid = as.integer(probe_pid),
        probe_elapsed_sec = as.numeric(probe_elapsed_sec),
        gamma = as.numeric(gamma_val),
        upper_cap_discovery_gamma = if (identical(probe_stage, "upper_cap_discovery")) as.numeric(gamma_val) else NA_real_,
        degenerate_high_gamma = isTRUE(degenerate_high_gamma),
        scheduled_probe_workers = as.integer(scheduled_probe_workers),
        coarse_probe_count = as.integer(coarse_probe_count),
        discovered_upper_gamma = NA_real_,
        upper_cap_stop_reason = NA_character_,
        refinement_interval_width = if ("refinement_interval_width" %in% names(metadata_row)) as.numeric(metadata_row$refinement_interval_width) else NA_real_,
        refinement_interval_id = if ("refinement_interval_id" %in% names(metadata_row)) as.integer(metadata_row$refinement_interval_id) else NA_integer_,
        refinement_points_per_interval = if ("refinement_points_per_interval" %in% names(metadata_row)) as.integer(metadata_row$refinement_points_per_interval) else NA_integer_,
        effective_cluster_count = as.numeric(effective_count),
        raw_cluster_count = as.numeric(raw_count),
        final_cluster_count = as.numeric(final_count),
        raw_class = as.character(state$raw_class),
        over_fragmented = isTRUE(state$over_fragmented),
        selected_for_refinement = identical(probe_stage, "refinement"),
        selected_for_target_interval = FALSE,
        plateau_round = 0L
      )
    },
    mc.cores = scheduled_probe_workers,
    mc.preschedule = FALSE
  )
  probe_results_dt <- data.table::rbindlist(probe_results, fill = TRUE)
  if (verbose && nrow(probe_results_dt) > 0L) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Sweep round", sweep_round,
        "- completed", nrow(probe_results_dt), "probe(s)",
        "| unique probe pid(s) =", length(unique(probe_results_dt$probe_pid)),
        "| elapsed sec [min/median/max] =",
        paste(
          round(min(probe_results_dt$probe_elapsed_sec), 3),
          round(stats::median(probe_results_dt$probe_elapsed_sec), 3),
          round(max(probe_results_dt$probe_elapsed_sec), 3),
          sep = "/"
        )
      )
    )
  }

  probe_results_dt
}

discover_cpm_upper_gamma <- function(igraph_obj, gamma_bounds, requested_max,
                                     n_iter_preliminary, beta_preliminary,
                                     min_cluster_size, snn_graph,
                                     active_probe_workers, verbose, seed) {
  lower_gamma <- max(as.numeric(gamma_bounds[[1]]), .Machine$double.xmin)
  hard_cap_gamma <- max(as.numeric(gamma_bounds[[2]]), lower_gamma)
  current_gamma <- lower_gamma
  batch_size <- min(max(1L, as.integer(active_probe_workers)), 6L)
  discovery_step_ratio <- 4
  discovery_results <- data.table::data.table()
  seen_non_degenerate <- FALSE
  consecutive_degenerate <- 0L
  last_non_degenerate_gamma <- NA_real_
  discovered_upper_gamma <- hard_cap_gamma
  upper_cap_stop_reason <- "hard_cap"
  discovery_round <- 0L

  repeat {
    discovery_round <- discovery_round + 1L
    batch_gamma_values <- build_cpm_discovery_batch_gamma_values(
      current_gamma = current_gamma,
      hard_cap_gamma = hard_cap_gamma,
      batch_size = batch_size,
      step_ratio = discovery_step_ratio
    )
    batch_gamma_values <- batch_gamma_values[is.finite(batch_gamma_values)]
    if (length(batch_gamma_values) == 0L) {
      break
    }
    if (verbose) {
      scice_message(
        paste(
          "RESOLUTION_SEARCH: Upper-cap discovery round", discovery_round,
          "- batch size =", length(batch_gamma_values),
          "| step ratio =", format(signif(discovery_step_ratio, 6)),
          "| gamma values =",
          paste(signif(batch_gamma_values, 6), collapse = ", ")
        )
      )
    }
    new_probe <- global_resolution_search_probe_batch(
      igraph_obj = igraph_obj,
      gamma_values = batch_gamma_values,
      sweep_round = 0L,
      objective_function = "CPM",
      n_iter_preliminary = n_iter_preliminary,
      beta_preliminary = beta_preliminary,
      requested_max = requested_max,
      min_cluster_size = min_cluster_size,
      snn_graph = snn_graph,
      active_probe_workers = active_probe_workers,
      verbose = verbose,
      seed = seed,
      probe_stage = "upper_cap_discovery",
      probe_metadata = data.table::data.table(
        gamma = as.numeric(batch_gamma_values),
        discovery_round = as.integer(discovery_round)
      )
    )
    discovery_results <- data.table::rbindlist(list(discovery_results, new_probe), fill = TRUE)
    discovery_results <- unique(discovery_results, by = "gamma")
    data.table::setorder(discovery_results, gamma)
    batch_rows <- new_probe[order(gamma)]
    stop_found <- FALSE
    for (idx in seq_len(nrow(batch_rows))) {
      probe_row <- batch_rows[idx]
      current_gamma_value <- as.numeric(probe_row$gamma)
      current_degenerate <- isTRUE(probe_row$degenerate_high_gamma)
      current_final <- as.numeric(probe_row$final_cluster_count)

      if (!current_degenerate) {
        seen_non_degenerate <- TRUE
        consecutive_degenerate <- 0L
        last_non_degenerate_gamma <- current_gamma_value
        if (is.finite(current_final) && current_final >= requested_max) {
          discovered_upper_gamma <- current_gamma_value
          upper_cap_stop_reason <- "target_covered"
          stop_found <- TRUE
          break
        }
      } else if (seen_non_degenerate) {
        consecutive_degenerate <- consecutive_degenerate + 1L
        if (consecutive_degenerate >= 2L) {
          discovered_upper_gamma <- last_non_degenerate_gamma
          upper_cap_stop_reason <- "high_gamma_degenerate"
          stop_found <- TRUE
          break
        }
      }
    }

    if (stop_found) {
      if (verbose) {
        scice_message(
          paste(
            "RESOLUTION_SEARCH: Upper-cap discovery stop - round", discovery_round,
            "| reason =", upper_cap_stop_reason,
            "| discovered upper gamma =", signif(discovered_upper_gamma, 6)
          )
        )
      }
      break
    }

    batch_max_gamma <- max(batch_gamma_values)
    if (batch_max_gamma >= hard_cap_gamma) {
      last_batch_row <- batch_rows[.N]
      discovered_upper_gamma <- if (seen_non_degenerate &&
                                    is.finite(last_non_degenerate_gamma) &&
                                    isTRUE(last_batch_row$degenerate_high_gamma)) {
        last_non_degenerate_gamma
      } else {
        batch_max_gamma
      }
      upper_cap_stop_reason <- "hard_cap"
      if (verbose) {
        scice_message(
          paste(
            "RESOLUTION_SEARCH: Upper-cap discovery stop - round", discovery_round,
            "| reason = hard_cap",
            "| discovered upper gamma =", signif(discovered_upper_gamma, 6)
          )
        )
      }
      break
    }

    next_gamma <- min(batch_max_gamma * discovery_step_ratio, hard_cap_gamma)
    if (!is.finite(next_gamma) || next_gamma <= batch_max_gamma) {
      discovered_upper_gamma <- if (seen_non_degenerate && is.finite(last_non_degenerate_gamma)) {
        last_non_degenerate_gamma
      } else {
        batch_max_gamma
      }
      upper_cap_stop_reason <- "hard_cap"
      break
    }
    current_gamma <- next_gamma
  }

  discovered_upper_gamma <- min(max(discovered_upper_gamma, lower_gamma), hard_cap_gamma)
  discovery_results[, discovered_upper_gamma := as.numeric(discovered_upper_gamma)]
  discovery_results[, upper_cap_stop_reason := as.character(upper_cap_stop_reason)]

  list(
    probe_results = discovery_results,
    discovered_upper_gamma = as.numeric(discovered_upper_gamma),
    upper_cap_stop_reason = as.character(upper_cap_stop_reason)
  )
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
    attr(gamma_dict, "discovered_upper_gamma") <- NA_real_
    attr(gamma_dict, "upper_cap_stop_reason") <- NA_character_
    attr(gamma_dict, "coarse_probe_count") <- NA_integer_
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
  discovered_upper_gamma <- as.numeric(gamma_bounds[[2]])
  upper_cap_stop_reason <- NA_character_
  discovery_probe_results <- data.table::data.table()
  coarse_probe_count <- as.integer(min(max(3L * active_probe_workers, 12L), 30L))

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

  if (identical(objective_function, "CPM")) {
    discovery_state <- discover_cpm_upper_gamma(
      igraph_obj = igraph_obj,
      gamma_bounds = gamma_bounds,
      requested_max = requested_max,
      n_iter_preliminary = n_iter_preliminary,
      beta_preliminary = beta_preliminary,
      min_cluster_size = min_cluster_size,
      snn_graph = snn_graph,
      active_probe_workers = active_probe_workers,
      verbose = verbose,
      seed = seed
    )
    discovery_probe_results <- discovery_state$probe_results
    discovered_upper_gamma <- discovery_state$discovered_upper_gamma
    upper_cap_stop_reason <- discovery_state$upper_cap_stop_reason
  }

  if (verbose) {
    scice_message(
      paste(
        "RESOLUTION_SEARCH: Initial coarse sweep - scheduled probe workers =", active_probe_workers,
        "| coarse probe count =", coarse_probe_count,
        "| discovered upper gamma =", signif(discovered_upper_gamma, 6),
        "| upper-cap stop reason =",
        ifelse(is.na(upper_cap_stop_reason), "not_applicable", upper_cap_stop_reason)
      )
    )
  }

  all_probe_results <- global_resolution_search_probe_batch(
    igraph_obj = igraph_obj,
    gamma_values = setdiff(
      build_gamma_sequence_for_range(
        gamma_range = c(gamma_bounds[[1]], discovered_upper_gamma),
        objective_function = objective_function,
        resolution_tolerance = resolution_tolerance,
        n_vertices = n_vertices,
        n_steps = coarse_probe_count
      ),
      discovery_probe_results$gamma
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
    seed = seed,
    probe_stage = "coarse",
    coarse_probe_count = coarse_probe_count
  )
  all_probe_results <- data.table::rbindlist(list(discovery_probe_results, all_probe_results), fill = TRUE)
  all_probe_results[, discovered_upper_gamma := as.numeric(discovered_upper_gamma)]
  all_probe_results[, upper_cap_stop_reason := as.character(upper_cap_stop_reason)]
  all_probe_results[, coarse_probe_count := as.integer(coarse_probe_count)]
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

    refinement_plan <- build_refinement_probe_plan(
      unresolved_intervals = interval_state$unresolved_intervals,
      objective_function = objective_function,
      resolution_tolerance = resolution_tolerance,
      active_probe_workers = active_probe_workers,
      existing_gamma_values = all_probe_results$gamma
    )
    next_probe_metadata <- refinement_plan$probe_metadata
    interval_summary <- refinement_plan$interval_summary
    next_probe_values <- as.numeric(next_probe_metadata$gamma)

    if (length(next_probe_values) == 0L) {
      break
    }

    if (verbose) {
      interval_points_label <- if (nrow(interval_summary) > 0L) {
        paste(
          paste0(
            interval_summary$refinement_interval_id,
            ":",
            interval_summary$refinement_points_per_interval
          ),
          collapse = ", "
        )
      } else {
        "none"
      }
      scice_message(
        paste(
          "RESOLUTION_SEARCH: Sweep round", sweep_round,
          "- unresolved intervals =", nrow(interval_summary),
          "| planned refinement probes =", length(next_probe_values),
          "| interval probe counts =", interval_points_label
        )
      )
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
      seed = seed,
      probe_stage = "refinement",
      coarse_probe_count = coarse_probe_count,
      probe_metadata = next_probe_metadata
    )
    all_probe_results <- data.table::rbindlist(list(all_probe_results, new_probe_results), fill = TRUE)
    all_probe_results <- unique(all_probe_results, by = "gamma")
    data.table::setorder(all_probe_results, gamma)
    all_probe_results[, discovered_upper_gamma := as.numeric(discovered_upper_gamma)]
    all_probe_results[, upper_cap_stop_reason := as.character(upper_cap_stop_reason)]
    all_probe_results[, coarse_probe_count := as.integer(coarse_probe_count)]

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
  attr(interval_state$gamma_dict, "discovered_upper_gamma") <- as.numeric(discovered_upper_gamma)
  attr(interval_state$gamma_dict, "upper_cap_stop_reason") <- upper_cap_stop_reason
  attr(interval_state$gamma_dict, "coarse_probe_count") <- as.integer(coarse_probe_count)
  interval_state$gamma_dict
}
