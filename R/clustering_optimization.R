# Apply raw-exact-first admission with progressively looser effective/raw fallbacks.
select_gamma_admission <- function(strict_flags, relaxed_flags,
                                   soft_guard_flags, hard_guard_flags,
                                   raw_strict_flags = NULL,
                                   raw_relaxed_flags = NULL) {
  if (is.null(raw_strict_flags)) {
    raw_strict_flags <- rep(FALSE, length(strict_flags))
  }
  if (is.null(raw_relaxed_flags)) {
    raw_relaxed_flags <- rep(FALSE, length(relaxed_flags))
  }

  candidate_sets <- list(
    raw_strict_soft = which(raw_strict_flags & soft_guard_flags),
    strict_soft = which(strict_flags & soft_guard_flags),
    relaxed_soft = which(relaxed_flags & soft_guard_flags),
    strict_hard = which(strict_flags & hard_guard_flags),
    relaxed_hard = which(relaxed_flags & hard_guard_flags),
    relaxed_unguarded = which(relaxed_flags),
    raw_relaxed_soft = which(raw_relaxed_flags & soft_guard_flags),
    raw_relaxed_hard = which(raw_relaxed_flags & hard_guard_flags),
    raw_relaxed_unguarded = which(raw_relaxed_flags)
  )

  for (mode in names(candidate_sets)) {
    indices <- candidate_sets[[mode]]
    if (length(indices) > 0L) {
      return(list(indices = indices, mode = mode))
    }
  }

  list(indices = integer(0), mode = "none")
}

extract_raw_median_gap <- function(result, target_clusters) {
  if (is.null(result$raw_median_gap)) {
    return(abs(as.numeric(result$mean_clusters) - target_clusters))
  }
  as.numeric(result$raw_median_gap)
}

refine_gamma_candidates_by_raw_gap <- function(valid_indices, admission_mode,
                                               gamma_results, target_clusters,
                                               min_cluster_size = 1L) {
  if (length(valid_indices) == 0L || min_cluster_size <= 1L) {
    return(list(
      indices = valid_indices,
      mode = admission_mode,
      raw_gaps = numeric(0),
      best_raw_gap = Inf
    ))
  }

  selected_raw_gaps <- vapply(
    gamma_results[valid_indices],
    extract_raw_median_gap,
    numeric(1),
    target_clusters = target_clusters
  )
  best_raw_gap <- Inf
  if (any(is.finite(selected_raw_gaps))) {
    best_raw_gap <- min(selected_raw_gaps[is.finite(selected_raw_gaps)])
  }

  if (length(valid_indices) > 1L && any(is.finite(selected_raw_gaps))) {
    keep_mask <- is.finite(selected_raw_gaps) & selected_raw_gaps == best_raw_gap
    if (any(keep_mask) && sum(keep_mask) < length(valid_indices)) {
      valid_indices <- valid_indices[keep_mask]
      selected_raw_gaps <- selected_raw_gaps[keep_mask]
    }
  }

  list(
    indices = valid_indices,
    mode = admission_mode,
    raw_gaps = selected_raw_gaps,
    best_raw_gap = best_raw_gap
  )
}

#' Sample one value using a fixed seed without leaking RNG state
#'
#' @param values Candidate values
#' @param seed Integer seed value
#' @return One sampled value
#' @keywords internal
sample_with_fixed_seed <- function(values, seed = 1L) {
  values <- as.integer(values)
  if (length(values) == 1L) {
    return(values)
  }

  has_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (has_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  as.integer(sample(values, size = 1))
}

#' Merge undersized clusters in one fast deterministic pass
#'
#' @description
#' Reassigns each undersized cluster to one eligible target cluster based on mean
#' SNN connectivity, without iterative re-merging loops.
#'
#' @param labels Integer cluster labels (0-based)
#' @param snn_graph Seurat SNN graph matrix
#' @param min_cluster_size Minimum required cells per cluster
#' @return Integer labels after one-pass merging (contiguous 0-based)
#' @keywords internal
merge_small_clusters_to_neighbors <- function(labels, snn_graph, min_cluster_size = 1L) {
  labels <- as.integer(labels)
  min_cluster_size <- max(1L, as.integer(min_cluster_size))

  if (min_cluster_size <= 1L || is.null(snn_graph) || length(labels) == 0L) {
    return(labels)
  }

  if (nrow(snn_graph) != length(labels) || ncol(snn_graph) != length(labels)) {
    stop("snn_graph dimensions must match label length when min_cluster_size > 1.")
  }

  base_labels <- reindex_cluster_labels(labels)
  cluster_sizes <- tabulate(base_labels + 1L)
  cluster_ids <- which(cluster_sizes > 0L) - 1L

  if (length(cluster_ids) <= 1L) {
    return(base_labels)
  }

  large_cluster_ids <- cluster_ids[cluster_sizes[cluster_ids + 1L] >= min_cluster_size]
  small_cluster_ids <- cluster_ids[cluster_sizes[cluster_ids + 1L] < min_cluster_size]

  if (length(small_cluster_ids) == 0L) {
    return(base_labels)
  }

  # If all clusters are undersized, collapse to the largest cluster (tie -> smallest id).
  if (length(large_cluster_ids) == 0L) {
    largest_size <- max(cluster_sizes[cluster_ids + 1L])
    largest_ids <- cluster_ids[cluster_sizes[cluster_ids + 1L] == largest_size]
    target_id <- min(largest_ids)
    merged <- base_labels
    merged[] <- as.integer(target_id)
    return(reindex_cluster_labels(merged))
  }

  target_map <- integer(length(small_cluster_ids))
  for (idx in seq_along(small_cluster_ids)) {
    small_id <- small_cluster_ids[idx]
    small_cells <- which(base_labels == small_id)
    connectivity <- vapply(large_cluster_ids, function(candidate_id) {
      candidate_cells <- which(base_labels == candidate_id)
      if (length(candidate_cells) == 0L) {
        return(-Inf)
      }
      sub_snn <- snn_graph[small_cells, candidate_cells, drop = FALSE]
      as.numeric(sum(sub_snn)) / (length(small_cells) * length(candidate_cells))
    }, numeric(1))

    if (length(connectivity) == 0L || all(!is.finite(connectivity))) {
      target_map[idx] <- min(large_cluster_ids)
      next
    }

    max_connectivity <- max(connectivity, na.rm = TRUE)
    tied_candidates <- large_cluster_ids[which(connectivity == max_connectivity)]
    target_map[idx] <- min(tied_candidates)
  }

  merged <- base_labels
  for (idx in seq_along(small_cluster_ids)) {
    small_id <- small_cluster_ids[idx]
    merged[base_labels == small_id] <- as.integer(target_map[idx])
  }

  reindex_cluster_labels(merged)
}

derive_manual_resolution_seed <- function(seed, resolution, offset = 0L) {
  if (is.null(seed)) {
    return(NULL)
  }

  gamma_component <- if (is.finite(resolution)) {
    as.integer(floor(abs(resolution * 100000) %% 100000000))
  } else {
    0L
  }

  derived_seed <- as.integer(
    (as.double(seed) + as.double(gamma_component) + as.double(offset)) %% .Machine$integer.max
  )
  if (is.na(derived_seed) || derived_seed <= 0L) {
    derived_seed <- 1L
  }
  derived_seed
}

finalize_selected_clustering <- function(matrix_ref, gamma, effective_cluster_median,
                                         raw_cluster_median, final_cluster_median = NA_real_,
                                         admission_mode,
                                         cluster_seed = NULL, n_bootstrap,
                                         n_workers, snn_graph = NULL,
                                         target_clusters = NULL,
                                         preferred_trial_indices = NULL,
                                         min_cluster_size = 1L, verbose = FALSE,
                                         worker_id = "OPTIMIZER",
                                         runtime_context = NULL) {
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  if (min_cluster_size > 1L && is.null(snn_graph)) {
    stop("snn_graph must be provided when min_cluster_size > 1.")
  }

  on.exit(release_cluster_matrix(matrix_ref), add = TRUE)
  best_clustering <- load_cluster_matrix(matrix_ref)
  heartbeat <- create_heartbeat_logger(verbose = verbose, context = worker_id)
  n_vertices <- nrow(best_clustering)
  n_trials <- ncol(best_clustering)
  preferred_trial_indices <- as.integer(preferred_trial_indices)
  preferred_trial_indices <- preferred_trial_indices[is.finite(preferred_trial_indices)]
  preferred_trial_indices <- preferred_trial_indices[
    preferred_trial_indices >= 1L & preferred_trial_indices <= n_trials
  ]
  preferred_trial_indices <- unique(preferred_trial_indices)
  if (length(preferred_trial_indices) == 0L && !is.null(target_clusters)) {
    target_clusters <- as.integer(target_clusters)[1]
    final_clusters_vec <- if (min_cluster_size > 1L) {
      vapply(
        seq_len(n_trials),
        function(trial_idx) {
          merged_labels <- merge_small_clusters_to_neighbors(
            labels = best_clustering[, trial_idx],
            snn_graph = snn_graph,
            min_cluster_size = min_cluster_size
          )
          as.integer(length(unique(merged_labels)))
        },
        integer(1)
      )
    } else {
      vapply(
        seq_len(n_trials),
        function(trial_idx) as.integer(length(unique(best_clustering[, trial_idx]))),
        integer(1)
      )
    }
    preferred_trial_indices <- which(final_clusters_vec == target_clusters)
  }

  if (verbose) {
    scice_message(paste(worker_id, ": Phase 5 - Bootstrap analysis with", n_bootstrap, "iterations"))
    bootstrap_start <- Sys.time()
  }

  bootstrap_workers <- cap_workers_by_memory(
    max(1L, as.integer(n_workers)),
    estimate_trial_matrix_bytes(n_vertices, n_trials, 1L),
    runtime_context
  )
  bootstrap_workers <- min(bootstrap_workers, max(1L, as.integer(n_bootstrap)))
  bootstrap_log_every <- max(1L, as.integer(floor(n_bootstrap / 5)))
  should_log_bootstrap_step <- function(idx) {
    idx == 1L || idx == n_bootstrap || (idx %% bootstrap_log_every) == 0L
  }

  if (bootstrap_workers == 1) {
    ic_bootstrap <- numeric(n_bootstrap)

    for (i in seq_len(n_bootstrap)) {
      if (!is.null(cluster_seed)) {
        bootstrap_seed <- cluster_seed + 10000 + i
        set.seed(bootstrap_seed)
      }

      sample_indices <- sample.int(ncol(best_clustering), ncol(best_clustering), replace = TRUE)
      bootstrap_matrix <- best_clustering[, sample_indices, drop = FALSE]
      extracted <- extract_clustering_array(bootstrap_matrix)
      ic_result <- calculate_ic_from_extracted(extracted)
      ic_bootstrap[i] <- 1 / ic_result
      heartbeat(function() {
        paste(
          "phase5 running - bootstrap", i, "/", n_bootstrap,
          "- latest IC =", round(ic_bootstrap[i], 4)
        )
      })

      if (verbose && should_log_bootstrap_step(i)) {
        scice_message(paste(worker_id, ": Phase 5 progress", i, "/", n_bootstrap, "- IC =", round(ic_bootstrap[i], 4)))
      }
    }
  } else {
    ic_bootstrap <- cross_platform_mclapply(seq_len(n_bootstrap), function(i) {
      if (!is.null(cluster_seed)) {
        bootstrap_seed <- cluster_seed + 10000 + i
        set.seed(bootstrap_seed)
      }

      sample_indices <- sample.int(ncol(best_clustering), ncol(best_clustering), replace = TRUE)
      bootstrap_matrix <- best_clustering[, sample_indices, drop = FALSE]
      extracted <- extract_clustering_array(bootstrap_matrix)
      ic_result <- calculate_ic_from_extracted(extracted)
      ic_value <- 1 / ic_result
      if (verbose && should_log_bootstrap_step(i)) {
        scice_message(paste(worker_id, ": Phase 5 progress", i, "/", n_bootstrap, "- IC =", round(ic_value, 4)))
      }
      ic_value
    }, mc.cores = bootstrap_workers)

    ic_bootstrap <- unlist(ic_bootstrap)
  }

  ic_median <- stats::median(ic_bootstrap)

  if (verbose) {
    bootstrap_time <- as.numeric(difftime(Sys.time(), bootstrap_start, units = "secs"))
    scice_message(paste(worker_id, ": Bootstrap analysis completed in", round(bootstrap_time, 3), "seconds"))
    scice_message(paste(worker_id, ": Bootstrap IC median:", round(ic_median, 4)))
    scice_message(paste(worker_id, ": Bootstrap IC range: [", round(min(ic_bootstrap), 4), ", ", round(max(ic_bootstrap), 4), "]", sep = ""))
    scice_message(paste(worker_id, ": OPTIMIZATION COMPLETE"))
  }

  extracted_all <- extract_clustering_array(best_clustering)
  selection_matrix <- best_clustering
  if (length(preferred_trial_indices) > 0L) {
    selection_matrix <- best_clustering[, preferred_trial_indices, drop = FALSE]
    if (verbose) {
      scice_message(
        paste(
          worker_id, ": Selecting representative best_labels from",
          length(preferred_trial_indices), "exact final-hit trial(s)",
          if (!is.null(target_clusters)) paste("for target", target_clusters) else ""
        )
      )
    }
  }
  extracted_best <- extract_clustering_array(selection_matrix)
  best_labels_raw <- get_best_clustering(extracted_best)
  best_labels <- if (min_cluster_size > 1L) {
    merge_small_clusters_to_neighbors(
      labels = best_labels_raw,
      snn_graph = snn_graph,
      min_cluster_size = min_cluster_size
    )
  } else {
    best_labels_raw
  }
  best_labels_raw_cluster_count <- as.integer(length(unique(best_labels_raw)))
  best_labels_final_cluster_count <- as.integer(length(unique(best_labels)))

  if (verbose && min_cluster_size > 1L) {
    scice_message(
      paste(
        worker_id, ": Final best_labels merged small clusters to satisfy min_cluster_size (value =",
        min_cluster_size,
        "; final clusters =",
        best_labels_final_cluster_count, ")"
      )
    )
  }
  if (verbose) {
    scice_message(
      paste(
        worker_id, ": Selected diagnostics - gamma =", signif(gamma, 6),
        "- effective_median =", signif(effective_cluster_median, 6),
        "- raw_median =", signif(raw_cluster_median, 6),
        "- final_median =", signif(final_cluster_median, 6),
        "- admission_mode =", admission_mode,
        "- best_labels_raw_clusters =", best_labels_raw_cluster_count,
        "- best_labels_final_clusters =", best_labels_final_cluster_count
      )
    )
  }

  rm(best_clustering)

  list(
    gamma = gamma,
    labels = extracted_all,
    ic_median = ic_median,
    ic_bootstrap = ic_bootstrap,
    best_labels = best_labels,
    effective_cluster_median = effective_cluster_median,
    raw_cluster_median = raw_cluster_median,
    final_cluster_median = final_cluster_median,
    admission_mode = admission_mode,
    best_labels_raw_cluster_count = best_labels_raw_cluster_count,
    best_labels_final_cluster_count = best_labels_final_cluster_count
  )
}

#' Optimize clustering within a resolution range
#'
#' @description
#' Refines one target cluster number inside a candidate gamma interval.
#'
#' @details
#' The optimizer samples gamma values within the interval, runs repeated Leiden
#' trials, retains gammas whose effective cluster count matches the target, then
#' selects the best gamma by IC score and performs bootstrap stability assessment.
#'
#' Returned values include best labels, IC bootstrap vector, iteration metadata,
#' and supporting diagnostics used by higher-level summary functions.
#'
#' @param igraph_obj igraph object to cluster
#' @param target_clusters Target cluster count
#' @param gamma_range Candidate resolution range
#' @param objective_function Objective function ("modularity" or "CPM")
#' @param n_trials Number of clustering trials per gamma
#' @param n_bootstrap Number of bootstrap iterations
#' @param seed Optional deterministic seed
#' @param beta Leiden beta parameter
#' @param n_iterations Initial Leiden iterations
#' @param max_iterations Maximum iterations for refinement
#' @param resolution_tolerance Resolution tolerance for gamma sequence logic
#' @param n_workers Requested worker count
#' @param snn_graph Seurat SNN graph matrix for final best-label merging
#' @param min_cluster_size Minimum required cells per effective cluster count
#' @param verbose Whether to print diagnostics
#' @param worker_id Worker label used in logs
#' @param in_parallel_context Whether called from an outer parallel worker
#' @param runtime_context Internal mutable context for spill/memory controls
#' @return List with best gamma, labels, IC summaries, and iteration metadata
#' @keywords internal
optimize_clustering <- function(igraph_obj, target_clusters, gamma_range, objective_function,
                               n_trials, n_bootstrap, seed = NULL, beta, n_iterations, max_iterations,
                               resolution_tolerance, n_workers, snn_graph = NULL,
                               gamma_seed_values = NULL,
                               min_cluster_size = 1L, verbose = FALSE, worker_id = "OPTIMIZER", 
                               in_parallel_context = FALSE, runtime_context = NULL) {
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  if (min_cluster_size > 1L && is.null(snn_graph)) {
    stop("snn_graph must be provided when min_cluster_size > 1.")
  }

  run_leiden_trial <- function(resolution, iterations, initial_membership = NULL) {
    leiden_clustering(
      igraph_obj, resolution, objective_function, iterations, beta, initial_membership
    )
  }
  
  # Set deterministic seeds for this cluster number if base seed provided
  cluster_seed <- NULL
  if (!is.null(seed)) {
    cluster_seed <- seed + target_clusters * 1000  # Different seed per cluster number
    set.seed(cluster_seed)
    if (verbose) {
      scice_message(paste(worker_id, ": Set deterministic seed:", cluster_seed))
    }
  }
  
  if (verbose) {
    raw_guard_limits <- raw_cluster_guard_limits(target_clusters)
    scice_message(paste(worker_id, ": Optimization parameters:"))
    scice_message(paste(worker_id, ":   Target clusters:", target_clusters))
    scice_message(paste(worker_id, ":   Trials per gamma:", n_trials))
    scice_message(paste(worker_id, ":   Bootstrap iterations:", n_bootstrap))
    scice_message(paste(worker_id, ":   Max iterations:", max_iterations))
    scice_message(paste(worker_id, ":   Beta:", beta))
    scice_message(paste(worker_id, ":   Leiden iterations:", n_iterations))
    scice_message(
      paste(
        worker_id,
        ":   Gamma admission: strict-first (as.integer(median)==target), fallback relaxed (any-hit and |median-target|<=1); IC uses all trials at admitted gamma"
      )
    )
    if (min_cluster_size > 1L) {
      scice_message(
        paste(
          worker_id, ":   Raw-cluster guard: soft <=", raw_guard_limits$soft,
          "hard <=", raw_guard_limits$hard
        )
      )
      scice_message(
        paste(
          worker_id, ":   Counting uses effective clusters (size >=", min_cluster_size,
          "); final merge applied only on best_labels"
        )
      )
    }
  }
  
  n_vertices <- igraph::vcount(igraph_obj)
  heartbeat <- create_heartbeat_logger(verbose = verbose, context = worker_id)
  range_width <- abs(gamma_range[2] - gamma_range[1])
  n_steps <- if (n_vertices >= 200000 && range_width <= max(resolution_tolerance * 10, .Machine$double.eps)) 5L else 11L
  delta_n <- 2
  
  if (verbose) {
    scice_message(paste(worker_id, ": Creating gamma sequence with", n_steps, "steps"))
    scice_message(paste(worker_id, ": Gamma range bounds: [", signif(gamma_range[1], 6), ", ", signif(gamma_range[2], 6), "]", sep = ""))
  }
  
  # Create gamma sequence using the same helper used by resolution-range plateau clamping.
  gamma_sequence <- build_gamma_sequence_for_range(
    gamma_range,
    objective_function = objective_function,
    resolution_tolerance = resolution_tolerance,
    n_vertices = n_vertices
  )
  gamma_seed_values <- as.numeric(gamma_seed_values)
  gamma_seed_values <- gamma_seed_values[is.finite(gamma_seed_values)]
  if (length(gamma_seed_values) > 0L) {
    gamma_lower <- min(gamma_range)
    gamma_upper <- max(gamma_range)
    tolerance <- max(.Machine$double.eps^0.5, abs(gamma_upper - gamma_lower) * 1e-8)
    gamma_seed_values <- gamma_seed_values[
      gamma_seed_values >= (gamma_lower - tolerance) &
        gamma_seed_values <= (gamma_upper + tolerance)
    ]
    if (length(gamma_seed_values) > 0L) {
      gamma_sequence <- sort(unique(c(gamma_sequence, gamma_seed_values)))
    }
  }

  if (verbose) {
    scice_message(paste(worker_id, ": Generated gamma sequence:"))
    for (i in 1:min(5, length(gamma_sequence))) {
      scice_message(paste(worker_id, ":   gamma[", i, "] = ", signif(gamma_sequence[i], 6), sep = ""))
    }
    if (length(gamma_sequence) > 5) {
      scice_message(paste(worker_id, ":   ... and", length(gamma_sequence) - 5, "more"))
    }
    if (length(gamma_seed_values) > 0L) {
      scice_message(
        paste(
          worker_id, ": Included", length(unique(gamma_seed_values)),
          "gamma seed values from shared search probes"
        )
      )
    }
  }
  if (verbose) {
    phase1_runs <- as.integer(length(gamma_sequence) * max(1L, as.integer(n_trials)))
    scice_message(paste(worker_id, ": Phase 1 expected Leiden runs:", format(phase1_runs, big.mark = ",")))
  }

  estimated_phase1_bytes <- estimate_trial_matrix_bytes(n_vertices, n_trials, length(gamma_sequence))
  if (should_enable_spill(runtime_context, estimated_phase1_bytes)) {
    activate_runtime_spill(runtime_context, estimated_bytes = estimated_phase1_bytes)
  }
  
  # Test initial clustering for each gamma and compute IC immediately for valid gammas
  if (verbose) {
    scice_message(paste(worker_id, ": Phase 1 - Testing", length(gamma_sequence), "gamma values with", n_trials, "trials each"))
    initial_phase_start <- Sys.time()
    if (in_parallel_context) {
      scice_message(
        paste(
          worker_id, ": Running in parallel context - deriving nested worker budget from",
          n_workers, "optimizer workers across", length(gamma_sequence), "gamma values"
        )
      )
    }
  }
  
  nested_workers <- if (in_parallel_context) {
    max(
      1L,
      as.integer(round(as.double(n_workers) / as.double(max(1L, length(gamma_sequence)))))
    )
  } else {
    max(1L, as.integer(n_workers))
  }
  if (!in_parallel_context) {
    nested_workers <- min(nested_workers, max(1L, as.integer(length(gamma_sequence))))
  }
  nested_workers <- cap_workers_by_memory(
    nested_workers,
    estimate_trial_matrix_bytes(n_vertices, n_trials, 1L),
    runtime_context
  )
  
  if (verbose) {
    scice_message(
      paste(
        worker_id, ": Worker budget for this cluster:", n_workers,
        "- using", nested_workers, "workers for Phase 1 gamma evaluation"
      )
    )
  }
  
  phase1_log_every <- max(1L, as.integer(floor(length(gamma_sequence) / 5)))
  should_log_phase1_step <- function(idx) {
    idx == 1L || idx == length(gamma_sequence) || (idx %% phase1_log_every) == 0L
  }
  
  evaluate_gamma <- function(gamma_val, gamma_idx) {
    log_this_gamma <- verbose && should_log_phase1_step(gamma_idx)
    gamma_start_time <- NULL
    if (log_this_gamma) {
      gamma_start_time <- Sys.time()
      scice_message(
        paste(
          worker_id, ": Phase 1 progress gamma", gamma_idx, "/", length(gamma_sequence),
          "started (gamma =", signif(gamma_val, 6), ")"
        )
      )
    }

    if (!is.null(seed)) {
      gamma_component <- if (is.finite(gamma_val)) {
        as.integer(floor(abs(gamma_val * 10000) %% 100000))
      } else {
        0L
      }
      gamma_seed <- as.integer((as.double(cluster_seed) + gamma_component) %% .Machine$integer.max)
      if (is.na(gamma_seed) || gamma_seed <= 0L) {
        gamma_seed <- 1L
      }
      set.seed(gamma_seed)
    }

    cluster_matrix <- matrix(0L, nrow = n_vertices, ncol = n_trials)
    for (trial_idx in seq_len(n_trials)) {
      cluster_matrix[, trial_idx] <- run_leiden_trial(gamma_val, n_iterations)
      heartbeat(function() {
        paste(
          "phase1 running - gamma", gamma_idx, "/", length(gamma_sequence),
          "- trial", trial_idx, "/", n_trials,
          "- target k =", target_clusters
        )
      })
    }

    effective_clusters_vec <- vapply(
      seq_len(n_trials),
      function(trial_idx) {
        count_effective_clusters(
          cluster_matrix[, trial_idx],
          min_cluster_size = min_cluster_size
        )
      },
      integer(1)
    )
    raw_clusters_vec <- vapply(
      seq_len(n_trials),
      function(trial_idx) {
        as.integer(length(unique(cluster_matrix[, trial_idx])))
      },
      integer(1)
    )
    final_clusters_vec <- if (min_cluster_size > 1L) {
      vapply(
        seq_len(n_trials),
        function(trial_idx) {
          merged_labels <- merge_small_clusters_to_neighbors(
            labels = cluster_matrix[, trial_idx],
            snn_graph = snn_graph,
            min_cluster_size = min_cluster_size
          )
          as.integer(length(unique(merged_labels)))
        },
        integer(1)
      )
    } else {
      raw_clusters_vec
    }
    median_effective_clusters <- stats::median(effective_clusters_vec)
    median_clusters_int <- as.integer(median_effective_clusters)
    raw_cluster_median <- stats::median(raw_clusters_vec)
    raw_cluster_median_int <- as.integer(raw_cluster_median)
    final_cluster_median <- stats::median(final_clusters_vec)
    final_cluster_median_int <- as.integer(final_cluster_median)
    hit_trials <- which(final_clusters_vec == target_clusters)
    hit_count <- as.integer(length(hit_trials))
    raw_hit_trials <- which(raw_clusters_vec == target_clusters)
    raw_hit_count <- as.integer(length(raw_hit_trials))
    hit_rate <- as.double(hit_count) / as.double(n_trials)
    median_gap <- abs(final_cluster_median - target_clusters)
    raw_median_gap <- abs(raw_cluster_median - target_clusters)
    within_median_window <- (median_gap <= 1)
    strict_valid <- (final_cluster_median_int == target_clusters)
    relaxed_valid <- (hit_count >= 1L) && within_median_window
    raw_within_median_window <- (raw_median_gap <= 1)
    raw_strict_valid <- (raw_cluster_median_int == target_clusters)
    raw_relaxed_valid <- (raw_hit_count >= 1L) && raw_within_median_window
    raw_guard_soft <- passes_raw_cluster_guard(
      raw_cluster_median,
      target_clusters,
      min_cluster_size = min_cluster_size,
      level = "soft"
    )
    raw_guard_hard <- passes_raw_cluster_guard(
      raw_cluster_median,
      target_clusters,
      min_cluster_size = min_cluster_size,
      level = "hard"
    )
    gamma_admitted <- strict_valid || relaxed_valid || raw_strict_valid || raw_relaxed_valid

    if (!gamma_admitted) {
      rm(cluster_matrix)
      if (log_this_gamma) {
        gamma_elapsed <- as.numeric(difftime(Sys.time(), gamma_start_time, units = "secs"))
        scice_message(
          paste(
            worker_id, ": Phase 1 progress gamma", gamma_idx, "/", length(gamma_sequence),
            "completed in", round(gamma_elapsed, 3), "seconds",
            "- median_effective =", signif(median_effective_clusters, 6),
            "- median_final =", signif(final_cluster_median, 6),
            "- median_raw =", signif(raw_cluster_median, 6),
            "- median gap =", round(median_gap, 3),
            "- final hit trials =", hit_count, "/", n_trials,
            "- raw hit trials =", raw_hit_count, "/", n_trials,
            "- strict_valid =", strict_valid,
            "- relaxed_valid =", relaxed_valid,
            "- raw_strict_valid =", raw_strict_valid,
            "- raw_relaxed_valid =", raw_relaxed_valid,
            "- raw_guard_soft =", raw_guard_soft,
            "- raw_guard_hard =", raw_guard_hard,
            "(target =", target_clusters, "; IC skipped)"
          )
        )
      }
      return(list(
        valid = FALSE,
        gamma = gamma_val,
        mean_clusters = final_cluster_median,
        median_clusters_raw = final_cluster_median,
        median_effective_clusters = median_effective_clusters,
        median_clusters_int = median_clusters_int,
        final_cluster_median = final_cluster_median,
        final_cluster_median_int = final_cluster_median_int,
        raw_cluster_median = raw_cluster_median,
        raw_cluster_median_int = raw_cluster_median_int,
        median_gap = median_gap,
        raw_median_gap = raw_median_gap,
        within_median_window = within_median_window,
        strict_valid = strict_valid,
        relaxed_valid = relaxed_valid,
        raw_within_median_window = raw_within_median_window,
        raw_strict_valid = raw_strict_valid,
        raw_relaxed_valid = raw_relaxed_valid,
        hit_count = hit_count,
        hit_rate = hit_rate,
        effective_hit_count = hit_count,
        raw_hit_count = raw_hit_count,
        raw_guard_soft = raw_guard_soft,
        raw_guard_hard = raw_guard_hard
      ))
    }

    extracted <- extract_clustering_array(cluster_matrix)
    ic_result <- calculate_ic_from_extracted(extracted)
    ic_score <- 1 / ic_result
    rm(extracted)

    matrix_ref <- store_cluster_matrix(
      cluster_matrix,
      runtime_context = runtime_context,
      prefix = sprintf("k%d_g%03d", target_clusters, gamma_idx)
    )
    rm(cluster_matrix)
    
    if (log_this_gamma) {
      gamma_elapsed <- as.numeric(difftime(Sys.time(), gamma_start_time, units = "secs"))
      scice_message(
        paste(
          worker_id, ": Phase 1 progress gamma", gamma_idx, "/", length(gamma_sequence),
          "completed in", round(gamma_elapsed, 3), "seconds",
          "- median_effective =", signif(median_effective_clusters, 6),
          "- median_final =", signif(final_cluster_median, 6),
          "- median_raw =", signif(raw_cluster_median, 6),
          "- median gap =", round(median_gap, 3),
          "- final hit trials =", hit_count, "/", n_trials,
          "- raw hit trials =", raw_hit_count, "/", n_trials,
          "- strict_valid =", strict_valid,
          "- relaxed_valid =", relaxed_valid,
          "- raw_strict_valid =", raw_strict_valid,
          "- raw_relaxed_valid =", raw_relaxed_valid,
          "- raw_guard_soft =", raw_guard_soft,
          "- raw_guard_hard =", raw_guard_hard,
          "- IC (all trials) =", round(ic_score, 4)
        )
      )
    }

    list(
      valid = TRUE,
      gamma = gamma_val,
      mean_clusters = final_cluster_median,
      median_clusters_raw = final_cluster_median,
      median_effective_clusters = median_effective_clusters,
      median_clusters_int = median_clusters_int,
      final_cluster_median = final_cluster_median,
      final_cluster_median_int = final_cluster_median_int,
      raw_cluster_median = raw_cluster_median,
      raw_cluster_median_int = raw_cluster_median_int,
      median_gap = median_gap,
      raw_median_gap = raw_median_gap,
      within_median_window = within_median_window,
      strict_valid = strict_valid,
      relaxed_valid = relaxed_valid,
      raw_within_median_window = raw_within_median_window,
      raw_strict_valid = raw_strict_valid,
      raw_relaxed_valid = raw_relaxed_valid,
      hit_count = hit_count,
      hit_trials = hit_trials,
      hit_rate = hit_rate,
      effective_hit_count = hit_count,
      raw_hit_count = raw_hit_count,
      raw_guard_soft = raw_guard_soft,
      raw_guard_hard = raw_guard_hard,
      ic = ic_score,
      matrix_ref = matrix_ref
    )
  }

  if (nested_workers == 1) {
    gamma_results <- vector("list", length(gamma_sequence))
    
    for (gamma_idx in seq_along(gamma_sequence)) {
      gamma_val <- gamma_sequence[gamma_idx]
      gamma_results[[gamma_idx]] <- evaluate_gamma(gamma_val, gamma_idx)
    }
  } else {
    gamma_results <- cross_platform_mclapply(seq_along(gamma_sequence), function(gamma_idx) {
      gamma_val <- gamma_sequence[gamma_idx]
      evaluate_gamma(gamma_val, gamma_idx)
    }, mc.cores = nested_workers)
  }
  
  if (verbose) {
    initial_phase_time <- as.numeric(difftime(Sys.time(), initial_phase_start, units = "secs"))
    scice_message(paste(worker_id, ": Phase 1 completed in", round(initial_phase_time, 3), "seconds"))
  }
  
  mean_clusters <- vapply(
    gamma_results,
    function(x) {
      if (!is.null(x$final_cluster_median)) {
        return(as.numeric(x$final_cluster_median))
      }
      if (!is.null(x$median_effective_clusters)) {
        return(as.numeric(x$median_effective_clusters))
      }
      if (!is.null(x$median_clusters_raw)) {
        return(as.numeric(x$median_clusters_raw))
      }
      as.numeric(x$mean_clusters)
    },
    numeric(1)
  )
  mean_clusters_int <- vapply(
    gamma_results,
    function(x) {
      if (!is.null(x$final_cluster_median_int)) {
        return(as.integer(x$final_cluster_median_int))
      }
      if (!is.null(x$median_clusters_int)) {
        return(as.integer(x$median_clusters_int))
      }
      as.integer(x$mean_clusters)
    },
    integer(1)
  )
  hit_counts <- vapply(
    gamma_results,
    function(x) {
      if (!is.null(x$effective_hit_count)) {
        return(as.integer(x$effective_hit_count))
      }
      if (is.null(x$hit_count)) {
        return(0L)
      }
      as.integer(x$hit_count)
    },
    integer(1)
  )
  within_median_window_flags <- vapply(
    gamma_results,
    function(x) isTRUE(x$within_median_window),
    logical(1)
  )
  strict_flags <- vapply(
    gamma_results,
    function(x) isTRUE(x$strict_valid),
    logical(1)
  )
  relaxed_flags <- vapply(
    gamma_results,
    function(x) isTRUE(x$relaxed_valid),
    logical(1)
  )
  raw_strict_flags <- vapply(
    gamma_results,
    function(x) isTRUE(x$raw_strict_valid),
    logical(1)
  )
  raw_relaxed_flags <- vapply(
    gamma_results,
    function(x) isTRUE(x$raw_relaxed_valid),
    logical(1)
  )
  soft_raw_guard_flags <- vapply(
    gamma_results,
    function(x) {
      if (is.null(x$raw_guard_soft)) {
        return(TRUE)
      }
      isTRUE(x$raw_guard_soft)
    },
    logical(1)
  )
  hard_raw_guard_flags <- vapply(
    gamma_results,
    function(x) {
      if (is.null(x$raw_guard_hard)) {
        return(TRUE)
      }
      isTRUE(x$raw_guard_hard)
    },
    logical(1)
  )
  
  if (verbose) {
    scice_message(paste(worker_id, ": Phase 2 - Filtering for target final merged cluster count (strict-first with relaxed fallback):", target_clusters))
    cluster_counts <- table(mean_clusters)
    for (i in 1:length(cluster_counts)) {
      count_val <- names(cluster_counts)[i]
      freq <- cluster_counts[i]
      scice_message(paste(worker_id, ":   ", freq, "gammas -> ", count_val, " final merged clusters", sep = ""))
    }
    cluster_counts_int <- table(mean_clusters_int)
    for (i in 1:length(cluster_counts_int)) {
      count_val <- names(cluster_counts_int)[i]
      freq <- cluster_counts_int[i]
      scice_message(paste(worker_id, ":   ", freq, "gammas -> median_int ", count_val, sep = ""))
    }
    scice_message(
      paste(
        worker_id, ":   gammas passing median window (|median-target|<=1):",
        sum(within_median_window_flags), "/", length(gamma_results)
      )
    )
    scice_message(
      paste(
        worker_id, ":   strict-valid gammas (as.integer(median)==target):",
        sum(strict_flags), "/", length(gamma_results)
      )
    )
    scice_message(
      paste(
        worker_id, ":   relaxed-valid gammas (any-hit + median-window):",
        sum(relaxed_flags), "/", length(gamma_results)
      )
    )
    scice_message(
      paste(
        worker_id, ":   raw-strict-valid gammas (as.integer(raw_median)==target):",
        sum(raw_strict_flags), "/", length(gamma_results)
      )
    )
    scice_message(
      paste(
        worker_id, ":   raw-relaxed-valid gammas (raw any-hit + |raw_median-target|<=1):",
        sum(raw_relaxed_flags), "/", length(gamma_results)
      )
    )
    scice_message(
      paste(
        worker_id, ":   gammas passing soft raw-cluster guard:",
        sum(soft_raw_guard_flags), "/", length(gamma_results)
      )
    )
    scice_message(
      paste(
        worker_id, ":   gammas passing hard raw-cluster guard:",
        sum(hard_raw_guard_flags), "/", length(gamma_results)
      )
    )
    admitted_count <- sum((strict_flags | relaxed_flags) & hard_raw_guard_flags)
    scice_message(
      paste(
        worker_id, ":   raw-guarded admitted gammas (strict or relaxed):",
        admitted_count, "/", length(gamma_results)
      )
    )
    if (admitted_count > 0) {
      admitted_hit_counts <- table(hit_counts[(strict_flags | relaxed_flags) & hard_raw_guard_flags])
      for (i in 1:length(admitted_hit_counts)) {
        hit_val <- names(admitted_hit_counts)[i]
        freq <- admitted_hit_counts[i]
        scice_message(paste(worker_id, ":   ", freq, "admitted gammas -> ", hit_val, " hit trials", sep = ""))
      }
    }
  }
  
  # Phase 2 admission decision:
  # 1) raw_strict_soft
  # 2) strict_soft
  # 3) relaxed_soft
  # 4) strict_hard
  # 5) relaxed_hard
  # 6) relaxed_unguarded
  # 7) raw_relaxed_soft
  # 8) raw_relaxed_hard
  # 9) raw_relaxed_unguarded
  admission_decision <- select_gamma_admission(
    strict_flags = strict_flags,
    relaxed_flags = relaxed_flags,
    soft_guard_flags = soft_raw_guard_flags,
    hard_guard_flags = hard_raw_guard_flags,
    raw_strict_flags = raw_strict_flags,
    raw_relaxed_flags = raw_relaxed_flags
  )
  valid_indices <- admission_decision$indices
  admission_mode <- admission_decision$mode
  pre_refine_candidate_count <- length(valid_indices)

  refined_candidates <- refine_gamma_candidates_by_raw_gap(
    valid_indices = valid_indices,
    admission_mode = admission_mode,
    gamma_results = gamma_results,
    target_clusters = target_clusters,
    min_cluster_size = min_cluster_size
  )
  valid_indices <- refined_candidates$indices
  admission_mode <- refined_candidates$mode
  selected_raw_gaps <- refined_candidates$raw_gaps
  best_raw_gap <- refined_candidates$best_raw_gap

  if (length(valid_indices) > 0L &&
      length(valid_indices) < pre_refine_candidate_count &&
      min_cluster_size > 1L &&
      any(is.finite(selected_raw_gaps))) {
    if (verbose) {
      scice_message(
        paste(
          worker_id, ": Phase 2 refinement - retaining",
          length(valid_indices), "candidate gammas with minimum raw-cluster gap",
          signif(best_raw_gap, 6), "to target", target_clusters
        )
      )
    }
  }
  
  if (length(valid_indices) == 0) {
    failure_order <- order(
      vapply(gamma_results, function(x) {
        if (is.null(x$median_gap) || is.na(x$median_gap)) {
          Inf
        } else {
          as.numeric(x$median_gap)
        }
      }, numeric(1)),
      vapply(gamma_results, function(x) {
        if (is.null(x$raw_median_gap) || is.na(x$raw_median_gap)) {
          Inf
        } else {
          as.numeric(x$raw_median_gap)
        }
      }, numeric(1)),
      vapply(gamma_results, function(x) as.numeric(x$gamma), numeric(1))
    )
    best_failure_idx <- if (length(failure_order) > 0L) failure_order[[1]] else NA_integer_
    failure_result <- if (is.na(best_failure_idx)) NULL else gamma_results[[best_failure_idx]]
    failure_refs <- lapply(gamma_results, function(x) {
      if (is.null(x$matrix_ref)) {
        return(NULL)
      }
      x$matrix_ref
    })
    failure_refs <- failure_refs[!vapply(failure_refs, is.null, logical(1))]
    if (length(failure_refs) > 0L) {
      release_cluster_matrix_refs(failure_refs)
    }
    if (verbose) {
      scice_message(
        paste(
          worker_id,
          ": ERROR - No gammas satisfied final-count admission or bounded raw-count fallback admission for target final merged cluster count",
          target_clusters
        )
      )
    }
    return(list(
      success = FALSE,
      failure_reason = "optimization_admission_failed",
      gamma = if (!is.null(failure_result)) as.numeric(failure_result$gamma) else NA_real_,
      labels = NULL,
      ic_median = NA_real_,
      ic_bootstrap = NULL,
      best_labels = NULL,
      effective_cluster_median = if (!is.null(failure_result)) as.numeric(failure_result$median_effective_clusters) else NA_real_,
      raw_cluster_median = if (!is.null(failure_result)) as.numeric(failure_result$raw_cluster_median) else NA_real_,
      final_cluster_median = if (!is.null(failure_result)) as.numeric(failure_result$final_cluster_median) else NA_real_,
      admission_mode = "optimization_admission_failed",
      best_labels_raw_cluster_count = NA_integer_,
      best_labels_final_cluster_count = NA_integer_,
      n_iterations = as.integer(n_iterations),
      k = as.integer(n_iterations)
    ))
  }
  
  if (verbose) {
    admission_messages <- c(
      raw_strict_soft = "soft-guarded raw-count exact family exists; using raw_strict_soft candidate set before effective families.",
      strict_soft = "no soft-guarded raw-count exact family; using bounded strict effective-match candidate set.",
      relaxed_soft = "no soft-guarded strict matches; using soft-guarded relaxed candidate set.",
      strict_hard = "no soft-guarded candidates; using hard-guarded strict candidate set.",
      relaxed_hard = "no soft-guarded strict candidates; using hard-guarded relaxed candidate set.",
      relaxed_unguarded = "no raw-guarded candidates; falling back to unguarded relaxed candidate set.",
      raw_relaxed_soft = "no effective candidates; using soft-guarded raw-count relaxed fallback candidate set.",
      raw_relaxed_hard = "no strict effective candidates; using hard-guarded raw-count relaxed fallback candidate set.",
      raw_relaxed_unguarded = "no raw-guarded effective candidates; falling back to unguarded raw-count relaxed candidate set."
    )
    scice_message(paste(worker_id, ": Phase 2 decision -", admission_messages[[admission_mode]]))
    scice_message(
      paste(
        worker_id, ": Found", length(valid_indices),
        "gammas selected under", admission_mode, "admission for", target_clusters, "final merged clusters"
      )
    )
    scice_message(paste(worker_id, ": Phase 3 - IC scores already computed during Phase 1 for valid gammas (all trials per admitted gamma)"))
  }
  
  valid_results <- gamma_results[valid_indices]
  exact_hit_gamma_flags <- vapply(
    valid_results,
    function(x) {
      !is.null(x$hit_count) && is.finite(x$hit_count) && as.integer(x$hit_count) > 0L
    },
    logical(1)
  )
  if (any(exact_hit_gamma_flags)) {
    if (verbose && sum(exact_hit_gamma_flags) < length(valid_results)) {
      scice_message(
        paste(
          worker_id, ": Phase 2 refinement - prioritizing",
          sum(exact_hit_gamma_flags), "gamma(s) with at least one exact final-hit trial over",
          length(valid_results) - sum(exact_hit_gamma_flags), "near-hit/raw-fallback gamma(s)"
        )
      )
    }
    valid_indices <- valid_indices[exact_hit_gamma_flags]
    valid_results <- valid_results[exact_hit_gamma_flags]
  }
  prefer_exact_hits <- any(exact_hit_gamma_flags)
  discarded_refs <- lapply(gamma_results[-valid_indices], function(x) {
    if (is.null(x$matrix_ref)) {
      return(NULL)
    }
    x$matrix_ref
  })
  discarded_refs <- discarded_refs[!vapply(discarded_refs, is.null, logical(1))]
  if (length(discarded_refs) > 0L) {
    release_cluster_matrix_refs(discarded_refs)
  }
  gamma_sequence <- vapply(valid_results, function(x) x$gamma, numeric(1))
  ic_scores <- vapply(valid_results, function(x) x$ic, numeric(1))
  clustering_refs <- lapply(valid_results, function(x) x$matrix_ref)
  effective_cluster_medians <- vapply(
    valid_results,
    function(x) as.numeric(x$median_effective_clusters),
    numeric(1)
  )
  final_cluster_medians <- vapply(
    valid_results,
    function(x) as.numeric(x$final_cluster_median),
    numeric(1)
  )
  raw_cluster_medians <- vapply(
    valid_results,
    function(x) as.numeric(x$raw_cluster_median),
    numeric(1)
  )
  rm(gamma_results, valid_results)
  
  if (verbose) {
    scice_message(paste(worker_id, ": IC score range: [", round(min(ic_scores), 4), ", ", round(max(ic_scores), 4), "]", sep = ""))
    perfect_scores <- sum(ic_scores == 1)
    if (perfect_scores > 0) {
      scice_message(paste(worker_id, ": Found", perfect_scores, "perfect IC scores (= 1.0)"))
    }
  }
  
  # Find the best gamma
  best_index <- which(ic_scores == 1)[1]
  if (is.na(best_index)) {
    best_index <- which.min(ic_scores)
  }
  
  best_gamma <- gamma_sequence[best_index]
  best_ref <- clustering_refs[[best_index]]
  k <- n_iterations
  
  if (verbose) {
    scice_message(paste(worker_id, ": Best gamma:", signif(best_gamma, 6), "with IC score:", round(ic_scores[best_index], 4)))
  }
  
  # If no perfect IC found, iterate to improve
  if (ic_scores[best_index] != 1 && length(gamma_sequence) > 1) {
    if (verbose) {
      scice_message(paste(worker_id, ": Phase 4 - Iterative improvement (IC =", round(ic_scores[best_index], 4), "< 1.0)"))
      scice_message(paste(worker_id, ": Starting with", length(gamma_sequence), "gamma values"))
      iterative_start <- Sys.time()
    }
    
    current_refs <- clustering_refs
    current_gammas <- gamma_sequence
    current_ic <- ic_scores
    
    # Track stability using a matrix for efficiency
    ic_history <- matrix(rep(current_ic, 10), nrow = length(current_ic))
    
    iteration_count <- 0
    converged <- FALSE
    iteration_workers <- cap_workers_by_memory(
      nested_workers,
      estimate_trial_matrix_bytes(n_vertices, n_trials, 1L) * 2,
      runtime_context
    )
    while (k < max_iterations) {
      k <- k + delta_n
      iteration_count <- iteration_count + 1
      
      if (verbose) {
        scice_message(paste(worker_id, ": Iteration", iteration_count, "(k =", k, ") - Refining", length(current_gammas), "gamma values"))
        iter_start <- Sys.time()
      }
      
      # Update clustering results in parallel
      new_results <- cross_platform_mclapply(seq_along(current_gammas), function(i) {
        gamma_val <- current_gammas[i]
        current_matrix <- load_cluster_matrix(current_refs[[i]])
        
        # Use previous results as initialization
        # Set deterministic seed for this iteration if base seed provided
        if (!is.null(seed)) {
          iter_seed <- cluster_seed + k * 100 + i
          set.seed(iter_seed)
        }
        
        new_clustering <- matrix(0L, nrow = n_vertices, ncol = n_trials)
        for (trial_idx in seq_len(n_trials)) {
          init_membership <- current_matrix[, sample.int(ncol(current_matrix), 1)]
          new_clustering[, trial_idx] <- run_leiden_trial(gamma_val, delta_n, init_membership)
          heartbeat(function() {
            paste(
              "phase4 running - k iteration", k,
              "- gamma", i, "/", length(current_gammas),
              "- trial", trial_idx, "/", n_trials
            )
          })
        }
        final_clusters_vec <- if (min_cluster_size > 1L) {
          vapply(
            seq_len(n_trials),
            function(trial_idx) {
              merged_labels <- merge_small_clusters_to_neighbors(
                labels = new_clustering[, trial_idx],
                snn_graph = snn_graph,
                min_cluster_size = min_cluster_size
              )
              as.integer(length(unique(merged_labels)))
            },
            integer(1)
          )
        } else {
          vapply(
            seq_len(n_trials),
            function(trial_idx) as.integer(length(unique(new_clustering[, trial_idx]))),
            integer(1)
          )
        }

        # Calculate IC score
        extracted <- extract_clustering_array(new_clustering)
        ic_result <- calculate_ic_from_extracted(extracted)
        matrix_ref <- store_cluster_matrix(
          new_clustering,
          runtime_context = runtime_context,
          prefix = sprintf("k%d_iter%d_g%03d", target_clusters, k, i)
        )
        rm(current_matrix, extracted, new_clustering)
        
        list(
          matrix_ref = matrix_ref,
          ic = 1 / ic_result,
          exact_hit_count = as.integer(sum(final_clusters_vec == target_clusters))
        )
      }, mc.cores = iteration_workers)
      
      # Extract results
      new_refs <- lapply(new_results, function(x) x$matrix_ref)
      new_ic <- vapply(new_results, function(x) x$ic, numeric(1))
      new_exact_hit_counts <- vapply(new_results, function(x) {
        if (is.null(x$exact_hit_count)) {
          return(0L)
        }
        as.integer(x$exact_hit_count)
      }, integer(1))
      candidate_refs <- new_refs
      candidate_ic <- new_ic
      candidate_gammas <- current_gammas
      candidate_ic_history <- cbind(ic_history[, -1], new_ic)
      candidate_exact_hit_flags <- new_exact_hit_counts > 0L

      if (prefer_exact_hits) {
        if (!any(candidate_exact_hit_flags)) {
          release_cluster_matrix_refs(new_refs)
          best_index <- which(current_ic == 1)[1]
          if (is.na(best_index)) {
            best_index <- which.min(current_ic)
          }
          best_gamma <- current_gammas[best_index]
          best_ref <- current_refs[[best_index]]
          if (length(current_refs) > 1) {
            release_cluster_matrix_refs(current_refs[-best_index])
          }
          if (verbose) {
            scice_message(
              paste(
                worker_id, ": CONVERGED - Iterative refinement lost all exact final-hit trials;",
                "falling back to previous exact-hit candidate at gamma =",
                signif(best_gamma, 6)
              )
            )
          }
          converged <- TRUE
          break
        }
        if (!all(candidate_exact_hit_flags)) {
          if (verbose) {
            scice_message(
              paste(
                worker_id, ": Iteration", iteration_count, "- retaining",
                sum(candidate_exact_hit_flags), "gamma(s) that still contain exact final-hit trials"
              )
            )
          }
          release_cluster_matrix_refs(candidate_refs[!candidate_exact_hit_flags])
          candidate_refs <- candidate_refs[candidate_exact_hit_flags]
          candidate_ic <- candidate_ic[candidate_exact_hit_flags]
          candidate_gammas <- candidate_gammas[candidate_exact_hit_flags]
          candidate_ic_history <- candidate_ic_history[candidate_exact_hit_flags, , drop = FALSE]
        }
      }
      release_cluster_matrix_refs(current_refs)
      
      if (verbose) {
        iter_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
        scice_message(paste(worker_id, ": Iteration", iteration_count, "completed in", round(iter_time, 3), "seconds"))
        scice_message(paste(worker_id, ": New IC range: [", round(min(candidate_ic), 4), ", ", round(max(candidate_ic), 4), "]", sep = ""))
      }
      
      # Update IC history
      ic_history <- candidate_ic_history
      
      # Check for stability and convergence
      stable_indices <- apply(ic_history, 1, function(row) length(unique(row)) == 1)
      perfect_indices <- which(candidate_ic == 1)
      
      if (verbose) {
        scice_message(paste(worker_id, ": Stable gammas:", sum(stable_indices), "/", length(stable_indices)))
        scice_message(paste(worker_id, ": Perfect IC scores:", length(perfect_indices)))
      }
      
      # Selection criteria
      if (length(perfect_indices) > 0) {
        best_index <- perfect_indices[1]
        best_gamma <- candidate_gammas[best_index]
        best_ref <- candidate_refs[[best_index]]
        if (length(candidate_refs) > 1) {
          release_cluster_matrix_refs(candidate_refs[-best_index])
        }
        if (verbose) {
          scice_message(paste(worker_id, ": CONVERGED - Found perfect IC score at gamma =", signif(best_gamma, 6)))
        }
        converged <- TRUE
        break
      } else if (all(stable_indices)) {
        best_index <- which.min(candidate_ic)
        best_gamma <- candidate_gammas[best_index]
        best_ref <- candidate_refs[[best_index]]
        if (length(candidate_refs) > 1) {
          release_cluster_matrix_refs(candidate_refs[-best_index])
        }
        if (verbose) {
          scice_message(paste(worker_id, ": CONVERGED - All gammas stable, best IC =", round(candidate_ic[best_index], 4)))
        }
        converged <- TRUE
        break
      } else {
        # Continue with best performing gammas
        keep_indices <- (candidate_ic <= stats::quantile(candidate_ic, 0.5)) | stable_indices
        keep_indices[which.min(candidate_ic)] <- TRUE
        
        if (sum(keep_indices) == 1) {
          best_index <- which(keep_indices)
          best_gamma <- candidate_gammas[best_index]
          best_ref <- candidate_refs[[best_index]]
          release_cluster_matrix_refs(candidate_refs[!keep_indices])
          if (verbose) {
            scice_message(paste(worker_id, ": CONVERGED - Single gamma remaining, IC =", round(candidate_ic[best_index], 4)))
          }
          converged <- TRUE
          break
        }
        
        if (verbose) {
          scice_message(paste(worker_id, ": Continuing with", sum(keep_indices), "best gammas"))
        }

        release_cluster_matrix_refs(candidate_refs[!keep_indices])
        
        current_gammas <- candidate_gammas[keep_indices]
        current_refs <- candidate_refs[keep_indices]
        ic_history <- ic_history[keep_indices, , drop = FALSE]
        current_ic <- candidate_ic[keep_indices]
      }
    }

    if (!converged) {
      best_index <- which.min(current_ic)
      best_gamma <- current_gammas[best_index]
      best_ref <- current_refs[[best_index]]
      if (length(current_refs) > 1) {
        release_cluster_matrix_refs(current_refs[-best_index])
      }
    }
    
    if (verbose) {
      iterative_time <- as.numeric(difftime(Sys.time(), iterative_start, units = "secs"))
      scice_message(paste(worker_id, ": Phase 4 completed in", round(iterative_time, 3), "seconds after", iteration_count, "iterations"))
    }
  } else {
    if (length(clustering_refs) > 1) {
      release_cluster_matrix_refs(clustering_refs[-best_index])
    }
  }

  best_gamma_diag_index <- which.min(abs(gamma_sequence - best_gamma))
  selected_effective_cluster_median <- effective_cluster_medians[[best_gamma_diag_index]]
  selected_final_cluster_median <- final_cluster_medians[[best_gamma_diag_index]]
  selected_raw_cluster_median <- raw_cluster_medians[[best_gamma_diag_index]]
  finalized_result <- finalize_selected_clustering(
    matrix_ref = best_ref,
    gamma = best_gamma,
    effective_cluster_median = selected_effective_cluster_median,
    raw_cluster_median = selected_raw_cluster_median,
    final_cluster_median = selected_final_cluster_median,
    admission_mode = admission_mode,
    cluster_seed = cluster_seed,
    n_bootstrap = n_bootstrap,
    n_workers = nested_workers,
    snn_graph = snn_graph,
    target_clusters = target_clusters,
    min_cluster_size = min_cluster_size,
    verbose = verbose,
    worker_id = worker_id,
    runtime_context = runtime_context
  )

  finalized_result$n_iterations <- k
  finalized_result$k <- k
  finalized_result$success <- TRUE
  finalized_result
}

evaluate_fixed_resolution <- function(igraph_obj, resolution, objective_function,
                                      n_trials, n_bootstrap, seed = NULL, beta,
                                      n_iterations, n_workers, snn_graph = NULL,
                                      min_cluster_size = 1L, verbose = FALSE,
                                      worker_id = "RESOLUTION",
                                      in_parallel_context = FALSE,
                                      runtime_context = NULL) {
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  if (min_cluster_size > 1L && is.null(snn_graph)) {
    stop("snn_graph must be provided when min_cluster_size > 1.")
  }

  run_leiden_trial <- function(iterations) {
    leiden_clustering(
      igraph_obj, resolution, objective_function, iterations, beta, initial_membership = NULL
    )
  }

  cluster_seed <- derive_manual_resolution_seed(seed, resolution)
  if (!is.null(cluster_seed)) {
    set.seed(cluster_seed)
    if (verbose) {
      scice_message(paste(worker_id, ": Set deterministic seed:", cluster_seed))
    }
  }

  if (verbose) {
    scice_message(paste(worker_id, ": Fixed-resolution parameters:"))
    scice_message(paste(worker_id, ":   Resolution:", signif(resolution, 6)))
    scice_message(paste(worker_id, ":   Trials:", n_trials))
    scice_message(paste(worker_id, ":   Bootstrap iterations:", n_bootstrap))
    scice_message(paste(worker_id, ":   Beta:", beta))
    scice_message(paste(worker_id, ":   Leiden iterations:", n_iterations))
    if (min_cluster_size > 1L) {
      scice_message(
        paste(
          worker_id, ":   Counting uses effective clusters (size >=", min_cluster_size,
          "); final merge applied only on best_labels"
        )
      )
    }
  }

  n_vertices <- igraph::vcount(igraph_obj)
  heartbeat <- create_heartbeat_logger(verbose = verbose, context = worker_id)
  estimated_phase1_bytes <- estimate_trial_matrix_bytes(n_vertices, n_trials, 1L)
  if (should_enable_spill(runtime_context, estimated_phase1_bytes)) {
    activate_runtime_spill(runtime_context, estimated_bytes = estimated_phase1_bytes)
  }

  if (verbose) {
    scice_message(paste(worker_id, ": Phase 1 - Evaluating fixed resolution with", n_trials, "trials"))
    phase1_start <- Sys.time()
    if (in_parallel_context) {
      scice_message(paste(worker_id, ": Running in parallel context with worker budget", n_workers))
    }
  }

  trial_workers <- max(1L, as.integer(n_workers))
  if (!in_parallel_context) {
    trial_workers <- min(trial_workers, max(1L, as.integer(n_trials)))
  }
  trial_workers <- cap_workers_by_memory(
    trial_workers,
    estimate_trial_matrix_bytes(n_vertices, 1L, 1L),
    runtime_context
  )

  phase1_log_every <- max(1L, as.integer(floor(n_trials / 5)))
  should_log_trial_step <- function(idx) {
    idx == 1L || idx == n_trials || (idx %% phase1_log_every) == 0L
  }
  run_single_trial <- function(trial_idx) {
    if (!is.null(cluster_seed)) {
      trial_seed <- cluster_seed + trial_idx
      set.seed(trial_seed)
    }
    labels <- run_leiden_trial(n_iterations)
    heartbeat(function() {
      paste(
        "phase1 running - fixed resolution", signif(resolution, 6),
        "- trial", trial_idx, "/", n_trials
      )
    })
    if (verbose && should_log_trial_step(trial_idx)) {
      scice_message(
        paste(
          worker_id, ": Phase 1 progress trial", trial_idx, "/", n_trials,
          "(resolution =", signif(resolution, 6), ")"
        )
      )
    }
    labels
  }

  if (trial_workers == 1) {
    trial_results <- vector("list", n_trials)
    for (trial_idx in seq_len(n_trials)) {
      trial_results[[trial_idx]] <- run_single_trial(trial_idx)
    }
  } else {
    trial_results <- cross_platform_mclapply(
      seq_len(n_trials),
      run_single_trial,
      mc.cores = trial_workers
    )
  }

  cluster_matrix <- do.call(cbind, trial_results)
  effective_clusters_vec <- vapply(
    seq_len(n_trials),
    function(trial_idx) {
      count_effective_clusters(cluster_matrix[, trial_idx], min_cluster_size = min_cluster_size)
    },
    integer(1)
  )
  raw_clusters_vec <- vapply(
    seq_len(n_trials),
    function(trial_idx) {
      as.integer(length(unique(cluster_matrix[, trial_idx])))
    },
    integer(1)
  )
  final_clusters_vec <- if (min_cluster_size > 1L) {
    vapply(
      seq_len(n_trials),
      function(trial_idx) {
        merged_labels <- merge_small_clusters_to_neighbors(
          labels = cluster_matrix[, trial_idx],
          snn_graph = snn_graph,
          min_cluster_size = min_cluster_size
        )
        as.integer(length(unique(merged_labels)))
      },
      integer(1)
    )
  } else {
    raw_clusters_vec
  }
  median_effective_clusters <- stats::median(effective_clusters_vec)
  raw_cluster_median <- stats::median(raw_clusters_vec)
  final_cluster_median <- stats::median(final_clusters_vec)

  extracted <- extract_clustering_array(cluster_matrix)
  ic_result <- calculate_ic_from_extracted(extracted)
  phase1_ic <- 1 / ic_result
  rm(extracted)

  matrix_ref <- store_cluster_matrix(
    cluster_matrix,
    runtime_context = runtime_context,
    prefix = sprintf("fixed_resolution_%s", gsub("[^0-9A-Za-z]+", "_", format(resolution, scientific = FALSE)))
  )
  rm(cluster_matrix)

  if (verbose) {
    phase1_time <- as.numeric(difftime(Sys.time(), phase1_start, units = "secs"))
    scice_message(paste(worker_id, ": Phase 1 completed in", round(phase1_time, 3), "seconds"))
    scice_message(
      paste(
        worker_id, ": Phase 1 diagnostics - resolution =", signif(resolution, 6),
        "- median_effective =", signif(median_effective_clusters, 6),
        "- median_raw =", signif(raw_cluster_median, 6),
        "- median_final =", signif(final_cluster_median, 6),
        "- IC (all trials) =", round(phase1_ic, 4)
      )
    )
  }

  finalized_result <- finalize_selected_clustering(
    matrix_ref = matrix_ref,
    gamma = resolution,
    effective_cluster_median = median_effective_clusters,
    raw_cluster_median = raw_cluster_median,
    final_cluster_median = final_cluster_median,
    admission_mode = "manual_resolution",
    cluster_seed = cluster_seed,
    n_bootstrap = n_bootstrap,
    n_workers = trial_workers,
    snn_graph = snn_graph,
    min_cluster_size = min_cluster_size,
    verbose = verbose,
    worker_id = worker_id,
    runtime_context = runtime_context
  )

  finalized_result$phase1_ic <- phase1_ic
  finalized_result$n_iterations <- n_iterations
  finalized_result$k <- n_iterations
  finalized_result
}
