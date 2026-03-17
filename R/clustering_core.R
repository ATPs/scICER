#' @importFrom stats median setNames
#' @importFrom parallel detectCores mclapply
#' @importFrom data.table data.table rbindlist
NULL

# Global clustering cache environment for optimization
clustering_cache_env <- new.env(parent = emptyenv())

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

#' Cached version of leiden clustering to avoid redundant computations
#'
#' @description
#' Reuses previously computed Leiden results for identical parameter sets.
#'
#' @details
#' Cache keys are built from resolution, objective function, iteration settings,
#' beta, and a caller-provided suffix. Results are stored in a process-local
#' environment, which speeds up repeated evaluations during resolution search and
#' filtering steps.
#'
#' \code{initial_membership} is intentionally excluded from the default cache key
#' so broad reuse remains possible unless callers differentiate with
#' \code{cache_key_suffix}.
#'
#' @param igraph_obj igraph object to cluster
#' @param resolution Resolution parameter for clustering
#' @param objective_function Objective function ("modularity" or "CPM")
#' @param n_iterations Number of iterations
#' @param beta Beta parameter
#' @param initial_membership Initial cluster membership (optional)
#' @param use_cache Whether to use caching (default: TRUE)
#' @param cache_key_suffix Additional suffix for cache key (default: "")
#' @param snn_graph Seurat SNN graph matrix (reserved for compatibility)
#' @param min_cluster_size Minimum required cells per cluster
#' @return Vector of cluster assignments (0-based)
#' @keywords internal
cached_leiden_clustering <- function(igraph_obj, resolution, objective_function, 
                                   n_iterations, beta, initial_membership = NULL,
                                   use_cache = TRUE, cache_key_suffix = "",
                                   snn_graph = NULL, min_cluster_size = 1L) {
  min_cluster_size <- max(1L, as.integer(min_cluster_size))
  
  if (!use_cache) {
    return(leiden_clustering(
      igraph_obj, resolution, objective_function, n_iterations, beta, initial_membership
    ))
  }
  
  # Create cache key based on parameters (excluding initial_membership for broader reuse)
  cache_key <- paste(
    "r", round(resolution, 8),
    "obj", objective_function, 
    "iter", n_iterations,
    "beta", round(beta, 4),
    "suffix", cache_key_suffix,
    sep = "_"
  )
  
  # Check if result exists in cache
  if (exists(cache_key, envir = clustering_cache_env)) {
    return(clustering_cache_env[[cache_key]])
  }
  
  # Compute result and cache it
  result <- leiden_clustering(
    igraph_obj, resolution, objective_function, n_iterations, beta, initial_membership
  )
  clustering_cache_env[[cache_key]] <- result
  
  return(result)
}

#' Clear the global clustering cache
#'
#' @description
#' Deletes all cached Leiden results from the current R process.
#'
#' @details
#' Called at the start of top-level scICER runs to avoid stale reuse across
#' unrelated analyses. Keeping this explicit makes memory behavior predictable.
#'
#' @keywords internal
clear_clustering_cache <- function() {
  rm(list = ls(envir = clustering_cache_env), envir = clustering_cache_env)
}

#' Get cache statistics
#'
#' @description
#' Reports simple cache usage metadata for diagnostics/logging.
#'
#' @details
#' Currently returns only the number of cache entries, which is sufficient for
#' lightweight runtime reporting without expensive memory introspection.
#'
#' @keywords internal  
get_cache_stats <- function() {
  cache_size <- length(ls(envir = clustering_cache_env))
  return(list(cache_entries = cache_size))
}

#' Cross-platform parallel lapply wrapper
#'
#' @description
#' Unified wrapper for parallel iteration on Unix and Windows.
#'
#' @details
#' Uses \code{parallel::mclapply()} on Unix-like systems and falls back to
#' sequential \code{lapply()} on Windows, where fork-based multiprocessing is
#' unavailable. This keeps calling code portable and avoids repeated platform
#' conditionals.
#'
#' @param X Vector/list to iterate over
#' @param FUN Function to apply
#' @param mc.cores Number of cores (ignored on Windows)
#' @param mc.preschedule Whether to preschedule chunks for fork workers (Unix only)
#' @param ... Additional arguments to FUN
#' @return List of results
#' @keywords internal
cross_platform_mclapply <- function(X, FUN, mc.cores = 1, mc.preschedule = TRUE, ...) {
  if (.Platform$OS.type == "windows") {
    # Use regular lapply on Windows
    return(lapply(X, FUN, ...))
  } else {
    # Use mclapply on Unix-like systems
    return(parallel::mclapply(X, FUN, mc.cores = mc.cores, mc.preschedule = mc.preschedule, ...))
  }
}

#' Resolve heartbeat interval for verbose progress logging
#'
#' @description
#' Returns heartbeat interval in seconds for periodic liveness logs.
#'
#' @param default_seconds Fallback interval in seconds
#' @return Positive numeric scalar
#' @keywords internal
get_heartbeat_interval_seconds <- function(default_seconds = 60) {
  configured <- getOption("scICER.internal_heartbeat_seconds", default_seconds)
  if (!is.numeric(configured) || length(configured) != 1 || !is.finite(configured) || configured <= 0) {
    return(as.double(default_seconds))
  }
  as.double(configured)
}

#' Create a throttled heartbeat logger
#'
#' @description
#' Builds a closure that emits at most one progress message per interval.
#'
#' @param verbose Whether logging is enabled
#' @param context Optional message prefix
#' @param interval_seconds Optional heartbeat interval override
#' @return Function accepting a message string or message-producing function
#' @keywords internal
create_heartbeat_logger <- function(verbose, context = "", interval_seconds = NULL) {
  enabled <- isTRUE(verbose)
  interval <- if (is.null(interval_seconds)) {
    get_heartbeat_interval_seconds()
  } else {
    as.double(interval_seconds)
  }
  if (!is.finite(interval) || interval <= 0) {
    interval <- 60
  }

  context <- if (is.null(context) || !nzchar(as.character(context)[1])) "" else as.character(context)[1]
  state <- new.env(parent = emptyenv())
  state$last_emit <- Sys.time()

  function(message_text) {
    if (!enabled) {
      return(invisible(FALSE))
    }

    now <- Sys.time()
    elapsed <- as.numeric(difftime(now, state$last_emit, units = "secs"))
    if (!is.finite(elapsed) || elapsed < interval) {
      return(invisible(FALSE))
    }
    state$last_emit <- now

    if (is.function(message_text)) {
      message_text <- message_text()
    }
    if (!is.character(message_text) || length(message_text) == 0 || is.na(message_text[1])) {
      message_text <- "progress update"
    }

    prefix <- if (nzchar(context)) paste0(context, " ") else ""
    scice_message(paste0(prefix, "HEARTBEAT (", as.integer(round(interval)), "s): ", message_text[1]))
    invisible(TRUE)
  }
}

#' Estimate clustering matrix memory footprint
#'
#' @description
#' Approximates in-memory bytes for integer clustering matrices used in trial loops.
#'
#' @details
#' The estimate follows the conservative formula used for spill and worker limits:
#' \code{n_cells * n_trials * n_gamma * 4 bytes}.
#'
#' @param n_cells Number of cells (rows)
#' @param n_trials Number of trials (columns)
#' @param n_gamma Number of gamma values represented
#' @return Estimated bytes as numeric scalar
#' @keywords internal
estimate_trial_matrix_bytes <- function(n_cells, n_trials, n_gamma = 1L) {
  as.double(n_cells) * as.double(max(1L, as.integer(n_trials))) * as.double(max(1L, as.integer(n_gamma))) * 4
}

#' Detect runtime memory budget for worker capping
#'
#' @description
#' Resolves an approximate memory budget in bytes used to cap nested workers.
#'
#' @details
#' Priority order:
#' \enumerate{
#'   \item \code{options(scICER.internal_memory_budget_bytes = ...)};
#'   \item Linux \code{/proc/meminfo} \code{MemAvailable} with 70\% safety factor;
#'   \item fallback default (4 GiB).
#' }
#'
#' @param default_bytes Fallback memory budget in bytes
#' @return Numeric byte budget
#' @keywords internal
detect_memory_budget_bytes <- function(default_bytes = 4 * 1024^3) {
  option_budget <- getOption("scICER.internal_memory_budget_bytes", NULL)
  if (!is.null(option_budget) && is.numeric(option_budget) && is.finite(option_budget) && option_budget > 0) {
    return(as.double(option_budget))
  }

  if (.Platform$OS.type != "windows" && file.exists("/proc/meminfo")) {
    meminfo <- tryCatch(readLines("/proc/meminfo", warn = FALSE), error = function(e) character(0))
    mem_line <- meminfo[grepl("^MemAvailable:\\s+[0-9]+\\s+kB$", meminfo)]
    if (length(mem_line) > 0) {
      mem_kb <- as.numeric(gsub("^MemAvailable:\\s+([0-9]+)\\s+kB$", "\\1", mem_line[1]))
      if (is.finite(mem_kb) && mem_kb > 0) {
        return(as.double(mem_kb) * 1024 * 0.7)
      }
    }
  }

  as.double(default_bytes)
}

#' Cap nested workers by estimated per-task memory
#'
#' @description
#' Reduces requested worker count when estimated concurrent memory exceeds budget.
#'
#' @param requested_workers Requested workers before memory cap
#' @param bytes_per_task Estimated bytes per parallel task
#' @param runtime_context Optional runtime context with memory budget metadata
#' @return Worker count after capping (integer >= 1)
#' @keywords internal
cap_workers_by_memory <- function(requested_workers, bytes_per_task, runtime_context = NULL) {
  workers <- max(1L, as.integer(requested_workers))
  if (!is.finite(bytes_per_task) || bytes_per_task <= 0) {
    return(workers)
  }

  memory_budget <- if (!is.null(runtime_context) &&
                       !is.null(runtime_context$memory_budget_bytes) &&
                       is.finite(runtime_context$memory_budget_bytes) &&
                       runtime_context$memory_budget_bytes > 0) {
    runtime_context$memory_budget_bytes
  } else {
    detect_memory_budget_bytes()
  }

  overhead_factor <- getOption("scICER.internal_worker_memory_overhead", 1.3)
  if (!is.numeric(overhead_factor) || !is.finite(overhead_factor) || overhead_factor <= 0) {
    overhead_factor <- 1.3
  }

  concurrent_cap <- floor(memory_budget / (as.double(bytes_per_task) * overhead_factor))
  if (!is.finite(concurrent_cap) || concurrent_cap < 1) {
    concurrent_cap <- 1
  }
  concurrent_cap <- min(concurrent_cap, .Machine$integer.max)
  concurrent_cap <- as.integer(concurrent_cap)
  min(workers, concurrent_cap)
}

#' Create runtime context for spill and memory controls
#'
#' @description
#' Builds mutable internal state shared across scICER internal functions.
#'
#' @details
#' The context is internal-only and does not modify public API surface.
#' Spill defaults:
#' \itemize{
#'   \item mode: \code{auto}
#'   \item format: \code{qs}
#'   \item temp directory root: current working directory
#' }
#'
#' Internal options for tests:
#' \itemize{
#'   \item \code{scICER.internal_force_spill}
#'   \item \code{scICER.internal_spill_threshold_bytes}
#'   \item \code{scICER.internal_spill_dir}
#'   \item \code{scICER.internal_memory_budget_bytes}
#' }
#'
#' @return Environment with runtime fields
#' @keywords internal
create_runtime_context <- function() {
  context <- new.env(parent = emptyenv())

  spill_threshold <- getOption("scICER.internal_spill_threshold_bytes", 2 * 1024^3)
  if (!is.numeric(spill_threshold) || !is.finite(spill_threshold) || spill_threshold <= 0) {
    spill_threshold <- 2 * 1024^3
  }

  configured_dir <- getOption("scICER.internal_spill_dir", NULL)
  if (!is.null(configured_dir) && (!is.character(configured_dir) || length(configured_dir) != 1 || !nzchar(configured_dir))) {
    configured_dir <- NULL
  }

  context$spill_mode <- "auto"
  context$spill_format <- "qs"
  context$spill_threshold_bytes <- as.double(spill_threshold)
  context$force_spill <- isTRUE(getOption("scICER.internal_force_spill", FALSE))
  context$spill_dir <- if (is.null(configured_dir)) {
    file.path(
      getwd(),
      sprintf("scicer_tmp_%d_%s", Sys.getpid(), format(Sys.time(), "%Y%m%d%H%M%S"))
    )
  } else {
    configured_dir
  }
  context$spill_active <- FALSE
  context$spill_announced <- FALSE
  context$memory_budget_bytes <- detect_memory_budget_bytes()
  context
}

#' Decide whether spill should be enabled
#'
#' @description
#' Evaluates auto/force spill policy for the current memory estimate.
#'
#' @param runtime_context Runtime context environment
#' @param estimated_bytes Estimated bytes for pending matrix allocations
#' @return Logical scalar
#' @keywords internal
should_enable_spill <- function(runtime_context, estimated_bytes) {
  if (is.null(runtime_context)) {
    return(FALSE)
  }
  if (isTRUE(runtime_context$force_spill)) {
    return(TRUE)
  }
  is.finite(estimated_bytes) && estimated_bytes >= runtime_context$spill_threshold_bytes
}

#' Activate spill mode and create temp directory
#'
#' @description
#' Lazily initializes spill directory and validates \pkg{qs} availability.
#'
#' @param runtime_context Runtime context environment
#' @param estimated_bytes Optional estimate used for diagnostics
#' @return Invisibly returns TRUE when spill is active
#' @keywords internal
activate_runtime_spill <- function(runtime_context, estimated_bytes = NA_real_) {
  if (is.null(runtime_context)) {
    return(invisible(FALSE))
  }

  if (isTRUE(runtime_context$spill_active)) {
    return(invisible(TRUE))
  }

  if (!requireNamespace("qs", quietly = TRUE)) {
    stop("Package 'qs' is required for scICER spill mode but is not available.")
  }

  dir.create(runtime_context$spill_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(runtime_context$spill_dir)) {
    stop("Failed to create scICER spill directory: ", runtime_context$spill_dir)
  }

  runtime_context$spill_active <- TRUE
  if (!isTRUE(runtime_context$spill_announced)) {
    message(sprintf("scICER temp spill dir: %s", runtime_context$spill_dir))
    runtime_context$spill_announced <- TRUE
  }

  if (is.finite(estimated_bytes)) {
    runtime_context$last_spill_estimate_bytes <- as.double(estimated_bytes)
  }

  invisible(TRUE)
}

#' Cleanup spill directory after run completion
#'
#' @description
#' Removes runtime spill directory recursively when spill mode was enabled.
#'
#' @param runtime_context Runtime context environment
#' @return Invisibly TRUE
#' @keywords internal
cleanup_runtime_spill <- function(runtime_context) {
  if (is.null(runtime_context)) {
    return(invisible(TRUE))
  }

  if (!(isTRUE(runtime_context$spill_active) || isTRUE(runtime_context$spill_announced))) {
    return(invisible(TRUE))
  }

  if (!is.character(runtime_context$spill_dir) || length(runtime_context$spill_dir) != 1 || !nzchar(runtime_context$spill_dir)) {
    warning("scICER spill cleanup skipped: invalid spill_dir.")
    return(invisible(FALSE))
  }

  if (dir.exists(runtime_context$spill_dir)) {
    unlink(runtime_context$spill_dir, recursive = TRUE, force = TRUE)
  }

  if (dir.exists(runtime_context$spill_dir)) {
    warning("scICER spill cleanup failed: ", runtime_context$spill_dir)
    return(invisible(FALSE))
  }

  message(sprintf("scICER temp spill dir cleaned: %s", runtime_context$spill_dir))
  invisible(TRUE)
}

#' Store clustering matrix in memory or spill file
#'
#' @description
#' Returns a lightweight reference object for later loading/release.
#'
#' @param cluster_matrix Matrix to store
#' @param runtime_context Runtime context environment
#' @param prefix File prefix used in spill mode
#' @return Reference list with \code{type} and payload
#' @keywords internal
store_cluster_matrix <- function(cluster_matrix, runtime_context, prefix = "cluster_matrix") {
  if (is.null(runtime_context) || !isTRUE(runtime_context$spill_active)) {
    return(list(type = "memory", value = cluster_matrix))
  }

  file_path <- tempfile(
    pattern = paste0(prefix, "_"),
    tmpdir = runtime_context$spill_dir,
    fileext = ".qs"
  )
  qs::qsave(cluster_matrix, file_path, preset = "fast")
  list(type = "spill", path = file_path)
}

#' Load clustering matrix from reference
#'
#' @description
#' Materializes a matrix previously returned by \code{store_cluster_matrix()}.
#'
#' @param matrix_ref Matrix reference list
#' @return Matrix
#' @keywords internal
load_cluster_matrix <- function(matrix_ref) {
  if (is.null(matrix_ref) || is.null(matrix_ref$type)) {
    stop("Invalid matrix reference.")
  }

  if (identical(matrix_ref$type, "memory")) {
    return(matrix_ref$value)
  }

  if (!identical(matrix_ref$type, "spill")) {
    stop("Unknown matrix reference type: ", matrix_ref$type)
  }

  if (is.null(matrix_ref$path) || !file.exists(matrix_ref$path)) {
    stop("Spill file does not exist: ", matrix_ref$path)
  }

  qs::qread(matrix_ref$path)
}

#' Release matrix reference resources
#'
#' @description
#' Deletes spill files when no longer needed. In-memory references are no-op.
#'
#' @param matrix_ref Matrix reference list
#' @return Invisibly TRUE
#' @keywords internal
release_cluster_matrix <- function(matrix_ref) {
  if (is.null(matrix_ref) || is.null(matrix_ref$type)) {
    return(invisible(TRUE))
  }

  if (identical(matrix_ref$type, "spill") &&
      is.character(matrix_ref$path) &&
      length(matrix_ref$path) == 1 &&
      file.exists(matrix_ref$path)) {
    unlink(matrix_ref$path, force = TRUE)
  }

  invisible(TRUE)
}

#' Release multiple matrix references
#'
#' @description
#' Convenience helper to release a list of references returned by
#' \code{store_cluster_matrix()}.
#'
#' @param matrix_refs List of matrix references
#' @return Invisibly TRUE
#' @keywords internal
release_cluster_matrix_refs <- function(matrix_refs) {
  if (length(matrix_refs) == 0) {
    return(invisible(TRUE))
  }
  invisible(lapply(matrix_refs, release_cluster_matrix))
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
                          in_parallel_context = FALSE, runtime_context = NULL) {
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
  results_dt <- data.table::data.table(
    cluster_number = integer(),
    gamma = numeric(),
    labels = list(),
    ic = numeric(),
    ic_vec = list(),
    best_labels = list(),
    effective_cluster_median = numeric(),
    raw_cluster_median = numeric(),
    admission_mode = character(),
    best_labels_raw_cluster_count = integer(),
    n_iter = integer(),
    mei = list(),
    k = integer(),
    excluded = logical(),
    exclusion_reason = character()
  )
  
  # Determine resolution search bounds
  if (objective_function == "modularity") {
    start_g <- -13
    end_g <- 20  # Increased for higher cluster numbers
  } else { # CPM
    start_g <- log(resolution_tolerance)
    if (start_g < -20) start_g <- -20
    end_g <- 20  # Increased for higher cluster numbers
  }
  
  # Binary search for resolution ranges
  if (verbose) {
    scice_message("CLUSTERING_MAIN: Starting binary search for resolution ranges...")
    scice_message(paste("CLUSTERING_MAIN: Objective function:", objective_function))
    scice_message(paste("CLUSTERING_MAIN: Search bounds: [", start_g, ", ", end_g, "]", sep = ""))
    scice_message(paste("CLUSTERING_MAIN: Target cluster range:", paste(cluster_range, collapse = ", ")))
    resolution_search_start <- Sys.time()
  }
  
  gamma_dict <- find_resolution_ranges(
    igraph_obj, cluster_range, start_g, end_g, objective_function,
    resolution_tolerance, n_workers, verbose, seed, snn_graph, min_cluster_size,
    in_parallel_context, runtime_context
  )
  
  if (verbose) {
    resolution_search_time <- as.numeric(difftime(Sys.time(), resolution_search_start, units = "secs"))
    scice_message(paste("CLUSTERING_MAIN: Resolution search completed in", round(resolution_search_time, 3), "seconds"))
    scice_message(paste("CLUSTERING_MAIN: Found resolution ranges for", length(gamma_dict), "cluster numbers"))
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
        admission_mode = NA_character_,
        best_labels_raw_cluster_count = NA_integer_,
        n_iter = NA_integer_,
        mei = list(NULL),
        k = NA_integer_,
        excluded = TRUE,
        exclusion_reason = x$reason
      ))
    }
    return(NULL)
  })
  
  # Handle excluded entries properly
  excluded_entries_list <- excluded_entries[!sapply(excluded_entries, is.null)]
  if (length(excluded_entries_list) > 0) {
    excluded_entries <- data.table::rbindlist(excluded_entries_list)
  } else {
    # Create empty data.table with proper structure
    excluded_entries <- data.table::data.table(
      cluster_number = integer(),
      gamma = numeric(),
      labels = list(),
      ic = numeric(),
      ic_vec = list(),
      best_labels = list(),
      effective_cluster_median = numeric(),
      raw_cluster_median = numeric(),
      admission_mode = character(),
      best_labels_raw_cluster_count = integer(),
      n_iter = integer(),
      mei = list(),
      k = integer(),
      excluded = logical(),
      exclusion_reason = character()
    )
  }
  
  if (length(valid_clusters) == 0) {
    if (verbose) {
      scice_message("CLUSTERING_MAIN: WARNING - No valid cluster numbers found!")
      scice_message("CLUSTERING_MAIN: All clusters were excluded during filtering")
      scice_message("CLUSTERING_MAIN: This will result in empty IC results")
      scice_message("CLUSTERING_MAIN: Returning results with only excluded cluster information")
    }
    # Return results with only excluded clusters
    return(list(
      gamma = excluded_entries$gamma,
      labels = excluded_entries$labels,
      ic = excluded_entries$ic,
      ic_vec = excluded_entries$ic_vec,
      n_cluster = excluded_entries$cluster_number,
      best_labels = excluded_entries$best_labels,
      effective_cluster_median = excluded_entries$effective_cluster_median,
      raw_cluster_median = excluded_entries$raw_cluster_median,
      admission_mode = excluded_entries$admission_mode,
      best_labels_raw_cluster_count = excluded_entries$best_labels_raw_cluster_count,
      n_iter = excluded_entries$n_iter,
      mei = excluded_entries$mei,
      k = excluded_entries$k,
      excluded = excluded_entries$excluded,
      exclusion_reason = excluded_entries$exclusion_reason
    ))
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
    result <- optimize_clustering(
      igraph_obj, cluster_num, gamma_range, objective_function,
      n_trials, n_bootstrap, seed, beta, n_iterations, max_iterations,
      resolution_tolerance, cluster_worker_budget, snn_graph, min_cluster_size, verbose,
      worker_id, in_parallel_context = TRUE,
      runtime_context = runtime_context
    )
    
    if (verbose) {
      opt_time <- as.numeric(difftime(Sys.time(), opt_start_time, units = "secs"))
      scice_message(paste(worker_id, ": Optimization completed in", round(opt_time, 3), "seconds"))
    }
    
    if (!is.null(result)) {
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
        admission_mode = result$admission_mode,
        best_labels_raw_cluster_count = result$best_labels_raw_cluster_count,
        n_iter = result$n_iterations,
        mei = list(mei_scores),
        k = result$k,
        excluded = FALSE,
        exclusion_reason = "none"
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
    # Create empty data.table with proper structure
    successful_results <- data.table::data.table(
      cluster_number = integer(),
      gamma = numeric(),
      labels = list(),
      ic = numeric(),
      ic_vec = list(),
      best_labels = list(),
      effective_cluster_median = numeric(),
      raw_cluster_median = numeric(),
      admission_mode = character(),
      best_labels_raw_cluster_count = integer(),
      n_iter = integer(),
      mei = list(),
      k = integer(),
      excluded = logical(),
      exclusion_reason = character()
    )
  }
  
  # Combine all results (excluded + successful) with robust error handling
  tryCatch({
    if (nrow(excluded_entries) > 0 && nrow(successful_results) > 0) {
      # Ensure both data.tables have compatible column structures
      if (all(colnames(excluded_entries) %in% colnames(successful_results)) &&
          all(colnames(successful_results) %in% colnames(excluded_entries))) {
        results_dt <- data.table::rbindlist(list(excluded_entries, successful_results), fill = TRUE)
      } else {
        # If columns don't match, combine manually
        all_cols <- union(colnames(excluded_entries), colnames(successful_results))
        for (col in all_cols) {
          if (!col %in% colnames(excluded_entries)) {
            excluded_entries[[col]] <- switch(col,
              "cluster_number" = integer(nrow(excluded_entries)),
              "gamma" = numeric(nrow(excluded_entries)),
              "labels" = rep(list(NULL), nrow(excluded_entries)),
              "ic" = numeric(nrow(excluded_entries)),
              "ic_vec" = rep(list(NULL), nrow(excluded_entries)),
              "best_labels" = rep(list(NULL), nrow(excluded_entries)),
              "effective_cluster_median" = rep(NA_real_, nrow(excluded_entries)),
              "raw_cluster_median" = rep(NA_real_, nrow(excluded_entries)),
              "admission_mode" = rep(NA_character_, nrow(excluded_entries)),
              "best_labels_raw_cluster_count" = rep(NA_integer_, nrow(excluded_entries)),
              "n_iter" = integer(nrow(excluded_entries)),
              "mei" = rep(list(NULL), nrow(excluded_entries)),
              "k" = integer(nrow(excluded_entries)),
              "excluded" = logical(nrow(excluded_entries)),
              "exclusion_reason" = character(nrow(excluded_entries)))
          }
          if (!col %in% colnames(successful_results)) {
            successful_results[[col]] <- switch(col,
              "cluster_number" = integer(nrow(successful_results)),
              "gamma" = numeric(nrow(successful_results)),
              "labels" = rep(list(NULL), nrow(successful_results)),
              "ic" = numeric(nrow(successful_results)),
              "ic_vec" = rep(list(NULL), nrow(successful_results)),
              "best_labels" = rep(list(NULL), nrow(successful_results)),
              "effective_cluster_median" = rep(NA_real_, nrow(successful_results)),
              "raw_cluster_median" = rep(NA_real_, nrow(successful_results)),
              "admission_mode" = rep(NA_character_, nrow(successful_results)),
              "best_labels_raw_cluster_count" = rep(NA_integer_, nrow(successful_results)),
              "n_iter" = integer(nrow(successful_results)),
              "mei" = rep(list(NULL), nrow(successful_results)),
              "k" = integer(nrow(successful_results)),
              "excluded" = logical(nrow(successful_results)),
              "exclusion_reason" = character(nrow(successful_results)))
          }
        }
        results_dt <- data.table::rbindlist(list(excluded_entries, successful_results), fill = TRUE)
      }
    } else if (nrow(excluded_entries) > 0) {
      results_dt <- excluded_entries
    } else if (nrow(successful_results) > 0) {
      results_dt <- successful_results
    } else {
      # Both are empty - create empty data.table with proper structure
      results_dt <- data.table::data.table(
        cluster_number = integer(),
        gamma = numeric(),
        labels = list(),
        ic = numeric(),
        ic_vec = list(),
        best_labels = list(),
        effective_cluster_median = numeric(),
        raw_cluster_median = numeric(),
        admission_mode = character(),
        best_labels_raw_cluster_count = integer(),
        n_iter = integer(),
        mei = list(),
        k = integer(),
        excluded = logical(),
        exclusion_reason = character()
      )
    }
  }, error = function(e) {
    warning("Error combining results, creating empty results: ", e$message)
    results_dt <<- data.table::data.table(
      cluster_number = integer(),
      gamma = numeric(),
      labels = list(),
      ic = numeric(),
      ic_vec = list(),
      best_labels = list(),
      effective_cluster_median = numeric(),
      raw_cluster_median = numeric(),
      admission_mode = character(),
      best_labels_raw_cluster_count = integer(),
      n_iter = integer(),
      mei = list(),
      k = integer(),
      excluded = logical(),
      exclusion_reason = character()
    )
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
        "best_labels_raw_cluster_count", "excluded"
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
          )
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
  
  # Convert back to list format for compatibility - with defensive programming
  return(list(
    gamma = if ("gamma" %in% colnames(results_dt)) results_dt$gamma else numeric(0),
    labels = if ("labels" %in% colnames(results_dt)) results_dt$labels else list(),
    ic = if ("ic" %in% colnames(results_dt)) results_dt$ic else numeric(0),
    ic_vec = if ("ic_vec" %in% colnames(results_dt)) results_dt$ic_vec else list(),
    n_cluster = if ("cluster_number" %in% colnames(results_dt)) results_dt$cluster_number else integer(0),
    best_labels = if ("best_labels" %in% colnames(results_dt)) results_dt$best_labels else list(),
    effective_cluster_median = if ("effective_cluster_median" %in% colnames(results_dt)) results_dt$effective_cluster_median else numeric(0),
    raw_cluster_median = if ("raw_cluster_median" %in% colnames(results_dt)) results_dt$raw_cluster_median else numeric(0),
    admission_mode = if ("admission_mode" %in% colnames(results_dt)) results_dt$admission_mode else character(0),
    best_labels_raw_cluster_count = if ("best_labels_raw_cluster_count" %in% colnames(results_dt)) results_dt$best_labels_raw_cluster_count else integer(0),
    n_iter = if ("n_iter" %in% colnames(results_dt)) results_dt$n_iter else integer(0),
    mei = if ("mei" %in% colnames(results_dt)) results_dt$mei else list(),
    k = if ("k" %in% colnames(results_dt)) results_dt$k else integer(0),
    excluded = if ("excluded" %in% colnames(results_dt)) results_dt$excluded else logical(0),
    exclusion_reason = if ("exclusion_reason" %in% colnames(results_dt)) results_dt$exclusion_reason else character(0)
  ))
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
  
  if (verbose) {
    scice_message(paste(worker_id, ": Generated gamma sequence:"))
    for (i in 1:min(5, length(gamma_sequence))) {
      scice_message(paste(worker_id, ":   gamma[", i, "] = ", signif(gamma_sequence[i], 6), sep = ""))
    }
    if (length(gamma_sequence) > 5) {
      scice_message(paste(worker_id, ":   ... and", length(gamma_sequence) - 5, "more"))
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
    median_effective_clusters <- stats::median(effective_clusters_vec)
    median_clusters_int <- as.integer(median_effective_clusters)
    raw_cluster_median <- stats::median(raw_clusters_vec)
    raw_cluster_median_int <- as.integer(raw_cluster_median)
    hit_trials <- which(effective_clusters_vec == target_clusters)
    hit_count <- as.integer(length(hit_trials))
    raw_hit_trials <- which(raw_clusters_vec == target_clusters)
    raw_hit_count <- as.integer(length(raw_hit_trials))
    hit_rate <- as.double(hit_count) / as.double(n_trials)
    median_gap <- abs(median_effective_clusters - target_clusters)
    raw_median_gap <- abs(raw_cluster_median - target_clusters)
    within_median_window <- (median_gap <= 1)
    strict_valid <- (median_clusters_int == target_clusters)
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
            "- median_int =", median_clusters_int,
            "- median_raw =", signif(raw_cluster_median, 6),
            "- median gap =", round(median_gap, 3),
            "- hit trials =", hit_count, "/", n_trials,
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
        mean_clusters = median_effective_clusters,
        median_clusters_raw = median_effective_clusters,
        median_effective_clusters = median_effective_clusters,
        median_clusters_int = median_clusters_int,
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
          "- median_int =", median_clusters_int,
          "- median_raw =", signif(raw_cluster_median, 6),
          "- median gap =", round(median_gap, 3),
          "- hit trials =", hit_count, "/", n_trials,
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
      mean_clusters = median_effective_clusters,
      median_clusters_raw = median_effective_clusters,
      median_effective_clusters = median_effective_clusters,
      median_clusters_int = median_clusters_int,
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
    scice_message(paste(worker_id, ": Phase 2 - Filtering for target effective cluster count (strict-first with relaxed fallback):", target_clusters))
    cluster_counts <- table(mean_clusters)
    for (i in 1:length(cluster_counts)) {
      count_val <- names(cluster_counts)[i]
      freq <- cluster_counts[i]
      scice_message(paste(worker_id, ":   ", freq, "gammas -> ", count_val, " effective clusters", sep = ""))
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
    if (verbose) {
      scice_message(
        paste(
          worker_id,
          ": ERROR - No gammas satisfied raw-guarded strict/relaxed admission, unguarded relaxed fallback, or bounded raw-count fallback admission for effective cluster count",
          target_clusters
        )
      )
    }
    return(NULL)
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
        "gammas selected under", admission_mode, "admission for", target_clusters, "effective clusters"
      )
    )
    scice_message(paste(worker_id, ": Phase 3 - IC scores already computed during Phase 1 for valid gammas (all trials per admitted gamma)"))
  }
  
  valid_results <- gamma_results[valid_indices]
  gamma_sequence <- vapply(valid_results, function(x) x$gamma, numeric(1))
  ic_scores <- vapply(valid_results, function(x) x$ic, numeric(1))
  clustering_refs <- lapply(valid_results, function(x) x$matrix_ref)
  effective_cluster_medians <- vapply(
    valid_results,
    function(x) as.numeric(x$median_effective_clusters),
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
          ic = 1 / ic_result
        )
      }, mc.cores = iteration_workers)
      
      # Extract results
      new_refs <- lapply(new_results, function(x) x$matrix_ref)
      new_ic <- vapply(new_results, function(x) x$ic, numeric(1))
      release_cluster_matrix_refs(current_refs)
      
      if (verbose) {
        iter_time <- as.numeric(difftime(Sys.time(), iter_start, units = "secs"))
        scice_message(paste(worker_id, ": Iteration", iteration_count, "completed in", round(iter_time, 3), "seconds"))
        scice_message(paste(worker_id, ": New IC range: [", round(min(new_ic), 4), ", ", round(max(new_ic), 4), "]", sep = ""))
      }
      
      # Update IC history
      ic_history <- cbind(ic_history[, -1], new_ic)
      
      # Check for stability and convergence
      stable_indices <- apply(ic_history, 1, function(row) length(unique(row)) == 1)
      perfect_indices <- which(new_ic == 1)
      
      if (verbose) {
        scice_message(paste(worker_id, ": Stable gammas:", sum(stable_indices), "/", length(stable_indices)))
        scice_message(paste(worker_id, ": Perfect IC scores:", length(perfect_indices)))
      }
      
      # Selection criteria
      if (length(perfect_indices) > 0) {
        best_index <- perfect_indices[1]
        best_gamma <- current_gammas[best_index]
        best_ref <- new_refs[[best_index]]
        if (length(new_refs) > 1) {
          release_cluster_matrix_refs(new_refs[-best_index])
        }
        if (verbose) {
          scice_message(paste(worker_id, ": CONVERGED - Found perfect IC score at gamma =", signif(best_gamma, 6)))
        }
        converged <- TRUE
        break
      } else if (all(stable_indices)) {
        best_index <- which.min(new_ic)
        best_gamma <- current_gammas[best_index]
        best_ref <- new_refs[[best_index]]
        if (length(new_refs) > 1) {
          release_cluster_matrix_refs(new_refs[-best_index])
        }
        if (verbose) {
          scice_message(paste(worker_id, ": CONVERGED - All gammas stable, best IC =", round(new_ic[best_index], 4)))
        }
        converged <- TRUE
        break
      } else {
        # Continue with best performing gammas
        keep_indices <- (new_ic <= stats::quantile(new_ic, 0.5)) | stable_indices
        keep_indices[which.min(new_ic)] <- TRUE
        
        if (sum(keep_indices) == 1) {
          best_index <- which(keep_indices)
          best_gamma <- current_gammas[best_index]
          best_ref <- new_refs[[best_index]]
          release_cluster_matrix_refs(new_refs[!keep_indices])
          if (verbose) {
            scice_message(paste(worker_id, ": CONVERGED - Single gamma remaining, IC =", round(new_ic[best_index], 4)))
          }
          converged <- TRUE
          break
        }
        
        if (verbose) {
          scice_message(paste(worker_id, ": Continuing with", sum(keep_indices), "best gammas"))
        }

        release_cluster_matrix_refs(new_refs[!keep_indices])
        
        current_gammas <- current_gammas[keep_indices]
        current_refs <- new_refs[keep_indices]
        ic_history <- ic_history[keep_indices, , drop = FALSE]
        current_ic <- new_ic[keep_indices]
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

  best_clustering <- load_cluster_matrix(best_ref)
  best_gamma_diag_index <- which.min(abs(gamma_sequence - best_gamma))
  selected_effective_cluster_median <- effective_cluster_medians[[best_gamma_diag_index]]
  selected_raw_cluster_median <- raw_cluster_medians[[best_gamma_diag_index]]
  
  # Bootstrap analysis
  if (verbose) {
    scice_message(paste(worker_id, ": Phase 5 - Bootstrap analysis with", n_bootstrap, "iterations"))
    bootstrap_start <- Sys.time()
  }

  bootstrap_workers <- cap_workers_by_memory(
    nested_workers,
    estimate_trial_matrix_bytes(n_vertices, n_trials, 1L),
    runtime_context
  )
  bootstrap_log_every <- max(1L, as.integer(floor(n_bootstrap / 5)))
  should_log_bootstrap_step <- function(idx) {
    idx == 1L || idx == n_bootstrap || (idx %% bootstrap_log_every) == 0L
  }
  
  if (bootstrap_workers == 1) {
    ic_bootstrap <- numeric(n_bootstrap)
    
    for (i in seq_len(n_bootstrap)) {
      # Set deterministic seed for this bootstrap iteration if base seed provided
      if (!is.null(seed)) {
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
      # Set deterministic seed for this bootstrap iteration if base seed provided
      if (!is.null(seed)) {
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
      return(ic_value)
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
  
  # Extract raw best labels for diagnostics, then merge the final user-facing
  # best_labels so min_cluster_size affects the returned clustering.
  extracted_best <- extract_clustering_array(best_clustering)
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
        worker_id, ": Selected diagnostics - gamma =", signif(best_gamma, 6),
        "- effective_median =", signif(selected_effective_cluster_median, 6),
        "- raw_median =", signif(selected_raw_cluster_median, 6),
        "- admission_mode =", admission_mode,
        "- best_labels_raw_clusters =", best_labels_raw_cluster_count,
        "- best_labels_final_clusters =", best_labels_final_cluster_count
      )
    )
  }
  release_cluster_matrix(best_ref)
  rm(best_clustering)
  
  return(list(
    gamma = best_gamma,
    labels = extracted_best,
    ic_median = ic_median,
    ic_bootstrap = ic_bootstrap,
    best_labels = best_labels,
    effective_cluster_median = selected_effective_cluster_median,
    raw_cluster_median = selected_raw_cluster_median,
    admission_mode = admission_mode,
    best_labels_raw_cluster_count = best_labels_raw_cluster_count,
    n_iterations = k,
    k = k
  ))
} 
