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

gamma_seed_role_priority <- function(seed_role) {
  priorities <- c(
    selected = 1L,
    left = 2L,
    right = 2L,
    exact = 3L,
    near = 4L,
    seed = 5L
  )
  matched <- priorities[as.character(seed_role)]
  matched[is.na(matched)] <- 99L
  as.integer(matched)
}

normalize_gamma_seed_table <- function(gamma_seed_values, gamma_range) {
  empty_table <- data.frame(
    gamma = numeric(0),
    seed_role = character(0),
    final_cluster_count = numeric(0),
    raw_cluster_count = numeric(0),
    stringsAsFactors = FALSE
  )

  if (is.null(gamma_seed_values)) {
    return(empty_table)
  }

  if (is.data.frame(gamma_seed_values)) {
    seed_table <- gamma_seed_values
  } else if (is.list(gamma_seed_values) &&
             !is.null(gamma_seed_values$gamma) &&
             !is.data.frame(gamma_seed_values)) {
    seed_table <- data.frame(
      gamma = as.numeric(gamma_seed_values$gamma),
      seed_role = if (!is.null(gamma_seed_values$seed_role)) {
        as.character(gamma_seed_values$seed_role)
      } else {
        rep("seed", length(gamma_seed_values$gamma))
      },
      final_cluster_count = if (!is.null(gamma_seed_values$final_cluster_count)) {
        as.numeric(gamma_seed_values$final_cluster_count)
      } else {
        rep(NA_real_, length(gamma_seed_values$gamma))
      },
      raw_cluster_count = if (!is.null(gamma_seed_values$raw_cluster_count)) {
        as.numeric(gamma_seed_values$raw_cluster_count)
      } else {
        rep(NA_real_, length(gamma_seed_values$gamma))
      },
      stringsAsFactors = FALSE
    )
  } else {
    gamma_values <- as.numeric(gamma_seed_values)
    seed_table <- data.frame(
      gamma = gamma_values,
      seed_role = rep("seed", length(gamma_values)),
      final_cluster_count = rep(NA_real_, length(gamma_values)),
      raw_cluster_count = rep(NA_real_, length(gamma_values)),
      stringsAsFactors = FALSE
    )
  }

  required_cols <- c("gamma", "seed_role", "final_cluster_count", "raw_cluster_count")
  missing_cols <- setdiff(required_cols, colnames(seed_table))
  for (col_name in missing_cols) {
    seed_table[[col_name]] <- if (col_name == "seed_role") {
      rep("seed", nrow(seed_table))
    } else {
      rep(NA_real_, nrow(seed_table))
    }
  }

  gamma_lower <- min(gamma_range)
  gamma_upper <- max(gamma_range)
  tolerance <- max(.Machine$double.eps^0.5, abs(gamma_upper - gamma_lower) * 1e-8)
  seed_table <- seed_table[
    is.finite(seed_table$gamma) &
      seed_table$gamma >= (gamma_lower - tolerance) &
      seed_table$gamma <= (gamma_upper + tolerance),
    required_cols,
    drop = FALSE
  ]
  if (nrow(seed_table) == 0L) {
    return(empty_table)
  }

  seed_table$seed_role[!nzchar(seed_table$seed_role)] <- "seed"
  seed_table$seed_role[is.na(seed_table$seed_role)] <- "seed"
  seed_table$role_priority <- gamma_seed_role_priority(seed_table$seed_role)
  seed_table <- seed_table[order(seed_table$gamma, seed_table$role_priority), , drop = FALSE]
  seed_table <- seed_table[!duplicated(seed_table[, c("gamma", "seed_role")]), , drop = FALSE]
  seed_table$role_priority <- NULL
  rownames(seed_table) <- NULL
  seed_table
}

thin_gamma_candidates_by_gap <- function(values, protected_values = numeric(0),
                                         objective_function = "CPM",
                                         min_log_gap = 0.08) {
  values <- sort(unique(as.numeric(values[is.finite(values)])))
  protected_values <- sort(unique(as.numeric(protected_values[is.finite(protected_values)])))
  if (length(values) <= 1L || objective_function != "CPM") {
    return(values)
  }

  keep <- numeric(0)
  compare_against <- sort(unique(c(protected_values, keep)))
  for (value in values) {
    compare_against <- sort(unique(c(protected_values, keep)))
    if (length(compare_against) == 0L) {
      keep <- c(keep, value)
      next
    }
    log_distance <- min(abs(log(value) - log(compare_against)))
    if (!is.finite(log_distance) || log_distance >= min_log_gap) {
      keep <- c(keep, value)
    }
  }
  sort(unique(keep))
}

select_evenly_spaced_gamma_values <- function(values, n_keep) {
  values <- sort(unique(as.numeric(values[is.finite(values)])))
  n_keep <- as.integer(n_keep)[1]
  if (length(values) <= n_keep || n_keep <= 0L) {
    return(values)
  }
  keep_positions <- unique(as.integer(round(seq(1, length(values), length.out = n_keep))))
  values[keep_positions]
}

build_even_interior_gamma_points <- function(gamma_range, n_points, objective_function) {
  n_points <- as.integer(n_points)[1]
  if (n_points <= 0L) {
    return(numeric(0))
  }

  lower <- min(gamma_range)
  upper <- max(gamma_range)
  if (!is.finite(lower) || !is.finite(upper) || identical(lower, upper)) {
    return(rep(lower, n_points))
  }

  if (objective_function == "CPM") {
    points <- exp(seq(log(lower), log(upper), length.out = n_points + 2L))
  } else {
    points <- seq(lower, upper, length.out = n_points + 2L)
  }
  as.numeric(points[-c(1L, length(points))])
}

fill_gamma_values_to_budget <- function(existing_values, budget, gamma_range,
                                        objective_function) {
  existing_values <- sort(unique(as.numeric(existing_values[is.finite(existing_values)])))
  budget <- as.integer(budget)[1]
  if (budget <= length(existing_values)) {
    return(existing_values)
  }

  candidates <- build_even_interior_gamma_points(
    gamma_range = gamma_range,
    n_points = max(0L, budget * 2L),
    objective_function = objective_function
  )
  if (objective_function == "CPM") {
    candidates <- thin_gamma_candidates_by_gap(
      values = candidates,
      protected_values = existing_values,
      objective_function = objective_function
    )
  } else {
    candidates <- sort(unique(candidates))
  }
  candidates <- select_evenly_spaced_gamma_values(
    candidates,
    max(0L, budget - length(existing_values))
  )

  for (candidate in candidates) {
    if (!any(abs(existing_values - candidate) <= .Machine$double.eps^0.5)) {
      existing_values <- sort(unique(c(existing_values, candidate)))
    }
  }
  existing_values
}

build_secondary_gamma_points <- function(primary_values, gamma_range,
                                         objective_function, n_points) {
  n_points <- as.integer(n_points)[1]
  if (n_points <= 0L) {
    return(numeric(0))
  }

  current_values <- sort(unique(as.numeric(primary_values[is.finite(primary_values)])))
  lower <- min(gamma_range)
  upper <- max(gamma_range)
  if (!any(abs(current_values - lower) <= .Machine$double.eps^0.5)) {
    current_values <- sort(unique(c(current_values, lower)))
  }
  if (!any(abs(current_values - upper) <= .Machine$double.eps^0.5)) {
    current_values <- sort(unique(c(current_values, upper)))
  }

  secondary_values <- numeric(0)
  transform_gamma <- if (objective_function == "CPM") log else identity
  inverse_transform <- if (objective_function == "CPM") exp else identity

  for (idx in seq_len(n_points)) {
    sorted_values <- sort(unique(c(current_values, secondary_values)))
    if (length(sorted_values) < 2L) {
      break
    }
    transformed <- transform_gamma(sorted_values)
    gap_widths <- diff(transformed)
    gap_idx <- which.max(gap_widths)
    if (length(gap_idx) == 0L || !is.finite(gap_widths[gap_idx]) || gap_widths[gap_idx] <= 0) {
      break
    }
    midpoint <- inverse_transform(mean(transformed[c(gap_idx, gap_idx + 1L)]))
    if (!is.finite(midpoint)) {
      break
    }
    secondary_values <- sort(unique(c(secondary_values, midpoint)))
    current_values <- sort(unique(c(current_values, midpoint)))
  }

  secondary_values[!(secondary_values %in% primary_values)]
}

build_optimization_gamma_batches <- function(gamma_range, gamma_seed_values,
                                             target_clusters, objective_function,
                                             resolution_tolerance, n_vertices,
                                             primary_budget = 8L,
                                             secondary_budget = 4L) {
  primary_budget <- max(1L, as.integer(primary_budget))
  secondary_budget <- max(0L, as.integer(secondary_budget))

  seed_table <- normalize_gamma_seed_table(gamma_seed_values, gamma_range)
  base_grid <- sort(unique(build_gamma_sequence_for_range(
    gamma_range,
    objective_function = objective_function,
    resolution_tolerance = resolution_tolerance,
    n_vertices = n_vertices
  )))

  anchors <- sort(unique(as.numeric(gamma_range)))
  exact_values <- numeric(0)
  near_values <- numeric(0)
  generic_seed_values <- numeric(0)

  if (nrow(seed_table) > 0L) {
    anchors <- sort(unique(c(
      anchors,
      seed_table$gamma[seed_table$seed_role %in% c("left", "right", "selected")]
    )))
    exact_values <- sort(unique(seed_table$gamma[seed_table$seed_role == "exact"]))
    near_values <- sort(unique(seed_table$gamma[seed_table$seed_role == "near"]))
    generic_seed_values <- sort(unique(seed_table$gamma[seed_table$seed_role %in% c("seed")]))
  }

  primary_values <- anchors

  remaining_slots <- max(0L, primary_budget - length(primary_values))
  if (remaining_slots > 0L && length(exact_values) > 0L) {
    exact_candidates <- exact_values[!(exact_values %in% primary_values)]
    exact_candidates <- select_evenly_spaced_gamma_values(exact_candidates, remaining_slots)
    primary_values <- sort(unique(c(primary_values, exact_candidates)))
  }

  remaining_slots <- max(0L, primary_budget - length(primary_values))
  if (remaining_slots > 0L && length(near_values) > 0L) {
    near_candidates <- near_values[!(near_values %in% primary_values)]
    near_candidates <- thin_gamma_candidates_by_gap(
      values = near_candidates,
      protected_values = primary_values,
      objective_function = objective_function
    )
    near_candidates <- select_evenly_spaced_gamma_values(near_candidates, remaining_slots)
    primary_values <- sort(unique(c(primary_values, near_candidates)))
  }

  remaining_slots <- max(0L, primary_budget - length(primary_values))
  if (remaining_slots > 0L && length(generic_seed_values) > 0L) {
    generic_candidates <- generic_seed_values[!(generic_seed_values %in% primary_values)]
    generic_candidates <- thin_gamma_candidates_by_gap(
      values = generic_candidates,
      protected_values = primary_values,
      objective_function = objective_function
    )
    generic_candidates <- select_evenly_spaced_gamma_values(generic_candidates, remaining_slots)
    primary_values <- sort(unique(c(primary_values, generic_candidates)))
  }

  primary_values <- fill_gamma_values_to_budget(
    existing_values = primary_values,
    budget = primary_budget,
    gamma_range = gamma_range,
    objective_function = objective_function
  )
  primary_values <- primary_values[primary_values >= min(gamma_range) & primary_values <= max(gamma_range)]
  primary_values <- sort(unique(primary_values))

  secondary_values <- build_secondary_gamma_points(
    primary_values = primary_values,
    gamma_range = gamma_range,
    objective_function = objective_function,
    n_points = secondary_budget
  )

  list(
    primary_gammas = sort(unique(primary_values)),
    secondary_gammas = sort(unique(secondary_values)),
    seed_table = seed_table
  )
}

derive_gamma_admission_state <- function(gamma_results, target_clusters,
                                         min_cluster_size = 1L,
                                         verbose = FALSE,
                                         worker_id = "OPTIMIZER") {
  if (length(gamma_results) == 0L) {
    return(list(
      valid_indices = integer(0),
      admission_mode = "none",
      mean_clusters = numeric(0),
      mean_clusters_int = integer(0),
      hit_counts = integer(0),
      strict_flags = logical(0),
      relaxed_flags = logical(0),
      raw_strict_flags = logical(0),
      raw_relaxed_flags = logical(0),
      soft_raw_guard_flags = logical(0),
      hard_raw_guard_flags = logical(0),
      exact_hit_gamma_count = 0L,
      selected_raw_gaps = numeric(0),
      best_raw_gap = Inf
    ))
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
    for (i in seq_along(cluster_counts)) {
      count_val <- names(cluster_counts)[i]
      freq <- cluster_counts[i]
      scice_message(paste(worker_id, ":   ", freq, "gammas -> ", count_val, " final merged clusters", sep = ""))
    }
    cluster_counts_int <- table(mean_clusters_int)
    for (i in seq_along(cluster_counts_int)) {
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
      for (i in seq_along(admitted_hit_counts)) {
        hit_val <- names(admitted_hit_counts)[i]
        freq <- admitted_hit_counts[i]
        scice_message(paste(worker_id, ":   ", freq, "admitted gammas -> ", hit_val, " hit trials", sep = ""))
      }
    }
  }

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

  if (verbose &&
      length(valid_indices) > 0L &&
      length(valid_indices) < pre_refine_candidate_count &&
      min_cluster_size > 1L &&
      any(is.finite(selected_raw_gaps))) {
    scice_message(
      paste(
        worker_id, ": Phase 2 refinement - retaining",
        length(valid_indices), "candidate gammas with minimum raw-cluster gap",
        signif(best_raw_gap, 6), "to target", target_clusters
      )
    )
  }

  list(
    valid_indices = valid_indices,
    admission_mode = admission_mode,
    mean_clusters = mean_clusters,
    mean_clusters_int = mean_clusters_int,
    hit_counts = hit_counts,
    strict_flags = strict_flags,
    relaxed_flags = relaxed_flags,
    raw_strict_flags = raw_strict_flags,
    raw_relaxed_flags = raw_relaxed_flags,
    soft_raw_guard_flags = soft_raw_guard_flags,
    hard_raw_guard_flags = hard_raw_guard_flags,
    exact_hit_gamma_count = sum(hit_counts > 0L),
    selected_raw_gaps = selected_raw_gaps,
    best_raw_gap = best_raw_gap
  )
}

should_expand_phase1_secondary <- function(valid_indices, admission_mode,
                                           exact_hit_gamma_count) {
  guarded_modes <- c(
    "raw_strict_soft", "strict_soft", "relaxed_soft",
    "strict_hard", "relaxed_hard"
  )
  if (length(valid_indices) == 0L) {
    return(TRUE)
  }
  if (isTRUE(exact_hit_gamma_count > 0L)) {
    return(FALSE)
  }
  if (admission_mode %in% guarded_modes) {
    return(FALSE)
  }
  admission_mode %in% c("relaxed_unguarded", "raw_relaxed_unguarded")
}

should_skip_phase4_refinement <- function(candidate_count, best_ic,
                                          exact_hit_gamma_count) {
  candidate_count <- as.integer(candidate_count)[1]
  best_ic <- as.numeric(best_ic)[1]
  exact_hit_gamma_count <- as.integer(exact_hit_gamma_count)[1]
  candidate_count <= 2L &&
    is.finite(best_ic) &&
    best_ic <= 1.005 &&
    exact_hit_gamma_count > 0L
}

phase4_iteration_cap_for_mode <- function(admission_mode) {
  if (admission_mode %in% c("relaxed_unguarded", "raw_relaxed_unguarded")) {
    return(2L)
  }
  3L
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

  bootstrap_start <- Sys.time()
  if (verbose) {
    scice_message(paste(worker_id, ": Phase 5 - Bootstrap analysis with", n_bootstrap, "iterations"))
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
  phase5_elapsed_sec <- if (exists("bootstrap_time")) {
    as.numeric(bootstrap_time)
  } else {
    as.numeric(difftime(Sys.time(), bootstrap_start, units = "secs"))
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
    best_labels_final_cluster_count = best_labels_final_cluster_count,
    phase5_elapsed_sec = phase5_elapsed_sec
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
  
  optimization_start_time <- Sys.time()
  n_vertices <- igraph::vcount(igraph_obj)
  heartbeat <- create_heartbeat_logger(verbose = verbose, context = worker_id)
  range_width <- abs(gamma_range[2] - gamma_range[1])
  n_steps <- if (n_vertices >= 200000 && range_width <= max(resolution_tolerance * 10, .Machine$double.eps)) 5L else 11L
  delta_n <- 2
  
  if (verbose) {
    scice_message(paste(worker_id, ": Creating gamma sequence with", n_steps, "steps"))
    scice_message(paste(worker_id, ": Gamma range bounds: [", signif(gamma_range[1], 6), ", ", signif(gamma_range[2], 6), "]", sep = ""))
  }

  gamma_batches <- build_optimization_gamma_batches(
    gamma_range = gamma_range,
    gamma_seed_values = gamma_seed_values,
    target_clusters = target_clusters,
    objective_function = objective_function,
    resolution_tolerance = resolution_tolerance,
    n_vertices = n_vertices,
    primary_budget = 8L,
    secondary_budget = 4L
  )
  primary_gamma_sequence <- gamma_batches$primary_gammas
  secondary_gamma_sequence <- gamma_batches$secondary_gammas
  gamma_seed_table <- gamma_batches$seed_table

  if (verbose) {
    scice_message(
      paste(
        worker_id, ": Primary Phase 1 gamma batch (", length(primary_gamma_sequence), "): ",
        paste(signif(primary_gamma_sequence, 6), collapse = ", "),
        sep = ""
      )
    )
    if (length(secondary_gamma_sequence) > 0L) {
      scice_message(
        paste(
          worker_id, ": Secondary Phase 1 gamma batch (", length(secondary_gamma_sequence), "): ",
          paste(signif(secondary_gamma_sequence, 6), collapse = ", "),
          sep = ""
        )
      )
    }
    if (nrow(gamma_seed_table) > 0L) {
      scice_message(
        paste(
          worker_id, ": Included", nrow(gamma_seed_table),
          "target-specific gamma seeds from shared search diagnostics"
        )
      )
    }
  }

  compute_phase1_nested_workers <- function(batch_gamma_count) {
    nested_workers <- if (in_parallel_context) {
      max(
        1L,
        as.integer(round(as.double(n_workers) / as.double(max(1L, batch_gamma_count))))
      )
    } else {
      max(1L, as.integer(n_workers))
    }
    if (!in_parallel_context) {
      nested_workers <- min(nested_workers, max(1L, as.integer(batch_gamma_count)))
    }
    cap_workers_by_memory(
      nested_workers,
      estimate_trial_matrix_bytes(n_vertices, n_trials, 1L),
      runtime_context
    )
  }

  phase1_log_every <- max(
    1L,
    as.integer(
      floor(max(length(primary_gamma_sequence), length(secondary_gamma_sequence), 1L) / 5)
    )
  )
  should_log_phase1_step <- function(idx) {
    idx == 1L || idx %% phase1_log_every == 0L
  }
  current_phase1_gamma_total <- 1L
  
  evaluate_gamma <- function(gamma_val, gamma_idx) {
    log_this_gamma <- verbose && should_log_phase1_step(gamma_idx)
    gamma_start_time <- NULL
    if (log_this_gamma) {
      gamma_start_time <- Sys.time()
      scice_message(
        paste(
          worker_id, ": Phase 1 progress gamma", gamma_idx,
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
          "phase1 running - gamma", gamma_idx, "/", current_phase1_gamma_total,
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
            worker_id, ": Phase 1 progress gamma", gamma_idx, "/", current_phase1_gamma_total,
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
          worker_id, ": Phase 1 progress gamma", gamma_idx, "/", current_phase1_gamma_total,
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
  evaluate_gamma_batch <- function(gamma_sequence, batch_label) {
    if (length(gamma_sequence) == 0L) {
      return(list(
        results = list(),
        elapsed_sec = 0,
        gamma_count = 0L,
        leiden_runs = 0L,
        nested_workers = 1L
      ))
    }

    estimated_phase1_bytes <- estimate_trial_matrix_bytes(n_vertices, n_trials, length(gamma_sequence))
    if (should_enable_spill(runtime_context, estimated_phase1_bytes)) {
      activate_runtime_spill(runtime_context, estimated_bytes = estimated_phase1_bytes)
    }

    nested_workers <- compute_phase1_nested_workers(length(gamma_sequence))
    if (verbose) {
      scice_message(
        paste(
          worker_id, ":", batch_label, "- Testing", length(gamma_sequence),
          "gamma values with", n_trials, "trials each"
        )
      )
      if (in_parallel_context) {
        scice_message(
          paste(
            worker_id, ": Running in parallel context - deriving nested worker budget from",
            n_workers, "optimizer workers across", length(gamma_sequence), "gamma values"
          )
        )
      }
      scice_message(
        paste(
          worker_id, ": Worker budget for this cluster:", n_workers,
          "- using", nested_workers, "workers for", batch_label, "gamma evaluation"
        )
      )
    }

    batch_start <- Sys.time()
    current_phase1_gamma_total <<- max(1L, as.integer(length(gamma_sequence)))
    if (nested_workers == 1L) {
      gamma_results <- vector("list", length(gamma_sequence))
      for (gamma_idx in seq_along(gamma_sequence)) {
        gamma_results[[gamma_idx]] <- evaluate_gamma(gamma_sequence[[gamma_idx]], gamma_idx)
      }
    } else {
      gamma_results <- cross_platform_mclapply(
        seq_along(gamma_sequence),
        function(gamma_idx) {
          evaluate_gamma(gamma_sequence[[gamma_idx]], gamma_idx)
        },
        mc.cores = nested_workers
      )
    }
    elapsed_sec <- as.numeric(difftime(Sys.time(), batch_start, units = "secs"))
    if (verbose) {
      scice_message(
        paste(worker_id, ":", batch_label, "completed in", round(elapsed_sec, 3), "seconds")
      )
    }
    list(
      results = gamma_results,
      elapsed_sec = elapsed_sec,
      gamma_count = length(gamma_sequence),
      leiden_runs = as.integer(length(gamma_sequence) * max(1L, as.integer(n_trials))),
      nested_workers = nested_workers
    )
  }

  if (verbose) {
    phase1_expected_runs <- as.integer(
      (length(primary_gamma_sequence) + length(secondary_gamma_sequence)) *
        max(1L, as.integer(n_trials))
    )
    scice_message(paste(worker_id, ": Phase 1 maximum expected Leiden runs:", format(phase1_expected_runs, big.mark = ",")))
  }

  primary_phase1 <- evaluate_gamma_batch(primary_gamma_sequence, "Primary Phase 1")
  primary_admission <- derive_gamma_admission_state(
    gamma_results = primary_phase1$results,
    target_clusters = target_clusters,
    min_cluster_size = min_cluster_size,
    verbose = FALSE,
    worker_id = worker_id
  )
  secondary_phase1_used <- should_expand_phase1_secondary(
    valid_indices = primary_admission$valid_indices,
    admission_mode = primary_admission$admission_mode,
    exact_hit_gamma_count = primary_admission$exact_hit_gamma_count
  )

  secondary_phase1 <- list(
    results = list(),
    elapsed_sec = 0,
    gamma_count = 0L,
    leiden_runs = 0L,
    nested_workers = 1L
  )
  secondary_trigger_reason <- NULL
  if (isTRUE(secondary_phase1_used) && length(secondary_gamma_sequence) > 0L) {
    if (length(primary_admission$valid_indices) == 0L) {
      secondary_trigger_reason <- "primary batch produced no admitted gamma candidates"
    } else {
      secondary_trigger_reason <- paste(
        "primary batch ended at", primary_admission$admission_mode,
        "without exact final-hit gamma support"
      )
    }
    if (verbose) {
      scice_message(paste(worker_id, ": Secondary batch triggered because", secondary_trigger_reason))
    }
    secondary_phase1 <- evaluate_gamma_batch(secondary_gamma_sequence, "Secondary Phase 1")
  }

  gamma_results <- c(primary_phase1$results, secondary_phase1$results)
  phase1_primary_gamma_count <- as.integer(primary_phase1$gamma_count)
  phase1_secondary_gamma_count <- as.integer(secondary_phase1$gamma_count)
  phase1_total_gamma_count <- as.integer(phase1_primary_gamma_count + phase1_secondary_gamma_count)
  phase1_elapsed_sec <- as.numeric(primary_phase1$elapsed_sec + secondary_phase1$elapsed_sec)
  phase1_leiden_runs <- as.integer(primary_phase1$leiden_runs + secondary_phase1$leiden_runs)
  phase1_nested_workers <- max(
    1L,
    as.integer(max(primary_phase1$nested_workers, secondary_phase1$nested_workers))
  )
  if (verbose) {
    scice_message(
      paste(
        worker_id, ": Phase 1 completed in", round(phase1_elapsed_sec, 3), "seconds",
        "- primary gammas =", phase1_primary_gamma_count,
        "- secondary gammas =", phase1_secondary_gamma_count,
        "- total runs =", format(phase1_leiden_runs, big.mark = ",")
      )
    )
  }

  admission_state <- derive_gamma_admission_state(
    gamma_results = gamma_results,
    target_clusters = target_clusters,
    min_cluster_size = min_cluster_size,
    verbose = verbose,
    worker_id = worker_id
  )
  valid_indices <- admission_state$valid_indices
  admission_mode <- admission_state$admission_mode
  selected_raw_gaps <- admission_state$selected_raw_gaps
  best_raw_gap <- admission_state$best_raw_gap
  exact_hit_gamma_count <- as.integer(admission_state$exact_hit_gamma_count)

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
      phase1_primary_gamma_count = phase1_primary_gamma_count,
      phase1_secondary_gamma_count = phase1_secondary_gamma_count,
      phase1_total_gamma_count = phase1_total_gamma_count,
      phase1_elapsed_sec = phase1_elapsed_sec,
      phase1_leiden_runs = phase1_leiden_runs,
      secondary_phase1_used = isTRUE(secondary_phase1_used),
      exact_hit_gamma_count = exact_hit_gamma_count,
      phase4_iterations = 0L,
      phase4_elapsed_sec = 0,
      phase5_elapsed_sec = 0,
      optimization_elapsed_sec = as.numeric(difftime(Sys.time(), optimization_start_time, units = "secs")),
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
  exact_hit_priority_enabled <- !startsWith(admission_mode, "raw_")
  if (exact_hit_priority_enabled && any(exact_hit_gamma_flags)) {
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
  prefer_exact_hits <- exact_hit_priority_enabled && any(exact_hit_gamma_flags)
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
  phase4_iterations <- 0L
  phase4_elapsed_sec <- 0
  phase4_skip_reason <- NULL
  
  if (ic_scores[best_index] != 1 && length(gamma_sequence) > 1) {
    if (should_skip_phase4_refinement(length(gamma_sequence), ic_scores[best_index], sum(exact_hit_gamma_flags))) {
      phase4_skip_reason <- "candidate set <= 2 with best IC <= 1.005 and exact-hit support"
      if (verbose) {
        scice_message(paste(worker_id, ": Phase 4 skipped because", phase4_skip_reason))
      }
      if (length(clustering_refs) > 1) {
        release_cluster_matrix_refs(clustering_refs[-best_index])
      }
    } else {
      phase4_iteration_limit <- phase4_iteration_cap_for_mode(admission_mode)
      if (verbose) {
        scice_message(paste(worker_id, ": Phase 4 - Iterative improvement (IC =", round(ic_scores[best_index], 4), "< 1.0)"))
        scice_message(paste(worker_id, ": Starting with", length(gamma_sequence), "gamma values"))
        scice_message(paste(worker_id, ": Phase 4 iteration cap:", phase4_iteration_limit))
        iterative_start <- Sys.time()
      } else {
        iterative_start <- Sys.time()
      }
      
      current_refs <- clustering_refs
      current_gammas <- gamma_sequence
      current_ic <- ic_scores
      ic_history <- matrix(rep(current_ic, 10), nrow = length(current_ic))
      converged <- FALSE
      iteration_workers <- cap_workers_by_memory(
        phase1_nested_workers,
        estimate_trial_matrix_bytes(n_vertices, n_trials, 1L) * 2,
        runtime_context
      )

      while (k < max_iterations && phase4_iterations < phase4_iteration_limit) {
        k <- k + delta_n
        phase4_iterations <- phase4_iterations + 1L
        
        if (verbose) {
          scice_message(paste(worker_id, ": Iteration", phase4_iterations, "(k =", k, ") - Refining", length(current_gammas), "gamma values"))
          iter_start <- Sys.time()
        }
        
        new_results <- cross_platform_mclapply(seq_along(current_gammas), function(i) {
          gamma_val <- current_gammas[i]
          current_matrix <- load_cluster_matrix(current_refs[[i]])
          
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
        candidate_ic_history <- cbind(ic_history[, -1, drop = FALSE], new_ic)
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
                  worker_id, ": Iteration", phase4_iterations, "- retaining",
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
          scice_message(paste(worker_id, ": Iteration", phase4_iterations, "completed in", round(iter_time, 3), "seconds"))
          scice_message(paste(worker_id, ": New IC range: [", round(min(candidate_ic), 4), ", ", round(max(candidate_ic), 4), "]", sep = ""))
        }
        
        ic_history <- candidate_ic_history
        stable_indices <- apply(ic_history, 1, function(row) length(unique(row)) == 1)
        perfect_indices <- which(candidate_ic == 1)
        
        if (verbose) {
          scice_message(paste(worker_id, ": Stable gammas:", sum(stable_indices), "/", length(stable_indices)))
          scice_message(paste(worker_id, ": Perfect IC scores:", length(perfect_indices)))
        }
        
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
          keep_indices <- (candidate_ic <= stats::quantile(candidate_ic, 0.5)) | stable_indices
          keep_indices[which.min(candidate_ic)] <- TRUE
          keep_pool <- which(keep_indices)
          keep_limit <- if (phase4_iterations == 1L) 4L else 2L
          if (length(keep_pool) > keep_limit) {
            best_global_idx <- which.min(candidate_ic)
            stable_pool <- which(stable_indices)
            stable_best_idx <- if (length(stable_pool) > 0L) {
              stable_pool[[which.min(candidate_ic[stable_pool])]]
            } else {
              integer(0)
            }
            ordered_pool <- keep_pool[order(candidate_ic[keep_pool], candidate_gammas[keep_pool])]
            keep_pool <- unique(c(best_global_idx, stable_best_idx, ordered_pool))
            keep_pool <- keep_pool[seq_len(min(keep_limit, length(keep_pool)))]
          }
          keep_indices <- seq_along(candidate_ic) %in% keep_pool
          
          if (sum(keep_indices) == 1L) {
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
            scice_message(
              paste(
                worker_id, ": Continuing with", sum(keep_indices),
                "best gammas (cap =", keep_limit, ")"
              )
            )
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
      
      phase4_elapsed_sec <- as.numeric(difftime(Sys.time(), iterative_start, units = "secs"))
      if (verbose) {
        scice_message(paste(worker_id, ": Phase 4 completed in", round(phase4_elapsed_sec, 3), "seconds after", phase4_iterations, "iterations"))
      }
    }
  } else {
    if (length(clustering_refs) > 1) {
      release_cluster_matrix_refs(clustering_refs[-best_index])
    }
    if (length(gamma_sequence) <= 1L) {
      phase4_skip_reason <- "single admitted gamma candidate"
    } else if (ic_scores[best_index] == 1) {
      phase4_skip_reason <- "perfect IC already achieved"
    }
    if (verbose && !is.null(phase4_skip_reason)) {
      scice_message(paste(worker_id, ": Phase 4 skipped because", phase4_skip_reason))
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
    n_workers = phase1_nested_workers,
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
  finalized_result$phase1_primary_gamma_count <- phase1_primary_gamma_count
  finalized_result$phase1_secondary_gamma_count <- phase1_secondary_gamma_count
  finalized_result$phase1_total_gamma_count <- phase1_total_gamma_count
  finalized_result$phase1_elapsed_sec <- phase1_elapsed_sec
  finalized_result$phase1_leiden_runs <- phase1_leiden_runs
  finalized_result$secondary_phase1_used <- isTRUE(secondary_phase1_used) && phase1_secondary_gamma_count > 0L
  finalized_result$exact_hit_gamma_count <- exact_hit_gamma_count
  finalized_result$phase4_iterations <- as.integer(phase4_iterations)
  finalized_result$phase4_elapsed_sec <- as.numeric(phase4_elapsed_sec)
  finalized_result$optimization_elapsed_sec <- as.numeric(difftime(Sys.time(), optimization_start_time, units = "secs"))
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
