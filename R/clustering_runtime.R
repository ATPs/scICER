# Global clustering cache environment for optimization
clustering_cache_env <- new.env(parent = emptyenv())

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
