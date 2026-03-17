# Test basic scICER functions

test_that("calculate_ecs works correctly", {
  # Create simple test data
  cluster_a <- c(1, 1, 2, 2, 3, 3)
  cluster_b <- c(1, 1, 2, 2, 3, 3)
  
  # Test identical clusterings
  expect_equal(calculate_ecs(cluster_a, cluster_b), 1.0, tolerance = 1e-6)
  
  # Test with return_vector
  scores <- calculate_ecs(cluster_a, cluster_b, return_vector = TRUE)
  expect_length(scores, length(cluster_a))
  expect_true(all(scores >= 0 & scores <= 1))
})

test_that("calculate_ecs handles different clusterings", {
  cluster_a <- c(1, 1, 2, 2, 3, 3)
  cluster_b <- c(1, 2, 1, 2, 1, 2)
  
  # Different clusterings should have similarity < 1
  similarity <- calculate_ecs(cluster_a, cluster_b)
  expect_true(similarity < 1.0)
  expect_true(similarity >= 0.0)
})

test_that("extract_clustering_array works", {
  # Create test clustering matrix
  n_cells <- 10
  n_trials <- 5
  clustering_matrix <- matrix(
    sample(1:3, n_cells * n_trials, replace = TRUE),
    nrow = n_cells, ncol = n_trials
  )
  
  result <- extract_clustering_array(clustering_matrix)
  
  expect_true(is.list(result))
  expect_true("arr" %in% names(result))
  expect_true("parr" %in% names(result))
  expect_true(length(result$arr) > 0)
  expect_true(length(result$parr) > 0)
  expect_equal(length(result$arr), length(result$parr))
  expect_true(all(result$parr >= 0 & result$parr <= 1))
  expect_equal(sum(result$parr), 1, tolerance = 1e-6)
})

test_that("calculate_ic works with clustering array", {
  # Create test data with identical clusterings (should give IC = 1.0)
  n_cells <- 10
  n_trials <- 3
  base_clustering <- sample(1:3, n_cells, replace = TRUE)
  clustering_matrix <- matrix(
    rep(base_clustering, n_trials),
    nrow = n_cells, ncol = n_trials
  )
  
  clustering_array <- extract_clustering_array(clustering_matrix)
  ic_score <- calculate_ic(clustering_array)
  
  expect_equal(ic_score, 1.0, tolerance = 1e-6)
})

test_that("calculate_mei_from_array works", {
  # Create test clustering array
  n_cells <- 10
  clustering_array <- list(
    arr = list(sample(1:3, n_cells, replace = TRUE)),
    parr = 1.0
  )
  
  mei_scores <- calculate_mei_from_array(clustering_array)
  
  expect_length(mei_scores, n_cells)  # One score per cell
  expect_true(all(mei_scores >= 0 & mei_scores <= 1))
})

test_that("scICE_clustering validates inputs", {
  # Test with invalid input
  expect_error(scICE_clustering("not_a_seurat_object"))
  expect_error(scICE_clustering(list(test = "data")))
})

test_that("utility functions handle edge cases", {
  # Test with single cluster
  single_cluster <- rep(1, 10)
  ecs_score <- calculate_ecs(single_cluster, single_cluster)
  expect_equal(ecs_score, 1.0, tolerance = 1e-6)
  
  # Test with empty input (should handle gracefully)
  expect_error(calculate_ecs(numeric(0), numeric(0)))
  expect_error(calculate_ecs(c(1, 2), c(1)))  # Different lengths
})

test_that("calculate_ecs scalar mode is stable for large inputs", {
  set.seed(1)
  n <- 250000
  cluster_a <- sample(1:10, n, replace = TRUE)
  cluster_b <- sample(1:10, n, replace = TRUE)

  scalar_score <- calculate_ecs(cluster_a, cluster_b, return_vector = FALSE)
  compact_labels <- function(x) as.integer(factor(x, levels = unique(x)))
  expected_scalar <- ClustAssess::element_sim(
    compact_labels(cluster_a),
    compact_labels(cluster_b),
    alpha = 0.9
  )
  vector_score <- mean(calculate_ecs(cluster_a, cluster_b, return_vector = TRUE))

  expect_equal(scalar_score, expected_scalar, tolerance = 1e-10)
  expect_true(is.finite(scalar_score))
  expect_true(scalar_score >= 0 && scalar_score <= 1)
  expect_true(is.finite(vector_score))
  expect_true(vector_score >= 0 && vector_score <= 1)
})

test_that("get_robust_labels requires explicit Seurat object for return_seurat", {
  data("pbmc_small", package = "SeuratObject")

  results <- scICE_clustering(
    object = pbmc_small,
    graph_name = "RNA_snn",
    cluster_range = 2:3,
    n_workers = 1,
    n_trials = 2,
    n_bootstrap = 2,
    seed = 123,
    remove_threshold = Inf,
    min_cluster_size = 1,
    verbose = FALSE
  )

  expect_error(
    get_robust_labels(results, return_seurat = TRUE),
    "please provide `object`"
  )

  seurat_with_labels <- get_robust_labels(
    results,
    threshold = Inf,
    return_seurat = TRUE,
    object = pbmc_small
  )
  expect_s4_class(seurat_with_labels, "Seurat")
  expect_true(any(grepl("^clusters_", colnames(seurat_with_labels@meta.data))))

  smaller_object <- subset(pbmc_small, cells = SeuratObject::Cells(pbmc_small)[1:40])
  expect_error(
    get_robust_labels(results, threshold = Inf, return_seurat = TRUE, object = smaller_object),
    "Cell count mismatch"
  )
})

test_that("get_robust_labels keeps backward compatibility for old results", {
  data("pbmc_small", package = "SeuratObject")

  results <- scICE_clustering(
    object = pbmc_small,
    graph_name = "RNA_snn",
    cluster_range = 2:3,
    n_workers = 1,
    n_trials = 2,
    n_bootstrap = 2,
    seed = 123,
    remove_threshold = Inf,
    min_cluster_size = 1,
    verbose = FALSE
  )
  results$seurat_object <- pbmc_small

  expect_warning(
    get_robust_labels(results, threshold = Inf, return_seurat = TRUE),
    "deprecated"
  )
})

test_that("count_effective_clusters handles thresholds and all-small fallback", {
  labels <- c(0L, 0L, 1L, 1L, 1L, 2L)

  expect_equal(scICER:::count_effective_clusters(labels, min_cluster_size = 1L), 3L)
  expect_equal(scICER:::count_effective_clusters(labels, min_cluster_size = 2L), 2L)
  expect_equal(scICER:::count_effective_clusters(labels, min_cluster_size = 4L), 0L)
  expect_equal(scICER:::count_effective_clusters(integer(0), min_cluster_size = 3L), 0L)
})

test_that("raw cluster guard thresholds match the default policy", {
  limits <- scICER:::raw_cluster_guard_limits(14L)

  expect_equal(limits$soft, 17L)
  expect_equal(limits$hard, 21L)
  expect_true(scICER:::passes_raw_cluster_guard(17, 14L, min_cluster_size = 2L, level = "soft"))
  expect_false(scICER:::passes_raw_cluster_guard(18, 14L, min_cluster_size = 2L, level = "soft"))
  expect_true(scICER:::passes_raw_cluster_guard(21, 14L, min_cluster_size = 2L, level = "hard"))
  expect_false(scICER:::passes_raw_cluster_guard(22, 14L, min_cluster_size = 2L, level = "hard"))
  expect_true(scICER:::passes_raw_cluster_guard(100, 14L, min_cluster_size = 1L, level = "soft"))
})

test_that("resolution search uses a tighter raw-cluster upper bound than phase-2 admission", {
  expect_equal(scICER:::raw_cluster_search_upper(8L), 9L)
  expect_equal(scICER:::raw_cluster_search_upper(14L), 16L)
  expect_equal(scICER:::raw_cluster_search_upper(20L), 22L)
})

test_that("resolution search classifier uses raw brackets", {
  raw_below <- scICER:::classify_resolution_search_state(
    raw_cluster_median = 7,
    effective_cluster_median = 7,
    target_clusters = 8L,
    min_cluster_size = 2L
  )
  expect_equal(raw_below$raw_class, "raw_below")
  expect_false(raw_below$over_fragmented)
  expect_equal(raw_below$lower_action, "increase_gamma")
  expect_equal(raw_below$upper_action, "increase_gamma")

  over_fragmented <- scICER:::classify_resolution_search_state(
    raw_cluster_median = 8,
    effective_cluster_median = 7,
    target_clusters = 8L,
    min_cluster_size = 2L
  )
  expect_equal(over_fragmented$raw_class, "raw_in_band")
  expect_true(over_fragmented$over_fragmented)
  expect_equal(over_fragmented$lower_action, "decrease_gamma")
  expect_equal(over_fragmented$upper_action, "decrease_gamma")

  raw_above_soft <- scICER:::classify_resolution_search_state(
    raw_cluster_median = 12,
    effective_cluster_median = 8,
    target_clusters = 8L,
    min_cluster_size = 2L
  )
  expect_equal(raw_above_soft$raw_class, "raw_above_soft")
  expect_false(raw_above_soft$raw_guard_soft)
  expect_false(raw_above_soft$over_fragmented)
  expect_equal(raw_above_soft$lower_action, "decrease_gamma")
  expect_equal(raw_above_soft$upper_action, "decrease_gamma")
})

test_that("clamp_gamma_range_to_raw_plateau narrows to exact raw plateaus and fallbacks", {
  noisy_exact <- scICER:::clamp_gamma_range_to_raw_plateau(
    gamma_sequence = c(0.1, 0.2, 0.3, 0.4),
    raw_cluster_medians = c(3, 2, 4, 4),
    target_clusters = 3L
  )
  expect_equal(noisy_exact$mode, "raw_exact")
  expect_equal(noisy_exact$bounds, c(0.1, 0.2))
  expect_equal(noisy_exact$indices, c(1L, 2L))

  exact <- scICER:::clamp_gamma_range_to_raw_plateau(
    gamma_sequence = c(0.1, 0.2, 0.3, 0.4, 0.5),
    raw_cluster_medians = c(2, 3, 3, 3, 4),
    target_clusters = 3L
  )
  expect_equal(exact$mode, "raw_exact")
  expect_equal(exact$bounds, c(0.2, 0.4))
  expect_equal(exact$indices, 2:4)

  bracket <- scICER:::clamp_gamma_range_to_raw_plateau(
    gamma_sequence = c(0.1, 0.2, 0.3, 0.4),
    raw_cluster_medians = c(1, 2, 4, 5),
    target_clusters = 3L
  )
  expect_equal(bracket$mode, "raw_bracket")
  expect_equal(bracket$bounds, c(0.2, 0.3))
  expect_equal(bracket$indices, c(2L, 3L))

  near_target <- scICER:::clamp_gamma_range_to_raw_plateau(
    gamma_sequence = c(0.1, 0.2, 0.3, 0.4),
    raw_cluster_medians = c(4, 5, 5, 5),
    target_clusters = 6L
  )
  expect_equal(near_target$mode, "raw_near_target")
  expect_equal(near_target$bounds, c(0.2, 0.4))
  expect_equal(near_target$indices, 2:4)

  coarse <- scICER:::clamp_gamma_range_to_raw_plateau(
    gamma_sequence = c(0.1, 0.2, 0.3, 0.4),
    raw_cluster_medians = c(1, 1, 1, 1),
    target_clusters = 6L
  )
  expect_equal(coarse$mode, "coarse")
  expect_equal(coarse$bounds, c(0.1, 0.4))
})

test_that("select_gamma_admission prefers soft-guarded raw exact family first", {
  decision <- scICER:::select_gamma_admission(
    strict_flags = c(TRUE, FALSE, FALSE),
    relaxed_flags = c(TRUE, TRUE, FALSE),
    soft_guard_flags = c(TRUE, TRUE, TRUE),
    hard_guard_flags = c(TRUE, TRUE, TRUE),
    raw_strict_flags = c(FALSE, TRUE, FALSE),
    raw_relaxed_flags = c(FALSE, TRUE, TRUE)
  )

  expect_equal(decision$mode, "raw_strict_soft")
  expect_equal(decision$indices, 2L)
})

test_that("refine_gamma_candidates_by_raw_gap only prunes within the chosen family", {
  refined <- scICER:::refine_gamma_candidates_by_raw_gap(
    valid_indices = c(1L, 2L, 3L),
    admission_mode = "strict_soft",
    gamma_results = list(
      list(mean_clusters = 5, raw_median_gap = 1),
      list(mean_clusters = 5, raw_median_gap = 0),
      list(mean_clusters = 5, raw_median_gap = 0)
    ),
    target_clusters = 5L,
    min_cluster_size = 2L
  )

  expect_equal(refined$mode, "strict_soft")
  expect_equal(refined$indices, c(2L, 3L))
  expect_equal(refined$best_raw_gap, 0)
})

test_that("merge_small_clusters_to_neighbors prioritizes eligible target clusters", {
  labels <- c(0L, 0L, 0L, 0L, 0L, 1L, 2L, 2L)
  snn <- matrix(0, nrow = 8, ncol = 8)

  # Cell 6 (singleton) is closer to cluster 2, but cluster 2 is undersized.
  snn[6, 1:5] <- 1
  snn[1:5, 6] <- 1
  snn[6, 7:8] <- 10
  snn[7:8, 6] <- 10

  merged <- scICER:::merge_small_clusters_to_neighbors(
    labels = labels,
    snn_graph = Matrix::Matrix(snn, sparse = TRUE),
    min_cluster_size = 3L
  )

  expect_equal(merged[6], merged[1])
  final_sizes <- table(merged)
  expect_true(length(final_sizes) == 1L || min(final_sizes) >= 3L)
})

test_that("merge_small_clusters_to_neighbors uses deterministic smallest-id tie breaking", {
  labels <- c(0L, 0L, 0L, 1L, 1L, 1L, 2L)
  snn <- matrix(0, nrow = 7, ncol = 7)

  # Singleton (cell 7) ties exactly between cluster 0 and 1.
  snn[7, 1:3] <- 2
  snn[1:3, 7] <- 2
  snn[7, 4:6] <- 2
  snn[4:6, 7] <- 2

  graph <- Matrix::Matrix(snn, sparse = TRUE)
  set.seed(999)
  merged_a <- scICER:::merge_small_clusters_to_neighbors(labels, graph, min_cluster_size = 2L)
  set.seed(1234)
  merged_b <- scICER:::merge_small_clusters_to_neighbors(labels, graph, min_cluster_size = 2L)

  expect_equal(merged_a, merged_b)
  expect_equal(merged_a[7], 0L)
})

test_that("merge_small_clusters_to_neighbors collapses all-small labels to one cluster", {
  labels <- c(0L, 0L, 1L, 1L)
  graph <- Matrix::Matrix(0, nrow = 4, ncol = 4, sparse = TRUE)

  merged <- scICER:::merge_small_clusters_to_neighbors(
    labels = labels,
    snn_graph = graph,
    min_cluster_size = 3L
  )

  expect_equal(length(unique(merged)), 1L)
  expect_true(all(merged == 0L))
})

test_that("scICE_clustering validates min_cluster_size and preserves legacy mode", {
  data("pbmc_small", package = "SeuratObject")

  expect_error(
    scICE_clustering(
      object = pbmc_small,
      graph_name = "RNA_snn",
      min_cluster_size = 0,
      n_workers = 1,
      n_trials = 1,
      n_bootstrap = 1,
      verbose = FALSE
    ),
    "min_cluster_size must be >= 1"
  )

  expect_error(
    scICE_clustering(
      object = pbmc_small,
      graph_name = "RNA_snn",
      min_cluster_size = 1.5,
      n_workers = 1,
      n_trials = 1,
      n_bootstrap = 1,
      verbose = FALSE
    ),
    "must be an integer"
  )

  result <- scICE_clustering(
    object = pbmc_small,
    graph_name = "RNA_snn",
    cluster_range = 2:3,
    n_workers = 1,
    n_trials = 2,
    n_bootstrap = 2,
    seed = 123,
    remove_threshold = Inf,
    min_cluster_size = 1,
    verbose = FALSE
  )

  expect_equal(result$min_cluster_size, 1L)
})

test_that("scICE_clustering keeps raw best labels while storing min_cluster_size", {
  data("pbmc_small", package = "SeuratObject")

  result <- scICE_clustering(
    object = pbmc_small,
    graph_name = "RNA_snn",
    cluster_range = 2:3,
    n_workers = 1,
    n_trials = 2,
    n_bootstrap = 2,
    seed = 123,
    remove_threshold = Inf,
    min_cluster_size = 2,
    verbose = FALSE
  )

  expect_equal(result$min_cluster_size, 2L)
  non_null_labels <- Filter(Negate(is.null), result$best_labels)
  expect_true(length(non_null_labels) > 0)
  expect_true(all(vapply(non_null_labels, is.integer, logical(1))))
})

test_that("scICE_clustering returns per-k selection diagnostics", {
  data("pbmc_small", package = "SeuratObject")

  result <- scICE_clustering(
    object = pbmc_small,
    graph_name = "RNA_snn",
    cluster_range = 2:3,
    n_workers = 1,
    n_trials = 2,
    n_bootstrap = 2,
    seed = 123,
    remove_threshold = Inf,
    min_cluster_size = 2,
    verbose = FALSE
  )

  expect_equal(length(result$effective_cluster_median), length(result$n_cluster))
  expect_equal(length(result$raw_cluster_median), length(result$n_cluster))
  expect_equal(length(result$admission_mode), length(result$n_cluster))
  expect_equal(length(result$best_labels_raw_cluster_count), length(result$n_cluster))
  expect_true(all(is.finite(result$effective_cluster_median)))
  expect_true(all(is.finite(result$raw_cluster_median)))
  expect_true(all(nzchar(result$admission_mode)))
  final_cluster_counts <- vapply(
    result$best_labels,
    function(labels) as.integer(length(unique(labels))),
    integer(1)
  )
  expect_true(all(result$best_labels_raw_cluster_count >= final_cluster_counts))
})

test_that("optimize_clustering preserves raw labels and returns merged best_labels", {
  ig <- igraph::make_ring(6)
  raw_labels <- c(0L, 0L, 0L, 1L, 2L, 2L)
  zero_graph <- Matrix::Matrix(0, nrow = 6, ncol = 6, sparse = TRUE)

  local_mocked_bindings(
    leiden_clustering = function(...) raw_labels,
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.2),
    objective_function = "modularity",
    n_trials = 2L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    snn_graph = zero_graph,
    min_cluster_size = 2L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  expect_false(is.null(result))
  expect_true(any(table(raw_labels) < 2L))
  expect_equal(length(result$labels$arr), 1L)
  expect_equal(result$labels$arr[[1]], raw_labels)
  expect_equal(result$effective_cluster_median, 2)
  expect_equal(result$raw_cluster_median, 3)
  expect_equal(result$admission_mode, "strict_soft")
  expect_equal(result$best_labels_raw_cluster_count, 3L)
  expect_equal(length(unique(result$best_labels)), 2L)
  expect_true(min(table(result$best_labels)) >= 2L)
})

test_that("optimize_clustering prefers strict admission when strict and relaxed both exist", {
  ig <- igraph::make_ring(4)
  label_count1 <- c(0L, 0L, 0L, 0L) # effective clusters (min=1): 1
  label_count2 <- c(0L, 0L, 1L, 1L) # effective clusters (min=1): 2
  label_count3 <- c(0L, 1L, 2L, 2L) # effective clusters (min=1): 3
  state <- new.env(parent = emptyenv())
  state$last_res_key <- NULL
  state$trial_idx <- 0L

  local_mocked_bindings(
    leiden_clustering = function(igraph_obj, resolution, objective_function,
                                 n_iterations, beta, initial_membership = NULL) {
      res_key <- sprintf("%.6f", resolution)
      if (is.null(state$last_res_key) || !identical(state$last_res_key, res_key)) {
        state$last_res_key <- res_key
        state$trial_idx <- 1L
      } else {
        state$trial_idx <- state$trial_idx + 1L
      }

      if (resolution < 0.15) {
        # strict-valid only for target=2: counts {1, 3} -> median_raw=2, median_int=2, no hit
        if (state$trial_idx == 1L) label_count1 else label_count3
      } else if (resolution < 0.25) {
        # relaxed-valid only: counts {2, 1} -> median_raw=1.5, median_int=1, has hit
        if (state$trial_idx == 1L) label_count2 else label_count1
      } else {
        # neither strict nor relaxed
        label_count1
      }
    },
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.3),
    objective_function = "modularity",
    n_trials = 2L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    min_cluster_size = 1L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  # strict-first should select strict-only gammas and ignore relaxed-only candidates.
  expect_false(is.null(result))
  expected_arr <- sort(c(
    paste(label_count1, collapse = ","),
    paste(label_count3, collapse = ",")
  ))
  observed_arr <- sort(vapply(result$labels$arr, paste, collapse = ",", character(1)))
  expect_equal(observed_arr, expected_arr)
})

test_that("optimize_clustering falls back to relaxed admission and computes IC from all trials", {
  ig <- igraph::make_ring(4)
  hit_labels <- c(0L, 0L, 1L, 1L)   # effective clusters (min=2): 2
  miss_labels <- c(0L, 1L, 2L, 3L)  # effective clusters (min=2): 0
  call_count <- 0L

  local_mocked_bindings(
    leiden_clustering = function(...) {
      call_count <<- call_count + 1L
      if ((call_count %% 2L) == 1L) {
        return(hit_labels)
      }
      miss_labels
    },
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.1),
    objective_function = "modularity",
    n_trials = 2L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    snn_graph = Matrix::Diagonal(4),
    min_cluster_size = 2L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  # Counts per gamma are {2, 0}: strict fails (median_int=1), relaxed passes.
  expect_false(is.null(result))
  expect_equal(length(result$labels$arr), 2L)
  expected_arr <- sort(c(
    paste(hit_labels, collapse = ","),
    paste(miss_labels, collapse = ",")
  ))
  observed_arr <- sort(vapply(result$labels$arr, paste, collapse = ",", character(1)))
  expect_equal(observed_arr, expected_arr)
})

test_that("optimize_clustering rejects pathological strict gamma and uses bounded relaxed fallback", {
  ig <- igraph::make_ring(40)
  pathological_labels <- c(0L, 0L, 1L, 1L, seq.int(2L, 37L))
  bounded_hit_labels <- c(rep(0L, 20L), rep(1L, 20L))
  bounded_miss_labels <- rep(0L, 40L)
  state <- new.env(parent = emptyenv())
  state$last_res_key <- NULL
  state$trial_idx <- 0L

  local_mocked_bindings(
    leiden_clustering = function(igraph_obj, resolution, objective_function,
                                 n_iterations, beta, initial_membership = NULL) {
      res_key <- sprintf("%.6f", resolution)
      if (is.null(state$last_res_key) || !identical(state$last_res_key, res_key)) {
        state$last_res_key <- res_key
        state$trial_idx <- 1L
      } else {
        state$trial_idx <- state$trial_idx + 1L
      }

      if (resolution < 0.15) {
        pathological_labels
      } else if (resolution < 0.25) {
        if (state$trial_idx == 1L) bounded_hit_labels else bounded_miss_labels
      } else {
        bounded_miss_labels
      }
    },
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.3),
    objective_function = "modularity",
    n_trials = 2L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    snn_graph = Matrix::Diagonal(40),
    min_cluster_size = 2L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  expect_false(is.null(result))
  expected_arr <- sort(c(
    paste(bounded_hit_labels, collapse = ","),
    paste(bounded_miss_labels, collapse = ",")
  ))
  observed_arr <- sort(vapply(result$labels$arr, paste, collapse = ",", character(1)))
  expect_equal(observed_arr, expected_arr)
})

test_that("optimize_clustering prefers raw exact family when it exists on the gamma grid", {
  ig <- igraph::make_ring(6)
  effective_strict_labels <- c(0L, 0L, 1L, 1L, 1L, 2L)
  raw_strict_labels <- c(0L, 0L, 0L, 0L, 0L, 1L)

  local_mocked_bindings(
    leiden_clustering = function(igraph_obj, resolution, objective_function,
                                 n_iterations, beta, initial_membership = NULL) {
      if (resolution < 0.15) {
        return(effective_strict_labels)
      }
      raw_strict_labels
    },
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.2),
    objective_function = "modularity",
    n_trials = 2L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    snn_graph = Matrix::Diagonal(6),
    min_cluster_size = 2L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  expect_false(is.null(result))
  expect_equal(length(result$labels$arr), 1L)
  expect_equal(result$labels$arr[[1]], raw_strict_labels)
  expect_equal(result$admission_mode, "raw_strict_soft")
  expect_equal(result$best_labels_raw_cluster_count, 2L)
  expect_true(min(table(result$best_labels)) >= 2L)
  expect_equal(length(unique(result$best_labels)), 1L)
})

test_that("optimize_clustering falls back to bounded raw-count admission when effective admission is empty", {
  ig <- igraph::make_ring(6)
  raw_strict_labels <- c(0L, 0L, 0L, 0L, 0L, 1L)

  local_mocked_bindings(
    leiden_clustering = function(...) raw_strict_labels,
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.1),
    objective_function = "modularity",
    n_trials = 2L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    snn_graph = Matrix::Diagonal(6),
    min_cluster_size = 2L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  expect_false(is.null(result))
  expect_equal(length(result$labels$arr), 1L)
  expect_equal(result$labels$arr[[1]], raw_strict_labels)
})

test_that("optimize_clustering accepts strict-only gamma even when no trial hits target", {
  ig <- igraph::make_ring(4)
  label_count1 <- c(0L, 0L, 0L, 0L) # effective clusters (min=1): 1
  label_count3 <- c(0L, 1L, 2L, 2L) # effective clusters (min=1): 3
  state <- new.env(parent = emptyenv())
  state$trial_idx <- 0L

  local_mocked_bindings(
    leiden_clustering = function(igraph_obj, resolution, objective_function,
                                 n_iterations, beta, initial_membership = NULL) {
      state$trial_idx <- state$trial_idx + 1L
      if ((state$trial_idx %% 2L) == 1L) {
        return(label_count1)
      }
      label_count3
    },
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.1),
    objective_function = "modularity",
    n_trials = 2L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    min_cluster_size = 1L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  # Counts are {1, 3}: strict passes (median_int=2), relaxed fails (no hit).
  expect_false(is.null(result))
  expected_arr <- sort(c(
    paste(label_count1, collapse = ","),
    paste(label_count3, collapse = ",")
  ))
  observed_arr <- sort(vapply(result$labels$arr, paste, collapse = ",", character(1)))
  expect_equal(observed_arr, expected_arr)
})

test_that("optimize_clustering returns NULL when strict and relaxed admissions both fail", {
  ig <- igraph::make_ring(4)
  hit_labels <- c(0L, 0L, 1L, 1L)   # effective clusters (min=2): 2
  miss_labels <- c(0L, 1L, 2L, 3L)  # effective clusters (min=2): 0
  call_count <- 0L

  local_mocked_bindings(
    leiden_clustering = function(...) {
      call_count <<- call_count + 1L
      if (call_count == 1L) {
        return(hit_labels)
      }
      miss_labels
    },
    .package = "scICER"
  )

  result <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 2L,
    gamma_range = c(0.1, 0.1),
    objective_function = "modularity",
    n_trials = 3L,
    n_bootstrap = 2L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    snn_graph = Matrix::Diagonal(4),
    min_cluster_size = 2L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )

  # Effective counts are {2, 0, 0}: strict fails (median_int=0), relaxed fails (gap>1).
  expect_null(result)
})

test_that("plot_ic show_gamma moves gamma labels onto rotated x-axis ticks", {
  scice_results <- list(
    gamma = c(0.5, 0.8),
    labels = list(NULL, NULL),
    ic = c(1.001, 1.02),
    ic_vec = list(c(1.001, 1.002), c(1.02, 1.03)),
    n_cluster = c(2L, 3L),
    best_labels = list(NULL, NULL),
    n_iter = c(10L, 10L),
    mei = list(NULL, NULL),
    k = c(10L, 10L),
    excluded = c(FALSE, FALSE),
    exclusion_reason = c("none", "none")
  )
  class(scice_results) <- "scICE"

  plot_with_gamma <- plot_ic(scice_results, show_gamma = TRUE)
  x_scale_with_gamma <- plot_with_gamma$scales$get_scales("x")
  expect_identical(x_scale_with_gamma$labels[["2"]], "2\n5.00e-01")
  expect_identical(x_scale_with_gamma$labels[["3"]], "3\n8.00e-01")
  expect_identical(plot_with_gamma$theme$axis.text.x$angle, 45)

  plot_without_gamma <- plot_ic(scice_results, show_gamma = FALSE)
  x_scale_without_gamma <- plot_without_gamma$scales$get_scales("x")
  expect_identical(x_scale_without_gamma$labels[["2"]], "2")
  expect_identical(x_scale_without_gamma$labels[["3"]], "3")
  expect_null(plot_without_gamma$theme$axis.text.x$angle)
})

test_that("find_resolution_ranges reuses one representative preliminary clustering per gamma step", {
  ig <- igraph::make_ring(4)
  call_count <- 0L
  
  local_mocked_bindings(
    cached_leiden_clustering = function(...) {
      call_count <<- call_count + 1L
      c(0L, 0L, 1L, 1L)
    },
    .package = "scICER"
  )
  
  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 2L,
    start_g = 0,
    end_g = 0.2,
    objective_function = "modularity",
    resolution_tolerance = 1.0,
    n_workers = 1,
    verbose = FALSE
  )
  
  expect_true("2" %in% names(ranges))
  expect_equal(call_count, 1L)
})

test_that("find_resolution_ranges target matching uses effective cluster count", {
  ig <- igraph::make_ring(4)
  call_count <- 0L

  local_mocked_bindings(
    cached_leiden_clustering = function(...) {
      call_count <<- call_count + 1L
      c(0L, 0L, 0L, 1L)
    },
    .package = "scICER"
  )

  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 1L,
    start_g = 0,
    end_g = 0.2,
    objective_function = "modularity",
    resolution_tolerance = 1.0,
    n_workers = 1,
    verbose = FALSE,
    snn_graph = Matrix::Diagonal(4),
    min_cluster_size = 2L
  )

  expect_true("1" %in% names(ranges))
  expect_equal(call_count, 12L)
})

test_that("find_resolution_ranges CPM fallback keeps degenerate ranges on the gamma scale", {
  ig <- igraph::make_ring(6)

  local_mocked_bindings(
    cached_leiden_clustering = function(...) c(0L, 0L, 1L, 1L, 2L, 2L),
    .package = "scICER"
  )

  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 3L,
    start_g = log(1e-6),
    end_g = log(1),
    objective_function = "CPM",
    resolution_tolerance = 1,
    n_workers = 1,
    verbose = FALSE
  )

  expect_true("3" %in% names(ranges))
  bounds <- ranges[["3"]]
  expect_true(all(is.finite(bounds)))
  expect_true(all(bounds > 0))
  expect_true(bounds[1] <= bounds[2])
  expect_lt(bounds[2], 1)
})

test_that("find_resolution_ranges avoids drifting into all-small high-gamma regime", {
  ig <- igraph::make_ring(8)

  local_mocked_bindings(
    cached_leiden_clustering = function(igraph_obj, resolution, ...) {
      if (resolution < 30) {
        return(c(rep(0L, 4), rep(1L, 4)))
      }
      if (resolution <= 60) {
        return(rep(0:3, each = 2))
      }
      0:7
    },
    .package = "scICER"
  )

  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 4L,
    start_g = 0,
    end_g = 100,
    objective_function = "modularity",
    resolution_tolerance = 1.0,
    n_workers = 1,
    verbose = FALSE,
    snn_graph = Matrix::Diagonal(8),
    min_cluster_size = 2L
  )

  expect_true("4" %in% names(ranges))
  bounds <- ranges[["4"]]
  expect_true(is.numeric(bounds))
  expect_true(length(bounds) == 2L)
  expect_true(bounds[1] < bounds[2])
  expect_lt(bounds[2], 80)
})

test_that("find_resolution_ranges keeps representative preliminary summaries when no hit exists", {
  ig <- igraph::make_ring(4)
  call_count <- 0L
  
  local_mocked_bindings(
    cached_leiden_clustering = function(...) {
      call_count <<- call_count + 1L
      c(0L, 1L, 2L, 3L)
    },
    .package = "scICER"
  )
  
  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 2L,
    start_g = 0,
    end_g = 0.2,
    objective_function = "modularity",
    resolution_tolerance = 1.0,
    n_workers = 1,
    verbose = FALSE
  )
  
  # The nominal preliminary trial set is summarized from one representative clustering.
  expect_true("2" %in% names(ranges))
  expect_equal(call_count, 1L)
})

test_that("parallel resolution search still reuses one representative preliminary clustering per gamma step", {
  skip_on_os("windows")
  
  ig <- igraph::make_ring(4)
  trial_log <- tempfile(pattern = "scicer_pretrial_")
  
  local_mocked_bindings(
    cached_leiden_clustering = function(..., cache_key_suffix = "") {
      cat(sprintf("start:%s\n", cache_key_suffix), file = trial_log, append = TRUE)
      Sys.sleep(0.01)
      cat(sprintf("finish:%s\n", cache_key_suffix), file = trial_log, append = TRUE)
      c(0L, 1L, 2L, 3L)
    },
    .package = "scICER"
  )
  
  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 2L,
    start_g = 0,
    end_g = 0.2,
    objective_function = "modularity",
    resolution_tolerance = 1.0,
    n_workers = 4,
    verbose = FALSE
  )
  
  lines <- if (file.exists(trial_log)) readLines(trial_log, warn = FALSE) else character(0)
  start_lines <- grep("^start:", lines, value = TRUE)
  finish_lines <- grep("^finish:", lines, value = TRUE)
  
  expect_true("2" %in% names(ranges))
  expect_equal(length(start_lines), 1L)
  expect_equal(length(finish_lines), 1L)
  expect_true(all(grepl("^start:res_search_lower_2$", start_lines)))
})

test_that("find_resolution_ranges uses serial preliminary trials when per-target worker capacity is low", {
  skip_on_os("windows")

  ig <- igraph::make_ring(8)
  captured_messages <- character()
  withCallingHandlers(
    scICER:::find_resolution_ranges(
      igraph_obj = ig,
      cluster_range = 1:4,
      start_g = 0,
      end_g = 0.2,
      objective_function = "modularity",
      resolution_tolerance = 1.0,
      n_workers = 8,
      verbose = TRUE
    ),
    message = function(m) {
      captured_messages <<- c(captured_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl("Preliminary trial mode: serial", captured_messages, fixed = TRUE)))
  expect_true(any(grepl("Preliminary trial workers per gamma: 1", captured_messages, fixed = TRUE)))
})
