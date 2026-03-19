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

test_that("calculate_ecs scalar mode equals the mean element-wise score", {
  set.seed(1)
  n <- 100000
  cluster_a <- sample(1:10, n, replace = TRUE)
  cluster_b <- sample(1:10, n, replace = TRUE)

  scalar_score <- calculate_ecs(cluster_a, cluster_b, return_vector = FALSE)
  expected_scalar <- mean(calculate_ecs(cluster_a, cluster_b, return_vector = TRUE))
  expect_equal(scalar_score, expected_scalar, tolerance = 1e-10)
  expect_true(is.finite(scalar_score))
  expect_true(scalar_score >= 0 && scalar_score <= 1)
})

test_that("calculate_ecs scalar mode stays stable for very large inputs", {
  set.seed(1)
  n <- 250000
  cluster_a <- sample(1:10, n, replace = TRUE)
  cluster_b <- sample(1:10, n, replace = TRUE)

  scalar_score <- calculate_ecs(cluster_a, cluster_b, return_vector = FALSE)
  vector_score <- mean(calculate_ecs(cluster_a, cluster_b, return_vector = TRUE))

  expect_equal(scalar_score, vector_score, tolerance = 1e-10)
  expect_true(is.finite(scalar_score))
  expect_true(scalar_score >= 0 && scalar_score <= 1)
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
  expect_equal(result$analysis_mode, "cluster_range")
  expect_null(result$resolution_input)
  expect_null(result$resolution_diagnostics)
  expect_true(is.finite(result$best_cluster))
  expect_true(is.finite(result$best_resolution))
})

test_that("scICE_clustering supports manual resolution mode", {
  data("pbmc_small", package = "SeuratObject")

  result <- scICE_clustering(
    object = pbmc_small,
    graph_name = "RNA_snn",
    resolution = c(0.2, 0.4),
    n_workers = 1,
    n_trials = 2,
    n_bootstrap = 2,
    seed = 123,
    min_cluster_size = 1,
    verbose = FALSE
  )

  expect_equal(result$analysis_mode, "resolution")
  expect_equal(result$resolution_input, c(0.2, 0.4))
  expect_true(is.data.frame(result$resolution_diagnostics))
  expect_true(all(c(
    "resolution", "cluster_number", "ic", "effective_cluster_median",
    "raw_cluster_median", "best_labels_raw_cluster_count",
    "best_labels_final_cluster_count", "n_iter", "selected"
  ) %in% names(result$resolution_diagnostics)))
  expect_equal(result$cluster_range_tested, result$n_cluster)
  expect_true(all(result$gamma %in% result$resolution_input))
  expect_true(is.finite(result$best_cluster))
  expect_true(is.finite(result$best_resolution))

  plot_obj <- plot_ic(result)
  expect_s3_class(plot_obj, "ggplot")

  labels_df <- get_robust_labels(result, threshold = Inf)
  expect_true(is.data.frame(labels_df))
  expect_gte(ncol(labels_df), 2L)
})

test_that("scICE_clustering manual resolution mode ignores cluster_range and deduplicates by best IC", {
  data("pbmc_small", package = "SeuratObject")
  n_cells <- ncol(pbmc_small)

  build_mock_resolution_result <- function(resolution, final_clusters, ic_value) {
    labels <- rep(seq_len(final_clusters) - 1L, length.out = n_cells)
    list(
      gamma = resolution,
      labels = list(arr = list(as.integer(labels)), parr = 1),
      ic_median = ic_value,
      ic_bootstrap = c(ic_value, ic_value),
      best_labels = as.integer(labels),
      effective_cluster_median = as.numeric(final_clusters),
      raw_cluster_median = as.numeric(final_clusters),
      admission_mode = "manual_resolution",
      best_labels_raw_cluster_count = as.integer(final_clusters),
      best_labels_final_cluster_count = as.integer(final_clusters),
      n_iterations = 10L,
      k = 10L
    )
  }

  local_mocked_bindings(
    evaluate_fixed_resolution = function(igraph_obj, resolution, ...) {
      if (isTRUE(all.equal(resolution, 0.1))) {
        return(build_mock_resolution_result(0.1, 2L, 1.05))
      }
      if (isTRUE(all.equal(resolution, 0.2))) {
        return(build_mock_resolution_result(0.2, 2L, 1.01))
      }
      build_mock_resolution_result(0.3, 3L, 1.02)
    },
    .package = "scICER"
  )

  expect_message(
    result <- scICE_clustering(
      object = pbmc_small,
      graph_name = "RNA_snn",
      cluster_range = 2:4,
      resolution = c(0.1, 0.2, 0.3),
      n_workers = 1,
      n_trials = 1,
      n_bootstrap = 2,
      seed = 123,
      min_cluster_size = 1,
      verbose = FALSE
    ),
    "ignoring `cluster_range`"
  )

  expect_equal(result$analysis_mode, "resolution")
  expect_equal(result$n_cluster, c(2L, 3L))
  expect_equal(result$gamma, c(0.2, 0.3))
  expect_equal(nrow(result$resolution_diagnostics), 3L)
  expect_identical(result$resolution_diagnostics$selected, c(FALSE, TRUE, TRUE))
  expect_equal(result$best_cluster, 2L)
  expect_equal(result$best_resolution, 0.2)
})

test_that("scICE_clustering keeps raw best labels while storing min_cluster_size", {
  data("pbmc_small", package = "SeuratObject")
  n_cells <- ncol(pbmc_small)
  mock_target_results <- function(searched_targets, final_clusters, ic_values,
                                  gamma_values, raw_clusters, excluded,
                                  exclusion_reason) {
    best_labels <- lapply(seq_along(searched_targets), function(i) {
      if (isTRUE(excluded[[i]])) {
        return(NULL)
      }
      as.integer(rep(seq_len(final_clusters[[i]]) - 1L, length.out = n_cells))
    })
    raw_labels <- lapply(best_labels, function(labels) {
      if (is.null(labels)) {
        return(NULL)
      }
      list(arr = list(labels), parr = 1)
    })
    data.table::data.table(
      cluster_number = as.integer(searched_targets),
      gamma = as.numeric(gamma_values),
      labels = raw_labels,
      ic = as.numeric(ic_values),
      ic_vec = lapply(ic_values, function(x) c(x, x + 0.001)),
      best_labels = best_labels,
      effective_cluster_median = ifelse(excluded, NA_real_, as.numeric(final_clusters)),
      raw_cluster_median = as.numeric(raw_clusters),
      admission_mode = rep("mock_admission", length(searched_targets)),
      best_labels_raw_cluster_count = as.integer(raw_clusters),
      best_labels_final_cluster_count = ifelse(excluded, NA_integer_, as.integer(final_clusters)),
      n_iter = rep(10L, length(searched_targets)),
      mei = lapply(seq_along(searched_targets), function(i) rep(0.5, n_cells)),
      k = as.integer(searched_targets),
      source_target_cluster = as.integer(searched_targets),
      excluded = as.logical(excluded),
      exclusion_reason = as.character(exclusion_reason),
      selected_main_result = rep(FALSE, length(searched_targets))
    )
  }

  local_mocked_bindings(
    find_resolution_ranges = function(...) {
      gamma_dict <- list("2" = c(0.2, 0.3))
      attr(gamma_dict, "resolution_search_diagnostics") <- data.frame(
        sweep_round = 1L,
        gamma = 0.25,
        effective_cluster_count = 2,
        raw_cluster_count = 2,
        final_cluster_count = 2,
        raw_class = "raw_in_band",
        over_fragmented = FALSE,
        selected_for_refinement = FALSE,
        selected_for_target_interval = TRUE,
        plateau_round = 0L,
        stringsAsFactors = FALSE
      )
      attr(gamma_dict, "coverage_complete") <- FALSE
      attr(gamma_dict, "plateau_stop") <- TRUE
      attr(gamma_dict, "uncovered_targets") <- 3L
      gamma_dict
    },
    clustering_main = function(...) {
      list(
        target_results = mock_target_results(
          searched_targets = c(2L, 3L),
          final_clusters = c(2L, NA_integer_),
          ic_values = c(1.001, NA_real_),
          gamma_values = c(0.25, NA_real_),
          raw_clusters = c(2L, NA_integer_),
          excluded = c(FALSE, TRUE),
          exclusion_reason = c("none", "resolution_search_failed")
        )
      )
    },
    .package = "scICER"
  )

  suppressWarnings(
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
  )

  expect_equal(result$min_cluster_size, 2L)
  non_null_labels <- Filter(Negate(is.null), result$best_labels)
  expect_true(length(non_null_labels) > 0)
  expect_true(all(vapply(non_null_labels, is.integer, logical(1))))
})

test_that("scICE_clustering returns per-k selection diagnostics", {
  data("pbmc_small", package = "SeuratObject")
  n_cells <- ncol(pbmc_small)
  mock_target_results <- function(searched_targets, final_clusters, ic_values,
                                  gamma_values, raw_clusters) {
    best_labels <- lapply(seq_along(searched_targets), function(i) {
      as.integer(rep(seq_len(final_clusters[[i]]) - 1L, length.out = n_cells))
    })
    raw_labels <- lapply(best_labels, function(labels) list(arr = list(labels), parr = 1))
    data.table::data.table(
      cluster_number = as.integer(searched_targets),
      gamma = as.numeric(gamma_values),
      labels = raw_labels,
      ic = as.numeric(ic_values),
      ic_vec = lapply(ic_values, function(x) c(x, x + 0.001)),
      best_labels = best_labels,
      effective_cluster_median = as.numeric(final_clusters),
      raw_cluster_median = as.numeric(raw_clusters),
      admission_mode = rep("mock_admission", length(searched_targets)),
      best_labels_raw_cluster_count = as.integer(raw_clusters),
      best_labels_final_cluster_count = as.integer(final_clusters),
      n_iter = rep(10L, length(searched_targets)),
      mei = lapply(seq_along(searched_targets), function(i) rep(0.5, n_cells)),
      k = as.integer(searched_targets),
      source_target_cluster = as.integer(searched_targets),
      excluded = rep(FALSE, length(searched_targets)),
      exclusion_reason = rep("none", length(searched_targets)),
      selected_main_result = rep(FALSE, length(searched_targets))
    )
  }

  local_mocked_bindings(
    find_resolution_ranges = function(...) {
      gamma_dict <- list(
        "2" = c(0.2, 0.3),
        "3" = c(0.3, 0.4)
      )
      attr(gamma_dict, "resolution_search_diagnostics") <- data.frame(
        sweep_round = c(1L, 1L),
        gamma = c(0.25, 0.35),
        effective_cluster_count = c(2, 3),
        raw_cluster_count = c(2, 3),
        final_cluster_count = c(2, 3),
        raw_class = c("raw_in_band", "raw_in_band"),
        over_fragmented = c(FALSE, FALSE),
        selected_for_refinement = c(FALSE, FALSE),
        selected_for_target_interval = c(TRUE, TRUE),
        plateau_round = c(0L, 0L),
        stringsAsFactors = FALSE
      )
      attr(gamma_dict, "coverage_complete") <- TRUE
      attr(gamma_dict, "plateau_stop") <- FALSE
      attr(gamma_dict, "uncovered_targets") <- integer(0)
      gamma_dict
    },
    clustering_main = function(...) {
      list(
        target_results = mock_target_results(
          searched_targets = c(2L, 3L),
          final_clusters = c(2L, 3L),
          ic_values = c(1.001, 1.002),
          gamma_values = c(0.25, 0.35),
          raw_clusters = c(2L, 3L)
        )
      )
    },
    .package = "scICER"
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
  expect_equal(result$target_diagnostics$gamma_left, c(0.2, 0.3))
  expect_equal(result$target_diagnostics$gamma_right, c(0.3, 0.4))
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

test_that("optimize_clustering can fall back to raw-count admission when final-count admission fails", {
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

  expect_type(result, "list")
  expect_true("best_labels_final_cluster_count" %in% names(result))
})

test_that("optimize_clustering reuses shared-search gamma seeds during phase 1", {
  ig <- igraph::make_ring(6)
  labels_2 <- c(0L, 0L, 0L, 1L, 1L, 1L)
  labels_3 <- c(0L, 0L, 1L, 1L, 2L, 2L)
  labels_4 <- c(0L, 0L, 1L, 2L, 3L, 3L)

  local_mocked_bindings(
    leiden_clustering = function(igraph_obj, resolution, objective_function,
                                 n_iterations, beta, initial_membership = NULL) {
      if (abs(resolution - 0.05) < 1e-8) {
        return(labels_3)
      }
      if (resolution < 0.05) {
        return(labels_2)
      }
      labels_4
    },
    .package = "scICER"
  )

  result_without_seed <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 3L,
    gamma_range = c(0, 0.2),
    objective_function = "modularity",
    n_trials = 1L,
    n_bootstrap = 1L,
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
  expect_false(isTRUE(result_without_seed$success))
  expect_identical(result_without_seed$failure_reason, "optimization_admission_failed")

  result_with_seed <- scICER:::optimize_clustering(
    igraph_obj = ig,
    target_clusters = 3L,
    gamma_range = c(0, 0.2),
    objective_function = "modularity",
    n_trials = 1L,
    n_bootstrap = 1L,
    seed = 123,
    beta = 0.1,
    n_iterations = 1L,
    max_iterations = 3L,
    resolution_tolerance = 1e-3,
    n_workers = 1L,
    gamma_seed_values = 0.05,
    min_cluster_size = 1L,
    verbose = FALSE,
    worker_id = "TEST",
    in_parallel_context = FALSE
  )
  expect_true(isTRUE(result_with_seed$success))
  expect_equal(result_with_seed$best_labels_final_cluster_count, 3L)
})

test_that("build_optimization_gamma_batches respects budget and preserves anchor seeds", {
  seed_table <- data.frame(
    gamma = c(1e-6, 2e-6, 3e-6, 4e-6, 5e-6, 6e-6, 7e-6, 8e-6, 9e-6, 1e-5),
    seed_role = c("left", "selected", "exact", "near", "near", "near", "seed", "seed", "seed", "right"),
    final_cluster_count = c(10, 11, 12, 11, 13, 12, 10, 11, 12, 14),
    raw_cluster_count = c(10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
  )

  batches <- scICER:::build_optimization_gamma_batches(
    gamma_range = c(1e-6, 1e-5),
    gamma_seed_values = seed_table,
    target_clusters = 12L,
    objective_function = "CPM",
    resolution_tolerance = 1e-8,
    n_vertices = 1000L
  )

  expect_lte(length(batches$primary_gammas), 8L)
  expect_lte(length(c(batches$primary_gammas, batches$secondary_gammas)), 12L)
  expect_true(1e-6 %in% batches$primary_gammas)
  expect_true(2e-6 %in% batches$primary_gammas)
  expect_true(1e-5 %in% batches$primary_gammas)
  expect_true(3e-6 %in% batches$primary_gammas)
})

test_that("phase1 secondary expansion only runs for unresolved unguarded paths", {
  expect_false(scICER:::should_expand_phase1_secondary(
    valid_indices = 1:2,
    admission_mode = "strict_soft",
    exact_hit_gamma_count = 0L
  ))
  expect_false(scICER:::should_expand_phase1_secondary(
    valid_indices = 1:2,
    admission_mode = "relaxed_unguarded",
    exact_hit_gamma_count = 1L
  ))
  expect_true(scICER:::should_expand_phase1_secondary(
    valid_indices = integer(0),
    admission_mode = "none",
    exact_hit_gamma_count = 0L
  ))
  expect_true(scICER:::should_expand_phase1_secondary(
    valid_indices = 1L,
    admission_mode = "relaxed_unguarded",
    exact_hit_gamma_count = 0L
  ))
})

test_that("phase4 skip and iteration cap helpers enforce bounded refinement", {
  expect_true(scICER:::should_skip_phase4_refinement(
    candidate_count = 2L,
    best_ic = 1.005,
    exact_hit_gamma_count = 1L
  ))
  expect_false(scICER:::should_skip_phase4_refinement(
    candidate_count = 3L,
    best_ic = 1.001,
    exact_hit_gamma_count = 2L
  ))
  expect_equal(scICER:::phase4_iteration_cap_for_mode("relaxed_unguarded"), 2L)
  expect_equal(scICER:::phase4_iteration_cap_for_mode("strict_soft"), 3L)
})

test_that("finalize_selected_clustering prefers exact final-hit trials for best_labels", {
  best_matrix <- cbind(
    c(0L, 0L, 1L, 2L, 3L, 3L),
    c(0L, 0L, 1L, 1L, 2L, 2L),
    c(0L, 1L, 1L, 2L, 3L, 3L)
  )

  local_mocked_bindings(
    load_cluster_matrix = function(matrix_ref) best_matrix,
    release_cluster_matrix = function(matrix_ref) invisible(NULL),
    extract_clustering_array = function(clustering_matrix) {
      arr <- lapply(seq_len(ncol(clustering_matrix)), function(idx) clustering_matrix[, idx])
      list(arr = arr, parr = rep(1 / length(arr), length(arr)))
    },
    calculate_ic_from_extracted = function(extracted) 1,
    get_best_clustering = function(clustering_array) clustering_array$arr[[1]],
    .package = "scICER"
  )

  result <- scICER:::finalize_selected_clustering(
    matrix_ref = "fake_ref",
    gamma = 0.1,
    effective_cluster_median = 3,
    raw_cluster_median = 3,
    final_cluster_median = 3,
    admission_mode = "strict_soft",
    cluster_seed = 123,
    n_bootstrap = 2L,
    n_workers = 1L,
    target_clusters = 3L,
    min_cluster_size = 1L,
    verbose = FALSE
  )

  expect_equal(result$best_labels_final_cluster_count, 3L)
  expect_equal(length(unique(result$best_labels)), 3L)
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
  diagnostics <- attr(ranges, "resolution_search_diagnostics")
  expect_true(is.data.frame(diagnostics))
  expect_equal(call_count, nrow(diagnostics))
})

test_that("find_resolution_ranges target matching uses final merged cluster count", {
  ig <- igraph::make_ring(4)
  call_count <- 0L

  local_mocked_bindings(
    cached_leiden_clustering = function(...) {
      call_count <<- call_count + 1L
      c(0L, 0L, 0L, 1L)
    },
    merge_small_clusters_to_neighbors = function(...) c(0L, 0L, 0L, 0L),
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
  expect_gte(call_count, 1L)
})

test_that("find_resolution_ranges refines beyond coarse brackets when targets are not optimization-ready", {
  ig <- igraph::make_ring(6)
  call_count <- 0L
  labels_2 <- c(0L, 0L, 0L, 1L, 1L, 1L)
  labels_3 <- c(0L, 0L, 1L, 1L, 2L, 2L)
  labels_4 <- c(0L, 0L, 1L, 2L, 3L, 3L)

  local_mocked_bindings(
    cached_leiden_clustering = function(igraph_obj, resolution, ...) {
      call_count <<- call_count + 1L
      if (abs(resolution - 0.05) < 1e-8) {
        return(labels_3)
      }
      if (resolution < 0.05) {
        return(labels_2)
      }
      labels_4
    },
    .package = "scICER"
  )

  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 3L,
    start_g = 0,
    end_g = 0.2,
    objective_function = "modularity",
    resolution_tolerance = 1e-3,
    n_workers = 1,
    verbose = FALSE
  )

  diagnostics <- attr(ranges, "resolution_search_diagnostics")
  expect_true("3" %in% names(ranges))
  expect_gt(nrow(diagnostics), 11L)
  expect_true(isTRUE(attr(ranges, "coverage_complete")))
  expect_true(any(abs(diagnostics$gamma - 0.05) < 1e-8))
  seed_values <- attr(ranges, "target_gamma_seeds")[["3"]]
  expect_true(any(abs(seed_values - 0.05) < 1e-8))
  expect_equal(call_count, nrow(diagnostics))
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
  expect_lte(bounds[2], 100)
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
  
  expect_false("2" %in% names(ranges))
  expect_gte(call_count, 1L)
  expect_true(is.data.frame(attr(ranges, "resolution_search_diagnostics")))
  expect_equal(attr(ranges, "uncovered_targets"), 2L)
})

test_that("parallel resolution search uses shared global gamma cache keys", {
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
  
  expect_false("2" %in% names(ranges))
  expect_gte(length(start_lines), 1L)
  expect_equal(length(start_lines), length(finish_lines))
  expect_true(all(grepl("^start:global_resolution_search$", start_lines)))
  expect_equal(attr(ranges, "uncovered_targets"), 2L)
})

test_that("find_resolution_ranges uses serial preliminary trials when per-target worker capacity is low", {
  skip_on_os("windows")

  ig <- igraph::make_ring(8)
  captured_messages <- character()
  withCallingHandlers(
    scICER:::find_resolution_ranges(
      igraph_obj = ig,
      cluster_range = 1:2,
      start_g = 0,
      end_g = 0.2,
      objective_function = "modularity",
      resolution_tolerance = 1.0,
      n_workers = 2,
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

test_that("find_resolution_ranges CPM upper-cap discovery stops at high-gamma degeneracy", {
  ig <- igraph::make_ring(8)

  local_mocked_bindings(
    cached_leiden_clustering = function(igraph_obj, resolution, ...) {
      if (resolution < 1e-5) {
        return(c(rep(0L, 4), rep(1L, 4)))
      }
      0:7
    },
    .package = "scICER"
  )

  ranges <- suppressWarnings(scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 6L,
    start_g = log(1e-8),
    end_g = 20,
    objective_function = "CPM",
    resolution_tolerance = 1e-8,
    n_workers = 2,
    verbose = FALSE,
    snn_graph = Matrix::Diagonal(8),
    min_cluster_size = 2L
  ))

  diagnostics <- attr(ranges, "resolution_search_diagnostics")
  expect_equal(attr(ranges, "upper_cap_stop_reason"), "high_gamma_degenerate")
  expect_equal(attr(ranges, "discovered_upper_gamma"), 2.56e-06, tolerance = 1e-12)
  expect_equal(unique(diagnostics$probe_stage[diagnostics$probe_stage == "upper_cap_discovery"]),
               "upper_cap_discovery")
  expect_true(any(diagnostics$degenerate_high_gamma))
  expect_true(all(table(diagnostics$discovery_round[diagnostics$probe_stage == "upper_cap_discovery"]) <= 2L))
  expect_true(all(
    diagnostics$gamma[diagnostics$probe_stage == "coarse"] <= attr(ranges, "discovered_upper_gamma") + 1e-12
  ))
})

test_that("find_resolution_ranges CPM upper-cap discovery stops once requested max is covered", {
  ig <- igraph::make_ring(12)

  local_mocked_bindings(
    cached_leiden_clustering = function(igraph_obj, resolution, ...) {
      if (resolution < 1e-5) {
        return(rep(0:1, each = 6))
      }
      if (resolution < 1e-4) {
        return(rep(0:5, each = 2))
      }
      0:11
    },
    .package = "scICER"
  )

  ranges <- suppressWarnings(scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 6L,
    start_g = log(1e-8),
    end_g = 20,
    objective_function = "CPM",
    resolution_tolerance = 1e-8,
    n_workers = 6,
    verbose = FALSE
  ))

  expect_equal(attr(ranges, "upper_cap_stop_reason"), "target_covered")
  expect_equal(attr(ranges, "discovered_upper_gamma"), 1.024e-05, tolerance = 1e-12)
  diagnostics <- attr(ranges, "resolution_search_diagnostics")
  expect_equal(sum(diagnostics$probe_stage == "upper_cap_discovery"), 6L)
  expect_true(all(table(diagnostics$discovery_round[diagnostics$probe_stage == "upper_cap_discovery"]) <= 6L))
})

test_that("build_cpm_discovery_batch_gamma_values narrows each discovery batch", {
  batch_values <- scICER:::build_cpm_discovery_batch_gamma_values(
    current_gamma = 1e-08,
    hard_cap_gamma = exp(20),
    batch_size = 8L,
    step_ratio = 4
  )

  expect_length(batch_values, 8L)
  expect_equal(batch_values[[1]], 1e-08, tolerance = 1e-16)
  expect_equal(max(batch_values), 1.6384e-04, tolerance = 1e-12)
  expect_true(all(diff(log(batch_values)) > 0))
})

test_that("derive_cpm_discovery_batch_plan narrows near the requested max", {
  default_plan <- scICER:::derive_cpm_discovery_batch_plan(
    active_probe_workers = 40,
    requested_max = 20,
    frontier_final_cluster_count = 4
  )
  expect_equal(default_plan$batch_size, 6L)
  expect_equal(default_plan$step_ratio, 4)

  mid_plan <- scICER:::derive_cpm_discovery_batch_plan(
    active_probe_workers = 40,
    requested_max = 20,
    frontier_final_cluster_count = 10
  )
  expect_equal(mid_plan$batch_size, 3L)
  expect_equal(mid_plan$step_ratio, 2)

  frontier_plan <- scICER:::derive_cpm_discovery_batch_plan(
    active_probe_workers = 40,
    requested_max = 20,
    frontier_final_cluster_count = 15
  )
  expect_equal(frontier_plan$batch_size, 2L)
  expect_equal(frontier_plan$step_ratio, 2)

  near_target_plan <- scICER:::derive_cpm_discovery_batch_plan(
    active_probe_workers = 40,
    requested_max = 20,
    frontier_final_cluster_count = 19
  )
  expect_equal(near_target_plan$batch_size, 2L)
  expect_equal(near_target_plan$step_ratio, 1.5)
})

test_that("find_resolution_ranges CPM narrows discovery batches near requested max", {
  ig <- igraph::make_ring(16)
  captured_messages <- character()

  local_mocked_bindings(
    cached_leiden_clustering = function(igraph_obj, resolution, ...) {
      if (resolution < 1e-5) {
        return(rep(0:1, each = 8))
      }
      if (resolution < 2.048e-05) {
        return(rep(0:5, c(3, 3, 3, 3, 2, 2)))
      }
      if (resolution < 4.096e-05) {
        return(rep(0:7, each = 2))
      }
      return(rep(0:15, each = 1))
    },
    .package = "scICER"
  )

  withCallingHandlers(
    suppressWarnings(scICER:::find_resolution_ranges(
      igraph_obj = ig,
      cluster_range = 8L,
      start_g = log(1e-8),
      end_g = 20,
      objective_function = "CPM",
      resolution_tolerance = 1e-8,
      n_workers = 40,
      verbose = TRUE,
      min_cluster_size = 1L
    )),
    message = function(m) {
      captured_messages <<- c(captured_messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )

  expect_true(any(grepl(
    "Upper-cap discovery round 2 - batch size = 2 | step ratio = 2 | gamma values = 2.048e-05, 4.096e-05",
    captured_messages,
    fixed = TRUE
  )))
})

test_that("find_resolution_ranges oversubscribes the first coarse sweep only", {
  ig <- igraph::make_ring(6)
  call_count <- 0L
  scICER:::clear_clustering_cache()

  local_mocked_bindings(
    cached_leiden_clustering = function(igraph_obj, resolution, ...) {
      call_count <<- call_count + 1L
      if (resolution < 0.05) {
        return(c(0L, 0L, 0L, 1L, 1L, 1L))
      }
      if (resolution > 0.05) {
        return(c(0L, 0L, 1L, 2L, 3L, 3L))
      }
      c(0L, 0L, 1L, 1L, 2L, 2L)
    },
    .package = "scICER"
  )

  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = 3L,
    start_g = 0,
    end_g = 0.2,
    objective_function = "modularity",
    resolution_tolerance = 1e-3,
    n_workers = 10,
    verbose = FALSE
  )

  diagnostics <- attr(ranges, "resolution_search_diagnostics")
  expect_equal(attr(ranges, "coarse_probe_count"), 30L)
  expect_equal(sum(diagnostics$probe_stage == "coarse"), 30L)
  expect_true(any(diagnostics$probe_stage == "refinement"))
  expect_true(all(diagnostics$scheduled_probe_workers[diagnostics$probe_stage == "coarse"] == 10L))
  expect_gte(nrow(diagnostics), 31L)
})

test_that("build_refinement_probe_plan densifies narrow intervals when workers exceed interval count", {
  plan <- scICER:::build_refinement_probe_plan(
    unresolved_intervals = list(
      "9" = c(1e-06, 2e-06),
      "10" = c(8e-06, 1.2e-05)
    ),
    objective_function = "CPM",
    resolution_tolerance = 1e-08,
    active_probe_workers = 10L,
    existing_gamma_values = numeric(0)
  )

  probe_metadata <- plan$probe_metadata
  interval_summary <- plan$interval_summary

  expect_equal(nrow(interval_summary), 2L)
  expect_gt(nrow(probe_metadata), 2L)
  expect_lte(nrow(probe_metadata), 10L)
  expect_true(all(interval_summary$refinement_points_per_interval == 5L))
  expect_true(all(probe_metadata$refinement_points_per_interval == 5L))
  expect_true(all(probe_metadata$gamma > min(c(1e-06, 8e-06))))
})

test_that("find_resolution_ranges refinement emits more than one probe per interval when gamma narrows", {
  ig <- igraph::make_ring(6)

  local_mocked_bindings(
    cached_leiden_clustering = function(igraph_obj, resolution, ...) {
      if (resolution < 0.04) {
        return(c(0L, 0L, 0L, 1L, 1L, 1L))
      }
      if (resolution < 0.08) {
        return(c(0L, 0L, 1L, 1L, 2L, 2L))
      }
      c(0L, 1L, 2L, 3L, 4L, 5L)
    },
    .package = "scICER"
  )

  ranges <- scICER:::find_resolution_ranges(
    igraph_obj = ig,
    cluster_range = c(4L, 5L),
    start_g = 0,
    end_g = 0.2,
    objective_function = "modularity",
    resolution_tolerance = 1e-3,
    n_workers = 10,
    verbose = FALSE
  )

  diagnostics <- attr(ranges, "resolution_search_diagnostics")
  refinement_round2 <- diagnostics[diagnostics$probe_stage == "refinement" & diagnostics$sweep_round == 2L, ]
  expect_gt(nrow(refinement_round2), 2L)
  expect_true(all(is.finite(refinement_round2$refinement_interval_width)))
  expect_true(all(refinement_round2$refinement_points_per_interval >= 1L))
  expect_true(any(refinement_round2$refinement_points_per_interval > 1L))
})
