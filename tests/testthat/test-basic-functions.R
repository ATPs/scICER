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

test_that("scICE_clustering enforces min_cluster_size on best labels", {
  data("pbmc_small", package = "SeuratObject")

  result <- scICE_clustering(
    object = pbmc_small,
    graph_name = "RNA_snn",
    cluster_range = 0L,
    n_workers = 1,
    n_trials = 2,
    n_bootstrap = 2,
    seed = 123,
    remove_threshold = Inf,
    min_cluster_size = 2,
    verbose = FALSE
  )

  non_null_labels <- Filter(Negate(is.null), result$best_labels)
  expect_true(length(non_null_labels) > 0)

  for (label_vec in non_null_labels) {
    cluster_sizes <- table(label_vec)
    expect_true(length(cluster_sizes) == 1L || min(cluster_sizes) >= 2L)
  }
})

test_that("optimize_clustering keeps raw labels and applies final-only merge", {
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

  expected_best <- scICER:::merge_small_clusters_to_neighbors(
    labels = raw_labels,
    snn_graph = zero_graph,
    min_cluster_size = 2L
  )
  expect_equal(result$best_labels, expected_best)
  expect_true(all(table(result$best_labels) >= 2L))
})

test_that("optimize_clustering admits gamma with any-hit plus median-window and computes IC from all trials", {
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

  # Median effective count per gamma is 1 (from {2, 0}); this satisfies
  # |median-target| <= 1 and has a hit, so gamma is admitted.
  expect_false(is.null(result))
  expect_equal(length(result$labels$arr), 2L)
  expected_arr <- sort(c(
    paste(hit_labels, collapse = ","),
    paste(miss_labels, collapse = ",")
  ))
  observed_arr <- sort(vapply(result$labels$arr, paste, collapse = ",", character(1)))
  expect_equal(observed_arr, expected_arr)
})

test_that("optimize_clustering rejects gamma when median-window fails even with a hit trial", {
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

  # Effective counts are {2, 0, 0}: has hit, but median is 0 so |0-2|=2 > 1.
  expect_null(result)
})

test_that("plot_ic show_gamma controls selected-gamma subtitle", {
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
  expect_true(grepl("Selected gamma:", plot_with_gamma$labels$subtitle, fixed = TRUE))
  expect_true(grepl("k=2:", plot_with_gamma$labels$subtitle, fixed = TRUE))

  plot_without_gamma <- plot_ic(scice_results, show_gamma = FALSE)
  expect_false(grepl("Selected gamma:", plot_without_gamma$labels$subtitle, fixed = TRUE))
})

test_that("find_resolution_ranges stops early when a preliminary trial hits target", {
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
  expect_equal(call_count, 1L)
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

test_that("find_resolution_ranges falls back to full preliminary trials when no hit exists", {
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
  
  # Small graphs use 15 preliminary trials per step.
  expect_true("2" %in% names(ranges))
  expect_equal(call_count, 15L)
})

test_that("parallel preliminary trials avoid launching all trials after early hit", {
  skip_on_os("windows")
  
  ig <- igraph::make_ring(4)
  trial_log <- tempfile(pattern = "scicer_pretrial_")
  
  local_mocked_bindings(
    cached_leiden_clustering = function(..., cache_key_suffix = "") {
      trial_idx <- suppressWarnings(as.integer(sub(".*_trial_([0-9]+)$", "\\1", cache_key_suffix)))
      if (is.na(trial_idx)) {
        trial_idx <- 1L
      }
      cat(sprintf("start:%d\n", trial_idx), file = trial_log, append = TRUE)
      if (trial_idx == 2L) {
        Sys.sleep(0.05)
        cat(sprintf("finish:%d\n", trial_idx), file = trial_log, append = TRUE)
        return(c(0L, 0L, 1L, 1L))
      }
      Sys.sleep(1)
      cat(sprintf("finish:%d\n", trial_idx), file = trial_log, append = TRUE)
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
  
  expect_true("2" %in% names(ranges))
  expect_true(any(grepl("^start:2$", start_lines)))
  expect_true(length(start_lines) >= 2L)
  expect_lte(length(start_lines), 4L)
})
