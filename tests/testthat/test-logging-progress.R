capture_verbose_messages <- function(expr) {
  messages <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  list(value = value, messages = messages)
}

test_that("scICE_clustering verbose emits startup, checkpoint, and completion logs", {
  data("pbmc_small", package = "SeuratObject")
  
  captured <- capture_verbose_messages(
    scICE_clustering(
      object = pbmc_small,
      graph_name = "RNA_snn",
      cluster_range = 2:3,
      n_workers = 1,
      n_trials = 2,
      n_bootstrap = 5,
      seed = 123,
      remove_threshold = Inf,
      min_cluster_size = 1,
      verbose = TRUE
    )
  )
  
  expect_s3_class(captured$value, "scICE")
  expect_true(any(grepl("Starting scICE clustering analysis", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("RESOLUTION_SEARCH: Sweep round", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Preliminary trial strategy - reuse one representative preliminary clustering per shared gamma probe", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Initial shared sweep covered final merged cluster counts up to", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("shared-sweep bounds", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Phase 1 progress", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Phase 1 completed", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Selected diagnostics - gamma =", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Per-k selection diagnostics", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Phase 5 progress", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Optimization scheduling - dynamic worker queue", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Long step started - extracting sparse entries from matrix slots", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Long step started - constructing igraph edge structure", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("ANALYSIS COMPLETE", captured$messages, fixed = TRUE)))
})

test_that("scICE_clustering verbose logs manual resolution mode and ignored cluster_range", {
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
      if (isTRUE(all.equal(resolution, 0.2))) {
        return(build_mock_resolution_result(0.2, 2L, 1.01))
      }
      build_mock_resolution_result(0.4, 3L, 1.02)
    },
    .package = "scICER"
  )

  captured <- capture_verbose_messages(
    scICE_clustering(
      object = pbmc_small,
      graph_name = "RNA_snn",
      cluster_range = 2:4,
      resolution = c(0.2, 0.4),
      n_workers = 1,
      n_trials = 1,
      n_bootstrap = 2,
      seed = 123,
      min_cluster_size = 1,
      verbose = TRUE
    )
  )

  expect_s3_class(captured$value, "scICE")
  expect_true(any(grepl("ignoring `cluster_range`", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Manual resolution mode selected - skipping cluster_range search", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("CLUSTERING_MAIN: Using manual resolution mode", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Manual resolutions retained after per-cluster IC selection", captured$messages, fixed = TRUE)))
})
