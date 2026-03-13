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
  expect_true(any(grepl("RESOLUTION_SEARCH: k =", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Preliminary trial strategy - early stop", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("lower-bound progress", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("trial summary - early stops:", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Phase 1 progress", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Phase 1 completed", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Phase 5 progress", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Optimization scheduling - dynamic worker queue", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Long step started - extracting sparse entries from matrix slots", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("Long step started - constructing igraph edge structure", captured$messages, fixed = TRUE)))
  expect_true(any(grepl("ANALYSIS COMPLETE", captured$messages, fixed = TRUE)))
})
