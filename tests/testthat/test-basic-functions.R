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
  
  expect_length(mei_scores, 1)  # One score per unique clustering
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