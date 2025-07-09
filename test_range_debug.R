library(scICER)
library(Seurat)

# Set up parallel processing
n_workers <- 1  # Use single worker to simplify debugging

# Load and prepare Seurat object
data(pbmc_small)
seurat_obj <- pbmc_small
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)  
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Test different cluster ranges to find where the error occurs
test_ranges <- list(
  c(3, 5),      # Small range that works
  c(3, 8),      # Medium range
  c(3, 12),     # Larger range
  c(3, 15),     # Even larger
  c(3, 20)      # Full range that fails
)

for (i in seq_along(test_ranges)) {
  range_test <- test_ranges[[i]]
  cat("\n=== Testing cluster range:", range_test[1], "to", range_test[2], "===\n")
  
  result <- try({
    scice_results <- scICE_clustering(
      object = seurat_obj,
      cluster_range = range_test[1]:range_test[2],
      n_workers = n_workers,
      remove_threshold = Inf,
      verbose = FALSE  # Reduce output
    )
    cat("SUCCESS: Range", range_test[1], "to", range_test[2], "completed\n")
    cat("Consistent clusters found:", length(scice_results$consistent_clusters), "\n")
  }, silent = TRUE)
  
  if (inherits(result, "try-error")) {
    cat("ERROR at range", range_test[1], "to", range_test[2], ":\n")
    cat(as.character(result))
    break
  }
} 