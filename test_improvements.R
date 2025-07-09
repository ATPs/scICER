# Test script for scICER improvements
# Load required packages
library(scICER)
library(Seurat)

# Load test data
data(pbmc_small)
seurat_obj <- pbmc_small

# Ensure preprocessing is complete
if (!"pca" %in% names(seurat_obj@reductions)) {
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
}

if (!"snn" %in% names(seurat_obj@graphs)) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
}

cat("Testing scICER improvements...\n")

# Run scICER analysis with a range that should include some excluded clusters
scice_results <- scICE_clustering(
  object = seurat_obj,
  cluster_range = 3:15,  # Same as user's example
  n_workers = 1,         # Use single worker for consistent testing
  n_trials = 5,          # Reduced for faster testing
  n_bootstrap = 10,      # Reduced for faster testing
  verbose = TRUE
)

cat("\n=== TESTING RESULTS ===\n")

# Test 1: Check that cluster numbering is correct
cat("1. Testing cluster numbering fix:\n")
cat("   Cluster range tested:", paste(scice_results$cluster_range_tested, collapse = ", "), "\n")
cat("   Actual cluster numbers in results:", paste(scice_results$n_cluster, collapse = ", "), "\n")
cat("   Consistent cluster numbers:", paste(scice_results$consistent_clusters, collapse = ", "), "\n")

# Test 2: Check exclusion information
if (!is.null(scice_results$excluded)) {
  cat("\n2. Testing exclusion information:\n")
  excluded_clusters <- scice_results$n_cluster[scice_results$excluded]
  if (length(excluded_clusters) > 0) {
    cat("   Excluded clusters:", paste(excluded_clusters, collapse = ", "), "\n")
    exclusion_reasons <- unique(scice_results$exclusion_reason[scice_results$excluded])
    cat("   Exclusion reasons:", paste(exclusion_reasons, collapse = ", "), "\n")
  } else {
    cat("   No clusters were excluded\n")
  }
} else {
  cat("\n2. No exclusion information available (old format)\n")
}

# Test 3: Test plot_ic with new features
cat("\n3. Testing enhanced plot_ic function:\n")
tryCatch({
  ic_plot <- plot_ic(scice_results, threshold = 1.005)
  cat("   ✓ plot_ic function works with new format\n")
}, error = function(e) {
  cat("   ✗ plot_ic function failed:", e$message, "\n")
})

# Test 4: Test get_robust_labels
cat("\n4. Testing get_robust_labels function:\n")
tryCatch({
  robust_labels <- get_robust_labels(scice_results, threshold = 1.005)
  if (!is.null(robust_labels)) {
    cat("   ✓ get_robust_labels extracted", ncol(robust_labels) - 1, "clusterings\n")
    cat("   Available clusterings:", paste(names(robust_labels)[-1], collapse = ", "), "\n")
  } else {
    cat("   ⚠ No robust labels found (may be normal for test data)\n")
  }
}, error = function(e) {
  cat("   ✗ get_robust_labels function failed:", e$message, "\n")
})

# Test 5: Test extract_consistent_clusters
cat("\n5. Testing extract_consistent_clusters function:\n")
tryCatch({
  consistent_summary <- extract_consistent_clusters(scice_results, threshold = 1.005)
  if (nrow(consistent_summary) > 0) {
    cat("   ✓ extract_consistent_clusters found", nrow(consistent_summary), "consistent clusterings\n")
    print(consistent_summary)
  } else {
    cat("   ⚠ No consistent clusters found (may be normal for test data)\n")
  }
}, error = function(e) {
  cat("   ✗ extract_consistent_clusters function failed:", e$message, "\n")
})

cat("\n=== TEST COMPLETE ===\n")
cat("All improvements have been tested. Check the output above for any issues.\n") 