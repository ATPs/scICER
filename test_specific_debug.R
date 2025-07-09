library(scICER)
library(Seurat)

# Set up parallel processing
n_workers <- 1

# Load and prepare Seurat object
data(pbmc_small)
seurat_obj <- pbmc_small
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)  
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Test individual cluster numbers that might be problematic
test_numbers <- c(13, 14, 15)

for (cluster_num in test_numbers) {
  cat("\n=== Testing single cluster number:", cluster_num, "===\n")
  
  result <- try({
    scice_results <- scICE_clustering(
      object = seurat_obj,
      cluster_range = cluster_num,  # Test just one cluster number
      n_workers = n_workers,
      remove_threshold = Inf,
      verbose = TRUE
    )
    cat("SUCCESS: Cluster number", cluster_num, "completed\n")
  }, silent = FALSE)
  
  if (inherits(result, "try-error")) {
    cat("ERROR at cluster number", cluster_num, ":\n")
    print(result)
    break
  }
}

cat("\nNow testing the combination that fails:\n")
result <- try({
  scice_results <- scICE_clustering(
    object = seurat_obj,
    cluster_range = 13:15,
    n_workers = n_workers,
    remove_threshold = Inf,
    verbose = TRUE
  )
  cat("SUCCESS: Range 13-15 completed\n")
}, silent = FALSE)

if (inherits(result, "try-error")) {
  cat("ERROR in range 13-15:\n")
  print(result)
} 