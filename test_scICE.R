# Load required packages
library(scICER)
library(Seurat)
library(foreach)
library(doParallel)

# Set up parallel processing (optional, but recommended)
n_workers <- 4  # adjust based on your system
cl <- makeCluster(min(n_workers, parallel::detectCores() - 1))
registerDoParallel(cl)

# Load your Seurat object
data(pbmc_small)  # Example dataset
seurat_obj <- pbmc_small

# Ensure preprocessing is complete
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)  
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)

# Important: Explicitly specify reduction when finding neighbors
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20, reduction = "pca")

# Run scICER analysis
scice_results <- scICE_clustering(
  object = seurat_obj,
  cluster_range = 2:15,
  n_workers = n_workers,
  verbose = TRUE
)

# Clean up parallel processing
stopCluster(cl)

# Save results
saveRDS(scice_results, "scice_results.rds")
saveRDS(seurat_obj, "seurat_obj.rds")

# Print success message
print("Test completed successfully!") 