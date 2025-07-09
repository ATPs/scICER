library(scICER)
library(Seurat)

# Set up parallel processing (optional, but recommended)
n_workers <- 1  # Use single worker to simplify debugging

# Load your Seurat object
data(pbmc_small)  # Example dataset
seurat_obj <- pbmc_small

# Ensure preprocessing is complete
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)  
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Run scICER analysis with reduced range for debugging
cat("Starting scICER analysis...\n")
try({
  scice_results <- scICE_clustering(
    object = seurat_obj,
    cluster_range = 3:5,  # Reduced range for faster debugging
    n_workers = n_workers,
    remove_threshold = Inf,
    verbose = TRUE
  )
  cat("scICER analysis completed successfully!\n")
  print(names(scice_results))
}, silent = FALSE) 