library(scICER)
library(Seurat)
library(foreach)
library(doParallel)
library(ggplot2)

# Set up parallel processing (optional, but recommended)
n_workers <- 4  # adjust based on your system

# Load your Seurat object
data(pbmc_small)  # Example dataset
seurat_obj <- pbmc_small

# Ensure preprocessing is complete
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)  
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Run scICER analysis
scice_results <- scICE_clustering(
  object = seurat_obj,
  cluster_range = 3:20,
  n_workers = n_workers,
  remove_threshold = Inf,
  verbose = TRUE
)



# Visualize results
plot_ic(scice_results, threshold=Inf)

# Extract consistent clustering labels
seurat_obj <- get_robust_labels(scice_results, return_seurat = TRUE)

# Visualize with consistent clusters, cluster ange from 3 to 15
xl.fig <- list()
for (n in 3:15) {
  if (paste0("clusters_", n) %in% colnames(seurat_obj@meta.data))
    xl.fig[[n]] <- DimPlot(seurat_obj, group.by = paste0("clusters_", n), label = TRUE) + ggtitle(paste("Clusters:", n))
}

cowplot::plot_grid(plotlist = xl.fig, ncol = 3)
