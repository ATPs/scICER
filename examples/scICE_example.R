#' scICE Example Analysis
#' 
#' This script demonstrates how to use scICE for clustering consistency evaluation
#' with a Seurat object. It covers the complete workflow from data preparation
#' to result interpretation.

library(Seurat)
library(scICER)
library(ggplot2)
library(dplyr)

# Load example data (replace with your own data)
# For this example, we'll use the built-in pbmc_small dataset
data("pbmc_small")
seurat_obj <- pbmc_small

# Display basic information about the dataset
cat("Dataset Information:\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")
cat("Available assays:", names(seurat_obj@assays), "\n")

# ============================================================================
# Step 1: Data Preprocessing (if not already done)
# ============================================================================

# Standard Seurat preprocessing workflow
if (!"pca" %in% names(seurat_obj@reductions)) {
  cat("Running standard preprocessing...\n")
  
  # Normalize data
  seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
  
  # Find variable features
  seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
  
  # Scale data
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  # Run PCA
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
}

# Find neighbors (required for scICE)
if (!"snn" %in% names(seurat_obj@graphs)) {
  cat("Computing SNN graph...\n")
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
}

# ============================================================================
# Step 2: Run scICE Analysis
# ============================================================================

cat("Running scICE clustering consistency evaluation...\n")

# Basic scICE analysis
scice_results <- scICE_clustering(
  object = seurat_obj,
  cluster_range = 2:10,          # Test 2-10 clusters
  n_workers = 2,                 # Use 2 workers (adjust based on your system)
  n_trials = 10,                 # Reduced for faster example
  n_bootstrap = 20,              # Reduced for faster example
  verbose = TRUE
)

# Display basic results
cat("\nscICE Analysis Results:\n")
cat("Tested cluster numbers:", paste(scice_results$n_cluster, collapse = ", "), "\n")
cat("Consistent clusters (IC < 1.005):", 
    paste(scice_results$consistent_clusters, collapse = ", "), "\n")

# ============================================================================
# Step 3: Visualize Results
# ============================================================================

cat("Creating visualizations...\n")

# Plot IC scores
ic_plot <- plot_ic(scice_results, 
                   title = "scICE Clustering Consistency Analysis",
                   threshold = 1.005)
print(ic_plot)

# Plot stability analysis
stability_plot <- plot_stability(scice_results)
print(stability_plot)

# ============================================================================
# Step 4: Extract and Use Consistent Clustering Results
# ============================================================================

# Get summary of consistent clusters
consistent_summary <- extract_consistent_clusters(scice_results, threshold = 1.005)
if (nrow(consistent_summary) > 0) {
  cat("\nConsistent Clusters Summary:\n")
  print(consistent_summary)
  
  # Extract clustering labels
  cluster_labels <- get_robust_labels(scice_results, threshold = 1.005)
  
  if (!is.null(cluster_labels)) {
    cat("\nExtracted clustering labels:\n")
    cat("Columns:", paste(names(cluster_labels), collapse = ", "), "\n")
    
    # Add results to Seurat object
    seurat_obj <- get_robust_labels(scice_results, return_seurat = TRUE)
    
    # Display which cluster assignments were added
    new_columns <- grep("clusters_", names(seurat_obj@meta.data), value = TRUE)
    cat("Added to Seurat metadata:", paste(new_columns, collapse = ", "), "\n")
  }
} else {
  cat("\nNo clusters meet the consistency threshold. Consider:\n")
  cat("- Adjusting the threshold (e.g., 1.01 instead of 1.005)\n")
  cat("- Expanding the cluster_range\n")
  cat("- Checking data preprocessing steps\n")
}

# ============================================================================
# Step 5: Visualization with Consistent Clusters (if available)
# ============================================================================

# If we have UMAP/t-SNE, visualize with consistent clusters
if (!"umap" %in% names(seurat_obj@reductions)) {
  cat("Computing UMAP for visualization...\n")
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)
}

# Plot with consistent clustering if available
cluster_columns <- grep("clusters_", names(seurat_obj@meta.data), value = TRUE)

if (length(cluster_columns) > 0) {
  cat("Creating UMAP plots with consistent clusters...\n")
  
  # Plot the first consistent clustering result
  first_clustering <- cluster_columns[1]
  umap_plot <- DimPlot(seurat_obj, 
                       group.by = first_clustering,
                       label = TRUE,
                       label.size = 4) +
    ggtitle(paste("UMAP with", first_clustering)) +
    theme(legend.position = "right")
  
  print(umap_plot)
  
  # If multiple consistent clusterings, create comparison plots
  if (length(cluster_columns) > 1) {
    comparison_plots <- list()
    for (i in seq_along(cluster_columns)) {
      plot_title <- paste("Clusters:", 
                         sub("clusters_", "", cluster_columns[i]))
      comparison_plots[[i]] <- DimPlot(seurat_obj, 
                                      group.by = cluster_columns[i],
                                      label = TRUE) +
        ggtitle(plot_title) +
        theme(legend.position = "none")
    }
    
    # Combine plots if multiple clusterings
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      combined_plot <- gridExtra::grid.arrange(grobs = comparison_plots, ncol = 2)
      print(combined_plot)
    }
  }
}

# ============================================================================
# Step 6: Advanced Analysis Examples
# ============================================================================

cat("\nAdvanced Analysis Examples:\n")

# Example 1: Custom parameter analysis
cat("Example 1: Custom parameter analysis for higher resolution...\n")
high_res_results <- scICE_clustering(
  object = seurat_obj,
  cluster_range = 8:15,          # Focus on higher cluster numbers
  objective_function = "CPM",     # Use CPM objective
  n_trials = 8,                  # Reduced for speed
  n_bootstrap = 15,              # Reduced for speed
  ic_threshold = 1.01,           # More lenient threshold
  verbose = FALSE
)

if (length(high_res_results$consistent_clusters) > 0) {
  cat("High resolution consistent clusters:", 
      paste(high_res_results$consistent_clusters, collapse = ", "), "\n")
}

# Example 2: Detailed stability analysis for a specific cluster number
if (7 %in% scice_results$n_cluster) {
  cat("Example 2: Detailed stability analysis for 7 clusters...\n")
  stability_detail <- plot_stability(scice_results, cluster_number = 7)
  print(stability_detail)
}

# ============================================================================
# Step 7: Save Results
# ============================================================================

cat("Saving results...\n")

# Save the scICE results object
saveRDS(scice_results, file = "scICE_results.rds")

# Save the updated Seurat object
saveRDS(seurat_obj, file = "seurat_with_scICE_clusters.rds")

# Export clustering labels as CSV
if (exists("cluster_labels") && !is.null(cluster_labels)) {
  write.csv(cluster_labels, file = "scICE_cluster_labels.csv", row.names = FALSE)
}

cat("\nAnalysis complete! Files saved:\n")
cat("- scICE_results.rds: Complete scICE analysis results\n")
cat("- seurat_with_scICE_clusters.rds: Seurat object with added cluster labels\n")
if (exists("cluster_labels") && !is.null(cluster_labels)) {
  cat("- scICE_cluster_labels.csv: Cluster labels as CSV\n")
}

# ============================================================================
# Step 8: Interpretation Guide
# ============================================================================

cat("\n" + paste(rep("=", 60), collapse = "") + "\n")
cat("INTERPRETATION GUIDE:\n")
cat(paste(rep("=", 60), collapse = "") + "\n")

cat("\n1. IC Scores:\n")
cat("   - Values close to 1.0 indicate high consistency\n")
cat("   - IC < 1.005 suggests ~0.25% inconsistency (very good)\n")
cat("   - IC < 1.01 suggests ~0.5% inconsistency (good)\n")
cat("   - IC > 1.05 suggests >2.5% inconsistency (poor)\n")

cat("\n2. Cluster Selection:\n")
cat("   - Choose cluster numbers with IC below your threshold\n")
cat("   - Consider biological relevance alongside statistical consistency\n")
cat("   - Multiple consistent options may reflect data structure\n")

cat("\n3. Next Steps:\n")
cat("   - Validate consistent clusters with known markers\n")
cat("   - Perform differential expression analysis\n")
cat("   - Consider functional enrichment analysis\n")
cat("   - Compare with other clustering methods\n")

cat("\n" + paste(rep("=", 60), collapse = "") + "\n")
cat("scICE analysis completed successfully!\n")
cat(paste(rep("=", 60), collapse = "") + "\n") 