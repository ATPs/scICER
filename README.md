# scICER: Single-cell Inconsistency-based Clustering Evaluation in R

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue)](https://cran.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-%3E%3D4.0.0-green)](https://satijalab.org/seurat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**scICER** is an R package that implements a systematic and efficient workflow to evaluate clustering consistency in single-cell RNA-seq data. This is the R port of the original [scICE Julia package](https://github.com/Mathbiomed/scICE), designed to be fully compatible with the [Seurat](https://satijalab.org/seurat/) ecosystem.

### Key Features

- 🔬 **Automated Cluster Evaluation**: Systematically tests multiple cluster numbers to find consistent results
- ⚡ **Element-Centric Similarity (ECS)**: Efficient algorithm for comparing clustering results
- 🧬 **Seurat Integration**: Seamless integration with Seurat objects and workflows  
- 🚀 **Parallel Processing**: Leverages multiple cores for faster computation
- 📊 **Comprehensive Visualization**: Built-in plotting functions for result interpretation
- 📈 **Robust Statistics**: Bootstrap-based confidence intervals and stability assessment
- 🎯 **Reproducible Analysis**: Optional seed parameter for deterministic results
- 🔍 **Smart Parameter Recommendations**: Automatically suggests optimal parameters based on dataset size
- ✅ **Preprocessing Validation**: Checks and guides users through proper Seurat preprocessing
- 📋 **Comprehensive Reporting**: Generates detailed analysis summaries with interpretation guidance

### Algorithm Overview

scICER evaluates clustering consistency by:

1. **Multiple Clustering**: Runs Leiden clustering multiple times with different parameters
2. **Binary Search**: Automatically finds resolution parameters for target cluster numbers  
3. **ECS Calculation**: Uses Element-Centric Similarity to efficiently compare clustering results
4. **Optimization**: Iteratively refines clustering to minimize inconsistency
5. **Bootstrap Validation**: Assesses stability through bootstrap sampling
6. **Threshold Application**: Identifies cluster numbers meeting consistency criteria
7. **Reproducibility Control**: Optional seed parameter for deterministic results

## Installation

### Prerequisites

Ensure you have R ≥ 4.0.0 and the required dependencies:

```r
# Install required CRAN packages
install.packages(c(
  "Seurat", "SeuratObject", "igraph", "Matrix", 
  "parallel", "foreach", "doParallel", "dplyr", 
  "ggplot2", "devtools"
))

# Optional: Install leiden package for better performance
install.packages("leiden")
```

### Install scICER from GitHub

```r
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install scICER from GitHub
devtools::install_github("ATPs/scICER", upgrade = "never")#avoid upgrade packages and use existing ones
```

### Load the package

```r
library(scICER)
library(Seurat)
```

## Quick Start

### Basic Workflow

```r
# Load required packages
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

# For reproducible results, use the seed parameter
scice_results <- scICE_clustering(
  object = seurat_obj,
  cluster_range = 3:20,
  n_workers = n_workers,
  remove_threshold = Inf,
  verbose = TRUE,
  seed = 123  # Set seed for reproducible results
)

# Visualize results
plot_ic(scice_results, threshold = 1.005)

# Extract consistent clustering labels
seurat_obj <- get_robust_labels(scice_results, return_seurat = TRUE)

# Visualize with consistent clusters, cluster range from 3 to 15
xl.fig <- list()
for (n in 3:15) {
  if (paste0("clusters_", n) %in% colnames(seurat_obj@meta.data)){
    xl.fig[[paste0("clusters_", n)]] <- DimPlot(seurat_obj, group.by = paste0("clusters_", n), label = TRUE) + ggtitle(paste("Clusters:", n))
  }
}

cowplot::plot_grid(plotlist = xl.fig, ncol = 3)
```

## Detailed Usage

### 1. Standard Analysis

```r
# Comprehensive scICER analysis
results <- scICE_clustering(
  object = your_seurat_object,
  cluster_range = 2:20,           # Range of cluster numbers to test
  n_workers = 8,                  # Number of parallel workers  
  n_trials = 15,                  # Clustering trials per resolution
  n_bootstrap = 100,              # Bootstrap iterations
  seed = 42,                      # Optional: Set for reproducible results
  ic_threshold = 1.005,           # Consistency threshold
  verbose = TRUE                  # Progress messages
)
```

### 2. Advanced Parameter Customization

```r
# Fine-tuned parameters for large datasets
results <- scICE_clustering(
  object = large_seurat_object,
  graph_name = "RNA_snn",         # Optional: specify graph name
  cluster_range = 5:12,           # Focused range for efficiency
  objective_function = "CPM",     # "CPM" or "modularity"
  beta = 0.1,                     # Leiden beta parameter
  n_iterations = 10,              # Initial Leiden iterations
  max_iterations = 150,           # Maximum optimization iterations
  remove_threshold = 1.15,        # Filter inconsistent results
  resolution_tolerance = 1e-8,    # Resolution search precision
  n_trials = 10,                  # Reduced for speed
  n_bootstrap = 50,               # Reduced for speed
  seed = 12345                    # For reproducible results
)
```

It's OK to set a large `remove_threshold` value, because `plot_ic` and `get_robust_labels` all have default `threshold = 1.005`.

### 3. Visualization and Results

```r
# Plot Inconsistency Coefficient (IC) scores
plot_ic(results, 
        threshold = 1.005,
        title = "Clustering Consistency Analysis")

# Plot stability analysis
plot_stability(results)

# Extract summary of consistent clusters
summary <- extract_consistent_clusters(results, threshold = 1.005)
print(summary)

# Get clustering labels as data frame
labels_df <- get_robust_labels(results)

# Add results directly to Seurat object
seurat_obj <- get_robust_labels(results, return_seurat = TRUE)
```

### 4. Integration with Seurat Workflow

```r
# Complete Seurat + scICER workflow
seurat_obj <- CreateSeuratObject(counts = your_counts)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Run scICER with reproducibility
scice_results <- scICE_clustering(
  seurat_obj, 
  cluster_range = 3:15,
  seed = 42  # For reproducible results
)

# Add consistent clusters to Seurat object
seurat_obj <- get_robust_labels(scice_results, return_seurat = TRUE)

# Continue with standard Seurat analysis
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)
DimPlot(seurat_obj, group.by = "clusters_8")  # Use consistent clustering

# Find markers for consistent clusters
Idents(seurat_obj) <- "clusters_8"
markers <- FindAllMarkers(seurat_obj)
```

## Understanding the Output

### Main Results Object

The `scICE_clustering()` function returns a list containing:

- `gamma`: Resolution parameters for each cluster number
- `labels`: Clustering results for each cluster number  
- `ic`: Inconsistency scores (lower = more consistent)
- `ic_vec`: Bootstrap IC distributions
- `n_cluster`: Number of clusters tested
- `best_labels`: Best clustering labels for each cluster number
- `consistent_clusters`: Cluster numbers meeting consistency threshold
- `mei`: Mutual Element-wise Information scores

### Interpretation Guidelines

- **IC < 1.005**: Highly consistent clustering (recommended)
- **IC 1.005-1.01**: Moderately consistent clustering
- **IC > 1.01**: Low consistency (not recommended)

### Visualization Outputs

- **IC Plot**: Shows consistency scores across cluster numbers
- **Stability Plot**: Displays bootstrap distributions
- **UMAP/t-SNE**: Visualize clusters with consistent labels

## Performance Tips

### For Large Datasets (>50k cells)

```r
# Get recommended parameters based on dataset size
recommended <- get_recommended_parameters(seurat_obj)

# Optimized settings for large datasets
results <- scICE_clustering(
  object = large_seurat_obj,
  cluster_range = 5:15,      # Focused range
  n_trials = 8,              # Fewer trials
  n_bootstrap = 50,          # Fewer bootstrap samples
  n_workers = detectCores() - 1,  # Maximum parallel processing
  seed = 42                  # For reproducibility
)
```

### Memory Management

- Reduce `n_trials` and `n_bootstrap` for memory-limited systems
- Use focused `cluster_range` based on biological expectations
- Install the `leiden` package for improved performance
- Use `get_recommended_parameters()` for optimal settings

## Troubleshooting

### Common Issues

1. **"Graph not found" error**: 
   - Run `FindNeighbors()` on your Seurat object first
   - Use `check_seurat_ready()` to validate preprocessing
   - By default, scICER uses the SNN graph from the active assay (e.g., "RNA_snn")
   - You can check available graphs with `names(your_seurat_object@graphs)`
   - Specify a different graph with the `graph_name` parameter if needed

2. **No consistent clusters found**: 
   - Try adjusting `ic_threshold` (e.g., 1.01)
   - Expand `cluster_range`
   - Check preprocessing quality
   - Use `create_results_summary()` for detailed analysis

3. **Memory issues**: 
   - Reduce `n_workers`, `n_trials`, or `cluster_range`
   - Use `get_recommended_parameters()` for guidance

4. **Slow performance**: 
   - Install the `leiden` package
   - Increase `n_workers`
   - Use recommended parameters for your dataset size

5. **Reproducibility issues**:
   - Always set the `seed` parameter
   - Keep all other parameters constant
   - Document parameter values used

### Performance Optimization

```r
# Check if leiden package is available
if (requireNamespace("leiden", quietly = TRUE)) {
  message("leiden package detected - enhanced performance available")
} else {
  message("Consider installing 'leiden' package for better performance")
}

# Validate Seurat object preprocessing
check_seurat_ready(seurat_obj)

# Get recommended parameters
params <- get_recommended_parameters(seurat_obj)
```

## Citation

If you use scICER in your research, please cite:

```
Cao, X. (2024). scICER: Single-cell Inconsistency-based Clustering Evaluation in R. 
R package version 1.0.0. https://github.com/ATPs/scICER
```

*The original scICE algorithm citation will be added upon publication*

## Q&A

### How does scICER select the final clustering solution?

For each number of clusters (e.g., when testing 3-10 clusters):
1. The algorithm runs `n_trials` independent clustering attempts (default: 15)
2. Each trial uses the Leiden algorithm with different random initializations
3. For each cluster number:
   - Calculates pairwise similarities between all trial results using Element-Centric Similarity (ECS)
   - Computes an Inconsistency (IC) score measuring how different the clustering solutions are
   - A lower IC score indicates more consistent clustering across trials
4. For cluster numbers with IC scores below your threshold (e.g., < 1.05):
   - Selects the "best" clustering solution from the trials
   - The "best" solution is the one with highest average similarity to other solutions (most representative)
   - This becomes the reported clustering for that number of clusters

So when you get results showing clusters 3-10 have IC scores < 1.05, it means:
- Each of these cluster numbers produced stable results across different random initializations
- The reported clustering solution is the most representative one from all trials
- You can be confident in using these clustering solutions for downstream analysis

### Will increasing n_trials or n_bootstrap improve results?

**n_trials**:
- Increasing `n_trials` (default: 15) can improve the robustness of results
- More trials = more chances to find truly stable clustering solutions
- For large datasets, you might want to reduce trials to improve speed
- Recommended ranges:
  - Small datasets (<10k cells): 15-20 trials
  - Medium datasets (10k-50k cells): 10-15 trials
  - Large datasets (>50k cells): 8-10 trials

**n_bootstrap**:
- Increasing `n_bootstrap` (default: 100) improves confidence in IC score estimates
- More bootstrap iterations = more accurate stability assessment
- The default 100 iterations is usually sufficient
- Consider increasing for:
  - Critical analyses requiring high confidence
  - Publications or benchmark studies
  - Unusual or complex datasets
- Recommended ranges:
  - Standard analysis: 100 iterations
  - High-confidence analysis: 200-500 iterations
  - Note: Higher values increase computation time

Remember: Quality preprocessing (normalization, feature selection, dimensionality reduction) often has more impact than increasing these parameters.

## Support and Contributing

- 📋 **Issues**: Report bugs or request features on [GitHub Issues](https://github.com/ATPs/scICER/issues)
- 📧 **Contact**: atps@outlook.com
- 🤝 **Contributing**: Pull requests welcome!

## License

This package is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

## Acknowledgments

- Original [scICE Julia package](https://github.com/Mathbiomed/scICE) developers
- [Seurat](https://satijalab.org/seurat/) team for the excellent single-cell framework
- R community for the robust package ecosystem

## Related Resources

- 📚 [Seurat Documentation](https://satijalab.org/seurat/)
- 🔬 [Original scICE (Julia)](https://github.com/Mathbiomed/scICE)
- 📖 [Single-cell Analysis Best Practices](https://www.sc-best-practices.org/) 