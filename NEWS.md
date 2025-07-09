# scICER 1.0.0

## Major Features

* **Initial Release**: First stable release of scICER, the R port of the Julia scICE package
* **Seurat Integration**: Full compatibility with Seurat objects and workflows
* **Element-Centric Similarity (ECS)**: Efficient algorithm for comparing clustering results
* **Parallel Processing**: Multi-core support for faster computation
* **Comprehensive Visualization**: Built-in plotting functions for result interpretation

## Core Functions

* `scICE_clustering()`: Main function for clustering consistency evaluation
* `plot_ic()`: Visualize inconsistency coefficients across cluster numbers
* `plot_stability()`: Display bootstrap stability analysis
* `get_robust_labels()`: Extract consistent clustering labels
* `extract_consistent_clusters()`: Get summary of consistent clusters
* `calculate_ecs()`: Calculate Element-Centric Similarity between clusterings
* `calculate_ic()`: Calculate inconsistency scores
* `calculate_mei_from_array()`: Calculate Mutual Element-wise Information

## Algorithm Features

* **Binary Search Optimization**: Automatically finds resolution parameters for target cluster numbers
* **Bootstrap Validation**: Robust statistical assessment of clustering stability
* **Leiden Clustering**: High-quality community detection with modularity and CPM objectives
* **Memory Efficient**: Optimized data structures for large single-cell datasets
* **Error Handling**: Comprehensive error checking and informative messages

## Integration Features

* **Seurat Workflow**: Seamless integration with standard Seurat preprocessing
* **Graph Support**: Works with SNN and KNN graphs from Seurat
* **Metadata Integration**: Automatic addition of results to Seurat metadata
* **Visualization Support**: Compatible with Seurat's DimPlot and other functions

## Performance Features

* **Multi-threading**: Parallel processing using foreach and doParallel
* **Leiden Package Support**: Optional use of the leiden package for enhanced performance
* **Memory Management**: Efficient handling of large datasets
* **Scalable**: Tested on datasets from hundreds to hundreds of thousands of cells

## Documentation

* **Comprehensive Vignettes**: Step-by-step tutorials for common use cases
* **Function Documentation**: Detailed documentation for all exported functions
* **Examples**: Real-world examples using standard datasets
* **Best Practices**: Guidelines for parameter selection and interpretation

## Dependencies

### Required
* R >= 4.0.0
* Seurat >= 4.0.0
* SeuratObject
* igraph
* Matrix
* parallel
* foreach
* doParallel
* dplyr
* ggplot2
* stats, methods, utils, grDevices, graphics

### Suggested
* leiden (for improved performance)
* testthat (for testing)
* knitr, rmarkdown (for vignettes)
* gridExtra (for advanced plotting)

## Installation

```r
# Install from GitHub
devtools::install_github("ATPs/scICER")

# Or use the installation script
source("install_and_test.R")
```

## Quick Start

```r
library(scICER)
library(Seurat)

# Standard Seurat preprocessing
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)

# Run scICER analysis
results <- scICE_clustering(seurat_obj, cluster_range = 2:15)

# Visualize and extract results
plot_ic(results)
seurat_obj <- get_robust_labels(results, return_seurat = TRUE)
```

## Citation

If you use scICER in your research, please cite:

> Cao, X. (2024). scICER: Single-cell Inconsistency-based Clustering Evaluation in R. 
> R package version 1.0.0. https://github.com/ATPs/scICER

## Acknowledgments

* Original scICE Julia package developers
* Seurat development team
* R single-cell genomics community

---

For detailed usage instructions, see `vignette("scICER-quickstart")` or visit the [GitHub repository](https://github.com/ATPs/scICER). 