\name{scICE_clustering}
\alias{scICE_clustering}
\title{Run scICE clustering analysis on a Seurat object}
\usage{
scICE_clustering(object, cluster_range = 2:10, n_workers = max(1, parallel::detectCores() - 1),
                n_trials = 100, n_bootstrap = 100, beta = 0.01, n_iterations = 10,
                max_iterations = 10, objective_function = "modularity",
                remove_threshold = 1.1, resolution_tolerance = 0.001,
                graph_name = NULL, verbose = TRUE)
}
\arguments{
  \item{object}{A Seurat object with a precomputed SNN graph.}
  
  \item{cluster_range}{Vector of cluster numbers to test (default: 2:10).}
  
  \item{n_workers}{Number of parallel workers for computation (default: max(1, parallel::detectCores() - 1)).}
  
  \item{n_trials}{Number of clustering trials per resolution (default: 100).}
  
  \item{n_bootstrap}{Number of bootstrap iterations (default: 100).}
  
  \item{beta}{Beta parameter for Leiden clustering (default: 0.01).}
  
  \item{n_iterations}{Number of Leiden iterations (default: 10).}
  
  \item{max_iterations}{Maximum iterations for optimization (default: 10).}
  
  \item{objective_function}{Objective function to use, either "modularity" or "CPM" (default: "modularity").}
  
  \item{remove_threshold}{Threshold for removing inconsistent results (default: 1.1).}
  
  \item{resolution_tolerance}{Tolerance for resolution parameter search (default: 0.001).}
  
  \item{graph_name}{Name of the graph to use for clustering. If NULL (default), will use the default SNN graph from the active assay (e.g., "RNA_snn" for RNA assay).}
  
  \item{verbose}{Whether to print progress messages (default: TRUE).}
}
\value{
Returns a list with class "scICE" containing:
\itemize{
  \item{gamma}{Resolution parameters used for each cluster number}
  \item{labels}{List of clustering label arrays}
  \item{ic}{Vector of inconsistency scores}
  \item{ic_vec}{List of bootstrap inconsistency scores}
  \item{n_cluster}{Vector of cluster numbers tested}
  \item{best_labels}{List of best clustering labels}
  \item{n_iter}{Vector of iterations needed for convergence}
  \item{mei}{List of Mutual Element-wise Information scores}
  \item{k}{Vector of final iteration counts}
  \item{seurat_object}{Original Seurat object}
  \item{cell_names}{Vector of cell names}
  \item{consistent_clusters}{Vector of cluster numbers meeting consistency threshold}
  \item{graph_name}{Name of the graph used for clustering}
}
}
\description{
Performs scICE (single-cell Iterative Clustering Ensemble) analysis on a Seurat object to identify stable clustering solutions. The function implements parallel processing and comprehensive stability assessment through bootstrap analysis.

The algorithm uses the Leiden community detection method with adaptive resolution parameter search to find optimal clustering solutions for each desired number of clusters.
}
\details{
The scICE algorithm proceeds through several steps:

1. Graph Preparation:
\itemize{
  \item Uses the SNN graph from the specified assay
  \item Automatically detects the appropriate graph if not specified
  \item Validates graph structure and connectivity
}

2. Parameter Search:
\itemize{
  \item Binary search for resolution parameters
  \item Parallel processing for efficiency
  \item Adaptive ranges based on objective function
}

3. Stability Assessment:
\itemize{
  \item Bootstrap analysis for each cluster number
  \item Inconsistency score calculation
  \item Filtering of unstable solutions
}

4. Result Optimization:
\itemize{
  \item Iterative refinement of clustering
  \item MEI score calculation
  \item Selection of best solutions
}

Performance can be tuned through several parameters:
\itemize{
  \item n_workers: Controls parallel processing
  \item n_trials and n_bootstrap: Balance accuracy vs. speed
  \item resolution_tolerance: Affects parameter search precision
}

The function provides comprehensive progress updates when verbose = TRUE:
\itemize{
  \item Graph validation status
  \item Search progress for each cluster number
  \item Optimization iterations
  \item Final results summary
}
}
\examples{
\dontrun{
# Load example data
data(pbmc_small)

# Basic usage
results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# With custom parameters
results <- scICE_clustering(
  pbmc_small,
  cluster_range = 2:15,
  n_trials = 200,
  n_bootstrap = 150,
  objective_function = "CPM",
  graph_name = "RNA_snn",
  verbose = TRUE
)

# Access results
print(results$consistent_clusters)
plot_ic(results)
plot_stability(results)
}
}
\seealso{
\code{\link{plot_ic}}, \code{\link{plot_stability}}, \code{\link{extract_consistent_clusters}}
}
\export 