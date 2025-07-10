\name{scICE_clustering}
\alias{scICE_clustering}
\title{Run scICER clustering analysis on a Seurat object}
\usage{
scICE_clustering(object, cluster_range = 2:10, n_workers = max(1, parallel::detectCores() - 1),
                n_trials = 100, n_bootstrap = 100, seed = NULL, beta = 0.01, n_iterations = 10,
                max_iterations = 10, objective_function = "CPM",
                remove_threshold = Inf, resolution_tolerance = 0.001,
                graph_name = NULL, verbose = TRUE)
}
\arguments{
  \item{object}{A Seurat object with a precomputed SNN graph.}
  
  \item{cluster_range}{Vector of cluster numbers to test (default: 2:10).}
  
  \item{n_workers}{Number of parallel workers for computation (default: max(1, parallel::detectCores() - 1)).}
  
  \item{n_trials}{Number of clustering trials per resolution (default: 100).}
  
  \item{n_bootstrap}{Number of bootstrap iterations (default: 100).}
  
  \item{seed}{Random seed for reproducibility (default: NULL). When set, ensures identical results across multiple runs with the same parameters.}
  
  \item{beta}{Beta parameter for Leiden clustering (default: 0.01).}
  
  \item{n_iterations}{Number of Leiden iterations (default: 10).}
  
  \item{max_iterations}{Maximum iterations for optimization (default: 10).}
  
  \item{objective_function}{Objective function to use, either "modularity" or "CPM" (default: "CPM"). 
    \itemize{
      \item CPM (Constant Potts Model): Recommended for most cases, especially larger datasets. CPM is more flexible 
            and can overcome the resolution limit by using a resolution parameter (gamma) that controls cluster sizes. 
            Higher resolution finds smaller clusters, lower resolution finds larger clusters.
      \item modularity: Traditional approach that measures connection density within clusters compared to random networks. 
            Better suited for smaller datasets but may struggle with identifying small, well-defined clusters 
            (known as the resolution limit).
    }}
  
  \item{remove_threshold}{Threshold for removing inconsistent results (default: Inf). Set to Inf to include all cluster numbers for plotting and analysis.}
  
  \item{resolution_tolerance}{Tolerance for resolution parameter search (default: 0.001).}
  
  \item{graph_name}{Name of the graph to use for clustering. If NULL (default), will use the default SNN graph from the active assay (e.g., "RNA_snn" for RNA assay).}
  
  \item{verbose}{Whether to print progress messages (default: TRUE). When TRUE, provides detailed information about:
    \itemize{
      \item Graph validation status and properties
      \item Progress for each cluster number being tested
      \item Number of stable solutions found
      \item Optimization iteration progress
      \item Summary of consistent cluster numbers identified
      \item List of excluded cluster numbers and reasons for exclusion
      \item List of inconsistent cluster numbers (tested but above threshold)
      \item Quality metrics for each clustering solution
    }}
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
  \item{cluster_range_tested}{Original cluster range that was tested}
  \item{excluded}{Logical vector indicating which clusters were excluded}
  \item{exclusion_reason}{Character vector explaining why clusters were excluded}
}

The function now includes information about all attempted cluster numbers, including those that were excluded due to instability or failed resolution parameter search. This allows users to understand why certain cluster numbers are missing from the results.
}
\description{
Performs scICER (single-cell Inconsistency-based Clustering Evaluation in R) analysis on a Seurat object to identify stable clustering solutions. The function implements parallel processing and comprehensive stability assessment through bootstrap analysis.

The algorithm uses the Leiden community detection method with adaptive resolution parameter search to find optimal clustering solutions for each desired number of clusters.
}
\details{
The scICER algorithm proceeds through several steps:

1. Graph Preparation:
\itemize{
  \item Uses the SNN graph from the specified assay
  \item Automatically detects the appropriate graph if not specified
  \item Validates graph structure and connectivity
  \item Reports graph properties (nodes, edges, density)
}

2. Parameter Search:
\itemize{
  \item Binary search for resolution parameters
  \item Parallel processing for efficiency
  \item Adaptive ranges based on objective function
  \item Reports progress for each cluster number
}

3. Stability Assessment:
\itemize{
  \item Bootstrap analysis for each cluster number
  \item Inconsistency score calculation
  \item Filtering of unstable solutions
  \item Reports number of stable solutions found
}

4. Result Optimization:
\itemize{
  \item Iterative refinement of clustering
  \item MEI score calculation
  \item Selection of best solutions
  \item Reports quality metrics for final solutions
}

5. Reproducibility Control:
\itemize{
  \item Optional seed parameter for deterministic results
  \item Hierarchical seeding strategy for parallel processing
  \item Consistent results across multiple runs with same seed
  \item Different seeds produce different valid solutions
}

Performance can be tuned through several parameters:
\itemize{
  \item n_workers: Controls parallel processing
  \item n_trials and n_bootstrap: Balance accuracy vs. speed
  \item resolution_tolerance: Affects parameter search precision
  \item seed: Controls reproducibility of results
}

The function provides comprehensive progress updates when verbose = TRUE:
\itemize{
  \item Initial graph validation and properties
  \item Search progress for each cluster number
  \item Number of stable solutions identified
  \item Optimization iteration progress
  \item Final results summary including:
    \itemize{
      \item Total number of consistent clusters found
      \item Range of cluster numbers that passed stability criteria
      \item Inconsistency scores for each solution
      \item MEI scores for final clusterings
    }
}
}
\examples{
\dontrun{
# Load example data
data(pbmc_small)

# Basic usage with CPM (recommended for most cases)
results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Using modularity with seed for reproducibility
results <- scICE_clustering(
  pbmc_small,
  cluster_range = 2:15,
  n_trials = 200,
  n_bootstrap = 150,
  seed = 12345,  # Set seed for reproducible results
  objective_function = "modularity",
  graph_name = "RNA_snn",
  verbose = TRUE
)

# Reproduce exact results by using same seed and parameters
results_reproduced <- scICE_clustering(
  pbmc_small,
  cluster_range = 2:15,
  n_trials = 200,
  n_bootstrap = 150,
  seed = 12345,  # Same seed = same results
  objective_function = "modularity",
  graph_name = "RNA_snn",
  verbose = TRUE
)

# Try different seeds to explore alternative solutions
seeds_to_try <- c(123, 456, 789)
results_list <- lapply(seeds_to_try, function(s) {
  scICE_clustering(pbmc_small, cluster_range = 2:10, seed = s)
})

# Access results
print(results$consistent_clusters)  # Shows which cluster numbers were stable
plot_ic(results)  # Visualize inconsistency scores
plot_stability(results)  # Visualize stability metrics
}
}
\seealso{
\code{\link{plot_ic}}, \code{\link{plot_stability}}, \code{\link{extract_consistent_clusters}}
}
%\export 