\name{scICE_clustering}
\alias{scICE_clustering}
\title{Single-cell Inconsistency-based Clustering Evaluation (scICE)}
\description{
scICE implements a systematic and efficient workflow to evaluate clustering consistency 
in single-cell RNA-seq data. It generates multiple cluster labels using the Leiden 
algorithm, calculates pairwise similarity with Element-Centric Similarity (ECS), 
and provides automated framework for finding consistent clusters.
}
\usage{
scICE_clustering(object, 
                graph_name = "snn",
                cluster_range = 1:20,
                n_workers = 10,
                n_trials = 15,
                n_bootstrap = 100,
                beta = 0.1,
                n_iterations = 10,
                max_iterations = 150,
                ic_threshold = 1.005,
                objective_function = "CPM",
                remove_threshold = 1.15,
                resolution_tolerance = 1e-8,
                verbose = TRUE)
}
\arguments{
\item{object}{A Seurat object containing single-cell data}

\item{graph_name}{Name of the graph to use for clustering (default: "snn")}

\item{cluster_range}{Vector of cluster numbers to test (default: 1:20)}

\item{n_workers}{Number of parallel workers to use (default: 10)}

\item{n_trials}{Number of clustering trials per resolution (default: 15)}

\item{n_bootstrap}{Number of bootstrap iterations (default: 100)}

\item{beta}{Beta parameter for Leiden clustering (default: 0.1)}

\item{n_iterations}{Number of Leiden iterations (default: 10)}

\item{max_iterations}{Maximum iterations for optimization (default: 150)}

\item{ic_threshold}{IC threshold for consistent clustering (default: 1.005)}

\item{objective_function}{Objective function for Leiden ("modularity" or "CPM", default: "CPM")}

\item{remove_threshold}{Threshold for removing inconsistent results (default: 1.15)}

\item{resolution_tolerance}{Tolerance for resolution parameter search (default: 1e-8)}

\item{verbose}{Whether to print progress messages (default: TRUE)}
}
\value{
A list containing:
\describe{
  \item{gamma}{Resolution parameters for each cluster number}
  \item{labels}{Clustering results for each cluster number}
  \item{ic}{Inconsistency scores for each cluster number}
  \item{ic_vec}{Bootstrap IC distributions}
  \item{n_cluster}{Number of clusters tested}
  \item{best_labels}{Best clustering labels for each cluster number}
  \item{n_iter}{Number of iterations used}
  \item{mei}{Mutual Element-wise Information scores}
  \item{consistent_clusters}{Cluster numbers meeting consistency threshold}
}
}
\details{
The scICE algorithm works by:
\enumerate{
  \item Testing multiple cluster numbers within the specified range
  \item Using binary search to find appropriate resolution parameters
  \item Running multiple clustering trials with different random seeds
  \item Calculating Element-Centric Similarity (ECS) between clustering results
  \item Computing Inconsistency Coefficients (IC) to measure clustering stability
  \item Using bootstrap sampling to assess statistical significance
  \item Identifying cluster numbers that meet the consistency threshold
}

Lower IC scores indicate more consistent clustering results. An IC score of 1.0 
represents perfect consistency, while higher values indicate increasing inconsistency.

The algorithm supports both CPM (Constant Potts Model) and modularity objective 
functions for Leiden clustering. CPM is generally recommended for larger datasets, 
while modularity may work better for smaller datasets.
}
\examples{
\dontrun{
# Load Seurat object
library(Seurat)
data(pbmc_small)

# Ensure preprocessing is complete
pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- RunPCA(pbmc_small)
pbmc_small <- FindNeighbors(pbmc_small, dims = 1:10)

# Run scICE analysis
scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Plot IC results
plot_ic(scice_results)

# Extract consistent clustering labels
consistent_labels <- get_robust_labels(scice_results, threshold = 1.005)

# Add results to Seurat object
pbmc_small <- get_robust_labels(scice_results, return_seurat = TRUE)
}
}
\references{
Original scICE algorithm: [Citation to be added upon publication]
}
\seealso{
\code{\link{plot_ic}}, \code{\link{get_robust_labels}}, \code{\link{extract_consistent_clusters}}
}
\author{
Xiaolong Cao
}
\keyword{clustering}
\keyword{single-cell}
\keyword{consistency} 