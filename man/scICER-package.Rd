\name{scICER-package}
\alias{scICER-package}
\alias{scICER}
\docType{package}
\title{
\packageTitle{scICER}
}
\description{
\packageDescription{scICER}
}
\details{
The scICER package provides a systematic and efficient workflow for evaluating 
clustering consistency in single-cell RNA-seq data. It is fully compatible with 
the Seurat ecosystem and designed for seamless integration into standard 
single-cell analysis workflows.

\subsection{Key Features}{
\itemize{
  \item \strong{Automated Cluster Evaluation}: Systematically tests multiple cluster numbers
  \item \strong{Element-Centric Similarity (ECS)}: Efficient clustering comparison algorithm
  \item \strong{Seurat Integration}: Seamless workflow with Seurat objects
  \item \strong{Parallel Processing}: Multi-core support for faster computation  
  \item \strong{Comprehensive Visualization}: Built-in plotting functions
  \item \strong{Bootstrap Validation}: Robust statistical assessment
}
}

\subsection{Algorithm Overview}{
scICER addresses the fundamental question: "How many clusters should I choose?" 
by evaluating clustering consistency across multiple runs:

\enumerate{
  \item \strong{Multiple Clustering}: Tests different cluster numbers with multiple trials
  \item \strong{Binary Search}: Finds appropriate resolution parameters automatically
  \item \strong{ECS Calculation}: Compares clustering results efficiently
  \item \strong{Optimization}: Iteratively improves clustering consistency
  \item \strong{Bootstrap Validation}: Assesses statistical significance
  \item \strong{Threshold Application}: Identifies reliable cluster numbers
}
}

\subsection{Main Functions}{
\describe{
  \item{\code{\link{scICE_clustering}}}{Main function for clustering consistency evaluation}
  \item{\code{\link{plot_ic}}}{Visualize inconsistency coefficients}
  \item{\code{\link{get_robust_labels}}}{Extract consistent clustering labels}
  \item{\code{\link{extract_consistent_clusters}}}{Summarize consistent clusters}
  \item{\code{\link{plot_stability}}}{Display bootstrap stability analysis}
  \item{\code{\link{calculate_ecs}}}{Calculate Element-Centric Similarity}
  \item{\code{\link{calculate_ic}}}{Calculate inconsistency scores}
  \item{\code{\link{calculate_mei_from_array}}}{Calculate Mutual Element-wise Information}
}
}

\subsection{Typical Workflow}{
\preformatted{
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

# Visualize results
plot_ic(results)
plot_stability(results)

# Extract consistent clusters
seurat_obj <- get_robust_labels(results, return_seurat = TRUE)

# Continue with Seurat workflow
DimPlot(seurat_obj, group.by = "clusters_8")
}
}

\subsection{Interpretation Guidelines}{
\strong{Inconsistency Coefficient (IC) values:}
\itemize{
  \item \strong{IC < 1.005}: Highly consistent (recommended)
  \item \strong{IC 1.005-1.01}: Moderately consistent  
  \item \strong{IC 1.01-1.02}: Low consistency (use with caution)
  \item \strong{IC > 1.02}: Poor consistency (likely overfitting)
}
}

\subsection{Performance Tips}{
\itemize{
  \item Use \code{n_workers = detectCores() - 1} for maximum speed
  \item Install the \code{leiden} package for improved performance
  \item Reduce \code{n_trials} and \code{n_bootstrap} for large datasets
  \item Focus \code{cluster_range} based on biological expectations
}
}
}
\author{
\packageAuthor{scICER}

Maintainer: \packageMaintainer{scICER}
}
\references{
Original scICE algorithm: [Citation to be added upon publication]

Related resources:
\itemize{
  \item Seurat: \url{https://satijalab.org/seurat/}
  \item Original scICE (Julia): \url{https://github.com/Mathbiomed/scICE}
}
}
\keyword{package}
\keyword{single-cell}
\keyword{clustering}
\seealso{
\code{\link{scICE_clustering}} for the main analysis function.

Useful links:
\itemize{
  \item \url{https://github.com/ATPs/scICER}
}
}
\examples{
\dontrun{
# Quick start example
library(scICER)
library(Seurat)

# Load example data
data("pbmc_small")

# Ensure preprocessing is complete  
pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- RunPCA(pbmc_small)
pbmc_small <- FindNeighbors(pbmc_small, dims = 1:10)

# Run scICER analysis
results <- scICE_clustering(pbmc_small, cluster_range = 2:8, n_workers = 2)

# Visualize results
plot_ic(results)

# Extract consistent clusters
consistent_clusters <- extract_consistent_clusters(results)
print(consistent_clusters)

# Add to Seurat object
pbmc_small <- get_robust_labels(results, return_seurat = TRUE)

# View help for main functions
?scICE_clustering
?plot_ic
?get_robust_labels
}
} 