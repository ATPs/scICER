\name{get_robust_labels}
\alias{get_robust_labels}
\title{Extract robust clustering labels from scICE results}
\description{
Extracts clustering labels for cluster numbers that meet the consistency threshold.
Returns a data frame with cell names and cluster assignments for each consistent
cluster number.
}
\usage{
get_robust_labels(scice_results, threshold = 1.005, return_seurat = FALSE)
}
\arguments{
\item{scice_results}{Results object from scICE_clustering function}

\item{threshold}{IC threshold for determining consistent clusters (default: 1.005)}

\item{return_seurat}{Whether to return results as Seurat metadata (default: FALSE)}
}
\value{
If \code{return_seurat = FALSE}: A data frame with cell names and cluster assignments 
for each consistent cluster number.

If \code{return_seurat = TRUE}: The original Seurat object with clustering results 
added to metadata as columns named "clusters_X" where X is the cluster number.

Returns \code{NULL} if no cluster numbers meet the consistency threshold.
}
\details{
This function identifies cluster numbers with IC scores below the specified threshold 
and extracts the corresponding clustering labels. The cluster labels are converted 
from 0-based (internal) to 1-based (R standard) indexing.

When \code{return_seurat = TRUE}, the function adds new metadata columns to the 
original Seurat object. These columns can then be used with standard Seurat 
functions like \code{DimPlot()} for visualization.

If no clusters meet the threshold, consider:
\itemize{
  \item Using a more lenient threshold (e.g., 1.01 instead of 1.005)
  \item Expanding the cluster_range in the original analysis
  \item Checking data preprocessing steps
}
}
\examples{
\dontrun{
# Run scICE analysis
scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Extract consistent clustering labels as data frame
consistent_labels <- get_robust_labels(scice_results, threshold = 1.005)
head(consistent_labels)

# Add results directly to Seurat object
pbmc_small <- get_robust_labels(scice_results, return_seurat = TRUE)

# View available clusterings
cluster_columns <- grep("clusters_", names(pbmc_small@meta.data), value = TRUE)
print(cluster_columns)

# Visualize with Seurat
library(Seurat)
DimPlot(pbmc_small, group.by = "clusters_7")

# Use for downstream analysis
Idents(pbmc_small) <- "clusters_7"
markers <- FindAllMarkers(pbmc_small)
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{extract_consistent_clusters}}, \code{\link{plot_ic}}
}
\author{
Xiaolong Cao
}
\keyword{clustering}
\keyword{results} 