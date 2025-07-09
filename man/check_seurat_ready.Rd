\name{check_seurat_ready}
\alias{check_seurat_ready}
\title{Check if Seurat object is ready for scICER analysis}
\description{
Validates that a Seurat object has the necessary preprocessing steps completed
for scICER analysis. Checks for normalization, variable features, scaling, PCA,
and neighbor graphs.
}
\usage{
check_seurat_ready(seurat_obj, graph_name = "snn")
}
\arguments{
\item{seurat_obj}{A Seurat object to validate}

\item{graph_name}{Name of the graph to check for (default: "snn")}
}
\value{
List with validation results:
\describe{
  \item{ready}{Logical indicating if object is ready for analysis}
  \item{issues}{Character vector of identified issues}
  \item{suggestions}{Character vector of suggested fixes}
  \item{summary}{Summary message about readiness status}
}
}
\details{
This function performs a comprehensive check of Seurat object preprocessing
to ensure compatibility with scICER analysis. It validates:

\itemize{
  \item \strong{Normalization}: Checks if data has been normalized
  \item \strong{Variable Features}: Verifies that variable features have been identified
  \item \strong{Scaling}: Ensures data has been scaled for PCA
  \item \strong{PCA}: Confirms that principal component analysis has been performed
  \item \strong{Neighbor Graph}: Validates that the specified graph exists
}

If any issues are found, the function provides specific suggestions for
resolving them using standard Seurat workflow commands.
}
\examples{
\dontrun{
library(Seurat)
data("pbmc_small")

# Check readiness before preprocessing
validation <- check_seurat_ready(pbmc_small)
print(validation$summary)

if (!validation$ready) {
  cat("Issues found:\n")
  cat(paste(validation$issues, collapse = "\n"))
  cat("\nSuggestions:\n")
  cat(paste(validation$suggestions, collapse = "\n"))
}

# Complete preprocessing
pbmc_small <- NormalizeData(pbmc_small)
pbmc_small <- FindVariableFeatures(pbmc_small)
pbmc_small <- ScaleData(pbmc_small)
pbmc_small <- RunPCA(pbmc_small)
pbmc_small <- FindNeighbors(pbmc_small, dims = 1:10)

# Check again
validation2 <- check_seurat_ready(pbmc_small)
print(validation2$summary)  # Should now be ready
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{get_recommended_parameters}}
}
\author{
Xiaolong Cao
}
\keyword{validation}
\keyword{preprocessing} 