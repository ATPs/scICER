\name{check_seurat_ready}
\alias{check_seurat_ready}
\title{Check if a Seurat object is ready for scICE analysis}
\usage{
check_seurat_ready(seurat_obj, graph_name = NULL)
}
\arguments{
  \item{seurat_obj}{A Seurat object to check. This should be a valid Seurat object with necessary preprocessing completed.}
  
  \item{graph_name}{Name of the graph to check for. If NULL (default), will check for the default SNN graph from the active assay (e.g., "RNA_snn" for RNA assay).}
}
\value{
Returns a list containing:
\itemize{
  \item{ready}{Logical indicating whether the object is ready for analysis}
  \item{issues}{Character vector of identified issues, if any}
}
}
\description{
Validates a Seurat object to ensure it has all necessary components for scICE analysis. This includes checking for required preprocessing steps, graph presence, and data validity.

The function performs comprehensive checks and provides detailed feedback about any missing or problematic components.
}
\details{
The function checks for several requirements:

1. Basic Object Validation:
\itemize{
  \item Valid Seurat object
  \item Non-empty data matrix
  \item Presence of required slots
}

2. Preprocessing Steps:
\itemize{
  \item Data normalization
  \item Feature selection
  \item Dimensionality reduction
  \item Neighbor graph computation
}

3. Graph Requirements:
\itemize{
  \item Presence of specified graph
  \item Graph connectivity
  \item Valid graph format
}

The function returns detailed messages about:
\itemize{
  \item Missing preprocessing steps
  \item Invalid or missing graphs
  \item Data quality issues
  \item Recommended actions
}
}
\examples{
\dontrun{
# Load example data
data(pbmc_small)

# Basic check with default graph
status <- check_seurat_ready(pbmc_small)
print(status$ready)
print(status$issues)

# Check with specific graph
status <- check_seurat_ready(pbmc_small, graph_name = "RNA_snn")

# Handle issues if found
if (!status$ready) {
  message("Issues found:")
  print(status$issues)
}
}
}
\seealso{
\code{\link{scICE_clustering}}
}
\export 