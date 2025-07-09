\name{get_recommended_parameters}
\alias{get_recommended_parameters}
\title{Get recommended parameters for scICER analysis based on dataset size}
\description{
Provides parameter recommendations for scICER analysis based on the number of cells
in the dataset. Helps optimize performance and accuracy for different dataset sizes.
}
\usage{
get_recommended_parameters(n_cells, analysis_type = "standard")
}
\arguments{
\item{n_cells}{Number of cells in the dataset}

\item{analysis_type}{Type of analysis: "quick", "standard", or "thorough" (default: "standard")}
}
\value{
List of recommended parameters including:
\describe{
  \item{dataset_category}{Dataset size category (small/medium/large/very_large)}
  \item{cluster_range}{Recommended range of cluster numbers to test}
  \item{n_workers}{Recommended number of parallel workers}
  \item{n_trials}{Recommended number of clustering trials}
  \item{n_bootstrap}{Recommended number of bootstrap samples}
  \item{objective_function}{Recommended objective function (modularity or CPM)}
  \item{ic_threshold}{Recommended IC threshold}
  \item{estimated_runtime}{Estimated analysis runtime}
  \item{recommendations}{Additional performance recommendations}
}
}
\details{
This function categorizes datasets into size groups and provides optimized 
parameters for each category:

\strong{Dataset Categories:}
\itemize{
  \item \strong{Small} (<1,000 cells): Uses modularity, fewer trials for speed
  \item \strong{Medium} (1,000-10,000 cells): Balanced parameters with CPM
  \item \strong{Large} (10,000-50,000 cells): Optimized for performance
  \item \strong{Very Large} (>50,000 cells): Minimal trials, maximum parallelization
}

\strong{Analysis Types:}
\itemize{
  \item \strong{Quick}: Reduced trials and bootstrap samples for fast results
  \item \strong{Standard}: Balanced accuracy and performance (recommended)
  \item \strong{Thorough}: Maximum trials and bootstrap samples for highest accuracy
}

The recommendations automatically adjust the number of parallel workers based on
available CPU cores and provide runtime estimates to help with planning.
}
\examples{
\dontrun{
# Get recommendations for different dataset sizes
small_params <- get_recommended_parameters(1000)
medium_params <- get_recommended_parameters(10000)
large_params <- get_recommended_parameters(100000)

# View recommendations
print(small_params$recommendations)
print(paste("Estimated runtime:", small_params$estimated_runtime))

# Use recommendations for analysis
results <- scICE_clustering(
  seurat_obj, 
  cluster_range = small_params$cluster_range,
  n_workers = small_params$n_workers,
  n_trials = small_params$n_trials,
  n_bootstrap = small_params$n_bootstrap,
  objective_function = small_params$objective_function
)

# Quick analysis for exploration
quick_params <- get_recommended_parameters(5000, "quick")
quick_results <- scICE_clustering(
  seurat_obj,
  cluster_range = quick_params$cluster_range,
  n_workers = quick_params$n_workers,
  n_trials = quick_params$n_trials,
  n_bootstrap = quick_params$n_bootstrap
)

# Thorough analysis for final results
thorough_params <- get_recommended_parameters(5000, "thorough")
final_results <- scICE_clustering(
  seurat_obj,
  cluster_range = thorough_params$cluster_range,
  n_workers = thorough_params$n_workers,
  n_trials = thorough_params$n_trials,
  n_bootstrap = thorough_params$n_bootstrap
)
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{check_seurat_ready}}
}
\author{
Xiaolong Cao
}
\keyword{parameters}
\keyword{optimization} 