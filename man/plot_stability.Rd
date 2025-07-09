\name{plot_stability}
\alias{plot_stability}
\title{Plot stability analysis for scICE clustering results}
\description{
Creates visualization of clustering stability across bootstrap samples, 
showing the distribution of IC scores for each cluster number or detailed 
analysis for a specific cluster number.
}
\usage{
plot_stability(scice_results, cluster_number = NULL)
}
\arguments{
\item{scice_results}{Results object from scICE_clustering function}

\item{cluster_number}{Specific cluster number to analyze in detail (optional)}
}
\value{
A ggplot object showing stability analysis. If \code{cluster_number} is specified,
returns a detailed plot for that specific cluster number. Otherwise, returns
an overview plot of all tested cluster numbers.
}
\details{
The stability plot helps assess the robustness of clustering results by showing
the distribution of IC scores across bootstrap samples. This provides insight
into the statistical significance of the consistency measurements.

When \code{cluster_number = NULL} (default):
\itemize{
  \item Shows boxplots of bootstrap IC distributions for all cluster numbers
  \item Helpful for comparing stability across different cluster numbers
  \item Identifies which cluster numbers have the most stable results
}

When a specific \code{cluster_number} is provided:
\itemize{
  \item Shows detailed distribution for that cluster number
  \item Includes histogram and density plots
  \item Shows confidence intervals and summary statistics
}

Interpretation guidelines:
\itemize{
  \item Tight distributions indicate stable clustering
  \item Wide distributions suggest high variability
  \item Lower median values indicate better consistency
}
}
\examples{
\dontrun{
# Run scICE analysis
scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Overview stability plot
plot_stability(scice_results)

# Detailed analysis for 7 clusters
plot_stability(scice_results, cluster_number = 7)

# Compare stability of different cluster numbers
stability_overview <- plot_stability(scice_results)
print(stability_overview)
}
}
\seealso{
\code{\link{plot_ic}}, \code{\link{scICE_clustering}}, \code{\link{extract_consistent_clusters}}
}
\author{
Xiaolong Cao
}
\keyword{visualization}
\keyword{stability} 