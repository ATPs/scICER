\name{plot_ic}
\alias{plot_ic}
\title{Plot Inconsistency (IC) scores across different cluster numbers}
\description{
Creates a boxplot showing the distribution of IC scores for each cluster number 
tested. Lower IC scores indicate more consistent clustering results.
}
\usage{
plot_ic(scice_results, threshold = 1.005, figure_size = c(10, 6),
       title = "Clustering Consistency Analysis", show_threshold = TRUE)
}
\arguments{
\item{scice_results}{Results object from scICE_clustering function}

\item{threshold}{IC threshold line to display (default: 1.005)}

\item{figure_size}{Figure size as c(width, height) (default: c(10, 6))}

\item{title}{Plot title (default: "Clustering Consistency Analysis")}

\item{show_threshold}{Whether to show threshold line (default: TRUE)}
}
\value{
A ggplot object showing IC score distributions across cluster numbers.
}
\details{
The plot displays boxplots of IC scores for each tested cluster number, with 
individual bootstrap samples shown as points. A horizontal dashed line indicates 
the consistency threshold. Cluster numbers with IC scores below the threshold 
are highlighted in green, indicating good consistency.

Key interpretation guidelines:
\itemize{
  \item IC < 1.005: Highly consistent clustering (recommended)
  \item IC 1.005-1.01: Moderately consistent clustering
  \item IC > 1.01: Low consistency (use with caution)
}
}
\examples{
\dontrun{
# Run scICE analysis
scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Basic IC plot
plot_ic(scice_results)

# Customize plot
plot_ic(scice_results, threshold = 1.01, title = "My Clustering Analysis")

# Hide threshold line
plot_ic(scice_results, show_threshold = FALSE)
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{plot_stability}}, \code{\link{extract_consistent_clusters}}
}
\author{
Xiaolong Cao
}
\keyword{visualization}
\keyword{clustering} 