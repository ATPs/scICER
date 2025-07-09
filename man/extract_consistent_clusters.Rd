\name{extract_consistent_clusters}
\alias{extract_consistent_clusters}
\title{Extract consistent clusters summary}
\description{
Provides a summary of consistent clusters found by scICE analysis, including
cluster numbers, IC scores, and other relevant statistics.
}
\usage{
extract_consistent_clusters(scice_results, threshold = 1.005)
}
\arguments{
\item{scice_results}{Results object from scICE_clustering function}

\item{threshold}{IC threshold for determining consistent clusters (default: 1.005)}
}
\value{
A data frame summarizing consistent clusters with columns:
\describe{
  \item{cluster_number}{Number of clusters}
  \item{ic_median}{Median IC score from bootstrap samples}
  \item{resolution_parameter}{Leiden resolution parameter used}
  \item{n_iterations}{Number of optimization iterations}
  \item{consistency_rank}{Rank by consistency (1 = most consistent)}
}

Returns an empty data frame if no cluster numbers meet the consistency threshold.
}
\details{
This function identifies and summarizes cluster numbers that meet the specified
consistency threshold. It provides key statistics to help users understand
which cluster numbers produced the most consistent results.

The summary includes:
\itemize{
  \item \strong{Cluster number}: The number of clusters tested
  \item \strong{IC median}: The median inconsistency score from bootstrap samples
  \item \strong{Resolution parameter}: The Leiden resolution used to achieve this clustering
  \item \strong{Iterations}: Number of optimization iterations required
  \item \strong{Consistency rank}: Ranking from most consistent (1) to least consistent
}

If no clusters meet the threshold, consider:
\itemize{
  \item Using a more lenient threshold (e.g., 1.01 instead of 1.005)
  \item Expanding the cluster range in the original analysis
  \item Reviewing data preprocessing steps
  \item Checking if the dataset has clear cluster structure
}
}
\examples{
\dontrun{
# Run scICE analysis
scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Extract consistent clusters with strict threshold
strict_clusters <- extract_consistent_clusters(scice_results, threshold = 1.005)
print(strict_clusters)

# Try with more lenient threshold if needed
if (nrow(strict_clusters) == 0) {
  lenient_clusters <- extract_consistent_clusters(scice_results, threshold = 1.01)
  print("With more lenient threshold:")
  print(lenient_clusters)
}

# View the most consistent clustering
if (nrow(strict_clusters) > 0) {
  best_cluster_num <- strict_clusters$cluster_number[1]
  cat("Most consistent clustering has", best_cluster_num, "clusters")
}
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{get_robust_labels}}, \code{\link{plot_ic}}
}
\author{
Xiaolong Cao
}
\keyword{clustering}
\keyword{summary} 