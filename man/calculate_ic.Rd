\name{calculate_ic}
\alias{calculate_ic}
\title{Calculate Inconsistency (IC) score from clustering results}
\description{
Calculate Inconsistency (IC) score from extracted clustering results.
Lower IC scores indicate more consistent clustering, with 1.0 representing 
perfect consistency.
}
\usage{
calculate_ic(clustering_array)
}
\arguments{
\item{clustering_array}{Result from extract_clustering_array function containing
unique clustering patterns and their probabilities}
}
\value{
IC score (numeric value). Lower values indicate better consistency:
\itemize{
  \item 1.0: Perfect consistency (all clustering trials identical)
  \item > 1.0: Increasing inconsistency
}
}
\details{
The Inconsistency Coefficient (IC) quantifies how consistent multiple clustering 
results are with each other. It is calculated using Element-Centric Similarity (ECS) 
between all pairs of unique clustering patterns, weighted by their frequencies.

The algorithm:
\enumerate{
  \item Extracts unique clustering patterns and their probabilities
  \item Calculates pairwise ECS similarity between all unique patterns
  \item Computes weighted average similarity based on pattern frequencies
  \item Returns the similarity score as the IC (1.0 = perfect consistency)
}

IC interpretation guidelines:
\itemize{
  \item \strong{IC = 1.0}: Perfect consistency - all trials produced identical clustering
  \item \strong{IC < 1.005}: Highly consistent - recommended for use
  \item \strong{IC 1.005-1.01}: Moderately consistent - reasonable choice
  \item \strong{IC 1.01-1.02}: Low consistency - use with caution
  \item \strong{IC > 1.02}: Poor consistency - likely indicates overfitting
}

This function is typically used internally by scICE but can be called directly
for custom clustering consistency analysis.
}
\examples{
\dontrun{
# This function is typically used internally, but can be used directly:

# Simulate multiple clustering results
n_cells <- 100
n_trials <- 10
clustering_matrix <- matrix(
  sample(1:5, n_cells * n_trials, replace = TRUE), 
  nrow = n_cells, ncol = n_trials
)

# Extract clustering array
clustering_array <- extract_clustering_array(clustering_matrix)

# Calculate IC score
ic_score <- calculate_ic(clustering_array)
print(paste("IC score:", round(ic_score, 4)))

# Perfect consistency example (all trials identical)
perfect_matrix <- matrix(
  rep(sample(1:5, n_cells, replace = TRUE), n_trials),
  nrow = n_cells, ncol = n_trials
)
perfect_array <- extract_clustering_array(perfect_matrix)
perfect_ic <- calculate_ic(perfect_array)
print(paste("Perfect consistency IC:", perfect_ic))  # Should be 1.0
}
}
\seealso{
\code{\link{calculate_ecs}}, \code{\link{scICE_clustering}}
}
\author{
Xiaolong Cao
}
\keyword{clustering}
\keyword{consistency} 