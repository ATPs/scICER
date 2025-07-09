\name{calculate_mei_from_array}
\alias{calculate_mei_from_array}
\title{Calculate Mutual Element-wise Information (MEI) from clustering array}
\description{
Calculate Mutual Element-wise Information (MEI) scores from extracted clustering 
results. MEI provides cell-level consistency scores showing how reliably each 
cell is clustered across multiple trials.
}
\usage{
calculate_mei_from_array(clustering_array)
}
\arguments{
\item{clustering_array}{Result from extract_clustering_array function containing
unique clustering patterns and their probabilities}
}
\value{
Vector of MEI scores, one for each cell. Higher scores indicate more consistent 
clustering for that cell across trials.
}
\details{
Mutual Element-wise Information (MEI) provides cell-level insight into clustering 
consistency. Unlike the global IC score, MEI shows which individual cells are 
consistently clustered together across multiple trials and which cells show 
variable cluster assignments.

The algorithm:
\enumerate{
  \item Calculates pairwise ECS similarities between all unique clustering patterns
  \item Weights similarities by pattern frequencies  
  \item Computes average similarity scores for each clustering pattern
  \item Returns pattern-wise consistency scores
}

MEI interpretation:
\itemize{
  \item \strong{High MEI}: Cell consistently assigned to similar clusters across trials
  \item \strong{Low MEI}: Cell shows variable cluster assignments (potential boundary cell)
  \item \strong{MEI = 1.0}: Perfect consistency for this clustering pattern
}

MEI scores are useful for:
\itemize{
  \item Identifying cells with uncertain cluster membership
  \item Finding potential boundary or transitional cells
  \item Assessing local clustering quality
  \item Filtering cells for downstream analysis based on consistency
}

If only one unique clustering pattern exists, all cells receive MEI = 1.0.
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

# Calculate MEI scores
mei_scores <- calculate_mei_from_array(clustering_array)
print(paste("Number of MEI scores:", length(mei_scores)))
print(paste("Mean MEI:", round(mean(mei_scores), 3)))

# Visualize MEI distribution
hist(mei_scores, 
     main = "Distribution of MEI Scores", 
     xlab = "MEI Score", 
     breaks = 20)

# Identify cells with low consistency
low_consistency_cells <- which(mei_scores < quantile(mei_scores, 0.1))
print(paste("Cells with lowest 10% MEI scores:", length(low_consistency_cells)))
}
}
\seealso{
\code{\link{calculate_ic}}, \code{\link{calculate_ecs}}, \code{\link{scICE_clustering}}
}
\author{
Xiaolong Cao
}
\keyword{clustering}
\keyword{consistency} 