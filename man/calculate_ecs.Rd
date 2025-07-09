\name{calculate_ecs}
\alias{calculate_ecs}
\title{Calculate Element-Centric Similarity (ECS) between two clustering results}
\description{
ECS calculates similarity between two clustering results based on how consistently 
each individual cell is grouped with other cells. This is more efficient than 
traditional methods that build consensus matrices.
}
\usage{
calculate_ecs(cluster_a, cluster_b, d = 0.9, return_vector = FALSE)
}
\arguments{
\item{cluster_a}{First clustering result (vector of cluster assignments)}

\item{cluster_b}{Second clustering result (vector of cluster assignments)}

\item{d}{Damping parameter (default: 0.9)}

\item{return_vector}{Whether to return similarity scores for each cell (default: FALSE)}
}
\value{
If \code{return_vector = FALSE}: Mean ECS score (single numeric value)

If \code{return_vector = TRUE}: Vector of ECS scores for each cell
}
\details{
Element-Centric Similarity (ECS) is an efficient algorithm for comparing clustering 
results that focuses on how consistently individual cells are grouped together. 
Unlike traditional methods that build full consensus matrices, ECS computes similarity 
scores directly for each cell based on its neighborhood consistency.

The algorithm works by:
\enumerate{
  \item Creating personalized PageRank vectors for each cell in both clusterings
  \item Computing the L1 distance between these vectors
  \item Converting distances to similarity scores
}

The damping parameter \code{d} controls the influence of cluster membership. 
Higher values (closer to 1) give more weight to direct cluster membership, 
while lower values consider broader neighborhood effects.

ECS scores range from 0 to 1, where:
\itemize{
  \item 1.0: Perfect agreement between clusterings
  \item 0.5: Random agreement
  \item 0.0: Complete disagreement
}
}
\examples{
\dontrun{
# Create example clustering results
n_cells <- 100
cluster_a <- sample(1:5, n_cells, replace = TRUE)
cluster_b <- sample(1:5, n_cells, replace = TRUE)

# Calculate mean ECS similarity
ecs_score <- calculate_ecs(cluster_a, cluster_b)
print(paste("ECS similarity:", round(ecs_score, 3)))

# Calculate cell-wise ECS scores
cell_scores <- calculate_ecs(cluster_a, cluster_b, return_vector = TRUE)
hist(cell_scores, main = "Distribution of Cell-wise ECS Scores")

# Test with identical clusterings (should give 1.0)
identical_score <- calculate_ecs(cluster_a, cluster_a)
print(paste("Identical clustering ECS:", identical_score))
}
}
\references{
Element-Centric Similarity algorithm from the scICE package
}
\seealso{
\code{\link{calculate_ic}}, \code{\link{scICE_clustering}}
}
\author{
Xiaolong Cao
}
\keyword{clustering}
\keyword{similarity} 