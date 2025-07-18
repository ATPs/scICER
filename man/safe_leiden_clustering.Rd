% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/leiden_wrapper.R
\name{safe_leiden_clustering}
\alias{safe_leiden_clustering}
\title{Run clustering with error handling and retries}
\usage{
safe_leiden_clustering(
  igraph_obj,
  resolution,
  objective_function,
  n_iterations,
  beta,
  initial_membership = NULL,
  max_retries = 3
)
}
\arguments{
\item{igraph_obj}{igraph object}

\item{resolution}{Resolution parameter}

\item{objective_function}{Objective function}

\item{n_iterations}{Number of iterations}

\item{beta}{Beta parameter}

\item{initial_membership}{Initial membership}

\item{max_retries}{Maximum number of retries}
}
\value{
Cluster assignments
}
\description{
Run clustering with error handling and retries
}
\keyword{internal}
