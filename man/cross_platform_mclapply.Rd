% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering_core.R
\name{cross_platform_mclapply}
\alias{cross_platform_mclapply}
\title{Cross-platform parallel lapply wrapper}
\usage{
cross_platform_mclapply(X, FUN, mc.cores = 1, ...)
}
\arguments{
\item{X}{Vector/list to iterate over}

\item{FUN}{Function to apply}

\item{mc.cores}{Number of cores (ignored on Windows)}

\item{...}{Additional arguments to FUN}
}
\value{
List of results
}
\description{
Cross-platform parallel lapply wrapper
}
\keyword{internal}
