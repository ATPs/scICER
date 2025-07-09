\name{create_results_summary}
\alias{create_results_summary}
\title{Create a summary report of scICER results}
\description{
Generates a comprehensive text summary of scICER analysis results,
including consistent clusters, parameter settings, and interpretation guidance.
}
\usage{
create_results_summary(scice_results, threshold = 1.005)
}
\arguments{
\item{scice_results}{Results object from scICE_clustering function}

\item{threshold}{IC threshold used for consistency (default: 1.005)}
}
\value{
Character vector with formatted summary report. Each element represents
a line of the report.
}
\details{
This function creates a comprehensive summary report that includes:

\strong{Report Sections:}
\itemize{
  \item \strong{Header}: Analysis date and package version
  \item \strong{Parameters}: Settings used for the analysis
  \item \strong{Results Overview}: Number of consistent clusters found
  \item \strong{IC Statistics}: Summary statistics for all tested cluster numbers
  \item \strong{Detailed Results}: Specific information for each consistent cluster
  \item \strong{Interpretation}: Guidelines for understanding IC scores
  \item \strong{Recommendations}: Specific advice based on the results
  \item \strong{Next Steps}: Commands for visualization and label extraction
}

The summary is designed to be both human-readable and suitable for saving
to a file for documentation purposes. It provides context-specific
recommendations based on whether consistent clusters were found and how
many were identified.

\strong{IC Score Interpretation:}
\itemize{
  \item IC < 1.005: Highly consistent (recommended)
  \item IC 1.005-1.01: Moderately consistent
  \item IC 1.01-1.02: Low consistency (use with caution)
  \item IC > 1.02: Poor consistency (likely overfitting)
}
}
\examples{
\dontrun{
# Run scICER analysis
results <- scICE_clustering(seurat_obj, cluster_range = 2:10)

# Generate and display summary
summary_report <- create_results_summary(results)
cat(paste(summary_report, collapse = "\n"))

# Save summary to file
writeLines(summary_report, "scICER_analysis_summary.txt")

# Generate summary with different threshold
lenient_summary <- create_results_summary(results, threshold = 1.01)
writeLines(lenient_summary, "scICER_lenient_summary.txt")

# Use in combination with other functions
validation <- check_seurat_ready(seurat_obj)
if (validation$ready) {
  params <- get_recommended_parameters(ncol(seurat_obj))
  results <- scICE_clustering(seurat_obj, 
                             cluster_range = params$cluster_range,
                             n_workers = params$n_workers)
  summary_report <- create_results_summary(results)
  
  # Print key findings
  cat("Analysis Summary:\n")
  cat(paste(summary_report[1:20], collapse = "\n"))
}
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{extract_consistent_clusters}}, \code{\link{plot_ic}}
}
\author{
Xiaolong Cao
}
\keyword{summary}
\keyword{reporting} 