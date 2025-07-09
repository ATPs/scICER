\name{plot_ic}
\alias{plot_ic}
\title{Plot Inconsistency (IC) scores across different cluster numbers}
\usage{
plot_ic(scice_results, threshold = 1.005, figure_size = c(10, 6),
        title = "Clustering Consistency Analysis", show_threshold = TRUE)
}
\arguments{
  \item{scice_results}{Results object from scICE_clustering function. Must be a valid scICE object containing clustering results.}
  
  \item{threshold}{Numeric value specifying the IC threshold line to display (default: 1.005). Lower values are more stringent.}
  
  \item{figure_size}{Numeric vector of length 2 specifying the figure dimensions as c(width, height) in inches (default: c(10, 6)).}
  
  \item{title}{Character string specifying the plot title (default: "Clustering Consistency Analysis").}
  
  \item{show_threshold}{Logical indicating whether to display the threshold line (default: TRUE).}
}
\value{
Returns a ggplot2 object containing:
\itemize{
  \item Boxplots showing IC score distributions for each cluster number
  \item Individual points showing all IC scores (with jitter)
  \item Optional threshold line in red
  \item Color-coded boxes: green for consistent clusters (IC < threshold), gray for inconsistent clusters
  \item Legend indicating consistency status
  \item Information about excluded clusters in subtitle (if applicable)
  \item Customized theme and labels
}

The plot can be further customized using standard ggplot2 functions.
}
\description{
Creates a comprehensive visualization of clustering consistency using boxplots and individual points. The plot shows the distribution of IC scores for each tested cluster number, with optional threshold line and highlighting of consistent clusters.

The function uses ggplot2 for high-quality graphics and includes careful attention to visual clarity and information density.
}
\details{
The plot includes several visual elements:
\enumerate{
  \item Boxplots showing the distribution of IC scores
  \item Individual points with jitter to show all data
  \item Optional threshold line with label
  \item Color-coding to distinguish consistent vs inconsistent clusters
  \item Clear labels and gridlines
  \item Legend explaining color coding
  \item Information about excluded clusters (if any)
}

Visual features:
\itemize{
  \item Green boxplots for consistent clusters (IC < threshold)
  \item Gray boxplots for inconsistent clusters (IC ≥ threshold)
  \item Red dashed line for threshold
  \item Legend at bottom explaining color coding
  \item Minimal theme with removed minor gridlines
  \item Properly scaled axes with padding
  \item Subtitle includes information about excluded clusters
}

The function automatically handles both old result formats (without exclusion information) and new result formats (with detailed exclusion information). Only clusters that were successfully tested and have IC scores are displayed in the plot.

Customization options:
\itemize{
  \item Adjust figure size for different display contexts
  \item Modify threshold value and visibility
  \item Change title and subtitle text
  \item Control point and line aesthetics
}

The function performs input validation and provides informative error messages for:
\itemize{
  \item Invalid scICE objects
  \item Missing IC results
  \item Invalid parameter values
}
}
\examples{
\dontrun{
# Run scICE analysis
scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Create basic plot
plot_ic(scice_results)

# Customize plot
plot_ic(scice_results, 
        threshold = 1.01,
        figure_size = c(12, 8),
        title = "My Custom Analysis",
        show_threshold = TRUE)

# Save plot
ggplot2::ggsave("clustering_results.pdf", width = 10, height = 6)
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{extract_consistent_clusters}}
}
%\export 