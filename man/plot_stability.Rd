\name{plot_stability}
\alias{plot_stability}
\title{Plot clustering stability across resolution parameters}
\usage{
plot_stability(scice_results, cluster_number = NULL)
}
\arguments{
  \item{scice_results}{Results object from scICE_clustering function. Must be a valid scICE object containing clustering results.}
  
  \item{cluster_number}{Optional integer specifying a specific cluster number to analyze. If NULL (default), plots stability across all cluster numbers.}
}
\value{
Returns a ggplot2 object containing either:
\itemize{
  \item When cluster_number is NULL:
    \itemize{
      \item Line plot showing IC scores vs. cluster numbers
      \item Points colored by resolution parameter
      \item Continuous color scale for resolution values
    }
  \item When cluster_number is specified:
    \itemize{
      \item Histogram of bootstrap IC scores
      \item Vertical line showing median IC score
      \item Summary statistics in subtitle
    }
}

The plot can be further customized using standard ggplot2 functions.
}
\description{
Creates visualizations to analyze clustering stability either across all cluster numbers or for a specific cluster number. The function provides two different types of plots depending on whether a specific cluster number is specified.

The function uses ggplot2 for high-quality graphics and includes careful attention to visual clarity and statistical summary.
}
\details{
The function provides two distinct visualization modes:

1. Overview Mode (cluster_number = NULL):
\itemize{
  \item Shows stability across all tested cluster numbers
  \item Color gradient indicates resolution parameter values
  \item Connects points to show trends
  \item Includes both IC scores and resolution information
}

2. Detailed Mode (specific cluster_number):
\itemize{
  \item Shows distribution of IC scores from bootstrap analysis
  \item Includes median line and summary statistics
  \item Provides detailed view of stability for one cluster number
}

Visual features:
\itemize{
  \item Clear color schemes for resolution parameters
  \item Appropriate binning for histograms
  \item Minimal theme with informative gridlines
  \item Properly scaled axes with padding
}

The function performs input validation for:
\itemize{
  \item Valid scICE object
  \item Presence of required data
  \item Valid cluster number if specified
  \item Availability of bootstrap results
}

Statistical summaries include:
\itemize{
  \item Median IC score
  \item Distribution shape
  \item Stability metrics
}
}
\examples{
\dontrun{
# Run scICE analysis
scice_results <- scICE_clustering(pbmc_small, cluster_range = 2:10)

# Plot stability across all cluster numbers
p1 <- plot_stability(scice_results)

# Plot stability for specific cluster number
p2 <- plot_stability(scice_results, cluster_number = 5)

# Save plots
ggplot2::ggsave("overall_stability.pdf", p1, width = 10, height = 6)
ggplot2::ggsave("cluster5_stability.pdf", p2, width = 8, height = 6)
}
}
\seealso{
\code{\link{scICE_clustering}}, \code{\link{plot_ic}}
}
%\export 