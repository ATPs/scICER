\name{sclens}
\alias{sclens}
\title{Single-cell Linear Embedding of Neighborhoods and Signals (scLENS)}
\usage{
sclens(
  seurat_obj,
  assay = "RNA",
  slot = "counts",
  th = 60,
  p_step = 0.001,
  n_perturb = 20,
  centering = "mean",
  reduction_name_all = "sclens_pca_all",
  reduction_name_filtered = "sclens_pca_filtered",
  n_threads = 1,
  verbose = TRUE
)
}
\arguments{
\item{seurat_obj}{Seurat object containing single-cell RNA-seq data}

\item{assay}{Which assay to use (default: "RNA")}

\item{slot}{Which slot to use (default: "counts")}

\item{th}{Threshold angle in degrees for signal robustness test (default: 60)}

\item{p_step}{Decrement level for sparsity in robustness test (default: 0.001)}

\item{n_perturb}{Number of perturbations for robustness test (default: 20)}

\item{centering}{Centering method ("mean" or "median", default: "mean")}

\item{reduction_name_all}{Name for all PCs reduction (default: "sclens_pca_all")}

\item{reduction_name_filtered}{Name for filtered PCs reduction (default: "sclens_pca_filtered")}

\item{n_threads}{Number of threads to use for parallel processing (default: 1)}

\item{verbose}{Whether to print detailed progress messages (default: TRUE)}
}
\value{
Modified Seurat object with scLENS reductions added. The function adds two new dimensional reductions:
\describe{
  \item{reduction_name_all}{Contains all signals identified after Random Matrix Theory filtering}
  \item{reduction_name_filtered}{Contains only robust signals that passed the perturbation testing}
}

Additional metadata is stored in \code{seurat_obj@misc$sclens_results} containing:
\describe{
  \item{signal_eigenvalues}{Eigenvalues of signal components}
  \item{all_eigenvalues}{All eigenvalues from the analysis}
  \item{mp_eigenvalues}{Marchenko-Pastur distributed eigenvalues}
  \item{lambda_c}{Tracy-Widom threshold for signal detection}
  \item{robustness_scores}{Robustness scores for each signal}
  \item{robust_signal_indices}{Indices of signals that passed robustness testing}
  \item{n_signals_total}{Total number of signals detected}
  \item{n_signals_robust}{Number of robust signals}
}
}
\description{
scLENS is a dimensionality reduction method for single-cell RNA-seq data that uses
Random Matrix Theory (RMT) to distinguish signal from noise and performs robustness
testing to identify stable signals. Converted from Julia implementation.
}
\details{
scLENS performs the following steps:
\enumerate{
  \item Data preprocessing with log1p transformation and normalization
  \item Random Matrix Theory (RMT) filtering to identify signal eigenvalues
  \item Robustness testing through multiple perturbations
  \item Construction of dimensionally-reduced representations
}

The function adds two reductions to the Seurat object:
\itemize{
  \item All signals after RMT filtering
  \item Only robust signals after perturbation testing
}

The method uses Random Matrix Theory to establish a statistical threshold for distinguishing
signal eigenvalues from noise. The Marchenko-Pastur distribution is used to model the null
hypothesis (noise-only) eigenvalue distribution, and the Tracy-Widom statistic provides
the threshold for signal detection.

Robustness testing involves multiple perturbations of the data matrix through addition of
random non-zero elements, followed by correlation analysis of the resulting eigenvectors
to identify stable signal directions.

For parallel processing, the function can utilize multiple CPU cores during the perturbation
testing phase, which is typically the most computationally intensive step.
}
\examples{
\dontrun{
# Load Seurat object
data(pbmc_small)

# Basic scLENS analysis
pbmc_small <- sclens(pbmc_small)

# Check available reductions
Reductions(pbmc_small)

# Plot robust signals
DimPlot(pbmc_small, reduction = "sclens_pca_filtered")

# Advanced usage with custom parameters
pbmc_small <- sclens(
  pbmc_small, 
  th = 45,                    # More stringent robustness threshold
  n_perturb = 30,            # More perturbations for better stability
  n_threads = 4,             # Use 4 CPU cores
  centering = "median"       # Use median centering
)

# Access detailed results
sclens_results <- pbmc_small@misc$sclens_results
print(paste("Total signals:", sclens_results$n_signals_total))
print(paste("Robust signals:", sclens_results$n_signals_robust))

# Plot eigenvalue distribution
plot(sclens_results$all_eigenvalues, 
     main = "Eigenvalue Distribution",
     xlab = "Component", 
     ylab = "Eigenvalue")
abline(h = sclens_results$lambda_c, col = "red", lty = 2)
}
}
\seealso{
\code{\link{RunPCA}}, \code{\link{RunUMAP}}, \code{\link{DimPlot}}
}
\author{
Converted from Julia implementation by the scLENS development team
} 