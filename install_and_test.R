#!/usr/bin/env Rscript

#' Installation and Testing Script for scICER
#' 
#' This script helps install and test the scICER package locally
#' It handles dependency installation and basic functionality testing

# Set working directory to package root
if (basename(getwd()) != "scICER") {
  if (dir.exists("scICER")) {
    setwd("scICER")
  } else {
    stop("Please run this script from the scICER directory or its parent directory")
  }
}

cat("=================================================\n")
cat("scICER Package Installation and Testing Script\n")
cat("=================================================\n\n")

# Check R version
if (getRversion() < "4.0.0") {
  stop("scICER requires R >= 4.0.0. Current version: ", getRversion())
}

cat("✓ R version check passed:", as.character(getRversion()), "\n")

# Function to install packages if not already installed
install_if_missing <- function(packages, from_cran = TRUE, from_bioc = FALSE) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      cat("Installing", pkg, "...\n")
      if (from_bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager")
        }
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
    } else {
      cat("✓", pkg, "already installed\n")
    }
  }
}

# Install required dependencies
cat("\n--- Installing Dependencies ---\n")

# Core dependencies
core_deps <- c(
  "devtools", "roxygen2", "knitr", "rmarkdown",
  "Seurat", "SeuratObject", "igraph", "Matrix",
  "parallel", "foreach", "doParallel", "dplyr",
  "ggplot2", "testthat"
)

install_if_missing(core_deps)

# Optional dependencies for better performance
optional_deps <- c("leiden", "gridExtra")

cat("\n--- Installing Optional Dependencies ---\n")
for (pkg in optional_deps) {
  tryCatch({
    install_if_missing(pkg)
  }, error = function(e) {
    cat("⚠ Failed to install", pkg, "- this is optional\n")
  })
}

# Load devtools
library(devtools)

cat("\n--- Building and Installing scICER ---\n")

# Generate documentation
cat("Generating documentation...\n")
document()

# Check package
cat("Checking package...\n")
check_results <- check(quiet = TRUE)

if (length(check_results$errors) > 0) {
  cat("❌ Package check failed with errors:\n")
  print(check_results$errors)
} else if (length(check_results$warnings) > 0) {
  cat("⚠ Package check completed with warnings:\n")
  print(check_results$warnings)
} else {
  cat("✓ Package check passed!\n")
}

# Install the package
cat("Installing scICER...\n")
install(quiet = TRUE)

cat("\n--- Testing Installation ---\n")

# Test basic functionality
tryCatch({
  library(scICER)
  cat("✓ scICER loaded successfully\n")
  
  # Test with example data
  library(Seurat)
  data("pbmc_small")
  
  # Ensure preprocessing
  if (!"pca" %in% names(pbmc_small@reductions)) {
    pbmc_small <- NormalizeData(pbmc_small, verbose = FALSE)
    pbmc_small <- FindVariableFeatures(pbmc_small, verbose = FALSE)
    pbmc_small <- ScaleData(pbmc_small, verbose = FALSE)
    pbmc_small <- RunPCA(pbmc_small, verbose = FALSE)
  }
  
  if (!"snn" %in% names(pbmc_small@graphs)) {
    pbmc_small <- FindNeighbors(pbmc_small, dims = 1:10, verbose = FALSE)
  }
  
  cat("✓ Example data prepared\n")
  
  # Run quick test
  cat("Running quick functionality test...\n")
  test_results <- scICE_clustering(
    object = pbmc_small,
    cluster_range = 2:4,
    n_workers = 1,
    n_trials = 3,
    n_bootstrap = 5,
    verbose = FALSE
  )
  
  cat("✓ Basic clustering test passed\n")
  
  # Test visualization
  p <- plot_ic(test_results)
  cat("✓ Visualization functions work\n")
  
  # Test label extraction
  labels <- get_robust_labels(test_results, threshold = 1.1)  # Lenient threshold for test
  if (!is.null(labels)) {
    cat("✓ Label extraction works\n")
  } else {
    cat("⚠ No consistent clusters found (this is normal for small test)\n")
  }
  
}, error = function(e) {
  cat("❌ Testing failed with error:\n")
  print(e)
  cat("\nThis might be due to missing dependencies or system configuration.\n")
  cat("Please check the error message above and install any missing packages.\n")
})

cat("\n--- Installation Summary ---\n")
cat("scICER installation and testing completed!\n\n")

cat("Quick start example:\n")
cat("library(scICER)\n")
cat("library(Seurat)\n")
cat("data('pbmc_small')\n")
cat("seurat_obj <- FindNeighbors(pbmc_small, dims = 1:10)\n")
cat("results <- scICE_clustering(seurat_obj, cluster_range = 2:8)\n")
cat("plot_ic(results)\n\n")

cat("For detailed examples, see:\n")
cat("- examples/scICE_example.R\n")
cat("- vignette('scICER-quickstart')\n\n")

cat("If you encounter issues:\n")
cat("1. Ensure all dependencies are installed\n")
cat("2. Check that your R version is >= 4.0.0\n")
cat("3. For performance, install: install.packages('leiden')\n")
cat("4. Report issues at: https://github.com/ATPs/scICER/issues\n")

cat("\n=================================================\n")
cat("Installation Complete!\n")
cat("=================================================\n") 