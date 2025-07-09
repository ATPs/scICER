#' Installation script for scICE R package
#'
#' This script installs all required dependencies and sets up the scICE package
#' for single-cell clustering consistency evaluation.

# Function to check and install packages
install_if_missing <- function(packages, repo = "CRAN") {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing", pkg, "from", repo))
      if (repo == "CRAN") {
        install.packages(pkg, dependencies = TRUE)
      } else if (repo == "Bioconductor") {
        BiocManager::install(pkg)
      }
    } else {
      message(paste(pkg, "is already installed"))
    }
  }
}

# Main installation function
install_scICE <- function(install_optional = TRUE) {
  message("Starting scICE installation...")
  
  # Install BiocManager if needed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  
  # Required packages from CRAN
  required_cran <- c(
    "Seurat",
    "SeuratObject", 
    "igraph",
    "Matrix",
    "parallel",
    "foreach",
    "doParallel",
    "dplyr",
    "ggplot2",
    "devtools",
    "roxygen2"
  )
  
  # Optional packages for better performance
  optional_cran <- c(
    "leiden",          # Better Leiden clustering
    "viridis",         # Color scales
    "RColorBrewer",    # Color palettes
    "gridExtra",       # Plot arrangements
    "knitr",           # Documentation
    "rmarkdown",       # Documentation
    "testthat"         # Testing
  )
  
  # Install required packages
  message("Installing required packages...")
  install_if_missing(required_cran, "CRAN")
  
  # Install optional packages if requested
  if (install_optional) {
    message("Installing optional packages...")
    install_if_missing(optional_cran, "CRAN")
  }
  
  # Check if all required packages are available
  missing_packages <- c()
  for (pkg in required_cran) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    stop(paste("Failed to install required packages:", paste(missing_packages, collapse = ", ")))
  }
  
  message("All dependencies installed successfully!")
  message("You can now load scICE using: devtools::load_all('path/to/scICE')")
  message("Or install as a package using: devtools::install_local('path/to/scICE')")
  
  return(TRUE)
}

# Run installation if script is executed directly
if (interactive()) {
  install_scICE(install_optional = TRUE)
} else {
  # When sourced, just define the function
  message("scICE installation function loaded. Run install_scICE() to install dependencies.")
} 