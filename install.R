# Install required packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

# Source the dependency installation script
source("R/install_scICE.R")

# Install dependencies (force non-interactive mode)
install_scICE(install_optional = TRUE)

# Set up devtools to work non-interactively
Sys.setenv(R_TESTS = "")
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Install the package
message("Installing scICER package...")

# Get the current directory
pkg_dir <- getwd()

# Install the package from the current directory
install.packages(pkg_dir, repos = NULL, type = "source")

message("scICER package installed successfully!") 