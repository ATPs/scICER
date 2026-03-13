# Install required packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

install_if_missing <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("Installing ", pkg)
      install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    } else {
      message(pkg, " is already installed")
    }
  }
}

required_cran <- c(
  "Seurat",
  "SeuratObject",
  "igraph",
  "Matrix",
  "parallel",
  "ggplot2",
  "devtools",
  "roxygen2"
)

optional_cran <- c(
  "leiden",
  "viridis",
  "RColorBrewer",
  "gridExtra",
  "knitr",
  "rmarkdown",
  "testthat"
)

message("Installing required packages...")
install_if_missing(required_cran)
message("Installing optional packages...")
install_if_missing(optional_cran)

Sys.setenv(R_TESTS = "")
options(repos = c(CRAN = "https://cloud.r-project.org"))

message("Installing scICER package...")
pkg_dir <- getwd()
install.packages(pkg_dir, repos = NULL, type = "source")
message("scICER package installed successfully!")
