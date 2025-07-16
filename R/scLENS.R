#' @import Seurat
#' @import Matrix
#' @import parallel
#' @import foreach
#' @import doParallel
#' @importFrom stats sd median mad quantile rnorm cor
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom methods as
#'
#' scLENS R Implementation
#' Converted from Julia implementation for single-cell RNA-seq data analysis
#' Uses Random Matrix Theory (RMT) for signal detection and robustness testing

#' Helper function: Convert Seurat object to sparse matrix
#' Julia equivalent: df2sparr(inp_df)
#' 
#' @param seurat_obj Seurat object
#' @param assay Which assay to use (default: "RNA")
#' @param slot Which slot to use (default: "counts")
#' @return Sparse matrix (genes x cells)
seurat_to_sparse <- function(seurat_obj, assay = "RNA", slot = "counts") {
  # Get the count matrix from Seurat object
  # Note: Seurat stores as genes x cells, Julia used cells x genes
  
  # Handle SeuratObject v5+ deprecation of slot parameter
  if (packageVersion("SeuratObject") >= "5.0.0") {
    # Use layer parameter for SeuratObject v5+
    count_matrix <- GetAssayData(seurat_obj, assay = assay, layer = slot)
  } else {
    # Use slot parameter for older versions
    count_matrix <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  }
  
  # Convert to dgCMatrix if not already
  if (!inherits(count_matrix, "dgCMatrix")) {
    count_matrix <- as(count_matrix, "dgCMatrix")
  }
  
  # Transpose to match Julia format (cells x genes)
  return(t(count_matrix))
}

#' Helper function: Project to simplex (L1 normalization)
#' Julia equivalent: proj_l = x -> issparse(x) ? spdiagm(1 ./sum(x,dims=2)[:]) * x : x ./ sum(x,dims=2)
#' 
#' @param x Matrix to normalize
#' @return L1 normalized matrix
proj_l <- function(x) {
  # Calculate row sums
  row_sums <- Matrix::rowSums(x)
  
  # Avoid division by zero
  row_sums[row_sums == 0] <- 1
  
  # Create diagonal matrix for normalization
  diag_norm <- Matrix::Diagonal(x = 1 / row_sums)
  
  return(diag_norm %*% x)
}

#' Helper function: Z-score with L2 normalization
#' Julia equivalent: zscore_with_l2(X)
#' 
#' @param X Matrix to normalize
#' @return Z-score normalized matrix with L2 correction
zscore_with_l2 <- function(X) {
  # Original Julia code:
  # std_ = std(X,dims=1)[:]
  # X_norm = X * spdiagm(1. ./ std_)
  # mu = mean(X_norm, dims=1)
  # l2X = sqrt.(sum(X_norm.^2,dims=2)[:])
  # l2mu = norm(mu)
  # l2norm_ = sqrt.(l2X.^2 .- 2 .* (X_norm * mu')[:] .+ l2mu^2)
  # (Matrix(X_norm) .- mu) ./ (l2norm_ / mean(l2norm_))
  
  # Calculate column standard deviations
  std_cols <- apply(X, 2, sd)
  std_cols[std_cols == 0] <- 1  # Avoid division by zero
  
  # Normalize by standard deviation
  X_norm <- X %*% Matrix::Diagonal(x = 1 / std_cols)
  
  # Calculate column means
  mu <- matrix(colMeans(X_norm), nrow = 1)
  
  # Calculate L2 norms
  l2X <- sqrt(rowSums(X_norm^2))
  l2mu <- norm(mu, type = "F")
  
  # Calculate L2 normalized distances
  mu_matrix <- matrix(rep(mu, nrow(X_norm)), nrow = nrow(X_norm), byrow = TRUE)
  l2norm_ <- sqrt(l2X^2 - 2 * rowSums(X_norm * mu_matrix) + l2mu^2)
  
  # Final normalization
  result <- (as.matrix(X_norm) - mu_matrix) / (l2norm_ / mean(l2norm_))
  
  return(result)
}

#' Helper function: Scale data with centering
#' Julia equivalent: scaled_gdata(X;dim=1, position_="mean")
#' 
#' @param X Matrix to scale
#' @param position_ Centering method ("mean" or "median")
#' @return Scaled matrix
scaled_gdata <- function(X, position_ = "mean") {
  # Original Julia code handles both mean and median centering
  # We implement mean centering here as it's the most common
  
  if (position_ == "mean") {
    # Calculate means and standard deviations by column
    col_means <- colMeans(X)
    col_sds <- apply(X, 2, sd)
    col_sds[col_sds == 0] <- 1  # Avoid division by zero
    
    # Center and scale
    X_centered <- sweep(X, 2, col_means, "-")
    X_scaled <- sweep(X_centered, 2, col_sds, "/")
    
    return(X_scaled)
  } else if (position_ == "median") {
    # Calculate medians and MADs by column
    col_medians <- apply(X, 2, median)
    col_mads <- apply(X, 2, mad)
    col_mads[col_mads == 0] <- 1  # Avoid division by zero
    
    # Center and scale
    X_centered <- sweep(X, 2, col_medians, "-")
    X_scaled <- sweep(X_centered, 2, col_mads, "/")
    
    return(X_scaled)
  } else {
    stop("position_ must be 'mean' or 'median'")
  }
}

#' Helper function: Generate random matrix for null model
#' Julia equivalent: random_nz(pre_df;rmix=true,mix_p=nothing)
#' 
#' @param X Sparse matrix
#' @param rmix Whether to randomize structure
#' @return Randomized sparse matrix
random_nz <- function(X, rmix = TRUE) {
  # Original Julia code:
  # tmp_X, return_mat = if typeof(pre_df) <: DataFrame
  #     df2sparr(pre_df), false
  # else
  #     sparse(pre_df), true
  # end
  # nz_row, nz_col, nz_val = findnz(tmp_X)
  # tmp_X.nzval .= shuffle(nz_val)
  
  # Get non-zero indices and values
  X_triplet <- summary(X)
  
  # Shuffle the values while keeping the sparsity pattern
  shuffled_values <- sample(X_triplet$x)
  
  # Create new sparse matrix with shuffled values
  X_random <- Matrix::sparseMatrix(
    i = X_triplet$i,
    j = X_triplet$j,
    x = shuffled_values,
    dims = dim(X)
  )
  
  if (rmix) {
    # Additional randomization of the matrix structure
    # This is a simplified version of the Julia _random_matrix function
    X_random <- X_random[sample(nrow(X_random)), ]
  }
  
  return(X_random)
}

#' Helper function: Compute Wishart matrix (covariance matrix)
#' Julia equivalent: _wishart_matrix(X;device="cpu")
#' 
#' @param X Input matrix
#' @return Wishart matrix (covariance matrix)
wishart_matrix <- function(X) {
  # Original Julia code:
  # X = if issparse(X)
  #     Matrix{eltype(X)}(X)
  # else
  #     X
  # end
  # out = Array{eltype(X),2}(undef, size(X,1),size(X,1))
  # mul!(out,X,X')
  # return out ./ size(X,2)
  
  # Convert sparse matrix to dense if needed for computation
  if (inherits(X, "sparseMatrix")) {
    X <- as.matrix(X)
  }
  
  # Compute X %*% t(X) / ncol(X)
  cov_matrix <- tcrossprod(X) / ncol(X)
  
  return(cov_matrix)
}

#' Helper function: Compute eigendecomposition
#' Julia equivalent: _get_eigen(Y;device="cpu")
#' 
#' @param Y Symmetric matrix
#' @return List with eigenvalues and eigenvectors
get_eigen <- function(Y) {
  # Original Julia code:
  # tmp_L, tmp_V = eigen(Y)
  # return tmp_L, tmp_V
  
  # Compute eigendecomposition
  eigen_result <- eigen(Y, symmetric = TRUE)
  
  return(list(
    values = eigen_result$values,
    vectors = eigen_result$vectors
  ))
}

#' Helper function: Marchenko-Pastur parameters
#' Julia equivalent: _mp_parameters(L)
#' 
#' @param L Vector of eigenvalues
#' @return List of Marchenko-Pastur parameters
mp_parameters <- function(L) {
  # Original Julia code:
  # moment_1 = mean(L)
  # moment_2 = mean(L.^2)
  # gamma = moment_2 / moment_1^2 - 1
  # s = moment_1
  # sigma = moment_2
  # b_plus = s * (1 + sqrt(gamma))^2
  # b_minus = s * (1 - sqrt(gamma))^2
  # x_peak = s * (1.0-gamma)^2.0 / (1.0+gamma)
  
  moment_1 <- mean(L)
  moment_2 <- mean(L^2)
  gamma <- moment_2 / moment_1^2 - 1
  
  # Simple fix: ensure gamma is non-negative for sqrt calculation
  if (!is.finite(gamma) || gamma < 0) gamma <- 0
  
  s <- moment_1
  sigma <- moment_2
  b_plus <- s * (1 + sqrt(gamma))^2
  b_minus <- s * (1 - sqrt(gamma))^2
  x_peak <- s * (1.0 - gamma)^2.0 / (1.0 + gamma)
  
  return(list(
    moment_1 = moment_1,
    moment_2 = moment_2,
    gamma = gamma,
    b_plus = b_plus,
    b_minus = b_minus,
    s = s,
    peak = x_peak,
    sigma = sigma
  ))
}

#' Helper function: Marchenko-Pastur calculation
#' Julia equivalent: _mp_calculation(L, Lr, eta=1, eps=1e-6, max_iter=10000)
#' 
#' @param L Original eigenvalues
#' @param Lr Random eigenvalues
#' @param eta Learning rate
#' @param eps Convergence threshold
#' @param max_iter Maximum iterations
#' @return List with filtered eigenvalues and bounds
mp_calculation <- function(L, Lr, eta = 1, eps = 1e-6, max_iter = 10000) {
  # Original Julia code implements an iterative algorithm to find the optimal
  # boundaries for the Marchenko-Pastur distribution
  
  converged <- FALSE
  iter <- 0
  mpp_Lr <- mp_parameters(Lr)
  b_plus <- mpp_Lr$b_plus
  b_minus <- mpp_Lr$b_minus
  
  L_updated <- L[L > b_minus & L < b_plus]

  # Add a check to handle cases with few or no eigenvalues in the initial MP bounds.
  if (length(L_updated) < 2) {
    warning("Initial filtering resulted in fewer than 2 eigenvalues to estimate the MP distribution. Skipping iterative calculation.")
    # Return the initial bounds from the random matrix.
    return(list(
      L_filtered = L_updated,
      b_plus = b_plus,
      b_minus = b_minus
    ))
  }

  new_mpp_L <- mp_parameters(L_updated)
  new_b_plus <- new_mpp_L$b_plus
  new_b_minus <- new_mpp_L$b_minus
  
  while (!converged) {
    loss <- (1 - new_b_plus / b_plus)^2
    iter <- iter + 1
    
    # Defensive Check: Is the loss value valid?
    if (is.na(loss) || loss <= eps) {
      converged <- TRUE
    } else if (iter >= max_iter) {
      cat("Max iterations exceeded!\n")
      converged <- TRUE
    } else {
      gradient <- new_b_plus - b_plus
      new_b_plus <- b_plus + eta * gradient
      L_updated <- L[L > new_b_minus & L < new_b_plus]
      b_plus <- new_b_plus
      b_minus <- new_b_minus
      
      # Defensive Check: Are there enough points to continue?
      if (length(L_updated) < 2) {
        warning("Iteration stopped because fewer than 2 eigenvalues remained in the MP bounds.")
        converged <- TRUE # Stop the loop gracefully
      } else {
        # Continue with calculations...
        up_mpp_L <- mp_parameters(L_updated)
        new_b_plus <- up_mpp_L$b_plus
        new_b_minus <- up_mpp_L$b_minus
      }
    }
  }
  
  final_L <- L[L > new_b_minus & L < new_b_plus]
  
  return(list(
    L_filtered = final_L,
    b_plus = new_b_plus,
    b_minus = new_b_minus
  ))
}

#' Helper function: Tracy-Widom calculation
#' Julia equivalent: _tw(L, L_mp)
#' 
#' @param L Original eigenvalues
#' @param L_mp Marchenko-Pastur eigenvalues
#' @return List with Tracy-Widom parameters
tw_calculation <- function(L, L_mp) {
  # Original Julia code:
  # gamma = _mp_parameters(L_mp)["gamma"]
  # p = length(L) / gamma
  # sigma = 1 / p^(2/3) * gamma^(5/6) * (1 + sqrt(gamma))^(4/3)
  # lambda_c = mean(L_mp) * (1 + sqrt(gamma))^2 + sigma
  
  mpp <- mp_parameters(L_mp)
  gamma <- mpp$gamma
  p <- length(L) / gamma
  sigma <- 1 / p^(2/3) * gamma^(5/6) * (1 + sqrt(gamma))^(4/3)
  lambda_c <- mean(L_mp) * (1 + sqrt(gamma))^2 + sigma
  
  return(list(
    lambda_c = lambda_c,
    gamma = gamma,
    p = p,
    sigma = sigma
  ))
}

#' Helper function: Get eigenvalues and eigenvectors
#' Julia equivalent: get_eigvec(X;device="gpu")
#' 
#' @param X Input matrix
#' @return List with eigenvalues and eigenvectors
get_eigvec <- function(X) {
  # Original Julia code:
  # N, M = size(X)
  # if N > M
  #     Y = _wishart_matrix(X',device=device)
  #     ...
  # else
  #     Y = _wishart_matrix(X,device=device)
  #     ...
  
  N <- nrow(X)
  M <- ncol(X)
  
  if (N > M) {
    # Use X' for computational efficiency
    Y <- wishart_matrix(t(X))
    eigen_result <- get_eigen(Y)
    L <- eigen_result$values
    V <- eigen_result$vectors
    
    # Filter positive eigenvalues
    positive_idx <- L > 0
    L <- L[positive_idx]
    V <- V[, positive_idx, drop = FALSE]
    
    # Sort by eigenvalue magnitude (descending)
    sort_idx <- order(L, decreasing = TRUE)
    nL <- L[sort_idx]
    nVs <- V[, sort_idx, drop = FALSE]
    
    # Project back to original space
    # This is a simplified version of the Julia projection
    mul_X <- nVs %*% diag(1 / sqrt(nL))
    new_nVs <- X %*% mul_X
    new_nVs <- apply(new_nVs, 2, function(x) x / norm(x, type = "2"))
    
    return(list(eigenvalues = nL, eigenvectors = new_nVs))
  } else {
    Y <- wishart_matrix(X)
    eigen_result <- get_eigen(Y)
    L <- eigen_result$values
    V <- eigen_result$vectors
    
    # Filter positive eigenvalues
    positive_idx <- L > 0
    L <- L[positive_idx]
    V <- V[, positive_idx, drop = FALSE]
    
    # Sort by eigenvalue magnitude (descending)
    sort_idx <- order(L, decreasing = TRUE)
    nL <- L[sort_idx]
    nVs <- V[, sort_idx, drop = FALSE]
    
    return(list(eigenvalues = nL, eigenvectors = nVs))
  }
}

#' Helper function: Get significant eigenvalues using RMT
#' Julia equivalent: get_sigev(X,Xr;device="gpu")
#' 
#' @param X Scaled data matrix
#' @param Xr Randomized data matrix
#' @return List with signal eigenvalues, eigenvectors, and RMT parameters
get_sigev <- function(X, Xr) {
  # Original Julia code implements Random Matrix Theory filtering
  # to distinguish signal from noise eigenvalues
  
  n <- nrow(X)
  m <- ncol(X)
  
  if (n > m) {
    # Compute covariance matrices
    Y <- wishart_matrix(t(X))
    Yr <- wishart_matrix(t(Xr))
    
    # Get eigendecompositions
    eigen_Y <- get_eigen(Y)
    eigen_Yr <- get_eigen(Yr)
    
    L <- eigen_Y$values
    V <- eigen_Y$vectors
    Lr <- eigen_Yr$values
    
    # Apply Marchenko-Pastur filtering
    mp_result <- mp_calculation(L, Lr[1:(length(Lr)-1)])
    L_mp <- mp_result$L_filtered
    b_min <- mp_result$b_minus
    
    # Calculate Tracy-Widom threshold
    tw_result <- tw_calculation(L, L_mp)
    lambda_c <- tw_result$lambda_c
    
    cat("Number of signal eigenvalues:", sum(L > lambda_c), "\n")
    
    # Extract signal eigenvalues and eigenvectors
    sel_L <- L[L > lambda_c]
    sel_Vs <- V[, L > lambda_c, drop = FALSE]
    
    # Extract noise eigenvalues and eigenvectors
    noise_idx <- L >= b_min & L <= lambda_c
    noise_L <- L[noise_idx]
    noise_V <- V[, noise_idx, drop = FALSE]
    
    # Sort by eigenvalue magnitude
    if (length(sel_L) > 0) {
      signal_sort_idx <- order(sel_L, decreasing = TRUE)
      nL <- sel_L[signal_sort_idx]
      nVs <- sel_Vs[, signal_sort_idx, drop = FALSE]
      
      # Project to original space
      mul_X <- nVs %*% diag(1 / sqrt(nL))
      new_nVs <- X %*% mul_X
      new_nVs <- apply(new_nVs, 2, function(x) x / norm(x, type = "2"))
    } else {
      nL <- numeric(0)
      new_nVs <- matrix(nrow = n, ncol = 0)
    }
    
    if (length(noise_L) > 0) {
      noise_sort_idx <- order(noise_L, decreasing = TRUE)
      snL <- noise_L[noise_sort_idx]
      noise_Vs_sorted <- noise_V[, noise_sort_idx, drop = FALSE]
      
      # Project noise eigenvectors
      mul_X2 <- noise_Vs_sorted %*% diag(sqrt(snL))
      new_noise_V <- X %*% mul_X2
      new_noise_V <- apply(new_noise_V, 2, function(x) x / norm(x, type = "2"))
    } else {
      snL <- numeric(0)
      new_noise_V <- matrix(nrow = n, ncol = 0)
    }
    
    return(list(
      signal_eigenvalues = nL,
      signal_eigenvectors = new_nVs,
      all_eigenvalues = L,
      mp_eigenvalues = L_mp,
      lambda_c = lambda_c,
      noise_eigenvalues = snL,
      noise_eigenvectors = new_noise_V
    ))
  } else {
    # Similar logic for n <= m case
    Y <- wishart_matrix(X)
    Yr <- wishart_matrix(Xr)
    
    eigen_Y <- get_eigen(Y)
    eigen_Yr <- get_eigen(Yr)
    
    L <- eigen_Y$values
    V <- eigen_Y$vectors
    Lr <- eigen_Yr$values
    
    mp_result <- mp_calculation(L, Lr[1:(length(Lr)-1)])
    L_mp <- mp_result$L_filtered
    b_min <- mp_result$b_minus
    
    tw_result <- tw_calculation(L, L_mp)
    lambda_c <- tw_result$lambda_c
    
    cat("Number of signal eigenvalues:", sum(L > lambda_c), "\n")
    
    sel_L <- L[L > lambda_c]
    sel_Vs <- V[, L > lambda_c, drop = FALSE]
    
    noise_idx <- L >= b_min & L <= lambda_c
    noise_L <- L[noise_idx]
    noise_V <- V[, noise_idx, drop = FALSE]
    
    if (length(sel_L) > 0) {
      signal_sort_idx <- order(sel_L, decreasing = TRUE)
      nL <- sel_L[signal_sort_idx]
      nVs <- sel_Vs[, signal_sort_idx, drop = FALSE]
    } else {
      nL <- numeric(0)
      nVs <- matrix(nrow = n, ncol = 0)
    }
    
    return(list(
      signal_eigenvalues = nL,
      signal_eigenvectors = nVs,
      all_eigenvalues = L,
      mp_eigenvalues = L_mp,
      lambda_c = lambda_c,
      noise_eigenvalues = noise_L[order(noise_L, decreasing = TRUE)],
      noise_eigenvectors = noise_V[, order(noise_L, decreasing = TRUE), drop = FALSE]
    ))
  }
}

#' Single-cell Linear Embedding of Neighborhoods and Signals (scLENS)
#' 
#' @description
#' scLENS is a dimensionality reduction method for single-cell RNA-seq data that uses
#' Random Matrix Theory (RMT) to distinguish signal from noise and performs robustness
#' testing to identify stable signals. Converted from Julia implementation.
#' 
#' @param seurat_obj Seurat object containing single-cell RNA-seq data
#' @param assay Which assay to use (default: "RNA")
#' @param slot Which slot to use (default: "counts")
#' @param th Threshold angle in degrees for signal robustness test (default: 60)
#' @param p_step Decrement level for sparsity in robustness test (default: 0.001)
#' @param n_perturb Number of perturbations for robustness test (default: 20)
#' @param centering Centering method ("mean" or "median", default: "mean")
#' @param is_normalized Whether the data is already normalized (default: FALSE). If TRUE, skips scLENS normalization steps (log1p transformation, L1 normalization, z-score with L2 normalization). Use TRUE when working with pre-normalized data from "data" or "scale.data" slots
#' @param reduction_name_all Name for all PCs reduction (default: "sclens_pca_all")
#' @param reduction_name_filtered Name for filtered PCs reduction (default: "sclens_pca_filtered")
#' @param n_threads Number of threads to use for parallel processing (default: 1)
#' @param verbose Whether to print detailed progress messages (default: TRUE)
#' @return Modified Seurat object with scLENS reductions added
#' 
#' @details
#' scLENS performs the following steps:
#' \enumerate{
#'   \item Data preprocessing with log1p transformation and normalization
#'   \item Random Matrix Theory (RMT) filtering to identify signal eigenvalues
#'   \item Robustness testing through multiple perturbations
#'   \item Construction of dimensionally-reduced representations
#' }
#' 
#' The function adds two reductions to the Seurat object:
#' \itemize{
#'   \item All signals after RMT filtering
#'   \item Only robust signals after perturbation testing
#' }
#' 
#' @examples
#' \dontrun{
#' # Load Seurat object
#' data(pbmc_small)
#' 
#' # Basic scLENS analysis with raw counts
#' pbmc_small <- sclens(pbmc_small, n_threads = 2)
#' 
#' # Using pre-normalized data
#' pbmc_normalized <- NormalizeData(pbmc_small)
#' pbmc_normalized <- sclens(pbmc_normalized, 
#'                          slot = "data", 
#'                          is_normalized = TRUE)
#' 
#' # Using scaled data
#' pbmc_scaled <- pbmc_small %>% NormalizeData() %>% ScaleData()
#' pbmc_scaled <- sclens(pbmc_scaled, 
#'                      slot = "scale.data", 
#'                      is_normalized = TRUE)
#' 
#' # Check reductions
#' Reductions(pbmc_small)
#' 
#' # Plot results
#' DimPlot(pbmc_small, reduction = "sclens_pca_filtered")
#' }
#' 
#' @export
sclens <- function(seurat_obj, 
                   assay = "RNA", 
                   slot = "counts",
                   th = 60, 
                   p_step = 0.001, 
                   n_perturb = 20, 
                   centering = "mean",
                   is_normalized = FALSE,
                   reduction_name_all = "sclens_pca_all",
                   reduction_name_filtered = "sclens_pca_filtered",
                   n_threads = 5,
                   verbose = TRUE) {
  
  # Validate inputs
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (verbose) {
    start_time <- Sys.time()
    message(paste(rep("=", 80), collapse = ""))
    message("Starting scLENS analysis...")
    message(paste("Timestamp:", format(start_time, "%Y-%m-%d %H:%M:%S")))
    message(paste("R Session Info: R", R.version.string))
    message(paste("Platform:", R.version$platform))
    message(paste("OS Type:", .Platform$OS.type))
    message(paste("Process ID:", Sys.getpid()))
    message(paste(rep("-", 80), collapse = ""))
    message("INPUT PARAMETERS:")
    message(paste("  Seurat object:", class(seurat_obj)[1]))
    message(paste("  Number of cells:", ncol(seurat_obj)))
    message(paste("  Number of features:", nrow(seurat_obj)))
    message(paste("  Using assay:", assay))
    message(paste("  Using slot:", slot))
    message(paste("  Threshold angle:", th, "degrees"))
    message(paste("  Perturbation step:", p_step))
    message(paste("  Number of perturbations:", n_perturb))
    message(paste("  Centering method:", centering))
    message(paste("  Data already normalized:", ifelse(is_normalized, "YES", "NO")))
    message(paste("  Output reduction (all):", reduction_name_all))
    message(paste("  Output reduction (filtered):", reduction_name_filtered))
    message(paste("  Number of threads:", n_threads))
    message(paste(rep("-", 80), collapse = ""))
  }
  
  # Set up parallel processing
  if (verbose) {
    message("PARALLEL PROCESSING SETUP:")
    message(paste("  Detected cores:", parallel::detectCores()))
    message(paste("  Requested threads:", n_threads))
    message(paste("  Already registered:", foreach::getDoParRegistered()))
  }
  
  if (n_threads > 1) {
    if (!foreach::getDoParRegistered()) {
      actual_threads <- min(n_threads, parallel::detectCores() - 1)
      if (verbose) {
        message(paste("  Creating cluster with", actual_threads, "workers"))
        message(paste("  Thread setup: Creating parallel cluster"))
        message(paste("  Backend: doParallel with", actual_threads, "cores"))
      }
      cl <- parallel::makeCluster(actual_threads)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl))
      if (verbose) {
        message("  Parallel cluster created successfully")
      }
    } else {
      if (verbose) {
        message("  Using existing parallel backend")
      }
    }
  } else {
    if (verbose) {
      message("  Running in sequential mode (n_threads = 1)")
    }
  }
  
  # Extract sparse matrix from Seurat object
  # Julia equivalent: X_ = df2sparr(inp_df)
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("DATA EXTRACTION:")
    extraction_start <- Sys.time()
    message(paste("  Starting data extraction at:", format(extraction_start, "%H:%M:%S")))
  }
  
  X_ <- seurat_to_sparse(seurat_obj, assay = assay, slot = slot)
  
  if (verbose) {
    extraction_end <- Sys.time()
    extraction_time <- as.numeric(difftime(extraction_end, extraction_start, units = "secs"))
    message(paste("  Data extraction completed in:", round(extraction_time, 3), "seconds"))
    message(paste("  Matrix dimensions:", nrow(X_), "x", ncol(X_), "(cells x genes)"))
    message(paste("  Matrix class:", paste(class(X_), collapse = ", ")))
    if (inherits(X_, "sparseMatrix")) {
      n_nonzero <- sum(X_ > 0)
      message(paste("  Non-zero entries:", n_nonzero))
      message(paste("  Sparsity:", round((1 - n_nonzero / (nrow(X_) * ncol(X_))) * 100, 2), "%"))
      if (n_nonzero > 0) {
        values <- as.vector(X_[X_ > 0])
        message(paste("  Value range: [", round(min(values), 4), ", ", round(max(values), 4), "]", sep = ""))
        message(paste("  Mean value:", round(mean(values), 4)))
      }
    }
  }
  
  # Get matrix dimensions
  N <- nrow(X_)  # Number of cells
  M <- ncol(X_)  # Number of genes
  
  # Define preprocessing functions
  # Julia equivalent: pre_scale = x -> log1p.(proj_l(x))
  pre_scale <- function(x) {
    log1p(proj_l(x))
  }
  
  # Define log normalization and scaling function based on centering method
  # Julia equivalent: logn_scale = if centering == "mean" ...
  if (centering == "mean") {
    logn_scale <- function(x) {
      # Apply z-score with L2 normalization, then center
      scaled_gdata(zscore_with_l2(x), position_ = "mean") - 
        matrix(colMeans(scaled_gdata(zscore_with_l2(x), position_ = "mean")), 
               nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    }
  } else if (centering == "median") {
    logn_scale <- function(x) {
      scaled_gdata(as.matrix(x), position_ = "median")
    }
  } else {
    warning("Specified centering method not supported. Using mean centering.")
    logn_scale <- function(x) {
      scaled_gdata(zscore_with_l2(x), position_ = "mean") - 
        matrix(colMeans(scaled_gdata(zscore_with_l2(x), position_ = "mean")), 
               nrow = nrow(x), ncol = ncol(x), byrow = TRUE)
    }
  }
  
  # Apply preprocessing
  # Julia equivalent: scaled_X = logn_scale(pre_scale(X_))
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("DATA PREPROCESSING:")
    preprocessing_start <- Sys.time()
    message(paste("  Starting preprocessing at:", format(preprocessing_start, "%H:%M:%S")))
    if (is_normalized) {
      message(paste("  Data normalization: SKIPPED (is_normalized = TRUE)"))
      message(paste("  Using data as-is from slot:", slot))
    } else {
      message(paste("  Data normalization: APPLYING (is_normalized = FALSE)"))
      message(paste("  Centering method:", centering))
    }
  }
  
  if (is_normalized) {
    # Use pre-normalized data but ensure it's centered for RMT analysis
    # This preserves existing normalization while making data suitable for covariance matrix calculations
    scaled_X <- scale(as.matrix(X_), center = TRUE, scale = FALSE)
    if (verbose) {
      message(paste("  Using pre-normalized data with centering applied"))
      message(paste("  Data was centered but not scaled to preserve existing normalization"))
    }
  } else {
    # Apply scLENS normalization pipeline
    scaled_X <- logn_scale(pre_scale(X_))
  }
  
  if (verbose) {
    preprocessing_end <- Sys.time()
    preprocessing_time <- as.numeric(difftime(preprocessing_end, preprocessing_start, units = "secs"))
    message(paste("  Preprocessing completed in:", round(preprocessing_time, 3), "seconds"))
    message(paste("  Matrix dimensions:", nrow(scaled_X), "x", ncol(scaled_X)))
    message(paste("  Matrix class:", paste(class(scaled_X), collapse = ", ")))
    if (is.matrix(scaled_X)) {
      message(paste("  Value range: [", round(min(scaled_X), 4), ", ", round(max(scaled_X), 4), "]", sep = ""))
      message(paste("  Mean value:", round(mean(scaled_X), 4)))
      message(paste("  Standard deviation:", round(sd(as.vector(scaled_X)), 4)))
    }
  }
  
  # Generate randomized matrix for null model
  # Julia equivalent: X_r = df2sparr(random_nz(inp_df,rmix=true))
  X_r <- random_nz(X_, rmix = TRUE)
  
  # Extract signals using Random Matrix Theory
  # Julia equivalent: nL, nV, L, L_mp, lambda_c, _, noiseV = get_sigev(scaled_X,logn_scale(pre_scale(X_r)),device=device_)
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("RANDOM MATRIX THEORY ANALYSIS:")
    rmt_start <- Sys.time()
    message(paste("  Starting RMT analysis at:", format(rmt_start, "%H:%M:%S")))
    message(paste("  Generating randomized null model..."))
  }
  
  # Apply same normalization logic to randomized matrix
  if (is_normalized) {
    scaled_X_r <- scale(as.matrix(X_r), center = TRUE, scale = FALSE)
  } else {
    scaled_X_r <- logn_scale(pre_scale(X_r))
  }
  
  rmt_result <- get_sigev(scaled_X, scaled_X_r)
  
  if (verbose) {
    rmt_end <- Sys.time()
    rmt_time <- as.numeric(difftime(rmt_end, rmt_start, units = "secs"))
    message(paste("  RMT analysis completed in:", round(rmt_time, 3), "seconds"))
  }
  
  nL <- rmt_result$signal_eigenvalues
  nV <- rmt_result$signal_eigenvectors
  L <- rmt_result$all_eigenvalues
  L_mp <- rmt_result$mp_eigenvalues
  lambda_c <- rmt_result$lambda_c
  noiseV <- rmt_result$noise_eigenvectors
  
  # Check if any signals were detected
  min_s <- ncol(nV)
  if (min_s == 0) {
    warning("No signals detected!")
    # Return Seurat object with empty reductions
    seurat_obj[[reduction_name_all]] <- CreateDimReducObject(
      embeddings = matrix(nrow = ncol(seurat_obj), ncol = 0),
      key = paste0(toupper(substr(reduction_name_all, 1, 1)), 
                   substr(reduction_name_all, 2, nchar(reduction_name_all)), "_"),
      assay = assay
    )
    seurat_obj[[reduction_name_filtered]] <- CreateDimReducObject(
      embeddings = matrix(nrow = ncol(seurat_obj), ncol = 0),
      key = paste0(toupper(substr(reduction_name_filtered, 1, 1)), 
                   substr(reduction_name_filtered, 2, nchar(reduction_name_filtered)), "_"),
      assay = assay
    )
    return(seurat_obj)
  }
  
  # Perform robustness testing
  # This section corresponds to the perturbation analysis in the Julia code
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("SIGNAL ROBUSTNESS TESTING:")
    robustness_start <- Sys.time()
    message(paste("  Starting robustness testing at:", format(robustness_start, "%H:%M:%S")))
    message(paste("  Number of perturbations:", n_perturb))
    message(paste("  Perturbation step size:", p_step))
    message(paste("  Robustness threshold angle:", th, "degrees"))
  }
  
  # Calculate noise baseline
  # Julia equivalent: nm = min(N,M); model_norm = Normal(0,sqrt(1/nm))
  nm <- min(N, M)
  p_tharr <- replicate(5000, max(abs(rnorm(nm, mean = 0, sd = sqrt(1/nm)))))
  p_th <- mean(p_tharr)
  if (verbose) {
    message(paste("  Calculated noise threshold:", round(p_th, 6)))
    message(paste("  Matrix dimensions for baseline: min(", N, ",", M, ") =", nm))
  }
  
  # Find optimal sparsity level for perturbation
  # This is a simplified version of the Julia sparsity calculation
  p_ <- 0.999
  if (verbose) cat("Calculating sparsity level for perturbation...\n")
  
  # Create zero indices for perturbation
  # Julia equivalent: z_idx1,z_idx2 = begin...
  X_summary <- summary(X_)
  nz_indices <- cbind(X_summary$i, X_summary$j)
  
  # Generate potential zero indices
  total_indices <- expand.grid(1:N, 1:M)
  zero_indices <- setdiff(
    paste(total_indices[,1], total_indices[,2]), 
    paste(nz_indices[,1], nz_indices[,2])
  )
  
  # Parse back to row,col format
  zero_coords <- do.call(rbind, strsplit(zero_indices, " "))
  z_idx1 <- as.numeric(zero_coords[,1])
  z_idx2 <- as.numeric(zero_coords[,2])
  
  # Simplified sparsity selection (using fixed value for efficiency)
  selected_sparsity <- 0.995
  if (verbose) cat("Selected perturbation sparsity:", selected_sparsity, "\n")
  
  # Perform perturbations
  # Julia equivalent: @showprogress "perturbing..." for _ in 1:n_perturb
  if (verbose) {
    message(paste("  Performing", n_perturb, "perturbations..."))
    message(paste("  Selected sparsity level:", selected_sparsity))
    message(paste("  Available zero indices:", length(z_idx1)))
    message(paste("  Perturbation indices per trial:", min(length(z_idx1), round((1 - selected_sparsity) * M * N))))
    perturb_start <- Sys.time()
    message(paste("  Starting perturbations at:", format(perturb_start, "%H:%M:%S")))
  }
  
  min_pc <- min(ceiling(min_s * 1.5), 50)  # Limit for computational efficiency
  
  if (n_threads > 1 && n_perturb > 1) {
    # Parallel execution
    if (verbose) {
      message(paste("  Running perturbations in parallel with", n_threads, "threads"))
    }
    
    perturbation_results <- foreach(i = 1:n_perturb, .combine = list, .multicombine = TRUE,
                                             .export = c("get_eigvec", "logn_scale", "pre_scale", "proj_l", 
                                                        "zscore_with_l2", "scaled_gdata", "wishart_matrix", 
                                                        "get_eigen", "mp_parameters", "mp_calculation", 
                                                        "tw_calculation", "get_sigev", "random_nz",
                                                        "X_summary", "z_idx1", "z_idx2", "selected_sparsity", 
                                                        "M", "N", "min_pc", "verbose", "n_perturb"),
                                             .packages = c("Matrix", "Seurat")) %dopar% {
      if (verbose && i %% max(1, floor(n_perturb/10)) == 0) {
        cat("Perturbation", i, "of", n_perturb, "(Thread", Sys.getpid(), ")\n")
      }
      
      # Sample indices for perturbation
      n_perturb_indices <- min(length(z_idx1), round((1 - selected_sparsity) * M * N))
      if (n_perturb_indices > 0 && length(z_idx1) > 0) {
        sple_idx <- sample(length(z_idx1), n_perturb_indices, replace = FALSE)
        
        # Create perturbed matrix
        perturb_i <- c(X_summary$i, z_idx1[sple_idx])
        perturb_j <- c(X_summary$j, z_idx2[sple_idx])
        perturb_x <- c(X_summary$x, rep(1, length(sple_idx)))
        
        tmp_X <- Matrix::sparseMatrix(i = perturb_i, j = perturb_j, x = perturb_x, dims = c(N, M))
        
        # Get eigenvectors from perturbed matrix (apply same normalization logic)
        if (is_normalized) {
          tmp_result <- get_eigvec(scale(as.matrix(tmp_X), center = TRUE, scale = FALSE))
        } else {
          tmp_result <- get_eigvec(logn_scale(pre_scale(tmp_X)))
        }
        tmp_nV <- tmp_result$eigenvectors
        tmp_nL <- tmp_result$eigenvalues
        
        # Store results (limit to min_pc components)
        max_components <- min(min_pc, ncol(tmp_nV))
        return(list(
          eigenvectors = tmp_nV[, 1:max_components, drop = FALSE],
          eigenvalues = tmp_nL[1:max_components]
        ))
      } else {
        return(NULL)
      }
    }
    
    # Extract results
    nV_set <- lapply(perturbation_results, function(x) if(!is.null(x)) x$eigenvectors else NULL)
    nL_set <- lapply(perturbation_results, function(x) if(!is.null(x)) x$eigenvalues else NULL)
    
  } else {
    # Sequential execution
    if (verbose) {
      message(paste("  Running perturbations sequentially"))
    }
    
    nV_set <- list()
    nL_set <- list()
    
    for (i in 1:n_perturb) {
      if (verbose) {
        cat("Perturbation", i, "of", n_perturb, "\n")
      }
      
      # Sample indices for perturbation
      n_perturb_indices <- min(length(z_idx1), round((1 - selected_sparsity) * M * N))
      if (n_perturb_indices > 0 && length(z_idx1) > 0) {
        sple_idx <- sample(length(z_idx1), n_perturb_indices, replace = FALSE)
        
        # Create perturbed matrix
        perturb_i <- c(X_summary$i, z_idx1[sple_idx])
        perturb_j <- c(X_summary$j, z_idx2[sple_idx])
        perturb_x <- c(X_summary$x, rep(1, length(sple_idx)))
        
        tmp_X <- Matrix::sparseMatrix(i = perturb_i, j = perturb_j, x = perturb_x, dims = c(N, M))
        
        # Get eigenvectors from perturbed matrix (apply same normalization logic)
        if (is_normalized) {
          tmp_result <- get_eigvec(scale(as.matrix(tmp_X), center = TRUE, scale = FALSE))
        } else {
          tmp_result <- get_eigvec(logn_scale(pre_scale(tmp_X)))
        }
        tmp_nV <- tmp_result$eigenvectors
        tmp_nL <- tmp_result$eigenvalues
        
        # Store results (limit to min_pc components)
        max_components <- min(min_pc, ncol(tmp_nV))
        nV_set[[i]] <- tmp_nV[, 1:max_components, drop = FALSE]
        nL_set[[i]] <- tmp_nL[1:max_components]
      }
    }
  }
  
  if (verbose) {
    perturb_end <- Sys.time()
    perturb_time <- as.numeric(difftime(perturb_end, perturb_start, units = "secs"))
    message(paste("  Perturbations completed in:", round(perturb_time, 3), "seconds"))
    message(paste("  Average time per perturbation:", round(perturb_time / n_perturb, 3), "seconds"))
  }
  
  # Calculate robustness scores
  # Julia equivalent: th_ = cos(deg2rad(th))
  th_ <- cos(th * pi / 180)
  if (verbose) {
    message(paste("  Calculating robustness scores..."))
    message(paste("  Cosine threshold:", round(th_, 4)))
    robustness_calc_start <- Sys.time()
    message(paste("  Starting robustness calculation at:", format(robustness_calc_start, "%H:%M:%S")))
  }
  
  # This is a simplified version of the robustness calculation
  # Original Julia code is quite complex with correlation analysis
  if (length(nV_set) > 1) {
    # Calculate pairwise correlations between perturbed eigenvectors
    correlations <- array(0, dim = c(min_s, length(nV_set), length(nV_set)))
    
    for (i in 1:(length(nV_set)-1)) {
      for (j in (i+1):length(nV_set)) {
        if (!is.null(nV_set[[i]]) && !is.null(nV_set[[j]])) {
          # Calculate correlation between signal eigenvectors
          for (k in 1:min(min_s, ncol(nV_set[[i]]), ncol(nV_set[[j]]))) {
            if (k <= ncol(nV) && k <= ncol(nV_set[[i]]) && k <= ncol(nV_set[[j]])) {
              corr_orig_i <- abs(cor(nV[,k], nV_set[[i]][,k]))
              corr_orig_j <- abs(cor(nV[,k], nV_set[[j]][,k]))
              corr_ij <- abs(cor(nV_set[[i]][,k], nV_set[[j]][,k]))
              correlations[k, i, j] <- min(corr_orig_i, corr_orig_j, corr_ij)
            }
          }
        }
      }
    }
    
    # Calculate robustness scores
    rob_scores <- numeric(min_s)
    for (k in 1:min_s) {
      valid_corrs <- correlations[k, , ][correlations[k, , ] > 0]
      if (length(valid_corrs) > 0) {
        rob_scores[k] <- median(valid_corrs, na.rm = TRUE)
      } else {
        rob_scores[k] <- 0
      }
    }
  } else {
    # If not enough perturbations, use all signals
    rob_scores <- rep(th_ + 0.1, min_s)
  }
  
  # Filter signals based on robustness threshold
  # Julia equivalent: sig_id = findall(rob_score .> th_)
  sig_id <- which(rob_scores > th_)
  
  if (verbose) {
    robustness_calc_end <- Sys.time()
    robustness_calc_time <- as.numeric(difftime(robustness_calc_end, robustness_calc_start, units = "secs"))
    message(paste("  Robustness calculation completed in:", round(robustness_calc_time, 3), "seconds"))
    message(paste("  Number of robust signals found:", length(sig_id)))
    if (length(sig_id) > 0) {
      message(paste("  Robust signal indices:", paste(sig_id, collapse = ", ")))
      message(paste("  Robust signal robustness scores:", paste(round(rob_scores[sig_id], 4), collapse = ", ")))
    }
  }
  
  # Reconstruct reduced data
  # Julia equivalent: Xout0 = nV.*(sqrt.(nL))'; Xout1 = nV[:,sig_id].*sqrt.(nL[sig_id])'
  if (verbose) {
    message(paste(rep("-", 80), collapse = ""))
    message("DATA RECONSTRUCTION:")
    reconstruction_start <- Sys.time()
    message(paste("  Starting data reconstruction at:", format(reconstruction_start, "%H:%M:%S")))
    message(paste("  Total signals after RMT:", length(nL)))
    message(paste("  Robust signals after testing:", length(sig_id)))
  }
  
  # All signals (after RMT filtering)
  if (length(nL) > 0) {
    Xout0 <- nV %*% diag(sqrt(nL))
  } else {
    Xout0 <- matrix(nrow = N, ncol = 0)
  }
  
  # Robust signals only
  if (length(sig_id) > 0) {
    Xout1 <- nV[, sig_id, drop = FALSE] %*% diag(sqrt(nL[sig_id]), nrow = length(sig_id))
  } else {
    Xout1 <- matrix(nrow = N, ncol = 0)
  }
  
  # Create Seurat DimReduc objects
  # The embeddings should be cells x components
  if (ncol(Xout0) > 0) {
    colnames(Xout0) <- paste0(reduction_name_all, "_", 1:ncol(Xout0))
    rownames(Xout0) <- colnames(seurat_obj)
  }
  
  if (ncol(Xout1) > 0) {
    colnames(Xout1) <- paste0(reduction_name_filtered, "_", 1:ncol(Xout1))
    rownames(Xout1) <- colnames(seurat_obj)
  }
  
  # Add reductions to Seurat object
  seurat_obj[[reduction_name_all]] <- CreateDimReducObject(
    embeddings = Xout0,
    key = paste0(toupper(substr(reduction_name_all, 1, 1)), 
                 substr(reduction_name_all, 2, nchar(reduction_name_all)), "_"),
    assay = assay
  )
  
  seurat_obj[[reduction_name_filtered]] <- CreateDimReducObject(
    embeddings = Xout1,
    key = paste0(toupper(substr(reduction_name_filtered, 1, 1)), 
                 substr(reduction_name_filtered, 2, nchar(reduction_name_filtered)), "_"),
    assay = assay
  )
  
  # Store additional metadata
  seurat_obj@misc$sclens_results <- list(
    signal_eigenvalues = nL,
    all_eigenvalues = L,
    mp_eigenvalues = L_mp,
    lambda_c = lambda_c,
    robustness_scores = rob_scores,
    robust_signal_indices = sig_id,
    n_signals_total = length(nL),
    n_signals_robust = length(sig_id)
  )
  
  if (verbose) {
    reconstruction_end <- Sys.time()
    reconstruction_time <- as.numeric(difftime(reconstruction_end, reconstruction_start, units = "secs"))
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    message(paste("  Data reconstruction completed in:", round(reconstruction_time, 3), "seconds"))
    message(paste("  Final embedding dimensions:"))
    message(paste("    - All signals:", nrow(Xout0), "cells x", ncol(Xout0), "components"))
    message(paste("    - Robust signals:", nrow(Xout1), "cells x", ncol(Xout1), "components"))
    
    message(paste(rep("=", 80), collapse = ""))
    message("ANALYSIS COMPLETE")
    message(paste("  Total execution time:", round(total_time, 3), "seconds"))
    message(paste("  Total signals detected:", length(nL)))
    message(paste("  Robust signals identified:", length(sig_id)))
    message(paste("  Signal retention rate:", round(length(sig_id) / max(1, length(nL)) * 100, 1), "%"))
    message(paste("  Reductions added to Seurat object:"))
    message(paste("    -", reduction_name_all, "(", ncol(Xout0), "components )"))
    message(paste("    -", reduction_name_filtered, "(", ncol(Xout1), "components )"))
    
    # Timing breakdown
    message(paste("  Timing breakdown:"))
    message(paste("    - Data extraction:", round(extraction_time, 3), "seconds"))
    message(paste("    - Preprocessing:", round(preprocessing_time, 3), "seconds"))
    message(paste("    - RMT analysis:", round(rmt_time, 3), "seconds"))
    message(paste("    - Robustness testing:", round(perturb_time + robustness_calc_time, 3), "seconds"))
    message(paste("    - Data reconstruction:", round(reconstruction_time, 3), "seconds"))
    
    message(paste("  Analysis completed successfully at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
    message(paste(rep("=", 80), collapse = ""))
  }
  
  return(seurat_obj)
}
