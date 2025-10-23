# ==============================================================================
# Linear Combination Test (LCT) for Gene Set Analysis
# ==============================================================================
# Description: Implementation of LCT with multiple covariance matrix estimators
#              (Shrinkage, Ridge, Graphical Lasso, Adaptive Lasso) for 
#              gene set analysis in high-dimensional data
#
# Authors: Dr. Sara Khademioureh
# Date: 23 October 2025
# ==============================================================================

# Required Libraries -----------------------------------------------------------
library(MASS)
library(corpcor)
library(qvalue)
library(dplyr, warn.conflicts = FALSE)
library(faux)
library(glasso)
library(cvCovEst)
library(rags2ridges)

# Helper Functions -------------------------------------------------------------

#' Convert Gene Set DataFrame to List Format
#' 
#' @param GS Gene set as dataframe (genes x sets) or list
#' @return List of gene sets with gene identifiers
GS.format.dataframe.to.list <- function(GS) {
  if (is.data.frame(GS)) {
    genes <- rownames(GS)
    L <- NULL
    for (ags in names(GS)) {
      w <- which(GS[, ags] == 1)
      if (length(w) > 0) {
        L <- c(L, list(genes[w]))
        names(L)[length(L)] <- ags
      }
    }
    return(L)
  } else {
    return(GS)
  }
}

#' T2-like Test Statistic for SAM-GS
#' 
#' @param DATA Expression data matrix (genes x samples)
#' @param cl Factor defining phenotype groups
#' @return Sum of squared mean differences
T2.like.SAMGS <- function(DATA, cl) {
  if (!is.null(cl)) {
    cl <- as.factor(cl)
    C1 <- which(cl == levels(cl)[1])
    C2 <- which(cl == levels(cl)[2])
  }
  
  if (is.null(C1) | is.null(C2)) {
    stop("Error - SAMGS-PC.TlikeStat: classes 1 and 2 are undefined.")
  }
  
  GS.size <- dim(DATA)[1]
  
  if (GS.size >= 2) {
    means.C1 <- rowMeans(DATA[, C1])
    means.C2 <- rowMeans(DATA[, C2])
    diffmean.C1C2 <- means.C1 - means.C2
  } else {
    means.C1 <- mean(DATA[C1])
    means.C2 <- mean(DATA[C2])
    diffmean.C1C2 <- means.C1 - means.C2
  }
  
  return(sum(diffmean.C1C2^2))
}

# LCT Core Functions -----------------------------------------------------------

#' Linear Combination Test with Shrinkage Covariance Estimator
#' 
#' @param DATA Expression data matrix (genes x samples)
#' @param cl Factor defining phenotype groups
#' @param nbPermutations Number of permutations for p-value estimation
#' @param s0 Small positive constant for stability
#' @return P-value from permutation test
LCT.shrinkage <- function(DATA, cl, nbPermutations = 1000, s0 = 0) {
  GS.size <- dim(DATA)[1]
  nb.Samples <- ncol(DATA)
  
  if (GS.size >= 2) {
    # Transpose and compute shrinkage covariance
    DATA <- t(DATA)
    Cov.Pooled <- cov.shrink(DATA, verbose = FALSE, lambda.var = 0)
    
    # Eigenvalue decomposition
    EIGEN.decom <- eigen(Cov.Pooled)
    D <- EIGEN.decom$values + s0
    U <- EIGEN.decom$vectors
    
    # Transform data
    DATA <- t(DATA %*% U) / sqrt(D)
  }
  
  if (GS.size == 1) {
    cov.pooled <- c(cov(t(DATA)))
    DATA <- DATA / sqrt(cov.pooled + s0)
    DATA <- as.matrix(DATA)
  }
  
  # Observed test statistic
  sam.sumsquareT.obs <- T2.like.SAMGS(DATA, cl)
  
  # Permutation test
  nb.Samples <- length(cl)
  sam.sumsquareT.permut <- rep(NA, nbPermutations)
  
  for (i in 1:nbPermutations) {
    ind <- sample(nb.Samples)
    data <- matrix(DATA[, ind], ncol = nb.Samples)
    sam.sumsquareT.permut[i] <- T2.like.SAMGS(data, cl)
  }
  
  # Calculate p-value
  p.value <- sum(sam.sumsquareT.permut >= sam.sumsquareT.obs) / nbPermutations
  
  return(p.value)
}

#' Linear Combination Test with Graphical Lasso Covariance Estimator
#' 
#' @param DATA Expression data matrix (genes x samples)
#' @param cl Factor defining phenotype groups
#' @param nbPermutations Number of permutations for p-value estimation
#' @param s0 Small positive constant for stability
#' @param rho Penalty parameter for glasso
#' @return P-value from permutation test
LCT.glasso <- function(DATA, cl, nbPermutations = 1000, s0 = 0, rho = 0.4) {
  GS.size <- dim(DATA)[1]
  nb.Samples <- ncol(DATA)
  
  if (GS.size >= 2) {
    DATA <- t(DATA)
    s <- var(DATA)
    
    # Graphical lasso estimation
    a <- glasso(s, rho = rho)
    Cov.Pooled <- a$w
    
    # Eigenvalue decomposition
    EIGEN.decom <- eigen(Cov.Pooled)
    D <- EIGEN.decom$values + s0
    U <- EIGEN.decom$vectors
    
    DATA <- t(DATA %*% U) / sqrt(D)
  }
  
  if (GS.size == 1) {
    cov.pooled <- c(cov(t(DATA)))
    DATA <- DATA / sqrt(cov.pooled + s0)
  }
  
  sam.sumsquareT.obs <- T2.like.SAMGS(DATA, cl)
  
  nb.Samples <- length(cl)
  sam.sumsquareT.permut <- rep(NA, nbPermutations)
  
  for (i in 1:nbPermutations) {
    ind <- sample(nb.Samples)
    data <- matrix(DATA[, ind], ncol = nb.Samples)
    sam.sumsquareT.permut[i] <- T2.like.SAMGS(data, cl)
  }
  
  p.value <- sum(sam.sumsquareT.permut >= sam.sumsquareT.obs) / nbPermutations
  
  return(p.value)
}

#' Linear Combination Test with Adaptive Lasso Covariance Estimator
#' 
#' @param DATA Expression data matrix (genes x samples)
#' @param cl Factor defining phenotype groups
#' @param nbPermutations Number of permutations for p-value estimation
#' @param s0 Small positive constant for stability
#' @param lambda Penalty parameter
#' @param n Adaptive weight exponent
#' @return P-value from permutation test
LCT.adaptive.lasso <- function(DATA, cl, nbPermutations = 1000, s0 = 0.5, 
                                lambda = 0.5, n = 0.8) {
  GS.size <- dim(DATA)[1]
  nb.Samples <- ncol(DATA)
  
  if (GS.size >= 2) {
    DATA <- t(DATA)
    s <- var(DATA)
    
    # Adaptive lasso estimation
    a <- adaptiveLassoEst(dat = s, lambda = lambda, n = n)
    Cov.Pooled <- a
    
    EIGEN.decom <- eigen(Cov.Pooled)
    D <- EIGEN.decom$values + s0
    U <- EIGEN.decom$vectors
    
    DATA <- t(DATA %*% U) / sqrt(D)
  }
  
  if (GS.size == 1) {
    cov.pooled <- c(cov(t(DATA)))
    DATA <- DATA / sqrt(cov.pooled + s0)
  }
  
  sam.sumsquareT.obs <- T2.like.SAMGS(DATA, cl)
  
  nb.Samples <- length(cl)
  sam.sumsquareT.permut <- rep(NA, nbPermutations)
  
  for (i in 1:nbPermutations) {
    ind <- sample(nb.Samples)
    data <- matrix(DATA[, ind], ncol = nb.Samples)
    sam.sumsquareT.permut[i] <- T2.like.SAMGS(data, cl)
  }
  
  p.value <- sum(sam.sumsquareT.permut >= sam.sumsquareT.obs) / nbPermutations
  
  return(p.value)
}

#' Linear Combination Test with Ridge Covariance Estimator
#' 
#' @param DATA Expression data matrix (genes x samples)
#' @param cl Factor defining phenotype groups
#' @param nbPermutations Number of permutations for p-value estimation
#' @param s0 Small positive constant for stability
#' @param lambda Penalty parameter for ridge
#' @return P-value from permutation test
LCT.ridge <- function(DATA, cl, nbPermutations = 1000, s0 = 0, lambda = 1) {
  GS.size <- dim(DATA)[1]
  nb.Samples <- ncol(DATA)
  
  if (GS.size >= 2) {
    DATA <- t(DATA)
    s <- var(DATA)
    
    # Ridge estimation
    Cov.Pooled <- ridgeP(s, lambda = lambda)
    
    EIGEN.decom <- eigen(Cov.Pooled)
    D <- EIGEN.decom$values + s0
    U <- EIGEN.decom$vectors
    
    DATA <- t(DATA %*% U) / sqrt(D)
  }
  
  if (GS.size == 1) {
    cov.pooled <- c(cov(t(DATA)))
    DATA <- DATA / sqrt(cov.pooled + s0)
  }
  
  sam.sumsquareT.obs <- T2.like.SAMGS(DATA, cl)
  
  nb.Samples <- length(cl)
  sam.sumsquareT.permut <- rep(NA, nbPermutations)
  
  for (i in 1:nbPermutations) {
    ind <- sample(nb.Samples)
    data <- matrix(DATA[, ind], ncol = nb.Samples)
    sam.sumsquareT.permut[i] <- T2.like.SAMGS(data, cl)
  }
  
  p.value <- sum(sam.sumsquareT.permut >= sam.sumsquareT.obs) / nbPermutations
  
  return(p.value)
}

# Wrapper Functions for Gene Set Analysis -------------------------------------

#' Perform LCT on Multiple Gene Sets
#' 
#' @param GS Gene sets (dataframe or list)
#' @param DATA Expression data matrix
#' @param cl Phenotype labels
#' @param nbPermutations Number of permutations
#' @param method Covariance estimation method ("shrinkage", "glasso", 
#'               "adaptive_lasso", or "ridge")
#' @return Dataframe with gene set statistics
perform_LCT <- function(GS, DATA, cl, nbPermutations = 1000, 
                        method = "shrinkage") {
  # Preprocess gene sets
  genes <- rownames(DATA)
  nb.Samples <- ncol(DATA)
  nb.GeneSets <- length(GS)
  
  GS <- GS.format.dataframe.to.list(GS)
  GS <- lapply(GS, function(z) as.numeric(which(genes %in% z)))
  GS.sizes <- sapply(GS, length)
  GS.data <- lapply(GS, function(z) as.matrix(DATA[z, ], ncol = nb.Samples))
  
  # Select appropriate LCT function
  LCT_function <- switch(method,
    "shrinkage" = LCT.shrinkage,
    "glasso" = LCT.glasso,
    "adaptive_lasso" = LCT.adaptive.lasso,
    "ridge" = LCT.ridge,
    stop("Unknown method. Choose from: shrinkage, glasso, adaptive_lasso, ridge")
  )
  
  # Calculate p-values
  GeneSets.pval <- sapply(GS.data, function(z) {
    LCT_function(z, cl, nbPermutations = nbPermutations)
  })
  
  # Calculate q-values if multiple gene sets
  if (nb.GeneSets >= 2) {
    GeneSets.qval <- qvalue(GeneSets.pval, pi0 = 1)$qvalues
    
    res <- data.frame(
      "GS_name" = names(GS),
      "GS_size" = GS.sizes,
      "p_value" = GeneSets.pval,
      "q_value" = GeneSets.qval,
      "method" = method,
      stringsAsFactors = FALSE
    )
  } else {
    res <- data.frame(
      "GS_name" = names(GS),
      "GS_size" = GS.sizes,
      "p_value" = GeneSets.pval,
      "method" = method,
      stringsAsFactors = FALSE
    )
  }
  
  return(res)
}

# Simulation Functions ---------------------------------------------------------

#' Generate Random Gene Names
#' 
#' @param n Number of gene names to generate
#' @return Character vector of gene names
generate_gene_names <- function(n) {
  a <- do.call(paste0, replicate(7, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(999999, n, TRUE)), sample(LETTERS, n, TRUE))
}

#' Simulate Gene Expression Data
#' 
#' @param rho Correlation between genes
#' @param n_samples Total number of samples
#' @param n_genes Number of genes
#' @param gamma Mean difference between groups
#' @param n_controls Number of control samples
#' @param n_cases Number of case samples
#' @return Matrix of simulated gene expression data
simulate_gene_expression <- function(rho, n_samples, n_genes, gamma, 
                                      n_controls, n_cases) {
  # Generate control group
  dat1 <- rnorm_multi(
    n = n_genes,
    mu = rep(0, n_controls),
    sd = rnorm(n_controls, 0, 1),
    r = rho,
    empirical = FALSE
  )
  
  # Generate case group with mean shift
  dat2 <- rnorm_multi(
    n = n_genes,
    mu = rep(0 + gamma, n_cases),
    sd = rnorm(n_cases, 0, 1),
    r = rho,
    empirical = FALSE
  )
  
  # Combine and format
  dat <- cbind(dat1, dat2)
  dat <- abs(dat)
  rownames(dat) <- generate_gene_names(n_genes)
  
  return(dat)
}

#' Run LCT Simulation Study
#' 
#' @param n_samples Sample size
#' @param n_genes Number of genes
#' @param rho Correlation between genes
#' @param n_datasets Number of simulated datasets
#' @param method Covariance estimation method
#' @param gene_sets Predefined gene sets
#' @param phenotype Binary phenotype vector
#' @param gamma Mean difference parameter
#' @return List of results from each simulation
run_LCT_simulation <- function(n_samples, n_genes, rho, n_datasets, method, 
                                gene_sets, phenotype, gamma) {
  results_list <- list()
  n_controls <- sum(phenotype == 0)
  n_cases <- sum(phenotype == 1)
  
  for (b in 1:n_datasets) {
    cat("Simulation iteration:", b, "/", n_datasets, "\n")
    
    tryCatch({
      # Generate data
      G <- simulate_gene_expression(
        rho = rho,
        n_samples = n_samples,
        n_genes = n_genes,
        gamma = gamma,
        n_controls = n_controls,
        n_cases = n_cases
      )
      
      # Run LCT
      LCT_result <- perform_LCT(
        GS = gene_sets,
        DATA = as.data.frame(G),
        cl = phenotype,
        nbPermutations = 1000,
        method = method
      )
      
      results_list[[b]] <- LCT_result
      
    }, error = function(e) {
      cat("Error in iteration", b, ":", conditionMessage(e), "\n")
    })
  }
  
  return(results_list)
}

# Main Simulation Script -------------------------------------------------------

#' Main function to run complete simulation study
#' 
#' @param output_dir Directory to save results
#' @export
run_complete_simulation <- function(output_dir = ".") {
  # Set random seed for reproducibility
  set.seed(999)
  
  # Simulation parameters
  sample_sizes <- c(20, 50)
  correlations <- c(0.1, 0.9)
  methods <- c("shrinkage", "adaptive_lasso", "ridge", "glasso")
  outcome_probs <- c(0.2, 0.8)
  gamma_values <- seq(0, 2, by = 0.1)
  
  n_genes <- 10000
  n_gene_sets <- 1
  n_iterations <- 1000
  
  # Create gene sets
  gene_names <- generate_gene_names(n_genes)
  gene_set_matrix <- matrix(0, n_genes, n_gene_sets)
  
  for (i in 1:n_gene_sets) {
    gene_set_matrix[, i] <- rbinom(n_genes, 1, 0.002)
  }
  
  gene_sets <- data.frame(gene_set_matrix)
  rownames(gene_sets) <- gene_names
  
  # Run simulations
  for (n_samples in sample_sizes) {
    for (prob in outcome_probs) {
      # Generate phenotype
      phenotype <- rbinom(n_samples, 1, prob)
      
      cat("\n===================================\n")
      cat("Sample size:", n_samples, "| Outcome probability:", prob, "\n")
      cat("===================================\n")
      
      for (rho in correlations) {
        cat("Correlation:", rho, "\n")
        
        for (method in methods) {
          cat("Method:", method, "\n")
          
          for (gamma in gamma_values) {
            start_time <- Sys.time()
            
            results <- run_LCT_simulation(
              n_samples = n_samples,
              n_genes = n_genes,
              rho = rho,
              n_datasets = n_iterations,
              method = method,
              gene_sets = gene_sets,
              phenotype = phenotype,
              gamma = gamma
            )
            
            # Save results
            filename <- sprintf(
              "%s/Simulation_n%d_p%.1f_rho%.1f_method_%s_gamma%.1f.csv",
              output_dir, n_samples, prob, rho, method, gamma
            )
            
            results_df <- bind_rows(results, .id = "iteration")
            write.csv(results_df, filename, row.names = FALSE)
            
            elapsed <- difftime(Sys.time(), start_time, units = "mins")
            cat("Completed in", round(elapsed, 2), "minutes\n")
          }
        }
      }
    }
  }
  
  cat("\nSimulation complete!\n")
}

# ==============================================================================
# Example Usage
# ==============================================================================

# Uncomment to run the complete simulation:
# run_complete_simulation(output_dir = "./simulation_results")
