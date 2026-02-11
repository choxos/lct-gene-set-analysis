#' Linear Combination Test with Ridge Covariance Estimator
#'
#' @description
#' Performs the Linear Combination Test using ridge regularization (L2 penalty)
#' for covariance estimation. Recommended for unbalanced study designs where
#' group sizes differ substantially. Ridge regularization provides better Type I
#' error control in unbalanced scenarios compared to shrinkage.
#'
#' @param DATA Expression data matrix with genes as rows and samples as columns.
#'   Row names should be gene identifiers. Matrix should be numeric.
#' @param cl Phenotype vector defining group membership. Must be binary (two
#'   levels). Can be numeric (0/1), factor, or character. Length must equal
#'   number of columns in DATA.
#' @param nbPermutations Number of permutations for p-value calculation.
#'   Default is 1000. See \code{\link{LCT.shrinkage}} for guidelines.
#' @param s0 Small positive constant added to eigenvalues for numerical
#'   stability. Default is 0.
#' @param lambda Ridge penalty parameter. Default is 1. Higher values impose
#'   stronger regularization. Typical range: 0.1-10. For very high-dimensional
#'   data, consider lambda = 5-10.
#'
#' @return Numeric p-value between 0 and 1 from the permutation test.
#'
#' @details
#' Ridge regularization adds a penalty term λ to the diagonal of the covariance
#' matrix, ensuring positive definiteness and numerical stability. This is
#' particularly important for:
#' \itemize{
#'   \item Unbalanced designs (e.g., 20 controls vs 5 cases)
#'   \item High-dimensional data (genes >> samples)
#'   \item Multicollinear gene expression patterns
#' }
#'
#' The ridge estimator is:
#' \deqn{\hat{\Sigma}_{\lambda} = \frac{1}{1+\lambda}S + \frac{\lambda}{1+\lambda}I}
#' where S is the sample covariance matrix and I is the identity matrix.
#'
#' **Performance**: ~0.5-1 second per gene set (with 1000 permutations)
#'
#' **When to use**: Unbalanced designs, highly correlated genes
#'
#' @note
#' The lambda parameter can be tuned based on data characteristics. Start with
#' default (lambda = 1) and increase if numerical issues occur.
#'
#' @seealso
#' \code{\link{LCT.shrinkage}} for balanced designs,
#' \code{\link{LCT.glasso}} for sparse networks,
#' \code{\link{perform_LCT}} for multiple gene sets
#'
#' @references
#' Van Wieringen, W. N., & Peeters, C. F. (2016). Ridge estimation of inverse
#' covariance matrices from high-dimensional data. \emph{Computational
#' Statistics & Data Analysis}, 103, 284-303. \doi{10.1016/j.csda.2016.05.012}
#'
#' Wang, X., Dinu, I., Liu, W., & Yasui, Y. (2011). Linear combination test
#' for hierarchical gene set analysis. \emph{Statistical Applications in
#' Genetics and Molecular Biology}, 10(1). \doi{10.2202/1544-6115.1618}
#'
#' @export
#' @examples
#' # Example with unbalanced design
#' set.seed(456)
#' n_genes <- 50
#'
#' # Create unbalanced groups (15 controls, 5 cases)
#' phenotype <- c(rep(0, 15), rep(1, 5))
#' n_samples <- length(phenotype)
#'
#' # Create expression data
#' expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
#' rownames(expression) <- paste0("Gene", 1:n_genes)
#'
#' # Add differential expression
#' expression[1:10, phenotype == 1] <- expression[1:10, phenotype == 1] + 1.5
#'
#' # Run LCT with ridge estimator (better for unbalanced designs)
#' p_value <- LCT.ridge(expression, phenotype, nbPermutations = 1000, lambda = 1)
#' print(p_value)
#'
#' \donttest{
#' # Compare different lambda values
#' p_lambda1 <- LCT.ridge(expression, phenotype, nbPermutations = 1000, lambda = 1)
#' p_lambda5 <- LCT.ridge(expression, phenotype, nbPermutations = 1000, lambda = 5)
#'
#' # Higher lambda = more regularization
#' cat("Lambda = 1:", p_lambda1, "\n")
#' cat("Lambda = 5:", p_lambda5, "\n")
#' }
LCT.ridge <- function(DATA, cl, nbPermutations = 1000, s0 = 0, lambda = 1) {
  GS.size <- dim(DATA)[1]
  nb.Samples <- ncol(DATA)

  if (GS.size >= 2) {
    DATA <- t(DATA)
    s <- stats::var(DATA)

    # Ridge estimation
    Cov.Pooled <- rags2ridges::ridgeP(s, lambda = lambda)

    EIGEN.decom <- eigen(Cov.Pooled)
    D <- EIGEN.decom$values + s0
    U <- EIGEN.decom$vectors

    DATA <- t(DATA %*% U) / sqrt(D)
  }

  if (GS.size == 1) {
    cov.pooled <- c(stats::cov(t(DATA)))
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
