#' Linear Combination Test with Adaptive Lasso Covariance Estimator
#'
#' @description
#' Performs the Linear Combination Test using the adaptive lasso for covariance
#' estimation. Adaptive lasso applies weighted L1 penalties that can achieve
#' oracle properties (consistent variable selection and asymptotically efficient
#' estimation). Useful when you want automatic feature selection along with
#' covariance estimation.
#'
#' @param DATA Expression data matrix with genes as rows and samples as columns.
#'   Row names should be gene identifiers. Matrix should be numeric.
#' @param cl Phenotype vector defining group membership. Must be binary (two
#'   levels). Can be numeric (0/1), factor, or character. Length must equal
#'   number of columns in DATA.
#' @param nbPermutations Number of permutations for p-value calculation.
#'   Default is 1000.
#' @param s0 Small positive constant added to eigenvalues for numerical
#'   stability. Default is 0.5 (higher than other methods due to potential for
#'   smaller eigenvalues).
#' @param lambda Penalty parameter for adaptive lasso. Default is 0.5. Controls
#'   overall regularization strength. Typical range: 0.1-1.0.
#' @param n Adaptive weight exponent. Default is 0.8. Controls how weights are
#'   computed: weights = (1/|initial_estimate|)^n. Higher n = more aggressive
#'   penalization of small coefficients. Typical range: 0.5-1.0.
#'
#' @return Numeric p-value between 0 and 1 from the permutation test.
#'
#' @details
#' The adaptive lasso uses data-driven weights to penalize different covariance
#' elements differently:
#' \deqn{\hat{\Sigma} = \arg\min \{ ||S - \Sigma||_F^2 + \lambda \sum w_j |\sigma_j| \}}
#' where weights \eqn{w_j = 1/|\hat{\sigma}_j|^n} are based on an initial
#' estimate.
#'
#' **Advantages**:
#' \itemize{
#'   \item Oracle properties (consistent selection + efficient estimation)
#'   \item Automatic feature selection
#'   \item Less biased than standard lasso for large coefficients
#' }
#'
#' **Disadvantages**:
#' \itemize{
#'   \item Two tuning parameters (lambda and n)
#'   \item Computationally slower than shrinkage (~1-2 sec per gene set)
#'   \item Requires initial estimator (can be unstable for small samples)
#' }
#'
#' **Performance**: ~1-2 seconds per gene set (with 1000 permutations)
#'
#' **When to use**: Feature selection, oracle properties desired, moderate sample sizes
#'
#' @note
#' The adaptive lasso may be sensitive to the initial covariance estimate,
#' especially with small sample sizes. For very small samples (n < 30),
#' consider LCT.shrinkage or LCT.ridge instead.
#'
#' @seealso
#' \code{\link{LCT.shrinkage}} for faster computation,
#' \code{\link{LCT.glasso}} for sparse networks,
#' \code{\link{perform_LCT}} for multiple gene sets
#'
#' @references
#' Rothman, A. J., Levina, E., & Zhu, J. (2009). Generalized thresholding of
#' large covariance matrices. \emph{Journal of the American Statistical
#' Association}, 104(485), 177-186. \doi{10.1198/jasa.2009.0101}
#'
#' Zou, H. (2006). The adaptive lasso and its oracle properties. \emph{Journal
#' of the American Statistical Association}, 101(476), 1418-1429.
#' \doi{10.1198/016214506000000735}
#'
#' Wang, X., Dinu, I., Liu, W., & Yasui, Y. (2011). Linear combination test
#' for hierarchical gene set analysis. \emph{Statistical Applications in
#' Genetics and Molecular Biology}, 10(1). \doi{10.2202/1544-6115.1618}
#'
#' @export
#' @examples
#' # Example with feature selection
#' set.seed(321)
#' n_genes <- 40
#' n_samples <- 30
#'
#' # Create expression data with sparse effects
#' expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
#' rownames(expression) <- paste0("Gene", 1:n_genes)
#'
#' # Only 5 genes are truly differentially expressed
#' phenotype <- c(rep(0, 15), rep(1, 15))
#' expression[1:5, phenotype == 1] <- expression[1:5, phenotype == 1] + 2
#'
#' # Run LCT with adaptive lasso
#' p_value <- LCT.adaptive.lasso(expression, phenotype, nbPermutations = 1000,
#'                                lambda = 0.5, n = 0.8)
#' print(p_value)
#'
#' \donttest{
#' # Compare different n values (weight exponent)
#' p_n05 <- LCT.adaptive.lasso(expression, phenotype, nbPermutations = 1000,
#'                              lambda = 0.5, n = 0.5)
#' p_n10 <- LCT.adaptive.lasso(expression, phenotype, nbPermutations = 1000,
#'                              lambda = 0.5, n = 1.0)
#'
#' cat("n = 0.5 (less aggressive):", p_n05, "\n")
#' cat("n = 1.0 (more aggressive):", p_n10, "\n")
#' }
LCT.adaptive.lasso <- function(DATA, cl, nbPermutations = 1000, s0 = 0.5,
                               lambda = 0.5, n = 0.8) {
  GS.size <- dim(DATA)[1]
  nb.Samples <- ncol(DATA)

  if (GS.size >= 2) {
    DATA <- t(DATA)
    s <- stats::var(DATA)

    # Adaptive lasso estimation
    a <- cvCovEst::adaptiveLassoEst(dat = s, lambda = lambda, n = n)
    Cov.Pooled <- a

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
