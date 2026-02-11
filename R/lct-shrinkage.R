#' Linear Combination Test with Shrinkage Covariance Estimator
#'
#' @description
#' Performs the Linear Combination Test using the Schäfer-Strimmer shrinkage
#' covariance estimator. This is the default and fastest method, recommended for
#' balanced study designs. The shrinkage estimator provides a well-conditioned
#' covariance matrix even when the number of genes exceeds the number of samples.
#'
#' @param DATA Expression data matrix with genes as rows and samples as columns.
#'   Row names should be gene identifiers. Matrix should be numeric.
#' @param cl Phenotype vector defining group membership. Must be binary (two
#'   levels). Can be numeric (0/1), factor, or character. Length must equal
#'   number of columns in DATA.
#' @param nbPermutations Number of permutations for p-value calculation.
#'   Default is 1000. Minimum recommended: 100. For more precise p-values near
#'   significance threshold, use ≥10000. Rule of thumb: use ≥100/α where α is
#'   your significance level.
#' @param s0 Small positive constant added to eigenvalues for numerical
#'   stability. Default is 0. Use small values (0.01-0.1) if eigenvalues are
#'   very small or negative due to numerical precision issues.
#'
#' @return Numeric p-value between 0 and 1 representing the proportion of
#'   permutations where the test statistic was as extreme or more extreme than
#'   the observed value. Small p-values (e.g., < 0.05) suggest the gene set is
#'   differentially expressed between groups.
#'
#' @details
#' The Linear Combination Test tests whether any linear combination of genes in
#' the set is associated with the phenotype. The algorithm:
#' \enumerate{
#'   \item Estimates the covariance matrix using shrinkage (corpcor::cov.shrink)
#'   \item Performs eigenvalue decomposition
#'   \item Transforms data to uncorrelated features
#'   \item Calculates T²-like test statistic (sum of squared mean differences)
#'   \item Permutes phenotype labels to generate null distribution
#'   \item Computes p-value as proportion of permuted statistics ≥ observed
#' }
#'
#' **Performance**: ~0.5-1 second per gene set (with 1000 permutations)
#'
#' **When to use**: Default choice for most analyses, especially balanced designs.
#'
#' @note
#' Set \code{set.seed()} before calling this function to ensure reproducible
#' p-values across runs.
#'
#' @seealso
#' \code{\link{LCT.ridge}} for unbalanced designs,
#' \code{\link{LCT.glasso}} for sparse networks,
#' \code{\link{LCT.adaptive.lasso}} for feature selection,
#' \code{\link{perform_LCT}} for analyzing multiple gene sets
#'
#' @references
#' Schäfer, J., & Strimmer, K. (2005). A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.
#' \emph{Statistical Applications in Genetics and Molecular Biology}, 4(1).
#' \doi{10.2202/1544-6115.1175}
#'
#' Wang, X., Dinu, I., Liu, W., & Yasui, Y. (2011). Linear combination test
#' for hierarchical gene set analysis. \emph{Statistical Applications in
#' Genetics and Molecular Biology}, 10(1). \doi{10.2202/1544-6115.1618}
#'
#' @export
#' @examples
#' # Example with simulated data
#' set.seed(123)
#' n_genes <- 50
#' n_samples <- 20
#'
#' # Create expression data (genes x samples)
#' expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
#' rownames(expression) <- paste0("Gene", 1:n_genes)
#'
#' # Add differential expression to first 10 genes in group 2
#' phenotype <- c(rep(0, 10), rep(1, 10))
#' expression[1:10, phenotype == 1] <- expression[1:10, phenotype == 1] + 2
#'
#' # Run LCT with shrinkage estimator
#' p_value <- LCT.shrinkage(expression, phenotype, nbPermutations = 1000)
#' print(p_value)  # Should be small (< 0.05) due to differential expression
#'
#' \donttest{
#' # Example with real gene set (subset of genes)
#' pathway_genes <- 1:20  # First 20 genes
#' subset_expression <- expression[pathway_genes, ]
#' p_pathway <- LCT.shrinkage(subset_expression, phenotype, nbPermutations = 5000)
#' }
LCT.shrinkage <- function(DATA, cl, nbPermutations = 1000, s0 = 0) {
  GS.size <- dim(DATA)[1]
  nb.Samples <- ncol(DATA)

  if (GS.size >= 2) {
    # Transpose and compute shrinkage covariance
    DATA <- t(DATA)
    Cov.Pooled <- corpcor::cov.shrink(DATA, verbose = FALSE, lambda.var = 0)

    # Eigenvalue decomposition
    EIGEN.decom <- eigen(Cov.Pooled)
    D <- EIGEN.decom$values + s0
    U <- EIGEN.decom$vectors

    # Transform data
    DATA <- t(DATA %*% U) / sqrt(D)
  }

  if (GS.size == 1) {
    cov.pooled <- c(stats::cov(t(DATA)))
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
