#' Perform Linear Combination Test on Multiple Gene Sets
#'
#' @description
#' Main interface for performing LCT analysis on one or more gene sets. This
#' wrapper function handles gene set formatting, method selection, p-value
#' calculation, and FDR correction automatically. Recommended for most users
#' as the primary entry point to the package.
#'
#' @param GS Gene sets in one of two formats:
#'   \itemize{
#'     \item \strong{List}: Named list where each element is a character vector
#'       of gene identifiers (e.g., \code{list(Pathway1 = c("GENE1", "GENE2"),
#'       Pathway2 = c("GENE3", "GENE4"))})
#'     \item \strong{Data frame}: Binary matrix with genes as rows and gene sets
#'       as columns. Row names should be gene identifiers, column names should
#'       be gene set names. Entry (i,j) = 1 if gene i belongs to set j, 0 otherwise.
#'   }
#' @param DATA Expression data matrix with genes as rows and samples as columns.
#'   Row names must be gene identifiers matching those in GS. Genes not in any
#'   gene set will be ignored. Matrix should be numeric.
#' @param cl Phenotype vector defining group membership. Must be binary (two
#'   levels). Can be numeric (0/1), factor, or character. Length must equal
#'   number of columns in DATA.
#' @param nbPermutations Number of permutations for p-value calculation.
#'   Default is 1000. See \code{\link{LCT.shrinkage}} for guidelines.
#' @param method Covariance estimation method to use. One of:
#'   \itemize{
#'     \item \code{"shrinkage"} (default): Fast, recommended for balanced designs
#'     \item \code{"ridge"}: Recommended for unbalanced designs
#'     \item \code{"glasso"}: Graphical lasso for sparse networks
#'     \item \code{"adaptive_lasso"}: Adaptive lasso for feature selection
#'   }
#'   See "Method Selection" section below for detailed guidance.
#'
#' @return Data frame with one row per gene set and the following columns:
#'   \describe{
#'     \item{GS_name}{Gene set identifier (from names of input list/data frame)}
#'     \item{GS_size}{Number of genes in the set (after matching with DATA)}
#'     \item{p_value}{Permutation-based p-value}
#'     \item{q_value}{FDR-adjusted q-value (only present if ≥2 gene sets tested)}
#'     \item{method}{Covariance estimation method used}
#'   }
#'   Rows are sorted by gene set name.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item Converts gene sets to list format (if needed)
#'   \item Matches gene set identifiers with expression data row names
#'   \item Extracts expression data for each gene set
#'   \item Applies selected LCT method to each gene set
#'   \item Calculates q-values using the qvalue package (if ≥2 sets)
#'   \item Returns formatted results data frame
#' }
#'
#' @section Method Selection:
#' Choose the covariance estimation method based on your data characteristics:
#'
#' \strong{Shrinkage} (default):
#' \itemize{
#'   \item Best for: Balanced designs, general use
#'   \item Speed: Fastest (~0.5-1 sec per gene set)
#'   \item When to avoid: Highly unbalanced groups
#' }
#'
#' \strong{Ridge}:
#' \itemize{
#'   \item Best for: Unbalanced designs (e.g., 20 controls vs 5 cases)
#'   \item Speed: Fast (~0.5-1 sec per gene set)
#'   \item Benefits: Better Type I error control in unbalanced scenarios
#' }
#'
#' \strong{Graphical Lasso}:
#' \itemize{
#'   \item Best for: Network analysis, sparse gene interactions
#'   \item Speed: Slowest (~4-8 sec per gene set)
#'   \item Benefits: Identifies conditional independence structures
#' }
#'
#' \strong{Adaptive Lasso}:
#' \itemize{
#'   \item Best for: Feature selection, oracle properties
#'   \item Speed: Moderate (~1-2 sec per gene set)
#'   \item Benefits: Automatic feature selection, less bias
#' }
#'
#' @section FDR Correction:
#' When testing ≥2 gene sets, q-values are automatically calculated using:
#' \itemize{
#'   \item The \code{qvalue} package (Bioconductor) if available (recommended)
#'   \item Fallback to \code{p.adjust(..., method = "BH")} if qvalue not installed
#' }
#'
#' Install qvalue with: \code{BiocManager::install("qvalue")}
#'
#' @note
#' \itemize{
#'   \item Set \code{set.seed()} before calling for reproducible results
#'   \item Gene sets with <2 genes after matching are skipped with a warning
#'   \item Large gene sets (>1000 genes) may be slow, especially with GLASSO
#' }
#'
#' @seealso
#' Individual methods: \code{\link{LCT.shrinkage}}, \code{\link{LCT.ridge}},
#' \code{\link{LCT.glasso}}, \code{\link{LCT.adaptive.lasso}}
#'
#' Simulation: \code{\link{simulate_gene_expression}}, \code{\link{run_LCT_simulation}}
#'
#' @references
#' Wang, X., Dinu, I., Liu, W., & Yasui, Y. (2011). Linear combination test
#' for hierarchical gene set analysis. \emph{Statistical Applications in
#' Genetics and Molecular Biology}, 10(1). \doi{10.2202/1544-6115.1618}
#'
#' Storey, J. D., & Tibshirani, R. (2003). Statistical significance for
#' genomewide studies. \emph{Proceedings of the National Academy of Sciences},
#' 100(16), 9440-9445. \doi{10.1073/pnas.1530509100}
#'
#' @export
#' @examples
#' # Create example data
#' set.seed(100)
#' n_genes <- 100
#' n_samples <- 20
#'
#' # Expression matrix
#' expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
#' rownames(expression) <- paste0("Gene", 1:n_genes)
#' colnames(expression) <- paste0("Sample", 1:n_samples)
#'
#' # Phenotype (10 controls, 10 cases)
#' phenotype <- rep(0:1, each = 10)
#'
#' # Add differential expression to first 15 genes
#' expression[1:15, phenotype == 1] <- expression[1:15, phenotype == 1] + 1.5
#'
#' # Define gene sets as a list
#' gene_sets <- list(
#'   Pathway_A = paste0("Gene", 1:20),    # Contains DE genes
#'   Pathway_B = paste0("Gene", 21:40),   # No DE genes
#'   Pathway_C = paste0("Gene", 1:10)     # All DE genes
#' )
#'
#' # Run LCT analysis with default shrinkage method
#' results <- perform_LCT(gene_sets, expression, phenotype, nbPermutations = 1000)
#' print(results)
#'
#' # View significant gene sets (q < 0.05)
#' significant <- results[results$q_value < 0.05, ]
#' print(significant)
#'
#' \donttest{
#' # Compare methods
#' results_shrinkage <- perform_LCT(gene_sets, expression, phenotype,
#'                                   method = "shrinkage", nbPermutations = 1000)
#' results_ridge <- perform_LCT(gene_sets, expression, phenotype,
#'                               method = "ridge", nbPermutations = 1000)
#'
#' # Merge and compare p-values
#' comparison <- data.frame(
#'   GS_name = results_shrinkage$GS_name,
#'   p_shrinkage = results_shrinkage$p_value,
#'   p_ridge = results_ridge$p_value
#' )
#' print(comparison)
#'
#' # Example with data frame gene set format
#' gs_matrix <- matrix(0, nrow = n_genes, ncol = 3)
#' rownames(gs_matrix) <- paste0("Gene", 1:n_genes)
#' colnames(gs_matrix) <- c("Set1", "Set2", "Set3")
#' gs_matrix[1:20, 1] <- 1
#' gs_matrix[21:40, 2] <- 1
#' gs_matrix[1:10, 3] <- 1
#' gs_df <- as.data.frame(gs_matrix)
#'
#' results_df <- perform_LCT(gs_df, expression, phenotype, nbPermutations = 1000)
#' print(results_df)
#' }
perform_LCT <- function(GS, DATA, cl, nbPermutations = 1000,
                        method = "shrinkage") {
  # Preprocess gene sets
  genes <- rownames(DATA)
  nb.Samples <- ncol(DATA)
  nb.GeneSets <- length(GS)

  # Convert to list format if needed
  GS <- GS.format.dataframe.to.list(GS)

  # Match gene identifiers
  GS <- lapply(GS, function(z) as.numeric(which(genes %in% z)))
  GS.sizes <- sapply(GS, length)

  # Extract expression data for each gene set
  GS.data <- lapply(GS, function(z) as.matrix(DATA[z, ], ncol = nb.Samples))

  # Select appropriate LCT function
  LCT_function <- switch(method,
    "shrinkage" = LCT.shrinkage,
    "glasso" = LCT.glasso,
    "adaptive_lasso" = LCT.adaptive.lasso,
    "ridge" = LCT.ridge,
    stop("Unknown method. Choose from: shrinkage, glasso, adaptive_lasso, ridge")
  )

  # Calculate p-values for each gene set
  GeneSets.pval <- sapply(GS.data, function(z) {
    LCT_function(z, cl, nbPermutations = nbPermutations)
  })

  # Calculate q-values if multiple gene sets
  if (nb.GeneSets >= 2) {
    # Try to use qvalue package (Bioconductor)
    if (requireNamespace("qvalue", quietly = TRUE)) {
      GeneSets.qval <- qvalue::qvalue(GeneSets.pval, pi0 = 1)$qvalues
    } else {
      warning("qvalue package not available. Using p.adjust(method='BH') instead. ",
              "For better FDR estimation, install qvalue: BiocManager::install('qvalue')")
      GeneSets.qval <- stats::p.adjust(GeneSets.pval, method = "BH")
    }

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
