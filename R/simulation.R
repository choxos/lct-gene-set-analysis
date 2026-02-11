#' Run LCT Simulation Study
#'
#' @description
#' Performs multiple iterations of LCT analysis on simulated data to evaluate
#' statistical power, Type I error control, and method performance under
#' controlled conditions. Useful for power analysis, sample size planning,
#' and method comparison.
#'
#' @param n_samples Sample size (total number of samples across both groups).
#' @param n_genes Number of genes to simulate in each dataset.
#' @param rho Correlation coefficient between genes (0-1).
#' @param n_datasets Number of simulation iterations to run.
#' @param method Covariance estimation method: "shrinkage", "ridge", "glasso",
#'   or "adaptive_lasso".
#' @param gene_sets Predefined gene sets to test. Can be list or data frame
#'   format (see \code{\link{perform_LCT}} for details).
#' @param phenotype Binary phenotype vector of length n_samples defining group
#'   membership (0 = control, 1 = case).
#' @param gamma Mean difference parameter (effect size). Added to case group.
#'
#' @return List of length n_datasets, where each element is a data frame of
#'   results from \code{\link{perform_LCT}} for that iteration. Each data frame
#'   contains columns: GS_name, GS_size, p_value, q_value (if applicable), method.
#'
#' @details
#' For each iteration:
#' \enumerate{
#'   \item Generate expression data using \code{\link{simulate_gene_expression}}
#'   \item Apply LCT to predefined gene sets using \code{\link{perform_LCT}}
#'   \item Store results
#'   \item Handle errors gracefully (skip iteration if error occurs)
#' }
#'
#' Progress is printed to console showing current iteration.
#'
#' **Typical use cases**:
#' \itemize{
#'   \item Power analysis: Run with gamma > 0, calculate proportion of p < 0.05
#'   \item Type I error: Run with gamma = 0, verify proportion of p < 0.05 ≈ 0.05
#'   \item Method comparison: Run multiple methods, compare power and error rates
#' }
#'
#' @note
#' \itemize{
#'   \item Set \code{set.seed()} before calling for reproducible results
#'   \item Large n_datasets (>1000) may take considerable time
#'   \item Failed iterations are skipped with a warning message
#' }
#'
#' @seealso
#' \code{\link{simulate_gene_expression}} for data generation,
#' \code{\link{perform_LCT}} for LCT analysis
#'
#' @export
#' @examples
#' # Small power analysis (few iterations for speed)
#' set.seed(999)
#'
#' # Set up parameters
#' n_genes <- 1000
#' n_samples <- 30
#' phenotype <- c(rep(0, 15), rep(1, 15))
#'
#' # Create gene set
#' gene_names <- generate_gene_names(n_genes)
#' gene_set <- list(TestSet = gene_names[1:50])
#'
#' \donttest{
#' # Run simulation with medium effect
#' results <- run_LCT_simulation(
#'   n_samples = n_samples,
#'   n_genes = n_genes,
#'   rho = 0.3,
#'   n_datasets = 10,  # Use 100+ for real analysis
#'   method = "shrinkage",
#'   gene_sets = gene_set,
#'   phenotype = phenotype,
#'   gamma = 1.0
#' )
#'
#' # Extract p-values
#' p_values <- sapply(results, function(x) x$p_value[1])
#'
#' # Calculate power (proportion of significant results)
#' power <- mean(p_values < 0.05, na.rm = TRUE)
#' cat("Estimated power:", power, "\n")
#'
#' # Type I error control simulation (gamma = 0)
#' null_results <- run_LCT_simulation(
#'   n_samples = n_samples,
#'   n_genes = n_genes,
#'   rho = 0.3,
#'   n_datasets = 10,
#'   method = "shrinkage",
#'   gene_sets = gene_set,
#'   phenotype = phenotype,
#'   gamma = 0  # No effect
#' )
#'
#' null_p_values <- sapply(null_results, function(x) x$p_value[1])
#' type1_error <- mean(null_p_values < 0.05, na.rm = TRUE)
#' cat("Type I error rate:", type1_error, "\n")  # Should be ≈ 0.05
#' }
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

#' Run Complete Factorial Simulation Study
#'
#' @description
#' Internal function to run a comprehensive factorial simulation study testing
#' all combinations of sample sizes, correlations, methods, and effect sizes.
#' This is a computationally intensive function intended for method development
#' and validation, not routine analysis.
#'
#' @param output_dir Directory to save CSV result files. Default is current
#'   directory. Results are saved as separate CSV files for each parameter
#'   combination.
#'
#' @return NULL (invisibly). Results are written to CSV files in output_dir.
#'   Each file contains results from n_iterations datasets for one parameter
#'   combination.
#'
#' @details
#' This function runs a full factorial simulation with:
#' \itemize{
#'   \item Sample sizes: 20, 50
#'   \item Correlations: 0.1 (weak), 0.9 (strong)
#'   \item Methods: shrinkage, ridge, glasso, adaptive_lasso
#'   \item Outcome probabilities: 0.2, 0.8 (balanced vs unbalanced)
#'   \item Effect sizes (gamma): 0 to 2.0 by 0.1 increments
#'   \item 1000 iterations per scenario
#'   \item 10000 genes per dataset
#' }
#'
#' **Total runs**: 2 × 2 × 4 × 2 × 21 × 1000 = 672,000 LCT analyses
#'
#' **Estimated runtime**: 50-100 hours (depending on hardware)
#'
#' CSV files are named: \code{Simulation_n[N]_p[P]_rho[R]_method_[M]_gamma[G].csv}
#'
#' @keywords internal
#' @noRd
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
    gene_set_matrix[, i] <- stats::rbinom(n_genes, 1, 0.002)
  }

  gene_sets <- data.frame(gene_set_matrix)
  rownames(gene_sets) <- gene_names

  # Run simulations
  for (n_samples in sample_sizes) {
    for (prob in outcome_probs) {
      # Generate phenotype
      phenotype <- stats::rbinom(n_samples, 1, prob)

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

            results_df <- dplyr::bind_rows(results, .id = "iteration")
            utils::write.csv(results_df, filename, row.names = FALSE)

            elapsed <- difftime(Sys.time(), start_time, units = "mins")
            cat("Completed in", round(elapsed, 2), "minutes\n")
          }
        }
      }
    }
  }

  cat("\nSimulation complete!\n")
  invisible(NULL)
}
