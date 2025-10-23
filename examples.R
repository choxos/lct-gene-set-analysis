# ==============================================================================
# Example: Using LCT for Gene Set Analysis
# ==============================================================================

# Load required libraries
source("lct_simulation.R")

# ==============================================================================
# Example 1: Basic LCT Analysis with Real Data
# ==============================================================================

# Simulate example data (replace with your own data)
set.seed(123)

n_genes <- 500
n_samples <- 40
n_controls <- 20
n_cases <- 20

# Create expression data matrix
expression_data <- matrix(
  rnorm(n_genes * n_samples, mean = 5, sd = 2),
  nrow = n_genes,
  ncol = n_samples
)

# Add differential expression to first 50 genes
expression_data[1:50, (n_controls + 1):n_samples] <- 
  expression_data[1:50, (n_controls + 1):n_samples] + 2

# Set gene and sample names
rownames(expression_data) <- paste0("Gene_", 1:n_genes)
colnames(expression_data) <- paste0("Sample_", 1:n_samples)

# Create phenotype vector
phenotype <- c(rep(0, n_controls), rep(1, n_cases))

# Define gene sets
gene_sets <- list(
  Pathway_A = paste0("Gene_", 1:60),
  Pathway_B = paste0("Gene_", 40:100),
  Pathway_C = paste0("Gene_", 200:250)
)

# ==============================================================================
# Run LCT with different methods
# ==============================================================================

cat("\n=== Running LCT with Shrinkage Estimator ===\n")
results_shrinkage <- perform_LCT(
  GS = gene_sets,
  DATA = as.data.frame(expression_data),
  cl = phenotype,
  nbPermutations = 1000,
  method = "shrinkage"
)
print(results_shrinkage)

cat("\n=== Running LCT with Ridge Estimator ===\n")
results_ridge <- perform_LCT(
  GS = gene_sets,
  DATA = as.data.frame(expression_data),
  cl = phenotype,
  nbPermutations = 1000,
  method = "ridge"
)
print(results_ridge)

# ==============================================================================
# Example 2: Compare All Methods
# ==============================================================================

cat("\n=== Comparing All Covariance Estimation Methods ===\n")

methods <- c("shrinkage", "ridge", "glasso", "adaptive_lasso")
results_list <- list()

for (method in methods) {
  cat("\nTesting method:", method, "\n")
  
  tryCatch({
    results_list[[method]] <- perform_LCT(
      GS = gene_sets,
      DATA = as.data.frame(expression_data),
      cl = phenotype,
      nbPermutations = 500,  # Reduced for faster example
      method = method
    )
  }, error = function(e) {
    cat("Error with", method, ":", conditionMessage(e), "\n")
  })
}

# Combine results
if (length(results_list) > 0) {
  combined_results <- bind_rows(results_list, .id = "method")
  
  # Save results
  write.csv(combined_results, "lct_comparison_results.csv", row.names = FALSE)
  cat("\nResults saved to: lct_comparison_results.csv\n")
  
  # Display summary
  cat("\n=== Summary of Results ===\n")
  print(combined_results)
}

# ==============================================================================
# Example 3: Visualize Results
# ==============================================================================

# If you have ggplot2 installed
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  # Plot p-values comparison
  p <- ggplot(combined_results, 
              aes(x = GS_name, y = -log10(as.numeric(p_value)), 
                  fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    theme_minimal() +
    labs(
      title = "Comparison of LCT Methods",
      subtitle = "Lower p-values indicate stronger association",
      x = "Gene Set",
      y = "-log10(p-value)",
      fill = "Method"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  print(p)
  
  # Save plot
  ggsave("lct_comparison_plot.png", plot = p, width = 10, height = 6, dpi = 300)
  cat("\nPlot saved to: lct_comparison_plot.png\n")
}

# ==============================================================================
# Example 4: Small Simulation Study
# ==============================================================================

cat("\n=== Running Small Simulation Study ===\n")

# Parameters
n_simulations <- 100
gamma_values <- c(0, 0.5, 1.0, 1.5)

simulation_results <- list()

for (gamma in gamma_values) {
  cat("\nTesting gamma =", gamma, "\n")
  
  sim_results <- run_LCT_simulation(
    n_samples = 30,
    n_genes = 200,
    rho = 0.5,
    n_datasets = n_simulations,
    method = "shrinkage",
    gene_sets = list(TestPathway = paste0("Gene_", 1:50)),
    phenotype = c(rep(0, 15), rep(1, 15)),
    gamma = gamma
  )
  
  # Calculate power (proportion of significant results)
  p_values <- sapply(sim_results, function(x) {
    if (!is.null(x)) as.numeric(x$p_value[1]) else NA
  })
  
  power <- mean(p_values < 0.05, na.rm = TRUE)
  
  simulation_results[[as.character(gamma)]] <- data.frame(
    gamma = gamma,
    power = power,
    mean_p = mean(p_values, na.rm = TRUE),
    se_p = sd(p_values, na.rm = TRUE) / sqrt(sum(!is.na(p_values)))
  )
}

# Combine and display
power_analysis <- bind_rows(simulation_results)
print(power_analysis)

# Save power analysis
write.csv(power_analysis, "lct_power_analysis.csv", row.names = FALSE)
cat("\nPower analysis saved to: lct_power_analysis.csv\n")

# ==============================================================================
# Example 5: Working with Gene Set Databases (KEGG, GO, etc.)
# ==============================================================================

cat("\n=== Example: Loading Gene Sets from File ===\n")

# Example function to read gene sets from GMT format (common for MSigDB, KEGG)
read_gmt_file <- function(gmt_file) {
  lines <- readLines(gmt_file)
  gene_sets <- list()
  
  for (line in lines) {
    parts <- strsplit(line, "\t")[[1]]
    set_name <- parts[1]
    genes <- parts[-(1:2)]  # Skip name and description
    gene_sets[[set_name]] <- genes
  }
  
  return(gene_sets)
}

# Example: Create a sample GMT file
sample_gmt <- c(
  "Pathway1\tDescription1\tGene_1\tGene_2\tGene_3\tGene_4",
  "Pathway2\tDescription2\tGene_5\tGene_6\tGene_7",
  "Pathway3\tDescription3\tGene_1\tGene_8\tGene_9\tGene_10"
)

writeLines(sample_gmt, "sample_genesets.gmt")

# Read gene sets
gene_sets_from_file <- read_gmt_file("sample_genesets.gmt")
cat("\nLoaded", length(gene_sets_from_file), "gene sets from file\n")

# Clean up
file.remove("sample_genesets.gmt")

# ==============================================================================
# Summary
# ==============================================================================

cat("\n" , rep("=", 70), "\n")
cat("Example analysis complete!\n")
cat("Files created:\n")
cat("  - lct_comparison_results.csv\n")
cat("  - lct_power_analysis.csv\n")
if (requireNamespace("ggplot2", quietly = TRUE)) {
  cat("  - lct_comparison_plot.png\n")
}
cat(rep("=", 70), "\n")
