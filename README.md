# lctGSA: Linear Combination Test for Gene Set Analysis

<!-- badges: start -->
[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/SKhademi/lctGSA/workflows/R-CMD-check/badge.svg)](https://github.com/SKhademi/lctGSA/actions)
<!-- badges: end -->

## Overview

**lctGSA** implements the Linear Combination Test (LCT) for identifying differentially expressed gene sets in high-dimensional gene expression data. The package provides multiple covariance matrix estimators optimized for different experimental designs, including unbalanced group sizes.

### Key Features

✨ **Four covariance estimation methods**
- Shrinkage (fast, default)
- Ridge (unbalanced designs)
- Graphical Lasso (sparse networks)
- Adaptive Lasso (feature selection)

🚀 **High-dimensional capability** - Handles cases where genes >> samples

📊 **Permutation-based testing** - Non-parametric significance evaluation

⚖️ **Unbalanced design support** - Robust methods for imbalanced group sizes

🧬 **FDR correction** - Automatic q-value calculation for multiple gene sets

## Installation

Install from GitHub:

```r
# Install devtools if needed
install.packages("devtools")

# Install lctGSA
devtools::install_github("SKhademi/lctGSA")
```

### Dependencies

The package requires several CRAN packages which will be installed automatically. The Bioconductor package `qvalue` is optional but recommended for better FDR estimation:

```r
# Install qvalue (optional, recommended)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
```

## Quick Start

```r
library(lctGSA)

# Load your data
# expression_data: genes × samples matrix with row names as gene IDs
# phenotype: binary vector (0/1) of length = number of samples
# gene_sets: list of gene sets, e.g., list(Pathway1 = c("Gene1", "Gene2", ...))

# Run LCT analysis
results <- perform_LCT(
  GS = gene_sets,
  DATA = as.data.frame(expression_data),
  cl = phenotype,
  nbPermutations = 1000,
  method = "shrinkage"  # or "ridge", "glasso", "adaptive_lasso"
)

# View results
print(results)
#   GS_name GS_size  p_value  q_value    method
# Pathway1      50    0.001    0.003 shrinkage
# Pathway2      75    0.045    0.068 shrinkage
# Pathway3      42    0.523    0.523 shrinkage

# Identify significant gene sets
significant <- results[results$q_value < 0.05, ]
```

## Example with Simulated Data

```r
set.seed(123)

# Simulate expression data (500 genes × 40 samples)
n_genes <- 500
n_samples <- 40
expression <- matrix(rnorm(n_genes * n_samples, mean = 5, sd = 2),
                      nrow = n_genes, ncol = n_samples)
rownames(expression) <- paste0("Gene_", 1:n_genes)

# Add differential expression to first 50 genes
phenotype <- c(rep(0, 20), rep(1, 20))
expression[1:50, phenotype == 1] <- expression[1:50, phenotype == 1] + 2

# Define gene sets
gene_sets <- list(
  Pathway_A = paste0("Gene_", 1:60),    # Contains DE genes
  Pathway_B = paste0("Gene_", 40:100),  # Partially DE
  Pathway_C = paste0("Gene_", 200:250)  # No DE
)

# Run analysis
results <- perform_LCT(gene_sets, as.data.frame(expression), phenotype,
                        nbPermutations = 1000)
print(results)
```

## Method Selection Guide

### When to Use Each Method

| Scenario | Recommended Method | Reason |
|----------|-------------------|---------|
| Balanced designs (equal groups) | `shrinkage` | Optimal power, fastest |
| Unbalanced designs (e.g., 20 vs 5) | `ridge` | Better Type I error control |
| Network/pathway analysis | `glasso` | Identifies sparse structures |
| Feature selection | `adaptive_lasso` | Oracle properties |

### Performance

| Method | Speed | Typical Use Case |
|--------|-------|------------------|
| Shrinkage | 0.5-1 sec/gene set | Default choice |
| Ridge | 0.5-1 sec/gene set | Unbalanced groups |
| Adaptive Lasso | 1-2 sec/gene set | Feature selection |
| Graphical Lasso | 4-8 sec/gene set | Network analysis |

## Documentation

- **Getting Started**: `vignette("lctGSA")`
- **Method Selection**: `vignette("method-selection")`
- **Function Help**: `?perform_LCT`, `?LCT.shrinkage`
- **Website**: https://SKhademi.github.io/lctGSA/

## Advanced Features

### Individual LCT Methods

For more control, use individual LCT functions:

```r
# Shrinkage estimator
p_value <- LCT.shrinkage(expression_subset, phenotype, nbPermutations = 1000)

# Ridge with custom lambda
p_value <- LCT.ridge(expression_subset, phenotype,
                      nbPermutations = 1000, lambda = 5)

# Graphical Lasso with custom rho
p_value <- LCT.glasso(expression_subset, phenotype,
                       nbPermutations = 1000, rho = 0.6)
```

### Simulation Studies

```r
# Simulate gene expression data
expr_data <- simulate_gene_expression(
  rho = 0.5,           # Correlation between genes
  n_samples = 40,
  n_genes = 1000,
  gamma = 1.0,         # Effect size
  n_controls = 20,
  n_cases = 20
)

# Run power analysis
results <- run_LCT_simulation(
  n_samples = 40,
  n_genes = 1000,
  rho = 0.5,
  n_datasets = 100,
  method = "shrinkage",
  gene_sets = my_gene_sets,
  phenotype = phenotype,
  gamma = 1.0
)

# Calculate power
p_values <- sapply(results, function(x) x$p_value[1])
power <- mean(p_values < 0.05)
```

## Citation

If you use lctGSA in your research, please cite:

Khademioureh, S., et al. (2025). Stability and Performance of Linear Combination Tests of Gene Set Enrichment for Multiple Covariance Estimators in Unbalanced Studies. *bioRxiv*.

Wang, X., Dinu, I., Liu, W., & Yasui, Y. (2011). Linear combination test for hierarchical gene set analysis. *Statistical Applications in Genetics and Molecular Biology*, 10(1). https://doi.org/10.2202/1544-6115.1618

## Getting Help

- 🐛 **Bug reports**: [GitHub Issues](https://github.com/SKhademi/lctGSA/issues)
- 💬 **Questions**: Use GitHub Discussions
- 📖 **Documentation**: https://SKhademi.github.io/lctGSA/

## References

**Covariance Estimation Methods:**

- Schäfer, J., & Strimmer, K. (2005). A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics. *Statistical Applications in Genetics and Molecular Biology*, 4(1). https://doi.org/10.2202/1544-6115.1175

- Van Wieringen, W. N., & Peeters, C. F. (2016). Ridge estimation of inverse covariance matrices from high-dimensional data. *Computational Statistics & Data Analysis*, 103, 284-303. https://doi.org/10.1016/j.csda.2016.05.012

- Friedman, J., Hastie, T., & Tibshirani, R. (2008). Sparse inverse covariance estimation with the graphical lasso. *Biostatistics*, 9(3), 432-441. https://doi.org/10.1093/biostatistics/kxm045

- Rothman, A. J., Levina, E., & Zhu, J. (2009). Generalized thresholding of large covariance matrices. *Journal of the American Statistical Association*, 104(485), 177-186. https://doi.org/10.1198/jasa.2009.0101

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Authors

- Dr. Sara Khademioureh ([@SKhademi](https://github.com/SKhademi))
- Dr. Payam Amini
- Dr. Erfan Ghasemi
- Paul Calistrate-Petre
- Dr. Saumyadipta Pyne
- Dr. Irina Dinu
