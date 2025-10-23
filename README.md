# Linear Combination Test (LCT) for Gene Set Analysis

[![R](https://img.shields.io/badge/R-%3E%3D%204.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)

## Overview

This repository implements the **Linear Combination Test (LCT)** for gene set analysis in high-dimensional gene expression data, with support for multiple covariance matrix estimators. LCT is designed to identify differentially expressed gene sets between phenotypic groups while accounting for the complex correlation structure among genes.

## Features

- **Multiple Covariance Estimators**: 
  - Shrinkage (Schäfer & Strimmer, 2005)
  - Ridge regularization (Van Wieringen and Peeters, 2016)
  - Graphical Lasso (sparse precision matrix) (Friedman et al., 2008)
  - Adaptive Lasso (Rothman et al., 2009)

- **Robust to High-Dimensional Data**: Handles cases where the number of genes exceeds the number of samples

- **Permutation-Based Testing**: Non-parametric p-value estimation via phenotype label permutation

- **Simulation Framework**: Complete simulation study implementation for method comparison

## Background

Gene Set Analysis (GSA) evaluates whether predefined groups of genes show coordinated differential expression between phenotypic conditions. LCT tests the null hypothesis:

> H₀: No linear combination of genes in the set is associated with the phenotype

The method:
1. Estimates the covariance matrix of gene expressions using regularization
2. Performs eigenvalue decomposition to obtain orthogonal bases
3. Transforms data to uncorrelated features
4. Computes test statistic as L2-norm of mean differences
5. Evaluates significance via permutation testing

## Installation

### Prerequisites

Required R packages:
```r
install.packages(c(
  "MASS",
  "corpcor",
  "qvalue",
  "dplyr",
  "faux",
  "glasso",
  "cvCovEst",
  "rags2ridges"
))
```

### Download

```bash
git clone https://github.com/yourusername/lct-gene-set-analysis.git
cd lct-gene-set-analysis
```

## Usage

### Basic Example

```r
# Source the functions
source("lct_simulation.R")

# Prepare your data
# DATA: genes (rows) x samples (columns)
# phenotype: binary vector indicating group membership (0/1)
# gene_sets: list of gene identifiers for each set

# Run LCT with shrinkage estimator
results <- perform_LCT(
  GS = gene_sets,
  DATA = expression_data,
  cl = phenotype,
  nbPermutations = 1000,
  method = "shrinkage"
)

# View results
head(results)
```

### Method Comparison

```r
# Compare different covariance estimators
methods <- c("shrinkage", "ridge", "glasso", "adaptive_lasso")

comparison_results <- lapply(methods, function(m) {
  perform_LCT(
    GS = gene_sets,
    DATA = expression_data,
    cl = phenotype,
    nbPermutations = 1000,
    method = m
  )
})

names(comparison_results) <- methods
```

### Running Simulations

```r
# Run complete simulation study
# Results will be saved to ./simulation_results/
run_complete_simulation(output_dir = "./simulation_results")
```

## Input Format

### Expression Data
- **Format**: Matrix or data frame
- **Dimensions**: Genes (rows) × Samples (columns)
- **Row names**: Gene identifiers
- **Values**: Normalized expression values

Example:
```
         Sample1  Sample2  Sample3  Sample4
Gene1    5.23     4.87     6.12     5.45
Gene2    3.45     3.21     3.67     3.89
Gene3    7.89     8.12     7.45     7.67
```

### Gene Sets
- **Format**: Data frame or list
- **Data frame**: Genes (rows) × Gene Sets (columns), binary (0/1)
- **List**: Each element is a vector of gene identifiers

Example (data frame):
```
       Pathway1  Pathway2  Pathway3
Gene1    1         0         1
Gene2    1         1         0
Gene3    0         1         1
```

Example (list):
```r
gene_sets <- list(
  Pathway1 = c("Gene1", "Gene2", "Gene5"),
  Pathway2 = c("Gene2", "Gene3", "Gene4"),
  Pathway3 = c("Gene1", "Gene3", "Gene6")
)
```

### Phenotype
- **Format**: Numeric or factor vector
- **Length**: Equal to number of samples
- **Values**: Binary (0/1 or two-level factor)

Example:
```r
phenotype <- c(0, 0, 0, 1, 1, 1, 1, 0, 0, 1)
```

## Output

The `perform_LCT` function returns a data frame with:

| Column    | Description                                      |
|-----------|--------------------------------------------------|
| GS_name   | Gene set name/identifier                         |
| GS_size   | Number of genes in the set                       |
| p_value   | Permutation-based p-value                        |
| q_value   | FDR-adjusted q-value (if multiple sets)          |
| method    | Covariance estimation method used                |

## Simulation Study Parameters

The simulation framework allows testing under various scenarios:

- **Sample sizes**: 20, 50
- **Correlations**: 0.1 (low), 0.9 (high)
- **Outcome probabilities**: 0.2, 0.8 (unbalanced designs)
- **Mean differences** (γ): 0 to 2 by 0.1
- **Iterations**: 1000 per scenario

## Performance Considerations

- **Computational time** varies by method:
  - Shrinkage: Fastest (~0.5-1 sec per gene set)
  - Ridge: Fast (~0.5-1 sec per gene set)
  - Adaptive Lasso: Moderate (~1-2 sec per gene set)
  - Graphical Lasso: Slower for large gene sets (~2-5 sec per gene set)

- **Memory**: Scales with gene set size; datasets with >5000 genes per set may require substantial RAM

## Method Selection Guide

| Scenario | Recommended Method | Rationale |
|----------|-------------------|-----------|
| Balanced design, high correlation | Shrinkage | Optimal power and speed |
| Unbalanced design, high correlation | Ridge | Better Type I error control |
| Low correlation | Shrinkage or Ridge | Both perform well |
| Small sample size (n<30) | Ridge | Most stable |
| Large gene sets (>200 genes) | Avoid Glasso | Computational issues |

## Citation

If you use this code in your research, please cite:

Khademioureh, Sara, et al. "Stability and Performance of Linear Combination Tests of Gene Set Enrichment for Multiple Covariance Estimators in Unbalanced Studies." bioRxiv (2025): 2025-01.

Original LCT method:
```
Wang, Xiaoming, et al. "Linear combination test for hierarchical gene set analysis." Statistical Applications in Genetics & Molecular Biology 10.1 (2011).
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes:

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## Issues

Please report any bugs or feature requests through the [GitHub Issues](../../issues) page.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- Dr. Xiaoming Wang for the original LCT methodology
- Reviewers whose feedback improved the implementation
- Contributors to the R packages used in this project

## Contact

Dr. Sara Khademioureh - [sarakhademi@ualberta.ca]

Project Link: [https://github.com/sara-khademi/lct-gene-set-analysis](https://github.com/sara-khademi/lct-gene-set-analysis)

---

**Keywords**: Gene Set Analysis, Linear Combination Test, High-Dimensional Data, Covariance Estimation, Genomics, Bioinformatics, RNA-seq, Microarray
