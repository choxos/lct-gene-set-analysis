# LCT Function Documentation

## Table of Contents
1. [Core Functions](#core-functions)
2. [Helper Functions](#helper-functions)
3. [Simulation Functions](#simulation-functions)
4. [Parameter Guidelines](#parameter-guidelines)
5. [Troubleshooting](#troubleshooting)

---

## Core Functions

### `perform_LCT()`

Main wrapper function for performing Linear Combination Test on multiple gene sets.

**Usage:**
```r
perform_LCT(GS, DATA, cl, nbPermutations = 1000, method = "shrinkage")
```

**Arguments:**
- `GS`: Gene sets in one of two formats:
  - Data frame: genes (rows) × gene sets (columns), binary matrix (0/1)
  - List: each element contains vector of gene identifiers
- `DATA`: Expression data matrix or data frame (genes × samples)
- `cl`: Phenotype vector (binary: 0/1 or factor with 2 levels)
- `nbPermutations`: Number of permutations for p-value estimation (default: 1000)
- `method`: Covariance estimation method
  - `"shrinkage"`: Schäfer-Strimmer shrinkage estimator (default)
  - `"ridge"`: Ridge regularization
  - `"glasso"`: Graphical Lasso (L1-penalized)
  - `"adaptive_lasso"`: Adaptive Lasso

**Returns:**
Data frame with columns:
- `GS_name`: Gene set identifier
- `GS_size`: Number of genes in set
- `p_value`: Permutation-based p-value
- `q_value`: FDR-adjusted q-value (if multiple sets tested)
- `method`: Method used

**Example:**
```r
results <- perform_LCT(
  GS = my_gene_sets,
  DATA = expression_matrix,
  cl = c(rep(0, 20), rep(1, 20)),
  nbPermutations = 1000,
  method = "shrinkage"
)
```

---

### `LCT.shrinkage()`

Linear Combination Test using shrinkage covariance estimator.

**Usage:**
```r
LCT.shrinkage(DATA, cl, nbPermutations = 1000, s0 = 0)
```

**Arguments:**
- `DATA`: Expression matrix for single gene set (genes × samples)
- `cl`: Phenotype vector
- `nbPermutations`: Number of permutations (default: 1000)
- `s0`: Small constant for numerical stability (default: 0)

**Returns:**
Numeric p-value

**Details:**
- Uses optimal shrinkage intensity automatically determined
- Variance shrinkage intensity set to 0 (empirical variances recovered)
- Correlation shrinkage optimized via formula from Schäfer & Strimmer (2005)

**Example:**
```r
p_value <- LCT.shrinkage(
  DATA = gene_set_data,
  cl = phenotype,
  nbPermutations = 1000,
  s0 = 0
)
```

---

### `LCT.ridge()`

Linear Combination Test using ridge covariance estimator.

**Usage:**
```r
LCT.ridge(DATA, cl, nbPermutations = 1000, s0 = 0, lambda = 1)
```

**Arguments:**
- `DATA`: Expression matrix for single gene set
- `cl`: Phenotype vector
- `nbPermutations`: Number of permutations (default: 1000)
- `s0`: Small constant for stability (default: 0)
- `lambda`: Ridge penalty parameter (default: 1)

**Returns:**
Numeric p-value

**Details:**
- Ridge estimation adds L2 penalty to covariance matrix
- `lambda = 1` provides good balance for most applications
- Larger λ → more regularization → more conservative
- Smaller λ → less regularization → closer to sample covariance

**Recommended λ values:**
- Small gene sets (<50): 0.1 - 0.5
- Medium gene sets (50-200): 0.5 - 1.0
- Large gene sets (>200): 1.0 - 5.0

**Example:**
```r
p_value <- LCT.ridge(
  DATA = gene_set_data,
  cl = phenotype,
  nbPermutations = 1000,
  lambda = 1
)
```

---

### `LCT.glasso()`

Linear Combination Test using graphical lasso (sparse precision matrix).

**Usage:**
```r
LCT.glasso(DATA, cl, nbPermutations = 1000, s0 = 0, rho = 0.4)
```

**Arguments:**
- `DATA`: Expression matrix for single gene set
- `cl`: Phenotype vector
- `nbPermutations`: Number of permutations (default: 1000)
- `s0`: Small constant for stability (default: 0)
- `rho`: L1 penalty parameter (default: 0.4)

**Returns:**
Numeric p-value

**Details:**
- Estimates sparse inverse covariance (precision) matrix
- Good for identifying conditional independence structures
- Computationally intensive for large gene sets
- Higher ρ → more sparsity

**Recommended ρ values:**
- Sparse networks: 0.1 - 0.3
- Moderate sparsity: 0.3 - 0.5
- Dense networks: 0.5 - 0.7

**Warning:** Very small ρ (<0.1) increases computation time significantly

**Example:**
```r
p_value <- LCT.glasso(
  DATA = gene_set_data,
  cl = phenotype,
  nbPermutations = 1000,
  rho = 0.4
)
```

---

### `LCT.adaptive.lasso()`

Linear Combination Test using adaptive lasso covariance estimator.

**Usage:**
```r
LCT.adaptive.lasso(DATA, cl, nbPermutations = 1000, s0 = 0.5, 
                   lambda = 0.5, n = 0.8)
```

**Arguments:**
- `DATA`: Expression matrix for single gene set
- `cl`: Phenotype vector
- `nbPermutations`: Number of permutations (default: 1000)
- `s0`: Small constant for stability (default: 0.5)
- `lambda`: Penalty parameter (default: 0.5)
- `n`: Adaptive weight exponent (default: 0.8)

**Returns:**
Numeric p-value

**Details:**
- Adaptive weights reduce bias for large coefficients
- Maintains sparsity while improving estimation
- `n` controls adaptive weight strength (typically 0.5-1.0)

**Example:**
```r
p_value <- LCT.adaptive.lasso(
  DATA = gene_set_data,
  cl = phenotype,
  nbPermutations = 1000,
  lambda = 0.5,
  n = 0.8
)
```

---

## Helper Functions

### `GS.format.dataframe.to.list()`

Converts gene set format from data frame to list.

**Usage:**
```r
GS.format.dataframe.to.list(GS)
```

**Arguments:**
- `GS`: Gene set as data frame or list

**Returns:**
List of gene sets

**Example:**
```r
# Input: data frame
gs_df <- data.frame(
  Set1 = c(1, 1, 0, 0),
  Set2 = c(0, 1, 1, 0)
)
rownames(gs_df) <- c("Gene1", "Gene2", "Gene3", "Gene4")

# Convert
gs_list <- GS.format.dataframe.to.list(gs_df)
# Result: list(Set1 = c("Gene1", "Gene2"), Set2 = c("Gene2", "Gene3"))
```

---

### `T2.like.SAMGS()`

Calculates T²-like test statistic (sum of squared mean differences).

**Usage:**
```r
T2.like.SAMGS(DATA, cl)
```

**Arguments:**
- `DATA`: Expression matrix (genes × samples)
- `cl`: Phenotype vector

**Returns:**
Numeric test statistic

---

### `generate_gene_names()`

Generates random gene identifiers.

**Usage:**
```r
generate_gene_names(n)
```

**Arguments:**
- `n`: Number of gene names to generate

**Returns:**
Character vector of unique gene names

**Example:**
```r
genes <- generate_gene_names(100)
# Returns: c("ABCDEFG0123A", "HIJKLMN4567B", ...)
```

---

## Simulation Functions

### `simulate_gene_expression()`

Generates simulated gene expression data with specified correlation structure.

**Usage:**
```r
simulate_gene_expression(rho, n_samples, n_genes, gamma, 
                         n_controls, n_cases)
```

**Arguments:**
- `rho`: Correlation between genes (0-1)
- `n_samples`: Total number of samples
- `n_genes`: Number of genes to simulate
- `gamma`: Mean difference between groups
- `n_controls`: Number of control samples
- `n_cases`: Number of case samples

**Returns:**
Matrix of simulated expression values (genes × samples)

**Example:**
```r
expr_data <- simulate_gene_expression(
  rho = 0.5,
  n_samples = 40,
  n_genes = 100,
  gamma = 1.0,
  n_controls = 20,
  n_cases = 20
)
```

---

### `run_LCT_simulation()`

Runs multiple iterations of LCT simulation study.

**Usage:**
```r
run_LCT_simulation(n_samples, n_genes, rho, n_datasets, method, 
                   gene_sets, phenotype, gamma)
```

**Arguments:**
- `n_samples`: Sample size
- `n_genes`: Number of genes
- `rho`: Gene correlation
- `n_datasets`: Number of simulation iterations
- `method`: Covariance estimation method
- `gene_sets`: Predefined gene sets to test
- `phenotype`: Binary phenotype vector
- `gamma`: Mean difference parameter

**Returns:**
List of results from each simulation iteration

---

### `run_complete_simulation()`

Executes full factorial simulation study across all parameter combinations.

**Usage:**
```r
run_complete_simulation(output_dir = ".")
```

**Arguments:**
- `output_dir`: Directory to save results (default: current directory)

**Details:**
Tests combinations of:
- Sample sizes: 20, 50
- Correlations: 0.1, 0.9
- Outcome probabilities: 0.2, 0.8
- Methods: all four covariance estimators
- Effect sizes: γ from 0 to 2 by 0.1

**Output:**
CSV files for each parameter combination

---

## Parameter Guidelines

### Number of Permutations

| Application | Recommended | Rationale |
|-------------|-------------|-----------|
| Exploratory analysis | 100-500 | Faster, rough estimates |
| Standard analysis | 1000 | Balance of accuracy/speed |
| Publication | 5000-10000 | High precision |
| Multiple testing | 10000+ | Precise tail probabilities |

**Rule of thumb:** Minimum permutations = 100 / α
- For α = 0.05: ≥ 2000 permutations
- For α = 0.01: ≥ 10000 permutations

### Penalty Parameter Selection

#### Ridge (λ)
- **Default**: 1.0
- **Cross-validation**: Use `cv.glmnet()` from glmnet package
- **Rule of thumb**: λ ≈ trace(Σ) / p where p = number of genes

#### Graphical Lasso (ρ)
- **Default**: 0.4
- **Low**: 0.1-0.3 (less sparse, more connections)
- **High**: 0.5-0.8 (more sparse, fewer connections)
- **BIC selection**: Use `glassopath()` with BIC criterion

#### Adaptive Lasso (λ, n)
- **λ default**: 0.5
- **n default**: 0.8
- **Recommendation**: n ∈ [0.5, 1.0] for oracle properties

---

## Troubleshooting

### Common Issues

#### 1. "Error: covariance matrix is singular"

**Cause:** Gene set size exceeds sample size with standard estimation

**Solutions:**
- Use regularized methods (recommended: ridge or shrinkage)
- Increase sample size if possible
- Reduce gene set size
- Increase penalty parameter

#### 2. "Error in glasso(): did not converge"

**Cause:** Graphical lasso failed to converge

**Solutions:**
- Increase `rho` parameter (more regularization)
- Check for constant genes (zero variance)
- Try different initialization
- Use alternative method (ridge or shrinkage)

#### 3. Computational time too long

**Cause:** Large gene sets with glasso, or too many permutations

**Solutions:**
- Reduce number of permutations (start with 500)
- Use faster method (shrinkage or ridge)
- Parallelize across gene sets
- Filter gene sets by size

#### 4. All p-values are 1.0

**Causes:**
- No differential expression
- Incorrect phenotype labels
- Data not properly normalized
- Wrong gene identifiers

**Solutions:**
- Check phenotype vector: `table(phenotype)`
- Verify gene identifiers match between data and gene sets
- Check for sufficient variation: `apply(DATA, 1, sd)`
- Ensure data is properly normalized

#### 5. Memory issues with large datasets

**Solutions:**
- Process gene sets in batches
- Increase system RAM
- Use sparse matrix representation
- Filter low-variance genes

---

## Best Practices

1. **Data Preprocessing:**
   - Normalize expression data before analysis
   - Filter low-expression genes
   - Check for batch effects
   - Verify phenotype labels

2. **Method Selection:**
   - Start with shrinkage (fastest, good default)
   - Use ridge for unbalanced designs
   - Reserve glasso for network analysis
   - Adaptive lasso for feature selection

3. **Multiple Testing:**
   - Always report q-values when testing multiple gene sets
   - Consider Bonferroni correction for small number of tests
   - Use FDR control for large-scale studies

4. **Validation:**
   - Compare results across methods
   - Check consistency with biological knowledge
   - Validate significant findings experimentally
   - Report method and parameters used

5. **Reproducibility:**
   - Set random seed: `set.seed(123)`
   - Document R version and package versions
   - Save session info: `sessionInfo()`
   - Archive code and parameters

---

## References

1. Wang, Xiaoming, et al. "Linear combination test for hierarchical gene set analysis." *Statistical Applications in Genetics & Molecular Biology* 10.1 (2011).

2. Schäfer, Juliane, and Korbinian Strimmer. "A shrinkage approach to large-scale covariance matrix estimation and implications for functional genomics." *Statistical applications in genetics and molecular biology* 4.1 (2005).

3. Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. "Sparse inverse covariance estimation with the graphical lasso." *Biostatistics* 9.3 (2008): 432-441.

4. Rothman, Adam J., Elizaveta Levina, and Ji Zhu. "Generalized thresholding of large covariance matrices." *Journal of the American Statistical Association* 104.485 (2009): 177-186.

5. Van Wieringen, Wessel N., and Carel FW Peeters. "Ridge estimation of inverse covariance matrices from high-dimensional data." *Computational Statistics & Data Analysis* 103 (2016): 284-303.
