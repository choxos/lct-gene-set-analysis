# Quick Start Guide

## 5-Minute Setup

### 1. Install Required Packages

```r
# Install CRAN packages
install.packages(c(
  "MASS", "corpcor", "dplyr", "faux", 
  "glasso", "rags2ridges", "cvCovEst"
))

# Install Bioconductor packages (for qvalue)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
```

### 2. Download the Code

```bash
git clone https://github.com/sara-khademi/lct-gene-set-analysis.git
cd lct-gene-set-analysis
```

### 3. Load the Functions

```r
source("lct_simulation.R")
```

---

## Your First Analysis (3 steps)

### Step 1: Prepare Your Data

```r
# Load your expression data
# Rows = genes, Columns = samples
expression_data <- read.csv("your_expression_data.csv", row.names = 1)

# Create phenotype vector (0 = control, 1 = case)
phenotype <- c(rep(0, 20), rep(1, 20))  # 20 controls, 20 cases

# Define your gene sets
gene_sets <- list(
  Pathway_A = c("GENE1", "GENE2", "GENE3", "GENE4"),
  Pathway_B = c("GENE5", "GENE6", "GENE7"),
  Pathway_C = c("GENE1", "GENE8", "GENE9", "GENE10")
)
```

### Step 2: Run LCT

```r
results <- perform_LCT(
  GS = gene_sets,
  DATA = expression_data,
  cl = phenotype,
  nbPermutations = 1000,
  method = "shrinkage"  # Fast and reliable default
)
```

### Step 3: View Results

```r
# See all results
print(results)

# Find significant gene sets (p < 0.05)
significant <- results[results$p_value < 0.05, ]
print(significant)

# Save results
write.csv(results, "lct_results.csv", row.names = FALSE)
```

**That's it!** You've completed your first LCT analysis.

---

## Common Use Cases

### Case 1: I have microarray data from GEO

```r
# Load GEOquery if needed
# BiocManager::install("GEOquery")
library(GEOquery)

# Download and prepare data
gse <- getGEO("GSE12345")
expression_data <- exprs(gse[[1]])
phenotype <- pData(gse[[1]])$characteristics_ch1

# Define gene sets from MSigDB or your own
# Then run LCT as shown above
```

### Case 2: I want to compare different methods

```r
methods <- c("shrinkage", "ridge", "glasso", "adaptive_lasso")
all_results <- list()

for (method in methods) {
  cat("Testing method:", method, "\n")
  all_results[[method]] <- perform_LCT(
    GS = gene_sets,
    DATA = expression_data,
    cl = phenotype,
    nbPermutations = 1000,
    method = method
  )
}

# Combine and compare
library(dplyr)
comparison <- bind_rows(all_results, .id = "method")
print(comparison)
```

### Case 3: I have unbalanced design (many more cases than controls)

```r
# Use ridge method for better Type I error control
results <- perform_LCT(
  GS = gene_sets,
  DATA = expression_data,
  cl = phenotype,
  nbPermutations = 1000,
  method = "ridge"  # Recommended for unbalanced designs
)
```

### Case 4: I want to test power through simulation

```r
# Quick simulation to estimate power
power_test <- run_LCT_simulation(
  n_samples = 40,
  n_genes = 200,
  rho = 0.5,
  n_datasets = 100,
  method = "shrinkage",
  gene_sets = list(TestSet = paste0("Gene_", 1:50)),
  phenotype = c(rep(0, 20), rep(1, 20)),
  gamma = 1.0  # Effect size
)

# Calculate power
p_values <- sapply(power_test, function(x) x$p_value[1])
power <- mean(p_values < 0.05)
cat("Estimated power:", power, "\n")
```

---

## Frequently Asked Questions

### Q: How many permutations should I use?

**A:** 
- Exploratory analysis: 500-1000
- Publication quality: 5000-10000
- Rule of thumb: At least 100/α (e.g., 2000 for α=0.05)

### Q: Which method should I choose?

**A:** Quick decision tree:
```
Is your design unbalanced (e.g., 80% cases, 20% controls)?
├─ Yes → Use "ridge"
└─ No → Are your genes highly correlated?
    ├─ Yes → Use "shrinkage"
    └─ No → Use "shrinkage" or "ridge"
```

### Q: My analysis is too slow. What can I do?

**A:** Try these in order:
1. Reduce permutations temporarily (use 500 instead of 1000)
2. Switch to faster method (shrinkage > ridge > adaptive_lasso > glasso)
3. Filter gene sets by size (exclude very large sets >500 genes)
4. Run analysis in batches

### Q: All my p-values are 1. What's wrong?

**A:** Check these:
1. Phenotype labels are correct: `table(phenotype)`
2. Gene names match between data and gene sets
3. Data has sufficient variation: `summary(apply(expression_data, 1, sd))`
4. Data is properly normalized

### Q: How do I cite this work?

**A:** Include both the original LCT paper and this implementation:

```
Wang, X., Dinu, I., Liu, W., & Yasui, Y. (2011). 
Linear combination test for hierarchical gene set analysis. 
Statistical Applications in Genetics and Molecular Biology, 10(1), Article 13.

Khademioureh, Sara, et al. "Stability and Performance of Linear Combination Tests of Gene Set Enrichment for Multiple Covariance Estimators in Unbalanced Studies." bioRxiv (2025): 2025-01.
GitHub repository: https://github.com/sara-khademi/lct-gene-set-analysis
```

---

## Troubleshooting Checklist

Before asking for help, please check:

- [ ] All required packages are installed
- [ ] Gene names match between expression data and gene sets
- [ ] Phenotype vector length equals number of samples
- [ ] Expression data has genes as rows, samples as columns
- [ ] No missing values in data: `any(is.na(expression_data))`
- [ ] Data has been normalized
- [ ] You're using latest version of the code

---

## Need More Help?

1. **Check the Documentation**: See `DOCUMENTATION.md` for detailed function descriptions
2. **Run Examples**: See `examples.R` for working code snippets
3. **Read the Paper**: Original LCT methodology explained in Wang et al. (2011)
4. **Open an Issue**: [GitHub Issues](https://github.com/yourusername/lct-gene-set-analysis/issues)

---

## Next Steps

Once you're comfortable with basic usage:

1. Read the full [README](README.md) for comprehensive overview
2. Explore [DOCUMENTATION.md](DOCUMENTATION.md) for advanced parameters
3. Run simulations to understand power and Type I error
4. Try different covariance estimation methods
5. Integrate with your existing analysis pipeline

---

**Happy analyzing!** 🧬📊
