# 📁 Repository Structure Guide

## File Organization

```
lct-gene-set-analysis/
│
├── 📄 README.md                      [7.3 KB] ⭐ START HERE
│   └── Main project overview, installation, and usage
│
├── 🚀 QUICKSTART.md                  [5.8 KB] ⭐ BEGINNERS
│   └── Get started in 5 minutes
│
├── 📚 DOCUMENTATION.md                [12 KB] ⭐ REFERENCE
│   └── Complete function reference and parameter guides
│
├── 📜 LICENSE                         [1.1 KB]
│   └── MIT License
│
├── 🔧 gitignore                       [900 B]
│   └── Git ignore patterns
│
├── 💻 lct_simulation.R                [15 KB] ⭐ MAIN CODE
│   │
│   ├── Helper Functions
│   │   ├── GS.format.dataframe.to.list()
│   │   ├── T2.like.SAMGS()
│   │   └── generate_gene_names()
│   │
│   ├── Core LCT Functions
│   │   ├── LCT.shrinkage()          # Schäfer-Strimmer estimator
│   │   ├── LCT.ridge()              # Ridge regularization
│   │   ├── LCT.glasso()             # Graphical Lasso
│   │   └── LCT.adaptive.lasso()     # Adaptive Lasso
│   │
│   ├── Wrapper Function
│   │   └── perform_LCT()            # Easy-to-use interface
│   │
│   └── Simulation Functions
│       ├── simulate_gene_expression()
│       ├── run_LCT_simulation()
│       └── run_complete_simulation()
│
└── 📝 examples.R                      [7.0 KB] ⭐ LEARN BY DOING
    │
    ├── Example 1: Basic LCT Analysis
    ├── Example 2: Compare All Methods
    ├── Example 3: Visualize Results
    ├── Example 4: Small Simulation Study
    └── Example 5: Working with Gene Set Databases
```

---



## Function Map

### Core Analysis Pipeline

```
Your Data
    ↓
[perform_LCT]  ← Main wrapper function
    ↓
├─ GS.format.dataframe.to.list()    # Format gene sets
├─ Preprocess data                   # Extract genes, sizes
└─ Loop through gene sets
        ↓
    [Choose Method]
        ↓
    ├─ LCT.shrinkage()       # Fast, default
    ├─ LCT.ridge()           # Unbalanced designs
    ├─ LCT.glasso()          # Sparse networks
    └─ LCT.adaptive.lasso()  # Feature selection
            ↓
        [Covariance Estimation]
            ↓
        [Eigenvalue Decomposition]
            ↓
        [Transform Data]
            ↓
        [T2.like.SAMGS]      # Test statistic
            ↓
        [Permutation Test]
            ↓
        P-value
            ↓
Results DataFrame
    ↓
Save / Analyze
```

---

## Workflow Examples

### Basic Workflow
```
1. Load data              (your_data.csv)
2. Define phenotype       (0s and 1s)
3. Define gene sets       (list of genes)
4. source("lct_simulation.R")
5. results <- perform_LCT(GS, DATA, cl)
6. View results
7. Save results
```

### Comparison Workflow
```
1. Load data
2. methods <- c("shrinkage", "ridge", "glasso", "adaptive_lasso")
3. results_list <- lapply(methods, perform_LCT)
4. Compare p-values
5. Select best method
6. Save comparison
```

### Simulation Workflow
```
1. Set parameters          (n, rho, gamma)
2. Define gene sets
3. run_LCT_simulation()
4. Calculate power
5. Adjust parameters
6. Repeat
```

---

## Key Functions Summary

| Function | Purpose | Speed | Best For |
|----------|---------|-------|----------|
| `perform_LCT()` | Main interface | Depends | All analyses |
| `LCT.shrinkage()` | Shrinkage estimator | ⚡⚡⚡ | Default choice |
| `LCT.ridge()` | Ridge estimator | ⚡⚡⚡ | Unbalanced designs |
| `LCT.glasso()` | Graphical lasso | ⚡ | Network structure |
| `LCT.adaptive.lasso()` | Adaptive lasso | ⚡⚡ | Feature selection |

Speed: ⚡ = Slow, ⚡⚡ = Medium, ⚡⚡⚡ = Fast

---

## Data Flow

### Input Requirements

```
Expression Data (DATA)
├── Format: Matrix or DataFrame
├── Dimensions: genes × samples
├── Row names: Gene IDs
└── Values: Normalized expression

Gene Sets (GS)
├── Format: List or DataFrame
├── List: list(Set1 = c("G1", "G2"), ...)
└── DataFrame: genes × sets, binary (0/1)

Phenotype (cl)
├── Format: Numeric or Factor
├── Length: Number of samples
└── Values: Binary (0/1 or two levels)
```

### Output Format

```
Results DataFrame
├── GS_name     : Gene set identifier
├── GS_size     : Number of genes
├── p_value     : Permutation p-value
├── q_value     : FDR-adjusted (if multiple sets)
└── method      : Covariance estimator used
```


---

## Recommended Reading Order

1. `README.md` (overview)
2. `DOCUMENTATION.md` (details)
3. `lct_simulation.R` (code)
4. `examples.R` (adapt to your needs)

---

## Color-Coded Priority

🔴 **Must Read** (Everyone)
- README.md
- QUICKSTART.md

🟡 **Should Read** (Most Users)
- DOCUMENTATION.md
- examples.R

🟢 **Optional** (Power Users)
- lct_simulation.R (full code)
- CHANGELOG.md
- PROJECT_SUMMARY.md

---

## File Dependencies

```
examples.R
    ↓ requires
lct_simulation.R
    ↓ requires
Packages: MASS, corpcor, qvalue, dplyr, faux, glasso, 
          cvCovEst, rags2ridges
```

---

## Version Information

```
Current Version: 1.0.0
Release Date: October 2025
Last Updated: October 23, 2025
R Version Required: ≥ 4.0
```

---

## Quick Access URLs

```
Main Page:
https://github.com/sara-khademi/lct-gene-set-analysis

Quick Start:
https://github.com/sara-khademi/lct-gene-set-analysis/blob/main/QUICKSTART.md

Documentation:
https://github.com/sara-khademi/lct-gene-set-analysis/blob/main/DOCUMENTATION.md

Examples:
https://github.com/sara-khademi/lct-gene-set-analysis/blob/main/examples.R

Issues:
https://github.com/sara-khademi/lct-gene-set-analysis/issues
```

---

