# lctGSA 0.1.0

## Initial Release

This is the first release of lctGSA, converting the original `lct-gene-set-analysis` scripts into a proper R package.

### Features

* Implemented four LCT covariance estimation methods:
  - `LCT.shrinkage()` - Schäfer-Strimmer shrinkage estimator (fast, default)
  - `LCT.ridge()` - Ridge regularization (recommended for unbalanced designs)
  - `LCT.glasso()` - Graphical lasso (sparse precision matrices)
  - `LCT.adaptive.lasso()` - Adaptive lasso (feature selection)

* Added `perform_LCT()` wrapper function for easy analysis of multiple gene sets

* Created comprehensive simulation framework:
  - `simulate_gene_expression()` - Generate synthetic data
  - `run_LCT_simulation()` - Power analysis and method evaluation

* Published vignettes and documentation:
  - Getting Started guide
  - Method selection recommendations
  - Comparison tutorials
  - Simulation examples

* Established pkgdown website at https://sara-khademi.github.io/lctGSA/

### Statistical Methods

* Permutation-based significance testing (default: 1000 permutations)
* FDR correction for multiple gene sets using qvalue or p.adjust
* High-dimensional capability (handles p > n scenarios)
* Support for both balanced and unbalanced study designs

### Data Formats

* Accepts gene sets as lists or binary matrices
* Compatible with standard gene expression matrices
* Flexible phenotype specification (numeric, factor, character)

### References

* Wang, X., Dinu, I., Liu, W., & Yasui, Y. (2011). Linear combination test for
  hierarchical gene set analysis. *Statistical Applications in Genetics and
  Molecular Biology*, 10(1). doi: 10.2202/1544-6115.1618

* Khademioureh, S., et al. (2025). Stability and Performance of Linear Combination
  Tests of Gene Set Enrichment for Multiple Covariance Estimators in Unbalanced
  Studies. *bioRxiv*.
