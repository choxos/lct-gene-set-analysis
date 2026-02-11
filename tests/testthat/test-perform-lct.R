test_that("perform_LCT returns correct structure", {
  set.seed(100)

  # Create test data
  n_genes <- 50
  n_samples <- 20
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(expression) <- paste0("G", 1:n_genes)
  phenotype <- c(rep(0, 10), rep(1, 10))

  # Create gene sets
  gene_sets <- list(
    Set1 = paste0("G", 1:20),
    Set2 = paste0("G", 21:40)
  )

  # Run LCT
  result <- perform_LCT(gene_sets, as.data.frame(expression), phenotype,
                         nbPermutations = 100, method = "shrinkage")

  # Check output structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true(all(c("GS_name", "GS_size", "p_value", "q_value", "method") %in% names(result)))
  expect_equal(result$method, rep("shrinkage", 2))
})

test_that("perform_LCT works with all methods", {
  set.seed(200)

  n_genes <- 40
  n_samples <- 16
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(expression) <- paste0("Gene", 1:n_genes)
  phenotype <- c(rep(0, 8), rep(1, 8))

  gene_sets <- list(TestSet = paste0("Gene", 1:20))

  methods <- c("shrinkage", "ridge", "glasso", "adaptive_lasso")

  for (method in methods) {
    result <- perform_LCT(gene_sets, as.data.frame(expression), phenotype,
                           nbPermutations = 50, method = method)

    expect_equal(result$method[1], method)
    expect_true(result$p_value[1] >= 0 && result$p_value[1] <= 1)
  }
})

test_that("perform_LCT handles single gene set", {
  set.seed(300)

  n_genes <- 30
  n_samples <- 14
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(expression) <- paste0("G", 1:n_genes)
  phenotype <- c(rep(0, 7), rep(1, 7))

  gene_sets <- list(SingleSet = paste0("G", 1:15))

  result <- perform_LCT(gene_sets, as.data.frame(expression), phenotype,
                         nbPermutations = 50)

  expect_equal(nrow(result), 1)
  expect_false("q_value" %in% names(result))  # No q-value for single set
})

test_that("perform_LCT calculates q-values for multiple sets", {
  set.seed(400)

  n_genes <- 60
  n_samples <- 20
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(expression) <- paste0("Gene", 1:n_genes)
  phenotype <- c(rep(0, 10), rep(1, 10))

  gene_sets <- list(
    Set1 = paste0("Gene", 1:20),
    Set2 = paste0("Gene", 21:40),
    Set3 = paste0("Gene", 41:60)
  )

  result <- perform_LCT(gene_sets, as.data.frame(expression), phenotype,
                         nbPermutations = 100)

  expect_true("q_value" %in% names(result))
  expect_true(all(result$q_value >= result$p_value))  # Q-values ≥ p-values
})

test_that("perform_LCT detects differential expression", {
  set.seed(500)

  n_genes <- 50
  n_samples <- 20
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(expression) <- paste0("Gene", 1:n_genes)
  phenotype <- c(rep(0, 10), rep(1, 10))

  # Add strong differential expression
  expression[1:20, phenotype == 1] <- expression[1:20, phenotype == 1] + 3

  gene_sets <- list(
    DE_Set = paste0("Gene", 1:20),     # Differentially expressed
    NonDE_Set = paste0("Gene", 30:50)  # Not differentially expressed
  )

  result <- perform_LCT(gene_sets, as.data.frame(expression), phenotype,
                         nbPermutations = 500)

  # DE set should have small p-value
  de_pval <- result$p_value[result$GS_name == "DE_Set"]
  non_de_pval <- result$p_value[result$GS_name == "NonDE_Set"]

  expect_true(de_pval < 0.1)  # Should be significant
  expect_true(de_pval < non_de_pval)  # DE set should be more significant
})

test_that("perform_LCT handles data frame gene set format", {
  set.seed(600)

  n_genes <- 40
  n_samples <- 16
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(expression) <- paste0("G", 1:n_genes)
  phenotype <- c(rep(0, 8), rep(1, 8))

  # Create gene set as data frame (binary matrix)
  gs_matrix <- matrix(0, nrow = n_genes, ncol = 2)
  rownames(gs_matrix) <- paste0("G", 1:n_genes)
  colnames(gs_matrix) <- c("GeneSet1", "GeneSet2")
  gs_matrix[1:20, 1] <- 1
  gs_matrix[21:40, 2] <- 1
  gene_sets_df <- as.data.frame(gs_matrix)

  result <- perform_LCT(gene_sets_df, as.data.frame(expression), phenotype,
                         nbPermutations = 50)

  expect_equal(nrow(result), 2)
  expect_true(all(c("GeneSet1", "GeneSet2") %in% result$GS_name))
})

test_that("perform_LCT errors on invalid method", {
  set.seed(700)

  n_genes <- 30
  n_samples <- 14
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  rownames(expression) <- paste0("Gene", 1:n_genes)
  phenotype <- c(rep(0, 7), rep(1, 7))
  gene_sets <- list(Set = paste0("Gene", 1:15))

  expect_error(
    perform_LCT(gene_sets, as.data.frame(expression), phenotype,
                 method = "invalid_method"),
    "Unknown method"
  )
})
