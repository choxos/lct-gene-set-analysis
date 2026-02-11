test_that("LCT.shrinkage returns valid p-value", {
  set.seed(123)

  n_genes <- 20
  n_samples <- 16
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  phenotype <- c(rep(0, 8), rep(1, 8))

  p_value <- LCT.shrinkage(expression, phenotype, nbPermutations = 100)

  expect_type(p_value, "double")
  expect_length(p_value, 1)
  expect_true(p_value >= 0 && p_value <= 1)
})

test_that("LCT.shrinkage is reproducible with set seed", {
  n_genes <- 25
  n_samples <- 18
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  phenotype <- c(rep(0, 9), rep(1, 9))

  set.seed(456)
  p1 <- LCT.shrinkage(expression, phenotype, nbPermutations = 100)

  set.seed(456)
  p2 <- LCT.shrinkage(expression, phenotype, nbPermutations = 100)

  expect_equal(p1, p2)
})

test_that("LCT.shrinkage detects differential expression", {
  set.seed(789)

  n_genes <- 30
  n_samples <- 20
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  phenotype <- c(rep(0, 10), rep(1, 10))

  # Add strong signal
  expression[, phenotype == 1] <- expression[, phenotype == 1] + 2.5

  p_value <- LCT.shrinkage(expression, phenotype, nbPermutations = 500)

  expect_true(p_value < 0.05)  # Should detect strong signal
})

test_that("LCT.ridge returns valid p-value", {
  set.seed(234)

  n_genes <- 20
  n_samples <- 16
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  phenotype <- c(rep(0, 8), rep(1, 8))

  p_value <- LCT.ridge(expression, phenotype, nbPermutations = 100, lambda = 1)

  expect_type(p_value, "double")
  expect_true(p_value >= 0 && p_value <= 1)
})

test_that("LCT.glasso returns valid p-value", {
  set.seed(345)

  n_genes <- 15
  n_samples <- 16
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  phenotype <- c(rep(0, 8), rep(1, 8))

  p_value <- LCT.glasso(expression, phenotype, nbPermutations = 50, rho = 0.4)

  expect_type(p_value, "double")
  expect_true(p_value >= 0 && p_value <= 1)
})

test_that("LCT.adaptive.lasso returns valid p-value", {
  set.seed(456)

  n_genes <- 15
  n_samples <- 16
  expression <- matrix(rnorm(n_genes * n_samples), nrow = n_genes)
  phenotype <- c(rep(0, 8), rep(1, 8))

  p_value <- LCT.adaptive.lasso(expression, phenotype, nbPermutations = 50)

  expect_type(p_value, "double")
  expect_true(p_value >= 0 && p_value <= 1)
})

test_that("All LCT methods handle single gene", {
  set.seed(567)

  n_samples <- 14
  expression <- matrix(rnorm(n_samples), nrow = 1)
  phenotype <- c(rep(0, 7), rep(1, 7))

  p_shrinkage <- LCT.shrinkage(expression, phenotype, nbPermutations = 50)
  p_ridge <- LCT.ridge(expression, phenotype, nbPermutations = 50)
  p_glasso <- LCT.glasso(expression, phenotype, nbPermutations = 50)
  p_adaptive <- LCT.adaptive.lasso(expression, phenotype, nbPermutations = 50)

  expect_true(all(c(p_shrinkage, p_ridge, p_glasso, p_adaptive) >= 0))
  expect_true(all(c(p_shrinkage, p_ridge, p_glasso, p_adaptive) <= 1))
})
