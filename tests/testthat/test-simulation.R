test_that("generate_gene_names creates unique names", {
  n <- 100
  gene_names <- generate_gene_names(n)

  expect_length(gene_names, n)
  expect_type(gene_names, "character")
  expect_equal(length(unique(gene_names)), n)  # All unique
})

test_that("simulate_gene_expression returns correct dimensions", {
  set.seed(123)

  n_genes <- 50
  n_controls <- 10
  n_cases <- 12
  n_samples <- n_controls + n_cases

  expr_data <- simulate_gene_expression(
    rho = 0.3,
    n_samples = n_samples,
    n_genes = n_genes,
    gamma = 1.0,
    n_controls = n_controls,
    n_cases = n_cases
  )

  expect_equal(nrow(expr_data), n_genes)
  expect_equal(ncol(expr_data), n_samples)
  expect_true(!is.null(rownames(expr_data)))  # Has gene names
  expect_true(all(expr_data >= 0))  # Non-negative (abs values)
})

test_that("simulate_gene_expression creates mean difference", {
  set.seed(456)

  n_genes <- 100
  n_controls <- 15
  n_cases <- 15
  gamma <- 1.5

  expr_data <- simulate_gene_expression(
    rho = 0.2,
    n_samples = n_controls + n_cases,
    n_genes = n_genes,
    gamma = gamma,
    n_controls = n_controls,
    n_cases = n_cases
  )

  # Check that cases have higher mean (approximately)
  # Get row means for each group
  control_means <- rowMeans(expr_data[, 1:n_controls, drop = FALSE])
  case_means <- rowMeans(expr_data[, (n_controls + 1):(n_controls + n_cases), drop = FALSE])

  # Overall mean should be higher for cases due to gamma > 0
  expect_true(mean(case_means) > mean(control_means))
})

test_that("run_LCT_simulation returns list of results", {
  set.seed(789)

  n_samples <- 20
  n_genes <- 100
  phenotype <- c(rep(0, 10), rep(1, 10))

  gene_names <- generate_gene_names(n_genes)
  gene_sets <- list(TestSet = gene_names[1:30])

  # Run small simulation
  results <- run_LCT_simulation(
    n_samples = n_samples,
    n_genes = n_genes,
    rho = 0.3,
    n_datasets = 3,  # Just 3 iterations for speed
    method = "shrinkage",
    gene_sets = gene_sets,
    phenotype = phenotype,
    gamma = 0.5
  )

  expect_type(results, "list")
  expect_true(length(results) <= 3)  # Up to 3 iterations
  expect_s3_class(results[[1]], "data.frame")
  expect_true("p_value" %in% names(results[[1]]))
})

test_that("simulate_gene_expression handles different correlations", {
  set.seed(321)

  n_genes <- 40
  n_samples <- 20

  # Low correlation
  data_low <- simulate_gene_expression(
    rho = 0.1,
    n_samples = n_samples,
    n_genes = n_genes,
    gamma = 0,
    n_controls = 10,
    n_cases = 10
  )

  # High correlation
  data_high <- simulate_gene_expression(
    rho = 0.9,
    n_samples = n_samples,
    n_genes = n_genes,
    gamma = 0,
    n_controls = 10,
    n_cases = 10
  )

  expect_equal(dim(data_low), dim(data_high))
  expect_true(all(data_low >= 0))
  expect_true(all(data_high >= 0))
})
