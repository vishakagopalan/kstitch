test_that("prioritize_cca_components returns Shape and Factor data frames", {

  set.seed(1)
  n <- 200; k <- 3
  Z         <- matrix(rnorm(n * k), n, k)
  shape_mat <- scale(Z %*% matrix(rnorm(k * 8), k, 8) + matrix(rnorm(n * 8, sd = 0.3), n, 8))
  expr_mat  <- scale(Z %*% matrix(rnorm(k * 6), k, 6) + matrix(rnorm(n * 6, sd = 0.3), n, 6))
  colnames(shape_mat) <- paste0("ShapePC", seq_len(8))
  colnames(expr_mat)  <- paste0("Factor",  seq_len(6))
  rownames(shape_mat) <- rownames(expr_mat) <- paste0("cell_", seq_len(n))

  res <- prioritize_cca_components(
    shape_mat             = shape_mat,
    expr_mat              = expr_mat,
    ccs_to_consider       = 1:2,
    num_replicates        = 20,       # low for test speed
    num_times_among_top_ranks = 10    # low for test speed
  )

  expect_named(res, c("Shape", "Factor"))
  expect_true(is.data.frame(res$Shape))
  expect_true(is.data.frame(res$Factor))

  # expected columns
  expect_true(all(c("feature_name", "CC", "cross_loading_mag", "n") %in% colnames(res$Shape)))
  expect_true(all(c("feature_name", "CC", "cross_loading_mag", "n") %in% colnames(res$Factor)))

  # feature names are from the right matrices
  if (nrow(res$Shape) > 0)  expect_true(all(res$Shape$feature_name  %in% colnames(shape_mat)))
  if (nrow(res$Factor) > 0) expect_true(all(res$Factor$feature_name %in% colnames(expr_mat)))
})

test_that("cell_subsample restricts the cell pool", {

  set.seed(2)
  n <- 200
  Z         <- matrix(rnorm(n * 3), n, 3)
  shape_mat <- scale(Z %*% matrix(rnorm(3 * 6), 3, 6) + matrix(rnorm(n * 6, sd = 0.3), n, 6))
  expr_mat  <- scale(Z %*% matrix(rnorm(3 * 5), 3, 5) + matrix(rnorm(n * 5, sd = 0.3), n, 5))
  colnames(shape_mat) <- paste0("ShapePC", seq_len(6))
  colnames(expr_mat)  <- paste0("Factor",  seq_len(5))
  rownames(shape_mat) <- rownames(expr_mat) <- paste0("cell_", seq_len(n))

  sub <- paste0("cell_", seq_len(120))

  # should run without error using only the subsampled cells
  expect_no_error(
    prioritize_cca_components(
      shape_mat             = shape_mat,
      expr_mat              = expr_mat,
      ccs_to_consider       = 1,
      num_replicates        = 10,
      num_times_among_top_ranks = 5,
      cell_subsample        = sub
    )
  )
})
